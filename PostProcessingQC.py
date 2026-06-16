#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Post-Processing QC Tool for Microglia 2D Preprocessing Macro
=============================================================
Verifies that output TIFF files from the Fiji preprocessing macro
meet expected specifications and flags anomalies.
 
Compatible with Python 2.7 and Python 3.x.
No external packages required -- uses only the standard library.
 
How to run:
    python PostProcessingQC.py
 
Checks performed per file:
    - Image dimensions (width x height)
    - Channel count
    - Z-slice count (should be 1 for 2D output)
    - Timepoint count (should be 1)
    - Bit depth
    - Per-channel intensity statistics (min, max, mean, std dev)
    - Detection of blank (all-zero) or saturated channels
    - Voxel calibration presence
    - Filename convention (-IBA1-DAPI.tiff suffix)
"""
 
from __future__ import print_function, division
 
import os
import sys
import struct
import csv
import math
import array
import threading
from datetime import datetime
 
try:
    import Queue as queue_mod        # Python 2
except ImportError:
    import queue as queue_mod        # Python 3
 
# ---- Cross-version Tkinter imports ----
try:
    import Tkinter as tk
    import tkFileDialog as filedialog
    import tkMessageBox as messagebox
    import ScrolledText as scrolledtext_mod
    ScrolledText = scrolledtext_mod.ScrolledText
    import ttk
except ImportError:
    import tkinter as tk
    import tkinter.filedialog as filedialog
    import tkinter.messagebox as messagebox
    import tkinter.scrolledtext as scrolledtext_mod
    ScrolledText = scrolledtext_mod.ScrolledText
    import tkinter.ttk as ttk
 
 
# ==========================================================================
#  TIFF PARSER (focused on QC needs: metadata + pixel data)
# ==========================================================================
 
# TIFF data-type definitions: id -> (name, byte_size, struct_format)
TIFF_TYPES = {
    1:  ("BYTE",      1, "B"),
    2:  ("ASCII",     1, "s"),
    3:  ("SHORT",     2, "H"),
    4:  ("LONG",      4, "L"),
    5:  ("RATIONAL",  8, None),
    6:  ("SBYTE",     1, "b"),
    7:  ("UNDEFINED", 1, "s"),
    8:  ("SSHORT",    2, "h"),
    9:  ("SLONG",     4, "l"),
    10: ("SRATIONAL", 8, None),
    11: ("FLOAT",     4, "f"),
    12: ("DOUBLE",    8, "d"),
}
 
# Minimal set of tag IDs needed for QC
TAG_IMAGE_WIDTH       = 256
TAG_IMAGE_LENGTH      = 257
TAG_BITS_PER_SAMPLE   = 258
TAG_COMPRESSION       = 259
TAG_IMAGE_DESCRIPTION = 270
TAG_STRIP_OFFSETS     = 273
TAG_SAMPLES_PER_PIXEL = 277
TAG_ROWS_PER_STRIP    = 278
TAG_STRIP_BYTE_COUNTS = 279
TAG_X_RESOLUTION      = 282
TAG_Y_RESOLUTION      = 283
TAG_RESOLUTION_UNIT   = 296
TAG_SAMPLE_FORMAT     = 339
TAG_TILE_OFFSETS      = 324
TAG_TILE_BYTE_COUNTS  = 325
 
COMPRESSION_NAMES = {
    1: "Uncompressed", 5: "LZW", 7: "JPEG",
    8: "Adobe Deflate", 32773: "PackBits",
}
 
RESOLUTION_UNIT_NAMES = {1: "No unit", 2: "Inch", 3: "Centimeter"}
UNIT_TO_UM = {2: 25400.0, 3: 10000.0}
 
 
class TiffQCParser(object):
    """
    Lightweight TIFF parser for QC.  Reads IFD metadata and raw pixel
    data (uncompressed only) using the standard library.
    """
 
    def __init__(self, filepath):
        self.filepath = filepath
        self.endian = "<"
        self.is_bigtiff = False
        self.ifds = []           # list of {tag_id: value}
        self._data = None
        self._parse()
 
    # ---- low-level helpers ----
 
    def _read(self, fmt, offset):
        full = self.endian + fmt
        size = struct.calcsize(full)
        return struct.unpack(full, self._data[offset:offset + size])
 
    def _read_tag_value(self, type_id, count, value_offset):
        tinfo = TIFF_TYPES.get(type_id)
        if tinfo is None:
            return None
 
        type_name, unit_size, fmt = tinfo
        total_size = unit_size * count
 
        if total_size <= 4:
            data_offset = value_offset
        else:
            data_offset = self._read("L", value_offset)[0]
 
        # ASCII
        if type_id == 2:
            raw = self._data[data_offset:data_offset + count]
            if isinstance(raw, bytes):
                return raw.rstrip(b"\x00").decode("latin-1", "replace")
            return raw.rstrip("\x00")
 
        # UNDEFINED
        if type_id == 7:
            return self._data[data_offset:data_offset + count]
 
        # RATIONAL / SRATIONAL
        if type_id in (5, 10):
            values = []
            long_fmt = "L" if type_id == 5 else "l"
            for i in range(count):
                off = data_offset + i * 8
                num = self._read(long_fmt, off)[0]
                den = self._read(long_fmt, off + 4)[0]
                values.append((num, den))
            return values[0] if count == 1 else values
 
        # Numeric types
        values = []
        for i in range(count):
            off = data_offset + i * unit_size
            val = self._read(fmt, off)[0]
            values.append(val)
        return values[0] if count == 1 else values
 
    # ---- main parse ----
 
    def _parse(self):
        with open(self.filepath, "rb") as f:
            self._data = f.read()
 
        if len(self._data) < 8:
            raise ValueError("File too small to be a TIFF")
 
        bom = self._data[0:2]
        if bom == b"II":
            self.endian = "<"
        elif bom == b"MM":
            self.endian = ">"
        else:
            raise ValueError("Not a TIFF file")
 
        magic = self._read("H", 2)[0]
        if magic == 42:
            self.is_bigtiff = False
            ifd_offset = self._read("L", 4)[0]
        elif magic == 43:
            self.is_bigtiff = True
            ifd_offset = self._read("Q", 8)[0]
        else:
            raise ValueError("Unknown TIFF magic number: %d" % magic)
 
        visited = set()
        while ifd_offset != 0 and ifd_offset < len(self._data):
            if ifd_offset in visited:
                break
            visited.add(ifd_offset)
 
            tags = {}
            if self.is_bigtiff:
                n_entries = self._read("Q", ifd_offset)[0]
                entry_start = ifd_offset + 8
                entry_size = 20
            else:
                n_entries = self._read("H", ifd_offset)[0]
                entry_start = ifd_offset + 2
                entry_size = 12
 
            for i in range(n_entries):
                eoff = entry_start + i * entry_size
                if self.is_bigtiff:
                    tag_id = self._read("H", eoff)[0]
                    type_id = self._read("H", eoff + 2)[0]
                    count = self._read("Q", eoff + 4)[0]
                    val_offset = eoff + 12
                else:
                    tag_id = self._read("H", eoff)[0]
                    type_id = self._read("H", eoff + 2)[0]
                    count = self._read("L", eoff + 4)[0]
                    val_offset = eoff + 8
 
                try:
                    value = self._read_tag_value(type_id, count, val_offset)
                except Exception:
                    value = None
                tags[tag_id] = value
 
            self.ifds.append(tags)
 
            next_off_pos = entry_start + n_entries * entry_size
            if self.is_bigtiff:
                ifd_offset = self._read("Q", next_off_pos)[0]
            else:
                ifd_offset = self._read("L", next_off_pos)[0]
 
    # ---- metadata accessors ----
 
    def get_tag(self, ifd_index, tag_id, default=None):
        if ifd_index < len(self.ifds):
            return self.ifds[ifd_index].get(tag_id, default)
        return default
 
    def get_imagej_metadata(self):
        """Parse ImageJ-style ImageDescription from IFD 0."""
        desc = self.get_tag(0, TAG_IMAGE_DESCRIPTION, "")
        if not desc or not desc.startswith("ImageJ"):
            return None
        result = {}
        for line in desc.strip().split("\n"):
            line = line.strip()
            if "=" in line:
                key, _, value = line.partition("=")
                result[key.strip()] = value.strip()
            elif line and not result:
                result["_header"] = line
        return result
 
    def get_dimensions(self, ifd_index=0):
        """Return (width, height, bits_per_sample, samples_per_pixel)."""
        width = self.get_tag(ifd_index, TAG_IMAGE_WIDTH, 0)
        height = self.get_tag(ifd_index, TAG_IMAGE_LENGTH, 0)
        bits = self.get_tag(ifd_index, TAG_BITS_PER_SAMPLE, 8)
        if isinstance(bits, (list, tuple)):
            bits = bits[0]
        spp = self.get_tag(ifd_index, TAG_SAMPLES_PER_PIXEL, 1)
        return width, height, bits, spp
 
    def get_voxel_size_um(self):
        """Compute XY voxel size in microns from TIFF resolution tags."""
        xres = self.get_tag(0, TAG_X_RESOLUTION)
        yres = self.get_tag(0, TAG_Y_RESOLUTION)
        res_unit = self.get_tag(0, TAG_RESOLUTION_UNIT, 2)
        scale = UNIT_TO_UM.get(res_unit)
 
        result = {}
        for axis, res in [("X", xres), ("Y", yres)]:
            if res is None or scale is None:
                result[axis] = None
                continue
            if isinstance(res, tuple) and len(res) == 2:
                num, den = res
                if den == 0 or num == 0:
                    result[axis] = None
                    continue
                pixels_per_unit = float(num) / den
            elif isinstance(res, (int, float)):
                pixels_per_unit = float(res)
            else:
                result[axis] = None
                continue
            result[axis] = scale / pixels_per_unit
        return result
 
    # ---- pixel data reading ----
 
    def read_channel_pixels(self, ifd_index):
        """
        Read pixel data from a specific IFD as a Python array.
 
        Returns an array.array of pixel values, or None if the data
        cannot be read (e.g. compressed TIFF, BigTIFF, tiled layout).
        """
        if ifd_index >= len(self.ifds):
            return None
 
        tags = self.ifds[ifd_index]
        compression = tags.get(TAG_COMPRESSION, 1)
        if compression != 1:
            return None   # only uncompressed supported
 
        if self.is_bigtiff:
            return None   # BigTIFF pixel reading not supported
 
        bits = tags.get(TAG_BITS_PER_SAMPLE, 8)
        if isinstance(bits, (list, tuple)):
            bits = bits[0]
 
        # Get strip locations
        strip_offsets = tags.get(TAG_STRIP_OFFSETS)
        strip_byte_counts = tags.get(TAG_STRIP_BYTE_COUNTS)
 
        if strip_offsets is None or strip_byte_counts is None:
            return None
 
        # Normalize to lists
        if not isinstance(strip_offsets, (list, tuple)):
            strip_offsets = [strip_offsets]
        if not isinstance(strip_byte_counts, (list, tuple)):
            strip_byte_counts = [strip_byte_counts]
 
        # Read all strips into a single bytes object
        raw_parts = []
        for offset, count in zip(strip_offsets, strip_byte_counts):
            end = offset + count
            if end > len(self._data):
                return None   # data truncated
            raw_parts.append(self._data[offset:end])
        raw = b"".join(raw_parts)
 
        # Build array of appropriate type
        if bits == 8:
            pixels = array.array("B")
        elif bits == 16:
            pixels = array.array("H")
        elif bits == 32:
            sample_fmt = tags.get(TAG_SAMPLE_FORMAT, 1)
            if isinstance(sample_fmt, (list, tuple)):
                sample_fmt = sample_fmt[0]
            if sample_fmt == 3:
                pixels = array.array("f")
            else:
                pixels = array.array("L")
        else:
            return None
 
        # Load bytes into array
        # Python 2: fromstring;  Python 3.2+: frombytes
        if sys.version_info[0] >= 3:
            pixels.frombytes(raw)
        else:
            pixels.fromstring(raw)
 
        # Fix endianness if file byte order differs from system
        if (self.endian == ">" and sys.byteorder == "little") or \
           (self.endian == "<" and sys.byteorder == "big"):
            pixels.byteswap()
 
        return pixels
 
 
# ==========================================================================
#  STATISTICS HELPERS
# ==========================================================================
 
def compute_channel_stats(pixels):
    """
    Compute intensity statistics for an array of pixel values.
 
    Returns a dict with min, max, mean, std, zero_fraction,
    and saturation_fraction.
    """
    n = len(pixels)
    if n == 0:
        return {
            "min": 0, "max": 0, "mean": 0.0, "std": 0.0,
            "zero_fraction": 1.0, "sat_fraction": 0.0,
        }
 
    # min, max, sum use C-optimized builtins on array objects
    min_val = min(pixels)
    max_val = max(pixels)
    total = sum(pixels)
    mean_val = total / n
 
    # Standard deviation (single-pass Welford for stability)
    # But since we already computed mean, two-pass is fine here
    sq_diff_sum = 0.0
    zero_count = 0
    sat_count = 0
 
    # Determine saturation ceiling from array type
    type_code = pixels.typecode
    if type_code == "B":
        sat_ceiling = 255
    elif type_code == "H":
        sat_ceiling = 65535
    else:
        sat_ceiling = max_val   # no saturation concept for float/32-bit
 
    for x in pixels:
        sq_diff_sum += (x - mean_val) ** 2
        if x == 0:
            zero_count += 1
        if x >= sat_ceiling:
            sat_count += 1
 
    std_val = math.sqrt(sq_diff_sum / n)
 
    return {
        "min": min_val,
        "max": max_val,
        "mean": round(mean_val, 2),
        "std": round(std_val, 2),
        "zero_fraction": round(zero_count / n, 4),
        "sat_fraction": round(sat_count / n, 4),
    }
 
 
# ==========================================================================
#  QC INSPECTION ENGINE
# ==========================================================================
 
# Status constants
PASS = "PASS"
WARN = "WARN"
FAIL = "FAIL"
 
 
def inspect_file(filepath, expected):
    """
    Run all QC checks on a single output file.
 
    Parameters
    ----------
    filepath : str
        Path to the TIFF file to inspect.
    expected : dict
        Expected values with keys: width, height, channels, slices,
        frames, bit_depth.
 
    Returns
    -------
    report_lines : list of (text, tag) tuples for GUI display
    csv_row : dict of values for CSV export
    status : str, one of PASS / WARN / FAIL
    """
    lines = []
    warnings = []
    failures = []
    filename = os.path.basename(filepath)
    file_size = os.path.getsize(filepath)
 
    csv_row = {"filename": filename, "filepath": filepath}
 
    lines.append(("=" * 72, "separator"))
    lines.append(("FILE: %s" % filename, "heading"))
    lines.append(("Size: %.2f MB" % (file_size / (1024.0 * 1024.0)), "normal"))
 
    # ---- Parse TIFF ----
    try:
        tif = TiffQCParser(filepath)
    except Exception as e:
        lines.append(("  ERROR: Cannot parse TIFF - %s" % str(e), "fail"))
        csv_row["status"] = FAIL
        csv_row["errors"] = str(e)
        return lines, csv_row, FAIL
 
    # ---- ImageJ metadata ----
    ij_meta = tif.get_imagej_metadata()
    if ij_meta:
        n_channels = int(ij_meta.get("channels", 1))
        n_slices = int(ij_meta.get("slices", 1))
        n_frames = int(ij_meta.get("frames", 1))
    else:
        # Fall back to IFD count
        n_channels = len(tif.ifds)
        n_slices = 1
        n_frames = 1
        warnings.append("No ImageJ metadata in ImageDescription")
 
    # ---- Dimensions from IFD 0 ----
    width, height, bits, spp = tif.get_dimensions(0)
 
    csv_row["width"] = width
    csv_row["height"] = height
    csv_row["channels"] = n_channels
    csv_row["slices"] = n_slices
    csv_row["frames"] = n_frames
    csv_row["bit_depth"] = bits
 
    lines.append(("", "normal"))
    lines.append(("-- Dimensions --", "heading"))
 
    # Width check
    if width != expected["width"]:
        lines.append(("  Width:    %d  (expected %d)" %
                       (width, expected["width"]), "fail"))
        failures.append("Width %d != expected %d" %
                        (width, expected["width"]))
    else:
        lines.append(("  Width:    %d" % width, "pass"))
 
    # Height check
    if height != expected["height"]:
        lines.append(("  Height:   %d  (expected %d)" %
                       (height, expected["height"]), "fail"))
        failures.append("Height %d != expected %d" %
                        (height, expected["height"]))
    else:
        lines.append(("  Height:   %d" % height, "pass"))
 
    # Channel count check
    if n_channels != expected["channels"]:
        lines.append(("  Channels: %d  (expected %d)" %
                       (n_channels, expected["channels"]), "fail"))
        failures.append("Channels %d != expected %d" %
                        (n_channels, expected["channels"]))
    else:
        lines.append(("  Channels: %d" % n_channels, "pass"))
 
    # Slice count check
    if n_slices != expected["slices"]:
        lines.append(("  Slices:   %d  (expected %d)" %
                       (n_slices, expected["slices"]), "fail"))
        failures.append("Slices %d != expected %d (not 2D)" %
                        (n_slices, expected["slices"]))
    else:
        lines.append(("  Slices:   %d" % n_slices, "pass"))
 
    # Frame count check
    if n_frames != expected["frames"]:
        lines.append(("  Frames:   %d  (expected %d)" %
                       (n_frames, expected["frames"]), "fail"))
        failures.append("Frames %d != expected %d" %
                        (n_frames, expected["frames"]))
    else:
        lines.append(("  Frames:   %d" % n_frames, "pass"))
 
    # Bit depth check
    if bits != expected["bit_depth"]:
        lines.append(("  Bit depth: %d  (expected %d)" %
                       (bits, expected["bit_depth"]), "warn"))
        warnings.append("Bit depth %d != expected %d" %
                        (bits, expected["bit_depth"]))
    else:
        lines.append(("  Bit depth: %d" % bits, "pass"))
 
    # ---- Voxel calibration ----
    lines.append(("", "normal"))
    lines.append(("-- Voxel Calibration --", "heading"))
    voxel = tif.get_voxel_size_um()
    vx = voxel.get("X")
    vy = voxel.get("Y")
 
    if vx is not None:
        lines.append(("  Voxel X: %.4f um" % vx, "pass"))
        csv_row["voxel_x_um"] = round(vx, 4)
    else:
        lines.append(("  Voxel X: not set", "warn"))
        warnings.append("No X voxel calibration")
        csv_row["voxel_x_um"] = ""
 
    if vy is not None:
        lines.append(("  Voxel Y: %.4f um" % vy, "pass"))
        csv_row["voxel_y_um"] = round(vy, 4)
    else:
        lines.append(("  Voxel Y: not set", "warn"))
        warnings.append("No Y voxel calibration")
        csv_row["voxel_y_um"] = ""
 
    # ---- Compression check ----
    compression = tif.get_tag(0, TAG_COMPRESSION, 1)
    comp_name = COMPRESSION_NAMES.get(compression, "Unknown (%d)" % compression)
    csv_row["compression"] = comp_name
 
    # ---- Per-channel intensity statistics ----
    lines.append(("", "normal"))
    lines.append(("-- Channel Intensity Statistics --", "heading"))
 
    can_read_pixels = (compression == 1 and not tif.is_bigtiff)
    if not can_read_pixels:
        reason = "BigTIFF" if tif.is_bigtiff else "compressed (%s)" % comp_name
        lines.append(("  Pixel stats skipped: %s" % reason, "warn"))
        warnings.append("Cannot read pixel data (%s)" % reason)
 
    channel_labels = ["Ch1 (IBA1/glia)", "Ch2 (DAPI)"]
 
    for ch in range(n_channels):
        ch_label = channel_labels[ch] if ch < len(channel_labels) else "Ch%d" % (ch + 1)
        prefix = "ch%d" % (ch + 1)
 
        if not can_read_pixels or ch >= len(tif.ifds):
            for key in ["min", "max", "mean", "std",
                        "zero_fraction", "sat_fraction"]:
                csv_row["%s_%s" % (prefix, key)] = ""
            continue
 
        pixels = tif.read_channel_pixels(ch)
        if pixels is None:
            lines.append(("  %s: could not read pixel data" % ch_label, "warn"))
            warnings.append("%s pixel data unreadable" % ch_label)
            for key in ["min", "max", "mean", "std",
                        "zero_fraction", "sat_fraction"]:
                csv_row["%s_%s" % (prefix, key)] = ""
            continue
 
        stats = compute_channel_stats(pixels)
        csv_row["%s_min" % prefix] = stats["min"]
        csv_row["%s_max" % prefix] = stats["max"]
        csv_row["%s_mean" % prefix] = stats["mean"]
        csv_row["%s_std" % prefix] = stats["std"]
        csv_row["%s_zero_fraction" % prefix] = stats["zero_fraction"]
        csv_row["%s_sat_fraction" % prefix] = stats["sat_fraction"]
 
        lines.append(("  %s:" % ch_label, "normal"))
        lines.append(("    Min: %-8d  Max: %-8d" %
                       (stats["min"], stats["max"]), "normal"))
        lines.append(("    Mean: %-10s  Std: %-10s" %
                       (stats["mean"], stats["std"]), "normal"))
 
        # Flag blank channels
        if stats["max"] == 0:
            lines.append(("    BLANK CHANNEL - all pixels are zero", "fail"))
            failures.append("%s is entirely blank" % ch_label)
        elif stats["zero_fraction"] > 0.50:
            lines.append(("    %.1f%% of pixels are zero" %
                           (stats["zero_fraction"] * 100), "warn"))
            warnings.append("%s has %.1f%% zero pixels" %
                            (ch_label, stats["zero_fraction"] * 100))
 
        # Flag saturation
        if stats["sat_fraction"] > 0.01:
            lines.append(("    %.1f%% of pixels at saturation" %
                           (stats["sat_fraction"] * 100), "warn"))
            warnings.append("%s has %.1f%% saturated pixels" %
                            (ch_label, stats["sat_fraction"] * 100))
 
    # ---- Filename convention ----
    lines.append(("", "normal"))
    lines.append(("-- Filename Convention --", "heading"))
    expected_suffix = "IBA1-DAPI.tiff"
    if filename.endswith(expected_suffix):
        lines.append(("  Suffix '%s': present" % expected_suffix, "pass"))
        csv_row["name_valid"] = "Yes"
    elif filename.lower().endswith(expected_suffix.lower()):
        lines.append(("  Suffix '%s': present (case differs)" %
                       expected_suffix, "warn"))
        warnings.append("Filename suffix has unexpected case")
        csv_row["name_valid"] = "Case"
    else:
        lines.append(("  Suffix '%s': missing" % expected_suffix, "warn"))
        warnings.append("Filename does not end with '%s'" % expected_suffix)
        csv_row["name_valid"] = "No"
 
    # ---- Determine overall status ----
    if failures:
        status = FAIL
    elif warnings:
        status = WARN
    else:
        status = PASS
 
    csv_row["status"] = status
    csv_row["warnings"] = "; ".join(warnings) if warnings else ""
    csv_row["failures"] = "; ".join(failures) if failures else ""
 
    # Status summary line
    lines.append(("", "normal"))
    if status == PASS:
        lines.append(("  STATUS: PASS", "pass"))
    elif status == WARN:
        lines.append(("  STATUS: WARN - %d warning(s)" % len(warnings), "warn"))
    else:
        lines.append(("  STATUS: FAIL - %d failure(s)" % len(failures), "fail"))
        for f in failures:
            lines.append(("    - %s" % f, "fail"))
 
    return lines, csv_row, status
 
 
# ==========================================================================
#  CSV EXPORT
# ==========================================================================
 
CSV_COLUMNS = [
    "filename", "status",
    "width", "height", "channels", "slices", "frames", "bit_depth",
    "compression",
    "ch1_min", "ch1_max", "ch1_mean", "ch1_std",
    "ch1_zero_fraction", "ch1_sat_fraction",
    "ch2_min", "ch2_max", "ch2_mean", "ch2_std",
    "ch2_zero_fraction", "ch2_sat_fraction",
    "voxel_x_um", "voxel_y_um",
    "name_valid", "warnings", "failures",
]
 
 
def write_csv(filepath, rows):
    """Write QC results to a CSV file."""
    # Python 2: csv needs binary mode;  Python 3: needs text + newline=''
    if sys.version_info[0] >= 3:
        f = open(filepath, "w", newline="", encoding="utf-8")
    else:
        f = open(filepath, "wb")
 
    try:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS,
                                extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    finally:
        f.close()
 
 
# ==========================================================================
#  GUI APPLICATION
# ==========================================================================
 
class QCApp(object):
    """Main application window."""
 
    # Color palette (matches existing diagnostic tool)
    BG           = "#1e1e2e"
    BG_FRAME     = "#262637"
    FG           = "#d4d4e8"
    ACCENT       = "#7aa2f7"
    ACCENT_HOVER = "#89b4fa"
    HEADING_FG   = "#7aa2f7"
    PASS_FG      = "#9ece6a"
    WARN_FG      = "#f7b955"
    FAIL_FG      = "#f7768e"
    SEP_FG       = "#565676"
 
    def __init__(self, root):
        self.root = root
        self.root.title("Preprocessing QC - Microglia 2D Macro")
        self.root.configure(bg=self.BG)
        self.root.minsize(860, 680)
 
        if sys.platform == "win32":
            self.FONT_FAMILY = "Consolas"
        else:
            self.FONT_FAMILY = "Menlo"
        self.FONT_MONO = (self.FONT_FAMILY, 10)
 
        # Centre on screen
        sw = self.root.winfo_screenwidth()
        sh = self.root.winfo_screenheight()
        w = min(int(sw * 0.6), 1100)
        h = min(int(sh * 0.75), 860)
        x = (sw - w) // 2
        y = (sh - h) // 2
        self.root.geometry("%dx%d+%d+%d" % (w, h, x, y))
 
        # State
        self.folder_path = tk.StringVar(value="")
        self.report_lines = []
        self.csv_rows = []
        self.is_running = False
 
        # Expected value variables (with defaults from SOP)
        self.exp_width    = tk.StringVar(value="1024")
        self.exp_height   = tk.StringVar(value="1024")
        self.exp_channels = tk.StringVar(value="2")
        self.exp_slices   = tk.StringVar(value="1")
        self.exp_frames   = tk.StringVar(value="1")
        self.exp_bits     = tk.StringVar(value="16")
 
        self._queue = queue_mod.Queue()
        self._build_ui()
        self._poll_queue()
 
    # ---- UI construction ----
 
    def _build_ui(self):
        # Title
        title_frame = tk.Frame(self.root, bg=self.BG)
        title_frame.pack(fill="x", padx=16, pady=(14, 2))
 
        tk.Label(
            title_frame, text="Preprocessing QC",
            font=(self.FONT_FAMILY, 16, "bold"),
            bg=self.BG, fg=self.ACCENT,
        ).pack(side="left")
 
        tk.Label(
            title_frame, text="for Microglia 2D Macro output",
            font=(self.FONT_FAMILY, 11),
            bg=self.BG, fg=self.FG,
        ).pack(side="left", padx=(8, 0), pady=(4, 0))
 
        # Instructions
        instr = (
            "Verify that output TIFF files from the Fiji preprocessing "
            "macro meet expected specifications.\n"
            "Select the output folder, adjust expected values if needed, "
            "then click Run QC."
        )
        tk.Label(
            self.root, text=instr, font=(self.FONT_FAMILY, 10),
            bg=self.BG, fg="#8888aa", justify="left",
            wraplength=780, anchor="w",
        ).pack(fill="x", padx=18, pady=(2, 8))
 
        # Folder selection
        sel_frame = tk.Frame(self.root, bg=self.BG)
        sel_frame.pack(fill="x", padx=16, pady=(0, 6))
 
        tk.Label(
            sel_frame, text="Output folder:", font=(self.FONT_FAMILY, 10, "bold"),
            bg=self.BG, fg=self.FG,
        ).pack(side="left")
 
        path_entry = tk.Entry(
            sel_frame, textvariable=self.folder_path,
            font=self.FONT_MONO, bg=self.BG_FRAME, fg=self.FG,
            insertbackground=self.FG, relief="flat", highlightthickness=1,
            highlightbackground="#44446a", highlightcolor=self.ACCENT,
        )
        path_entry.pack(side="left", fill="x", expand=True,
                        padx=(8, 8), ipady=3)
 
        self.browse_btn = tk.Button(
            sel_frame, text="  Browse...  ", command=self._browse,
            font=(self.FONT_FAMILY, 10), bg=self.ACCENT, fg="#1e1e2e",
            activebackground=self.ACCENT_HOVER, activeforeground="#1e1e2e",
            relief="flat", cursor="hand2", bd=0, padx=10, pady=2,
        )
        self.browse_btn.pack(side="left")
 
        # ---- Expected values panel ----
        exp_outer = tk.Frame(self.root, bg=self.BG)
        exp_outer.pack(fill="x", padx=16, pady=(0, 6))
 
        tk.Label(
            exp_outer, text="Expected values:",
            font=(self.FONT_FAMILY, 10, "bold"),
            bg=self.BG, fg=self.FG,
        ).pack(side="left")
 
        exp_fields = [
            ("W", self.exp_width),
            ("H", self.exp_height),
            ("Ch", self.exp_channels),
            ("Slices", self.exp_slices),
            ("Frames", self.exp_frames),
            ("Bits", self.exp_bits),
        ]
        for label_text, var in exp_fields:
            tk.Label(
                exp_outer, text="  %s:" % label_text,
                font=(self.FONT_FAMILY, 9),
                bg=self.BG, fg="#8888aa",
            ).pack(side="left")
 
            entry = tk.Entry(
                exp_outer, textvariable=var, width=5,
                font=(self.FONT_FAMILY, 10), bg=self.BG_FRAME, fg=self.FG,
                insertbackground=self.FG, relief="flat", highlightthickness=1,
                highlightbackground="#44446a", highlightcolor=self.ACCENT,
                justify="center",
            )
            entry.pack(side="left", padx=(2, 0), ipady=1)
 
        # Buttons
        btn_frame = tk.Frame(self.root, bg=self.BG)
        btn_frame.pack(fill="x", padx=16, pady=(0, 6))
 
        self.run_btn = tk.Button(
            btn_frame, text="  Run QC  ", command=self._run,
            font=(self.FONT_FAMILY, 11, "bold"), bg=self.ACCENT, fg="#1e1e2e",
            activebackground=self.ACCENT_HOVER, activeforeground="#1e1e2e",
            relief="flat", cursor="hand2", bd=0, padx=14, pady=4,
        )
        self.run_btn.pack(side="left")
 
        self.save_csv_btn = tk.Button(
            btn_frame, text="  Save CSV...  ", command=self._save_csv,
            font=(self.FONT_FAMILY, 10), bg="#44446a", fg=self.FG,
            activebackground="#565676", activeforeground=self.FG,
            relief="flat", cursor="hand2", bd=0, padx=10, pady=4,
            state="disabled",
        )
        self.save_csv_btn.pack(side="left", padx=(10, 0))
 
        self.save_report_btn = tk.Button(
            btn_frame, text="  Save Report...  ", command=self._save_report,
            font=(self.FONT_FAMILY, 10), bg="#44446a", fg=self.FG,
            activebackground="#565676", activeforeground=self.FG,
            relief="flat", cursor="hand2", bd=0, padx=10, pady=4,
            state="disabled",
        )
        self.save_report_btn.pack(side="left", padx=(10, 0))
 
        self.clear_btn = tk.Button(
            btn_frame, text="  Clear  ", command=self._clear,
            font=(self.FONT_FAMILY, 10), bg="#44446a", fg=self.FG,
            activebackground="#565676", activeforeground=self.FG,
            relief="flat", cursor="hand2", bd=0, padx=10, pady=4,
        )
        self.clear_btn.pack(side="left", padx=(10, 0))
 
        # Progress bar
        style = ttk.Style()
        style.theme_use("default")
        style.configure(
            "Custom.Horizontal.TProgressbar",
            troughcolor=self.BG_FRAME, background=self.ACCENT,
            darkcolor=self.ACCENT, lightcolor=self.ACCENT, borderwidth=0,
        )
        self.progress = ttk.Progressbar(
            self.root, orient="horizontal", mode="determinate",
            style="Custom.Horizontal.TProgressbar",
        )
        self.progress.pack(fill="x", padx=16, pady=(0, 4))
 
        # Results text area
        text_frame = tk.Frame(self.root, bg=self.BG)
        text_frame.pack(fill="both", expand=True, padx=16, pady=(0, 6))
 
        self.text = ScrolledText(
            text_frame, wrap="none", font=self.FONT_MONO,
            bg=self.BG_FRAME, fg=self.FG, insertbackground=self.FG,
            relief="flat", highlightthickness=1,
            highlightbackground="#44446a", highlightcolor=self.ACCENT,
            state="disabled",
        )
        self.text.pack(fill="both", expand=True)
 
        # Color tags for results
        self.text.tag_configure(
            "heading", foreground=self.HEADING_FG,
            font=(self.FONT_FAMILY, 10, "bold"),
        )
        self.text.tag_configure("pass",      foreground=self.PASS_FG)
        self.text.tag_configure("warn",      foreground=self.WARN_FG)
        self.text.tag_configure("fail",      foreground=self.FAIL_FG)
        self.text.tag_configure("separator", foreground=self.SEP_FG)
        self.text.tag_configure("normal",    foreground=self.FG)
 
        # Horizontal scrollbar
        hscroll = tk.Scrollbar(text_frame, orient="horizontal",
                               command=self.text.xview)
        hscroll.pack(fill="x")
        self.text.configure(xscrollcommand=hscroll.set)
 
        # Status bar
        self.status_var = tk.StringVar(
            value="Ready - select a folder to begin.")
        tk.Label(
            self.root, textvariable=self.status_var,
            font=(self.FONT_FAMILY, 9), bg=self.BG, fg="#666680", anchor="w",
        ).pack(fill="x", padx=18, pady=(0, 8))
 
    # ---- Actions ----
 
    def _browse(self):
        folder = filedialog.askdirectory(
            title="Select the folder containing your preprocessed .tiff files"
        )
        if folder:
            self.folder_path.set(folder)
 
    def _clear(self):
        self.text.configure(state="normal")
        self.text.delete("1.0", "end")
        self.text.configure(state="disabled")
        self.report_lines = []
        self.csv_rows = []
        self.save_csv_btn.configure(state="disabled")
        self.save_report_btn.configure(state="disabled")
        self.progress["value"] = 0
        self.status_var.set("Cleared.")
 
    def _get_expected(self):
        """Read expected values from the GUI fields."""
        try:
            return {
                "width":    int(self.exp_width.get()),
                "height":   int(self.exp_height.get()),
                "channels": int(self.exp_channels.get()),
                "slices":   int(self.exp_slices.get()),
                "frames":   int(self.exp_frames.get()),
                "bit_depth": int(self.exp_bits.get()),
            }
        except ValueError:
            messagebox.showerror(
                "Invalid expected values",
                "All expected values must be integers."
            )
            return None
 
    def _collect_files(self, folder):
        """Find .tiff and .tif files in folder (non-recursive)."""
        found = []
        try:
            entries = os.listdir(folder)
        except OSError as e:
            messagebox.showerror("Cannot read folder", str(e))
            return found
 
        for fn in sorted(entries):
            ext = os.path.splitext(fn)[1].lower()
            if ext in (".tiff", ".tif"):
                found.append(os.path.join(folder, fn))
        return found
 
    def _run(self):
        folder = self.folder_path.get().strip()
        if not folder:
            messagebox.showwarning("No folder selected",
                                   "Please select a folder first.")
            return
        if not os.path.isdir(folder):
            messagebox.showerror("Invalid path",
                                  "'%s' is not a valid folder." % folder)
            return
 
        expected = self._get_expected()
        if expected is None:
            return
 
        files = self._collect_files(folder)
        if not files:
            messagebox.showinfo(
                "No files found",
                "No .tiff or .tif files found in:\n%s" % folder
            )
            return
 
        self.is_running = True
        self.run_btn.configure(state="disabled")
        self.browse_btn.configure(state="disabled")
        self.save_csv_btn.configure(state="disabled")
        self.save_report_btn.configure(state="disabled")
        self._clear()
 
        thread = threading.Thread(
            target=self._run_qc, args=(files, expected))
        thread.daemon = True
        thread.start()
 
    def _run_qc(self, files, expected):
        total = len(files)
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
 
        self._append_line(("PREPROCESSING QC REPORT", "heading"))
        self._append_line(("Generated: %s" % timestamp, "normal"))
        self._append_line(("Files found: %d" % total, "normal"))
        self._append_line((
            "Expected: %dx%d, %dch, %dz, %dt, %d-bit" % (
                expected["width"], expected["height"],
                expected["channels"], expected["slices"],
                expected["frames"], expected["bit_depth"],
            ), "normal"
        ))
        self._append_line(("", "normal"))
 
        counts = {PASS: 0, WARN: 0, FAIL: 0}
 
        for idx, filepath in enumerate(files, 1):
            filename = os.path.basename(filepath)
            self._set_status(
                "Checking file %d/%d: %s" % (idx, total, filename))
            self._set_progress(idx, total)
 
            try:
                file_lines, csv_row, status = inspect_file(
                    filepath, expected)
            except Exception as e:
                file_lines = [
                    ("=" * 72, "separator"),
                    ("FILE: %s" % filename, "heading"),
                    ("  ERROR: %s" % str(e), "fail"),
                ]
                csv_row = {
                    "filename": filename, "status": FAIL,
                    "errors": str(e),
                }
                status = FAIL
 
            counts[status] = counts.get(status, 0) + 1
            self.csv_rows.append(csv_row)
 
            for line in file_lines:
                self._append_line(line)
            self._append_line(("", "normal"))
 
        # ---- Batch summary ----
        self._append_line(("=" * 72, "separator"))
        self._append_line(("BATCH SUMMARY", "heading"))
        self._append_line(("  Total files:  %d" % total, "normal"))
        self._append_line(("  Passed:       %d" % counts[PASS], "pass"))
        if counts[WARN] > 0:
            self._append_line(("  Warnings:     %d" % counts[WARN], "warn"))
        else:
            self._append_line(("  Warnings:     %d" % counts[WARN], "normal"))
        if counts[FAIL] > 0:
            self._append_line(("  Failed:       %d" % counts[FAIL], "fail"))
        else:
            self._append_line(("  Failed:       %d" % counts[FAIL], "normal"))
        self._append_line(("=" * 72, "separator"))
 
        self._set_status("Done - %d file(s) checked.  %d passed, "
                         "%d warnings, %d failed." %
                         (total, counts[PASS], counts[WARN], counts[FAIL]))
        self._set_progress(total, total)
        self._queue.put(("finished", None))
 
    # ---- Thread-safe UI helpers (queue-based) ----
 
    def _append_line(self, line_tuple):
        self.report_lines.append(line_tuple)
        self._queue.put(("line", line_tuple))
 
    def _set_status(self, msg):
        self._queue.put(("status", msg))
 
    def _set_progress(self, current, total):
        self._queue.put(("progress", (current, total)))
 
    def _poll_queue(self):
        batch_limit = 80
        lines_this_tick = 0
 
        while True:
            try:
                msg_type, data = self._queue.get_nowait()
            except queue_mod.Empty:
                break
 
            if msg_type == "line":
                text, tag = data
                self.text.configure(state="normal")
                self.text.insert("end", text + "\n", tag)
                lines_this_tick += 1
                if lines_this_tick >= batch_limit:
                    self.text.see("end")
                    self.text.configure(state="disabled")
                    break
 
            elif msg_type == "status":
                self.status_var.set(data)
 
            elif msg_type == "progress":
                current, total = data
                self.progress.configure(maximum=total, value=current)
 
            elif msg_type == "finished":
                self._on_finished()
 
        if lines_this_tick > 0:
            self.text.see("end")
            self.text.configure(state="disabled")
 
        self.root.after(50, self._poll_queue)
 
    def _on_finished(self):
        self.is_running = False
        self.run_btn.configure(state="normal")
        self.browse_btn.configure(state="normal")
        self.save_csv_btn.configure(state="normal")
        self.save_report_btn.configure(state="normal")
 
    # ---- Save ----
 
    def _save_csv(self):
        if not self.csv_rows:
            messagebox.showinfo("Nothing to save", "Run QC first.")
            return
 
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filepath = filedialog.asksaveasfilename(
            title="Save QC Results as CSV",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
            initialfile="preprocessing_qc_%s.csv" % timestamp,
        )
        if not filepath:
            return
 
        try:
            write_csv(filepath, self.csv_rows)
            self.status_var.set("CSV saved to: %s" % filepath)
        except Exception as e:
            messagebox.showerror("Save failed", str(e))
 
    def _save_report(self):
        if not self.report_lines:
            messagebox.showinfo("Nothing to save", "Run QC first.")
            return
 
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filepath = filedialog.asksaveasfilename(
            title="Save QC Report",
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")],
            initialfile="preprocessing_qc_%s.txt" % timestamp,
        )
        if not filepath:
            return
 
        try:
            with open(filepath, "w") as f:
                for text, _tag in self.report_lines:
                    f.write(text + "\n")
            self.status_var.set("Report saved to: %s" % filepath)
        except Exception as e:
            messagebox.showerror("Save failed", str(e))
 
 
# ==========================================================================
#  ENTRY POINT
# ==========================================================================
 
def main():
    root = tk.Tk()
    QCApp(root)
    root.mainloop()
 
 
if __name__ == "__main__":
    main()