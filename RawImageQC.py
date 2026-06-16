#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Raw Image QC Tool for CZI Microscope Files
===========================================
Inspects raw .czi files from Zeiss confocal microscopes, extracts
acquisition metadata, and flags inconsistencies across a batch.

No external packages required -- uses only the Python standard library.

Compatible with Python 3.6+.

How to run:
    python RawImageQC.py

Checks performed per file:
    - File size on disk
    - XY frame dimensions (pixels)
    - Z-slice count
    - Channel count and names
    - Timepoint count
    - Scene count
    - Bit depth (pixel type)
    - Compression type
    - Voxel calibration (X, Y, Z in um)
    - Objective lens info
    - Scan speed, zoom, averaging
    - Acquisition date/time

Cross-file consistency checks:
    - Flags any property that differs from the batch majority
    - Computes theoretical file size to explain disk size differences
"""

import os
import sys
import csv
import struct
import threading
import traceback
import xml.etree.ElementTree as ET
from datetime import datetime
from collections import Counter

import tkinter as tk
import tkinter.filedialog as filedialog
import tkinter.messagebox as messagebox
import tkinter.scrolledtext as scrolledtext_mod
import tkinter.ttk as ttk

ScrolledText = scrolledtext_mod.ScrolledText
import queue as queue_mod


# =========================================================================
#  LIGHTWEIGHT CZI PARSER (standard library only)
# =========================================================================

# CZI pixel type codes -> (name, numpy-like dtype string, bits per pixel)
CZI_PIXEL_TYPES = {
    0:  ("Gray8",        "uint8",    8),
    1:  ("Gray16",       "uint16",   16),
    2:  ("Gray32Float",  "float32",  32),
    3:  ("Bgr24",        "uint8",    8),
    4:  ("Bgr48",        "uint16",   16),
    8:  ("Bgr96Float",   "float32",  32),
    9:  ("Bgra32",       "uint8",    8),
    10: ("Gray64ComplexFloat", "complex64", 64),
    11: ("Bgr192ComplexFloat", "complex64", 64),
    12: ("Gray32",       "uint32",   32),
    13: ("Gray64Float",  "float64",  64),
}

CZI_COMPRESSION_TYPES = {
    0: "Uncompressed",
    1: "JPEG",
    2: "LZW",
    4: "JPEGXR",
    5: "Camera specific",
    100: "zstd (100)",
    1000: "zstd (1000)",
}


class CziQCParser:
    """
    Lightweight CZI file parser for metadata extraction.

    Reads only the file header, metadata XML segment, and subblock
    directory -- no pixel data is loaded. Uses only the standard
    library (struct + xml.etree).
    """

    def __init__(self, filepath):
        self.filepath = filepath
        self.metadata_xml = None
        self.subblocks = []
        self._parse()

    def _parse(self):
        with open(self.filepath, "rb") as f:
            self._parse_file_header(f)
            self._parse_metadata_segment(f)
            self._parse_directory_segment(f)

    def _parse_file_header(self, f):
        """Read the ZISRAWFILE segment header and file header data."""
        # Segment header: 16 bytes SID + 8 bytes alloc + 8 bytes used
        sid = f.read(16)
        if not sid.startswith(b"ZISRAWFILE"):
            raise ValueError(
                "Not a CZI file (magic: %r)" % sid[:10])

        f.read(8)   # allocated_size
        f.read(8)   # used_size

        # File header data (starts at offset 32)
        self.version = struct.unpack("<ii", f.read(8))
        f.read(8)   # reserved
        f.read(16)  # primary_file_guid
        f.read(16)  # file_guid
        self.file_part = struct.unpack("<i", f.read(4))[0]
        self.directory_position = struct.unpack("<q", f.read(8))[0]
        self.metadata_position = struct.unpack("<q", f.read(8))[0]
        f.read(4)   # update_pending
        self.attachment_directory_position = struct.unpack("<q", f.read(8))[0]

    def _parse_metadata_segment(self, f):
        """Read the ZISRAWMETADATA segment and extract XML."""
        if self.metadata_position <= 0:
            return

        f.seek(self.metadata_position)
        sid = f.read(16)
        if not sid.startswith(b"ZISRAWMETADATA"):
            return

        f.read(8)  # allocated_size
        f.read(8)  # used_size

        xml_size = struct.unpack("<i", f.read(4))[0]
        f.read(4)    # binary_size
        f.read(248)  # spare

        xml_bytes = f.read(xml_size)
        self.metadata_xml = xml_bytes.decode("utf-8", errors="replace")

    def _parse_directory_segment(self, f):
        """Read the ZISRAWDIRECTORY segment and parse subblock entries."""
        if self.directory_position <= 0:
            return

        f.seek(self.directory_position)
        sid = f.read(16)
        if not sid.startswith(b"ZISRAWDIRECTORY"):
            return

        f.read(8)  # allocated_size
        f.read(8)  # used_size

        entry_count = struct.unpack("<i", f.read(4))[0]
        f.read(124)  # reserved

        for _ in range(entry_count):
            entry = self._parse_directory_entry(f)
            if entry is not None:
                self.subblocks.append(entry)

    def _parse_directory_entry(self, f):
        """Parse a single SubBlockDirectoryEntryDV."""
        schema = f.read(2)
        if schema != b"DV":
            # Unknown schema; we cannot determine the entry size,
            # so stop parsing further entries.
            return None

        pixel_type_code = struct.unpack("<i", f.read(4))[0]
        file_position = struct.unpack("<q", f.read(8))[0]
        file_part = struct.unpack("<i", f.read(4))[0]
        compression_code = struct.unpack("<i", f.read(4))[0]
        pyramid_type = struct.unpack("<B", f.read(1))[0]
        f.read(1)  # spare1
        f.read(4)  # spare2
        dimension_count = struct.unpack("<i", f.read(4))[0]

        # Parse dimension entries (20 bytes each)
        dimensions = {}
        for _ in range(dimension_count):
            dim_id = f.read(4).decode("ascii", errors="replace").rstrip("\x00")
            start = struct.unpack("<i", f.read(4))[0]
            size = struct.unpack("<i", f.read(4))[0]
            f.read(4)  # start_coordinate (float)
            stored_size = struct.unpack("<i", f.read(4))[0]
            dimensions[dim_id] = {
                "start": start,
                "size": size,
                "stored_size": stored_size,
            }

        # Resolve pixel type
        pt = CZI_PIXEL_TYPES.get(pixel_type_code)
        if pt:
            dtype_name, dtype_str, bits = pt
        else:
            dtype_name = "unknown(%d)" % pixel_type_code
            dtype_str = "unknown"
            bits = 0

        comp_name = CZI_COMPRESSION_TYPES.get(
            compression_code, "unknown(%d)" % compression_code)

        return {
            "pixel_type_code": pixel_type_code,
            "dtype_name": dtype_name,
            "dtype": dtype_str,
            "bits_per_component": bits,
            "compression_code": compression_code,
            "compression": comp_name,
            "pyramid_type": pyramid_type,
            "dimensions": dimensions,
        }


# =========================================================================
#  CZI METADATA EXTRACTION
# =========================================================================

def extract_czi_metadata(filepath):
    """
    Extract all QC-relevant metadata from a single .czi file.

    Returns a dict with standardized keys.
    """
    info = {
        "filepath": filepath,
        "filename": os.path.basename(filepath),
        "file_size_bytes": os.path.getsize(filepath),
        "file_size_mb": round(
            os.path.getsize(filepath) / (1024.0 * 1024.0), 2),
    }

    czi = CziQCParser(filepath)

    # -- Dimensions from subblock directory ----------------------------
    sbs = czi.subblocks
    info["subblock_count"] = len(sbs)

    z_vals = set()
    c_vals = set()
    t_vals = set()
    s_vals = set()
    xy_shapes = set()
    compressions = set()
    dtypes = set()
    bit_depths = set()

    for sb in sbs:
        dims = sb["dimensions"]
        for dim_id, dinfo in dims.items():
            if dim_id == "Z":
                z_vals.add(dinfo["start"])
            elif dim_id == "C":
                c_vals.add(dinfo["start"])
            elif dim_id == "T":
                t_vals.add(dinfo["start"])
            elif dim_id == "S":
                s_vals.add(dinfo["start"])

        x_size = dims.get("X", {}).get("size", 0)
        y_size = dims.get("Y", {}).get("size", 0)
        if x_size > 0 and y_size > 0:
            xy_shapes.add((y_size, x_size))

        compressions.add(sb["compression"])
        dtypes.add(sb["dtype"])
        bit_depths.add(sb["bits_per_component"])

    info["size_z"] = len(z_vals)
    info["size_c"] = len(c_vals)
    info["size_t"] = len(t_vals)
    info["size_s"] = len(s_vals) if s_vals else 1

    if xy_shapes:
        main_shape = max(xy_shapes, key=lambda s: s[0] * s[1])
        info["size_y"] = main_shape[0]
        info["size_x"] = main_shape[1]
    else:
        info["size_y"] = 0
        info["size_x"] = 0

    info["compression"] = ", ".join(sorted(compressions))
    info["dtype"] = ", ".join(sorted(dtypes))

    bd = max(bit_depths) if bit_depths else 0
    info["bit_depth"] = bd

    # Theoretical uncompressed data size
    bytes_per_pixel = bd // 8 if bd > 0 else 1
    info["theoretical_data_mb"] = round(
        info["size_x"] * info["size_y"] * info["size_z"]
        * info["size_c"] * info["size_t"] * info["size_s"]
        * bytes_per_pixel / (1024.0 * 1024.0), 2
    )

    # -- Parse metadata XML --------------------------------------------
    info["voxel_x_um"] = None
    info["voxel_y_um"] = None
    info["voxel_z_um"] = None
    info["fov_x_um"] = None
    info["fov_y_um"] = None
    info["channel_names"] = []
    info["channel_names_str"] = ""
    info["objective"] = ""
    info["magnification"] = ""
    info["na"] = ""
    info["immersion"] = ""
    info["frame_setting"] = ""
    info["scan_speed"] = ""
    info["zoom"] = ""
    info["bits_per_pixel_setting"] = ""
    info["scan_direction"] = ""
    info["averaging"] = ""
    info["acquisition_date"] = ""

    if czi.metadata_xml:
        try:
            root = ET.fromstring(czi.metadata_xml)
            _parse_xml_metadata(root, info)
        except ET.ParseError:
            info.setdefault("_warnings", []).append(
                "Could not parse metadata XML")

    return info


def _parse_xml_metadata(root, info):
    """Extract structured fields from the CZI metadata XML tree."""

    # -- Voxel calibration --
    for dist in root.iter("Distance"):
        axis_id = dist.get("Id")
        val_elem = dist.find("Value")
        if axis_id and val_elem is not None and val_elem.text:
            try:
                um = float(val_elem.text) * 1e6
                if axis_id == "X":
                    info["voxel_x_um"] = round(um, 4)
                elif axis_id == "Y":
                    info["voxel_y_um"] = round(um, 4)
                elif axis_id == "Z":
                    info["voxel_z_um"] = round(um, 4)
            except ValueError:
                pass

    # Physical field of view
    if info["voxel_x_um"] and info["size_x"]:
        info["fov_x_um"] = round(info["voxel_x_um"] * info["size_x"], 1)
    if info["voxel_y_um"] and info["size_y"]:
        info["fov_y_um"] = round(info["voxel_y_um"] * info["size_y"], 1)

    # -- Channel names --
    channel_names = []
    seen = set()
    for ch_elem in root.iter("Channel"):
        ch_id = ch_elem.get("Id")
        ch_name = ch_elem.get("Name")
        if ch_id and ch_id.startswith("Channel:") and ch_name:
            idx = ch_id.replace("Channel:", "")
            if idx not in seen:
                seen.add(idx)
                channel_names.append(ch_name)
    info["channel_names"] = channel_names
    info["channel_names_str"] = ", ".join(channel_names)

    # -- Objective (active one has Manufacturer/Model) --
    for obj in root.iter("Objective"):
        model_elem = obj.find("Manufacturer")
        nom_mag = obj.find("NominalMagnification")
        if model_elem is not None and nom_mag is not None:
            m = model_elem.find("Model")
            if m is not None and m.text:
                info["objective"] = m.text.strip()
            info["magnification"] = nom_mag.text.strip()
            na = obj.find("LensNA")
            if na is not None and na.text:
                info["na"] = na.text.strip()
            imm = obj.find("Immersion")
            if imm is not None and imm.text:
                info["immersion"] = imm.text.strip()
            break

    # -- Acquisition settings from Camera block --
    for cam in root.iter("Camera"):
        _get_xml_text(cam, "Frame", info, "frame_setting")
        _get_xml_text(cam, "BitsPerPixel", info, "bits_per_pixel_setting")
        _get_xml_text(cam, "Zoom", info, "zoom")
        _get_xml_text(cam, "ScanSpeed", info, "scan_speed")
        _get_xml_text(cam, "ScanDirection", info, "scan_direction")
        _get_xml_text(cam, "AveragingNumber", info, "averaging")
        break

    # -- Acquisition date --
    for tag_name in ("AcquisitionDateAndTime", "CreationDate"):
        for elem in root.iter(tag_name):
            if elem.text and elem.text.strip():
                info["acquisition_date"] = elem.text.strip()
                break
        if info["acquisition_date"]:
            break


def _get_xml_text(parent, tag, info, key):
    """Helper: set info[key] from parent.find(tag).text if present."""
    elem = parent.find(tag)
    if elem is not None and elem.text:
        info[key] = elem.text.strip()


# =========================================================================
#  CROSS-FILE COMPARISON
# =========================================================================

# Properties to check for consistency across the batch.
# (key, display_name, is_critical)
#   is_critical = True  -> FAIL if inconsistent
#   is_critical = False -> WARN if inconsistent
CONSISTENCY_CHECKS = [
    ("size_x",                  "Frame Width (px)",      True),
    ("size_y",                  "Frame Height (px)",     True),
    ("size_z",                  "Z-Slices",              True),
    ("size_c",                  "Channel Count",         True),
    ("channel_names_str",       "Channel Names",         True),
    ("bit_depth",               "Bit Depth",             True),
    ("voxel_x_um",              "Voxel X (um)",          True),
    ("voxel_y_um",              "Voxel Y (um)",          True),
    ("voxel_z_um",              "Voxel Z (um)",          True),
    ("size_t",                  "Timepoints",            False),
    ("size_s",                  "Scenes",                False),
    ("compression",             "Compression",           False),
    ("objective",               "Objective",             False),
    ("magnification",           "Magnification",         False),
    ("scan_speed",              "Scan Speed",            False),
    ("zoom",                    "Zoom",                  False),
    ("bits_per_pixel_setting",  "BitsPerPixel Setting",  False),
    ("averaging",               "Averaging",             False),
]


def compute_consensus(all_info):
    """
    Determine the majority (consensus) value for each checked property
    and flag files that deviate.

    Returns
    -------
    consensus : dict  {key: majority_value}
    deviations : list of (filename, key, display_name, file_value,
                          consensus_value, is_critical)
    """
    consensus = {}
    for key, display_name, is_critical in CONSISTENCY_CHECKS:
        values = [info.get(key) for info in all_info]
        counter = Counter(values)
        most_common_val, _ = counter.most_common(1)[0]
        consensus[key] = most_common_val

    deviations = []
    for info in all_info:
        for key, display_name, is_critical in CONSISTENCY_CHECKS:
            file_val = info.get(key)
            cons_val = consensus[key]
            if file_val != cons_val:
                deviations.append((
                    info["filename"], key, display_name,
                    file_val, cons_val, is_critical,
                ))

    return consensus, deviations


# =========================================================================
#  REPORT GENERATION (text + tag tuples for GUI)
# =========================================================================

PASS = "PASS"
WARN = "WARN"
FAIL = "FAIL"


def format_file_report(info):
    """Generate report lines for a single file."""
    lines = []
    fname = info["filename"]

    lines.append(("=" * 72, "separator"))
    lines.append(("FILE: %s" % fname, "heading"))
    lines.append((
        "  Disk size: %.2f MB" % info["file_size_mb"], "normal"))

    if "_error" in info:
        lines.append(("  ERROR: %s" % info["_error"], "fail"))
        return lines, FAIL

    # -- Dimensions --
    lines.append(("", "normal"))
    lines.append(("  -- Dimensions --", "heading"))
    lines.append(("  Frame (XY):   %d x %d px" % (
        info["size_x"], info["size_y"]), "normal"))
    lines.append(("  Z-Slices:     %d" % info["size_z"], "normal"))
    lines.append(("  Channels:     %d" % info["size_c"], "normal"))
    lines.append(("  Timepoints:   %d" % info["size_t"], "normal"))
    if info["size_s"] > 1:
        lines.append(("  Scenes:       %d" % info["size_s"], "normal"))
    lines.append(("  Bit Depth:    %d-bit (%s)" % (
        info["bit_depth"], info["dtype"]), "normal"))
    lines.append(("  Compression:  %s" % info["compression"], "normal"))
    lines.append(("  Subblocks:    %d" % info["subblock_count"], "normal"))
    lines.append(("  Theoretical uncompressed data: %.2f MB" % (
        info["theoretical_data_mb"]), "normal"))

    # -- Voxel calibration --
    lines.append(("", "normal"))
    lines.append(("  -- Voxel Calibration --", "heading"))

    for axis in ("X", "Y", "Z"):
        key = "voxel_%s_um" % axis.lower()
        val = info.get(key)
        if val is not None:
            lines.append(("  Voxel %s:  %.4f um" % (axis, val), "pass"))
        else:
            lines.append(("  Voxel %s:  NOT SET" % axis, "warn"))

    if info.get("fov_x_um") and info.get("fov_y_um"):
        lines.append(("  FOV (XY):  %.1f x %.1f um" % (
            info["fov_x_um"], info["fov_y_um"]), "normal"))

    # -- Channels --
    lines.append(("", "normal"))
    lines.append(("  -- Channels --", "heading"))
    for i, name in enumerate(info.get("channel_names", [])):
        lines.append(("  Ch%d: %s" % (i, name), "normal"))

    # -- Acquisition settings --
    lines.append(("", "normal"))
    lines.append(("  -- Acquisition Settings --", "heading"))
    if info.get("objective"):
        lines.append(("  Objective:    %s" % info["objective"], "normal"))
    if info.get("magnification"):
        mag_str = "%sx" % info["magnification"]
        if info.get("na"):
            mag_str += " / NA %s" % info["na"]
        if info.get("immersion"):
            mag_str += " (%s)" % info["immersion"]
        lines.append(("  Magnification: %s" % mag_str, "normal"))
    if info.get("frame_setting"):
        lines.append((
            "  Frame setting: %s" % info["frame_setting"], "normal"))
    if info.get("bits_per_pixel_setting"):
        lines.append((
            "  BitsPerPixel:  %s" % info["bits_per_pixel_setting"], "normal"))
    if info.get("scan_speed"):
        lines.append(("  Scan speed:    %s" % info["scan_speed"], "normal"))
    if info.get("zoom"):
        lines.append(("  Zoom:          %s" % info["zoom"], "normal"))
    if info.get("scan_direction"):
        lines.append((
            "  Scan direction: %s" % info["scan_direction"], "normal"))
    if info.get("averaging"):
        lines.append(("  Averaging:     %s" % info["averaging"], "normal"))
    if info.get("acquisition_date"):
        lines.append((
            "  Acquired:      %s" % info["acquisition_date"], "normal"))

    return lines, PASS


def format_consistency_report(consensus, deviations, all_info):
    """Generate the cross-file consistency comparison report."""
    lines = []
    lines.append(("=" * 72, "separator"))
    lines.append(("CROSS-FILE CONSISTENCY CHECK", "heading"))
    lines.append(("", "normal"))

    if not deviations:
        lines.append((
            "  All %d files are consistent across all checked properties."
            % len(all_info), "pass"))
        return lines

    # Group deviations by property
    by_prop = {}
    for fname, key, display_name, file_val, cons_val, is_crit in deviations:
        if key not in by_prop:
            by_prop[key] = {
                "display_name": display_name,
                "consensus": cons_val,
                "is_critical": is_crit,
                "files": [],
            }
        by_prop[key]["files"].append((fname, file_val))

    lines.append(("  Consensus values (majority across %d files):" % (
        len(all_info)), "normal"))
    for key, display_name, _ in CONSISTENCY_CHECKS:
        val = consensus.get(key)
        if val is not None and val != "":
            lines.append(("    %-22s = %s" % (display_name, val), "normal"))
    lines.append(("", "normal"))

    n_crit = sum(1 for v in by_prop.values() if v["is_critical"])
    n_warn = sum(1 for v in by_prop.values() if not v["is_critical"])

    if n_crit > 0:
        lines.append((
            "  CRITICAL inconsistencies (%d):" % n_crit, "fail"))
        for key, data in by_prop.items():
            if not data["is_critical"]:
                continue
            lines.append(("", "normal"))
            lines.append(("    %s  (consensus: %s)" % (
                data["display_name"], data["consensus"]), "fail"))
            for fname, fval in data["files"]:
                lines.append(("      %s -> %s" % (fname, fval), "fail"))

    if n_warn > 0:
        lines.append(("", "normal"))
        lines.append((
            "  Minor inconsistencies (%d):" % n_warn, "warn"))
        for key, data in by_prop.items():
            if data["is_critical"]:
                continue
            lines.append(("", "normal"))
            lines.append(("    %s  (consensus: %s)" % (
                data["display_name"], data["consensus"]), "warn"))
            for fname, fval in data["files"]:
                lines.append(("      %s -> %s" % (fname, fval), "warn"))

    return lines


# =========================================================================
#  CSV EXPORT
# =========================================================================

CSV_COLUMNS = [
    "filename", "status",
    "file_size_mb", "theoretical_data_mb",
    "size_x", "size_y", "size_z", "size_c", "size_t", "size_s",
    "bit_depth", "dtype", "compression",
    "voxel_x_um", "voxel_y_um", "voxel_z_um",
    "fov_x_um", "fov_y_um",
    "channel_names_str",
    "objective", "magnification", "na", "immersion",
    "frame_setting", "bits_per_pixel_setting",
    "scan_speed", "zoom", "scan_direction", "averaging",
    "acquisition_date",
    "subblock_count",
    "deviations",
]


def write_csv(filepath, rows):
    """Write QC results to a CSV file."""
    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f, fieldnames=CSV_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


# =========================================================================
#  GUI APPLICATION
# =========================================================================

class CZIQCApp:
    """Main application window."""

    # Color palette (matches existing QC diagnostic tool)
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
        self.root.title("Raw Image QC - CZI Microscope File Inspector")
        self.root.configure(bg=self.BG)
        self.root.minsize(900, 700)

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

        self._queue = queue_mod.Queue()
        self._build_ui()
        self._poll_queue()

    # ---- UI construction ----

    def _build_ui(self):
        # Title
        title_frame = tk.Frame(self.root, bg=self.BG)
        title_frame.pack(fill="x", padx=16, pady=(14, 2))

        tk.Label(
            title_frame, text="Raw Image QC",
            font=(self.FONT_FAMILY, 16, "bold"),
            bg=self.BG, fg=self.ACCENT,
        ).pack(side="left")

        tk.Label(
            title_frame, text="for CZI microscope files",
            font=(self.FONT_FAMILY, 11),
            bg=self.BG, fg=self.FG,
        ).pack(side="left", padx=(8, 0), pady=(4, 0))

        # Instructions
        instr = (
            "Inspect raw .czi files and flag acquisition-setting "
            "inconsistencies across a batch.\n"
            "Select a folder containing .czi files, then click Run QC. "
            "Files are compared against each other (majority = consensus)."
        )
        tk.Label(
            self.root, text=instr, font=(self.FONT_FAMILY, 10),
            bg=self.BG, fg="#8888aa", justify="left",
            wraplength=820, anchor="w",
        ).pack(fill="x", padx=18, pady=(2, 8))

        # Folder selection
        sel_frame = tk.Frame(self.root, bg=self.BG)
        sel_frame.pack(fill="x", padx=16, pady=(0, 6))

        tk.Label(
            sel_frame, text="CZI folder:",
            font=(self.FONT_FAMILY, 10, "bold"),
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

        # Recursive subfolder option
        opt_frame = tk.Frame(self.root, bg=self.BG)
        opt_frame.pack(fill="x", padx=16, pady=(0, 6))

        self.recursive_var = tk.BooleanVar(value=False)
        tk.Checkbutton(
            opt_frame, text="Include subfolders",
            variable=self.recursive_var,
            font=(self.FONT_FAMILY, 10), bg=self.BG, fg=self.FG,
            selectcolor=self.BG_FRAME, activebackground=self.BG,
            activeforeground=self.FG,
        ).pack(side="left")

        # Buttons
        btn_frame = tk.Frame(self.root, bg=self.BG)
        btn_frame.pack(fill="x", padx=16, pady=(0, 6))

        self.run_btn = tk.Button(
            btn_frame, text="  Run QC  ", command=self._run,
            font=(self.FONT_FAMILY, 11, "bold"),
            bg=self.ACCENT, fg="#1e1e2e",
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
            title="Select the folder containing your .czi files"
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

    def _collect_files(self, folder):
        """Find .czi files in folder."""
        found = []
        if self.recursive_var.get():
            for dirpath, _, filenames in os.walk(folder):
                for fn in sorted(filenames):
                    if fn.lower().endswith(".czi"):
                        found.append(os.path.join(dirpath, fn))
        else:
            try:
                entries = os.listdir(folder)
            except OSError as e:
                messagebox.showerror("Cannot read folder", str(e))
                return found
            for fn in sorted(entries):
                if fn.lower().endswith(".czi"):
                    found.append(os.path.join(folder, fn))
        return found

    def _run(self):
        folder = self.folder_path.get().strip()
        if not folder:
            messagebox.showwarning(
                "No folder selected", "Please select a folder first.")
            return
        if not os.path.isdir(folder):
            messagebox.showerror(
                "Invalid path", "'%s' is not a valid folder." % folder)
            return

        files = self._collect_files(folder)
        if not files:
            messagebox.showinfo(
                "No files found",
                "No .czi files found in:\n%s" % folder)
            return

        self.is_running = True
        self.run_btn.configure(state="disabled")
        self.browse_btn.configure(state="disabled")
        self.save_csv_btn.configure(state="disabled")
        self.save_report_btn.configure(state="disabled")
        self._clear()

        thread = threading.Thread(target=self._run_qc, args=(files,))
        thread.daemon = True
        thread.start()

    def _run_qc(self, files):
        total = len(files)
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        self._append_line(("RAW IMAGE QC REPORT", "heading"))
        self._append_line(("Generated: %s" % timestamp, "normal"))
        self._append_line(("Files found: %d" % total, "normal"))
        self._append_line(("", "normal"))

        # ---- Phase 1: Extract metadata from each file ----
        all_info = []
        for idx, filepath in enumerate(files, 1):
            filename = os.path.basename(filepath)
            self._set_status(
                "Reading metadata %d/%d: %s" % (idx, total, filename))
            self._set_progress(idx, total * 2)

            try:
                info = extract_czi_metadata(filepath)
            except Exception as e:
                info = {
                    "filename": filename,
                    "filepath": filepath,
                    "file_size_mb": round(
                        os.path.getsize(filepath) / (1024.0 * 1024.0), 2),
                    "_error": str(e),
                }
            all_info.append(info)

            file_lines, file_status = format_file_report(info)
            for line in file_lines:
                self._append_line(line)
            self._append_line(("", "normal"))

        # ---- Phase 2: Cross-file consistency ----
        self._set_status("Comparing files...")
        self._set_progress(total + 1, total * 2)

        valid_info = [i for i in all_info if "_error" not in i]

        if len(valid_info) >= 2:
            consensus, deviations = compute_consensus(valid_info)
            consistency_lines = format_consistency_report(
                consensus, deviations, valid_info)
            for line in consistency_lines:
                self._append_line(line)
            self._append_line(("", "normal"))
        elif len(valid_info) == 1:
            self._append_line(("=" * 72, "separator"))
            self._append_line(("CROSS-FILE CONSISTENCY CHECK", "heading"))
            self._append_line((
                "  Only 1 valid file - no cross-file comparison possible.",
                "warn"))
            self._append_line(("", "normal"))
            deviations = []
        else:
            deviations = []

        # ---- Phase 3: Build CSV rows and summary ----
        counts = {PASS: 0, WARN: 0, FAIL: 0}

        dev_by_file = {}
        for fname, key, display_name, fval, cval, is_crit in deviations:
            dev_by_file.setdefault(fname, []).append(
                "%s: %s (expected %s)" % (display_name, fval, cval))

        for info in all_info:
            fname = info["filename"]
            has_error = "_error" in info
            file_devs = dev_by_file.get(fname, [])
            has_critical = any(
                is_crit for fn, k, dn, fv, cv, is_crit in deviations
                if fn == fname
            )

            if has_error:
                status = FAIL
            elif has_critical:
                status = FAIL
            elif file_devs:
                status = WARN
            else:
                status = PASS
            counts[status] += 1

            csv_row = dict(info)
            csv_row["status"] = status
            csv_row["deviations"] = "; ".join(file_devs) if file_devs else ""
            self.csv_rows.append(csv_row)

        # ---- Batch summary ----
        self._append_line(("=" * 72, "separator"))
        self._append_line(("BATCH SUMMARY", "heading"))
        self._append_line(("  Total files:  %d" % total, "normal"))
        self._append_line(("  Passed:       %d" % counts[PASS], "pass"))

        tag = "warn" if counts[WARN] > 0 else "normal"
        self._append_line(("  Warnings:     %d" % counts[WARN], tag))

        tag = "fail" if counts[FAIL] > 0 else "normal"
        self._append_line(("  Failed:       %d" % counts[FAIL], tag))

        self._append_line(("=" * 72, "separator"))

        self._set_status(
            "Done - %d file(s) checked.  %d passed, "
            "%d warnings, %d failed." % (
                total, counts[PASS], counts[WARN], counts[FAIL]))
        self._set_progress(total * 2, total * 2)
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
            initialfile="raw_image_qc_%s.csv" % timestamp,
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
            initialfile="raw_image_qc_%s.txt" % timestamp,
        )
        if not filepath:
            return
        try:
            with open(filepath, "w", encoding="utf-8") as f:
                for text, _tag in self.report_lines:
                    f.write(text + "\n")
            self.status_var.set("Report saved to: %s" % filepath)
        except Exception as e:
            messagebox.showerror("Save failed", str(e))


# =========================================================================
#  ENTRY POINT
# =========================================================================

def main():
    root = tk.Tk()
    CZIQCApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
