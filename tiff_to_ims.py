#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
TIFF to IMS Converter
=====================
Converts ImageJ composite TIFFs to Imaris .ims (HDF5) format,
bypassing Imaris 9.0.1's buggy TIFF reader.

Requirements:
    h5py   (Python 2.7: pip install h5py==2.10.0)
    numpy  (Python 2.7: pip install numpy==1.16.6)
    For Python 3: pip install h5py numpy  (any recent version)

Usage:
    python tiff_to_ims.py
"""

from __future__ import print_function, division

import os
import sys
import struct
import array
import math
from datetime import datetime

# ---- External dependencies ----
# These are the only non-standard imports. If either is missing
# the script will exit with a clear message rather than a traceback.
try:
    import numpy as np
except ImportError:
    print("ERROR: numpy is required.  Install with:")
    print("  pip install numpy")
    sys.exit(1)

try:
    import h5py
except ImportError:
    print("ERROR: h5py is required.  Install with:")
    print("  pip install h5py")
    sys.exit(1)

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
#  TIFF READER
# ==========================================================================
# This is the same struct-based parser used in the QC tool.
# It reads metadata and raw pixel data from uncompressed TIFFs
# using only the standard library (no Pillow/tifffile needed).

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

# Tag IDs we need
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


class TiffReader(object):
    """
    Reads metadata and pixel data from standard (non-BigTIFF)
    uncompressed TIFFs produced by ImageJ.
    """

    def __init__(self, filepath):
        self.filepath = filepath
        self.endian = "<"
        self.ifds = []           # list of {tag_id: value}
        self._data = None
        self._parse()

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

        # ASCII string
        if type_id == 2:
            raw = self._data[data_offset:data_offset + count]
            if isinstance(raw, bytes):
                return raw.rstrip(b"\x00").decode("latin-1", "replace")
            return raw.rstrip("\x00")

        # UNDEFINED bytes
        if type_id == 7:
            return self._data[data_offset:data_offset + count]

        # RATIONAL / SRATIONAL  (two LONGs: numerator, denominator)
        if type_id in (5, 10):
            long_fmt = "L" if type_id == 5 else "l"
            values = []
            for i in range(count):
                off = data_offset + i * 8
                num = self._read(long_fmt, off)[0]
                den = self._read(long_fmt, off + 4)[0]
                values.append((num, den))
            return values[0] if count == 1 else values

        # All other numeric types
        values = []
        for i in range(count):
            off = data_offset + i * unit_size
            val = self._read(fmt, off)[0]
            values.append(val)
        return values[0] if count == 1 else values

    def _parse(self):
        with open(self.filepath, "rb") as f:
            self._data = f.read()

        # Byte-order marker
        bom = self._data[0:2]
        if bom == b"II":
            self.endian = "<"
        elif bom == b"MM":
            self.endian = ">"
        else:
            raise ValueError("Not a TIFF file")

        magic = self._read("H", 2)[0]
        if magic != 42:
            raise ValueError("Not a standard TIFF (magic=%d)" % magic)

        # Walk the IFD chain to read all pages
        ifd_offset = self._read("L", 4)[0]
        visited = set()
        while ifd_offset != 0 and ifd_offset < len(self._data):
            if ifd_offset in visited:
                break
            visited.add(ifd_offset)

            tags = {}
            n_entries = self._read("H", ifd_offset)[0]
            entry_start = ifd_offset + 2

            for i in range(n_entries):
                eoff = entry_start + i * 12
                tag_id  = self._read("H", eoff)[0]
                type_id = self._read("H", eoff + 2)[0]
                count   = self._read("L", eoff + 4)[0]
                val_off = eoff + 8
                try:
                    value = self._read_tag_value(type_id, count, val_off)
                except Exception:
                    value = None
                tags[tag_id] = value

            self.ifds.append(tags)
            next_off_pos = entry_start + n_entries * 12
            ifd_offset = self._read("L", next_off_pos)[0]

    # ---- Convenience accessors ----

    def get_imagej_metadata(self):
        """Parse the ImageJ-style ImageDescription from IFD 0."""
        desc = self.ifds[0].get(TAG_IMAGE_DESCRIPTION, "")
        if not desc or "ImageJ" not in desc:
            return {}
        result = {}
        for line in desc.strip().split("\n"):
            if "=" in line:
                k, _, v = line.partition("=")
                result[k.strip()] = v.strip()
        return result

    def get_pixel_size_um(self):
        """
        Compute XY pixel size in microns.

        ImageJ stores resolution as pixels-per-unit in the standard
        TIFF XResolution tag, with the unit name in its own
        ImageDescription metadata.  The TIFF ResolutionUnit tag is
        typically set to 1 (dimensionless) by ImageJ, so we rely on
        the ImageDescription 'unit' field instead.
        """
        ij = self.get_imagej_metadata()
        unit = ij.get("unit", "").lower()

        xres_tag = self.ifds[0].get(TAG_X_RESOLUTION)
        yres_tag = self.ifds[0].get(TAG_Y_RESOLUTION)

        def rational_to_float(r):
            if isinstance(r, tuple) and len(r) == 2:
                num, den = r
                return float(num) / den if den != 0 else 0.0
            return float(r) if r else 0.0

        xres = rational_to_float(xres_tag)  # pixels per unit
        yres = rational_to_float(yres_tag)

        # Convert to microns based on the unit string
        if "micron" in unit or "um" in unit or "µm" in unit:
            px = 1.0 / xres if xres > 0 else 0.0
            py = 1.0 / yres if yres > 0 else 0.0
        elif "mm" in unit:
            px = 1000.0 / xres if xres > 0 else 0.0
            py = 1000.0 / yres if yres > 0 else 0.0
        else:
            # Unknown or missing unit; return raw 1/resolution
            px = 1.0 / xres if xres > 0 else 0.0
            py = 1.0 / yres if yres > 0 else 0.0

        return px, py

    def read_channel_pixels(self, ifd_index):
        """
        Read pixel data from one IFD as a 2D numpy array
        with shape (height, width).
        """
        tags = self.ifds[ifd_index]

        compression = tags.get(TAG_COMPRESSION, 1)
        if compression != 1:
            raise ValueError("IFD %d is compressed (type %d); "
                             "only uncompressed TIFFs are supported"
                             % (ifd_index, compression))

        width  = tags.get(TAG_IMAGE_WIDTH, 0)
        height = tags.get(TAG_IMAGE_LENGTH, 0)
        bits   = tags.get(TAG_BITS_PER_SAMPLE, 8)
        if isinstance(bits, (list, tuple)):
            bits = bits[0]

        strip_offsets     = tags.get(TAG_STRIP_OFFSETS)
        strip_byte_counts = tags.get(TAG_STRIP_BYTE_COUNTS)
        if strip_offsets is None or strip_byte_counts is None:
            raise ValueError("IFD %d has no strip data" % ifd_index)

        # Normalize single values to lists
        if not isinstance(strip_offsets, (list, tuple)):
            strip_offsets = [strip_offsets]
        if not isinstance(strip_byte_counts, (list, tuple)):
            strip_byte_counts = [strip_byte_counts]

        # Concatenate all strips into one bytes object
        raw_parts = []
        for offset, count in zip(strip_offsets, strip_byte_counts):
            raw_parts.append(self._data[offset:offset + count])
        raw = b"".join(raw_parts)

        # Convert to numpy array
        if bits == 16:
            dt = np.dtype(np.uint16)
        elif bits == 8:
            dt = np.dtype(np.uint8)
        else:
            raise ValueError("Unsupported bit depth: %d" % bits)

        # Set byte order to match the TIFF
        if self.endian == "<":
            dt = dt.newbyteorder("<")
        else:
            dt = dt.newbyteorder(">")

        pixels = np.frombuffer(raw, dtype=dt).reshape((height, width))

        # Return in native byte order for h5py
        return pixels.astype(dt.newbyteorder("="))


# ==========================================================================
#  IMS WRITER
# ==========================================================================
# The .ims format is HDF5 with a specific group/attribute layout.
# This writer creates the structure Imaris 9.x expects, based on
# the known-good BLK-191 file we analyzed.

def ims_str(text):
    """
    Convert a Python string to the format Imaris uses for HDF5
    attributes: a numpy array of single-byte characters.

    Imaris stores every string attribute -- even numbers like "1024" --
    as an array where each character is a separate element of dtype S1
    (a fixed-length 1-byte string).  For example, "1024" is stored as:
        array([b'1', b'0', b'2', b'4'], dtype='|S1')
    """
    return np.array(list(text), dtype="S1")


def write_ims(output_path, channels, width, height,
              pixel_size_x, pixel_size_y, z_spacing,
              channel_colors, channel_names,
              recording_date=None):
    """
    Write a single-timepoint, single-Z-plane, multi-channel image
    to Imaris .ims (HDF5) format.

    Parameters
    ----------
    output_path : str
        Where to write the .ims file.
    channels : list of numpy.ndarray
        One 2D array per channel, each shape (height, width), uint16.
    width, height : int
        Image dimensions in pixels.
    pixel_size_x, pixel_size_y : float
        Pixel size in microns.
    z_spacing : float
        Z voxel depth in microns (from ImageJ 'spacing' metadata).
    channel_colors : list of str
        One "R G B" string per channel (e.g. "1.000 0.000 0.000").
    channel_names : list of str
        Human-readable name per channel (e.g. "IBA1").
    recording_date : str or None
        ISO-like timestamp.  Defaults to now.
    """
    if recording_date is None:
        recording_date = datetime.now().strftime(
            "%Y-%m-%d %H:%M:%S.000")

    n_channels = len(channels)

    # Physical extents in microns.
    # ExtMin is always 0; ExtMax = number_of_pixels * pixel_size.
    ext_x = width  * pixel_size_x
    ext_y = height * pixel_size_y
    ext_z = z_spacing   # single Z plane

    with h5py.File(output_path, "w", libver="earliest") as f:

        # ============================================================
        # Root attributes  --  required for Imaris to recognise the
        # file as a valid .ims image.  Without these, Imaris reports
        # "no supported image found".
        # ============================================================
        f.attrs.create("DataSetDirectoryName",
                       ims_str("DataSet"))
        f.attrs.create("DataSetInfoDirectoryName",
                       ims_str("DataSetInfo"))
        f.attrs.create("ImarisDataSet",
                       ims_str("ImarisDataSet"))
        f.attrs.create("ImarisVersion",
                       ims_str("5.5.0"))
        f.attrs.create("NumberOfDataSets",
                       np.array([1], dtype=np.int32))
        f.attrs.create("ThumbnailDirectoryName",
                       ims_str("Thumbnail"))

        # ============================================================
        # /Scene and /Scene8  --  empty scene groups.  Imaris uses
        # these for Surpass objects; present but empty for raw images.
        # ============================================================
        for scene_name in ["Scene", "Scene8"]:
            sg = f.create_group(scene_name)
            sg.create_dataset("Content",
                              data=np.void(b""),
                              dtype=np.dtype("V1"))
            sg.create_dataset("Data",
                              data=np.void(b""),
                              dtype=np.dtype("V1"))

        # ============================================================
        # /DataSet  --  the actual pixel data
        # ============================================================
        # Imaris organises data as:
        #   ResolutionLevel / TimePoint / Channel / Data
        # We write only one resolution level and one timepoint.

        for ch_idx in range(n_channels):
            grp_path = ("DataSet/ResolutionLevel 0/TimePoint 0/"
                        "Channel %d" % ch_idx)
            ch_grp = f.create_group(grp_path)

            # Pixel data: shape (Z, Y, X).  Z=1 for our 2D images.
            pixel_data = channels[ch_idx].reshape(
                (1, height, width))

            # Chunking and compression match what Imaris produces:
            #   chunks = (1, height, width//2)
            #   gzip compression level 2
            chunk_x = min(512, width)
            ch_grp.create_dataset(
                "Data", data=pixel_data, dtype=np.uint16,
                chunks=(1, height, chunk_x),
                compression="gzip", compression_opts=2)

            # Per-channel attributes that Imaris expects.
            # These are all stored as character arrays.
            ch_grp.attrs.create(
                "HistogramMin",
                ims_str("%.3f" % float(np.min(channels[ch_idx]))))
            ch_grp.attrs.create(
                "HistogramMax",
                ims_str("%.3f" % float(np.max(channels[ch_idx]))))
            ch_grp.attrs.create(
                "ImageSizeX", ims_str(str(width)))
            ch_grp.attrs.create(
                "ImageSizeY", ims_str(str(height)))
            ch_grp.attrs.create(
                "ImageSizeZ", ims_str("1"))

        # ============================================================
        # /DataSetInfo  --  metadata describing the image
        # ============================================================

        # ---- Per-channel metadata ----
        for ch_idx in range(n_channels):
            ci = f.create_group("DataSetInfo/Channel %d" % ch_idx)

            ci.attrs.create("Color",
                            ims_str(channel_colors[ch_idx]))
            ci.attrs.create("ColorMode",
                            ims_str("BaseColor"))
            ci.attrs.create("ColorOpacity",
                            ims_str("1.000"))

            ch_min = float(np.min(channels[ch_idx]))
            ch_max = float(np.max(channels[ch_idx]))
            ci.attrs.create("ColorRange",
                            ims_str("%.3f %.3f" % (ch_min, ch_max)))
            ci.attrs.create("Description",
                            ims_str(channel_names[ch_idx]))
            ci.attrs.create("GammaCorrection",
                            ims_str("1.000"))
            ci.attrs.create("Name",
                            ims_str(channel_names[ch_idx]))

        # ---- Image-level metadata ----
        img = f.create_group("DataSetInfo/Image")
        img.attrs.create("Description",
                         ims_str("Converted from ImageJ TIFF"))
        img.attrs.create("ExtMax0", ims_str("%.2f" % ext_x))
        img.attrs.create("ExtMax1", ims_str("%.2f" % ext_y))
        img.attrs.create("ExtMax2", ims_str("%.2f" % ext_z))
        img.attrs.create("ExtMin0", ims_str("0"))
        img.attrs.create("ExtMin1", ims_str("0"))
        img.attrs.create("ExtMin2", ims_str("0"))
        img.attrs.create("Name",
                         ims_str("(name not specified)"))
        img.attrs.create("RecordingDate",
                         ims_str(recording_date))
        img.attrs.create("ResampleDimensionX",
                         ims_str("true"))
        img.attrs.create("ResampleDimensionY",
                         ims_str("true"))
        img.attrs.create("ResampleDimensionZ",
                         ims_str("true"))
        img.attrs.create("Unit", ims_str("um"))
        img.attrs.create("X", ims_str(str(width)))
        img.attrs.create("Y", ims_str(str(height)))
        img.attrs.create("Z", ims_str("1"))

        # ---- Imaris version info ----
        iver = f.create_group("DataSetInfo/Imaris")
        iver.attrs.create("ThumbnailMode",
                          ims_str("thumbnailMIP"))
        iver.attrs.create("ThumbnailSize", ims_str("256"))
        iver.attrs.create("Version", ims_str("9.0"))

        ids = f.create_group("DataSetInfo/ImarisDataSet")
        ids.attrs.create("Creator",
                         ims_str("TIFF to IMS Converter"))
        ids.attrs.create("NumberOfImages", ims_str("1"))
        ids.attrs.create("Version", ims_str("9.0"))

        log = f.create_group("DataSetInfo/Log")
        log.attrs.create("Entries", ims_str("0"))

        # ---- TimeInfo: exactly 1 timepoint ----
        ti = f.create_group("DataSetInfo/TimeInfo")
        ti.attrs.create("DatasetTimePoints", ims_str("1"))
        ti.attrs.create("FileTimePoints",    ims_str("1"))
        ti.attrs.create("TimePoint1",
                        ims_str(recording_date))

        # ============================================================
        # /DataSetTimes  --  structured time table
        # ============================================================
        # Imaris expects two compound datasets here.  'Time' stores
        # the ID and birth/death times of each timepoint.  'TimeBegin'
        # stores the timestamp string.

        dst = f.create_group("DataSetTimes")

        time_dtype = np.dtype([
            ("ID",            "<i8"),
            ("Birth",         "<i8"),
            ("Death",         "<i8"),
            ("IDTimeBegin",   "<i8"),
        ])
        time_data = np.array(
            [(0, 0, 1000000000, 0)], dtype=time_dtype)
        dst.create_dataset("Time", data=time_data)

        tb_dtype = np.dtype([
            ("ID",              "<i8"),
            ("ObjectTimeBegin", "S256"),
        ])
        ts_bytes = recording_date.encode("utf-8")
        tb_data = np.array(
            [(0, ts_bytes)], dtype=tb_dtype)
        dst.create_dataset("TimeBegin", data=tb_data)

        # ============================================================
        # /Thumbnail  --  small preview image
        # ============================================================
        # Generate a simple MIP thumbnail from Channel 0 and scale
        # to 256 pixels tall.  Imaris will regenerate this on open,
        # but having one prevents a blank thumbnail in Arena.

        ch0 = channels[0].astype(np.float32)
        # Normalize to 0-255
        ch_min = ch0.min()
        ch_max = ch0.max()
        if ch_max > ch_min:
            ch0 = (ch0 - ch_min) / (ch_max - ch_min) * 255.0
        else:
            ch0 = np.zeros_like(ch0)
        ch0 = ch0.astype(np.uint8)

        # Resize to 256 tall using nearest-neighbor (no scipy needed)
        scale = 256.0 / height
        new_w = max(1, int(width * scale))
        row_idx = (np.arange(256) / scale).astype(int)
        col_idx = (np.arange(new_w) / scale).astype(int)
        row_idx = np.clip(row_idx, 0, height - 1)
        col_idx = np.clip(col_idx, 0, width - 1)
        thumb = ch0[np.ix_(row_idx, col_idx)]

        # Pad/crop to match the width Imaris expects (1024)
        thumb_final = np.zeros((256, 1024), dtype=np.uint8)
        paste_w = min(thumb.shape[1], 1024)
        thumb_final[:, :paste_w] = thumb[:, :paste_w]

        th_grp = f.create_group("Thumbnail")
        th_grp.create_dataset("Data", data=thumb_final)


# ==========================================================================
#  CONVERSION LOGIC
# ==========================================================================

def convert_one_file(tiff_path, output_dir):
    """
    Convert a single ImageJ TIFF to .ims format.

    Returns (output_path, status_message).
    """
    filename = os.path.basename(tiff_path)
    basename = os.path.splitext(filename)[0]
    output_path = os.path.join(output_dir, basename + ".ims")

    # ---- Read the TIFF ----
    tif = TiffReader(tiff_path)

    # ---- Extract metadata ----
    ij_meta = tif.get_imagej_metadata()
    n_channels = int(ij_meta.get("channels", len(tif.ifds)))
    z_spacing  = float(ij_meta.get("spacing", "1.0"))

    width  = tif.ifds[0].get(TAG_IMAGE_WIDTH, 0)
    height = tif.ifds[0].get(TAG_IMAGE_LENGTH, 0)
    bits   = tif.ifds[0].get(TAG_BITS_PER_SAMPLE, 16)
    if isinstance(bits, (list, tuple)):
        bits = bits[0]

    pixel_size_x, pixel_size_y = tif.get_pixel_size_um()

    # Warn if calibration is missing
    warnings = []
    if pixel_size_x == 0 or pixel_size_y == 0:
        warnings.append("No pixel size found; using 1.0 um default")
        pixel_size_x = pixel_size_x or 1.0
        pixel_size_y = pixel_size_y or 1.0

    # ---- Read pixel data for each channel ----
    channels = []
    for ch in range(n_channels):
        if ch < len(tif.ifds):
            pixels = tif.read_channel_pixels(ch)
            channels.append(pixels)
        else:
            warnings.append("Channel %d: no IFD found, filling zeros"
                            % ch)
            channels.append(np.zeros((height, width), dtype=np.uint16))

    # ---- Channel display settings ----
    # Default colors matching the SOP: Ch0 = red (IBA1), Ch1 = blue (DAPI)
    default_colors = [
        "1.000 0.000 0.000",   # red
        "0.000 0.000 1.000",   # blue
        "0.000 1.000 0.000",   # green
        "1.000 1.000 0.000",   # yellow
    ]
    default_names = [
        "IBA1 (microglia)",
        "DAPI",
        "Channel 3",
        "Channel 4",
    ]
    colors = [default_colors[i % len(default_colors)]
              for i in range(n_channels)]
    names  = [default_names[i]  if i < len(default_names)
              else "Channel %d" % (i + 1)
              for i in range(n_channels)]

    # ---- Write the .ims file ----
    write_ims(
        output_path=output_path,
        channels=channels,
        width=width,
        height=height,
        pixel_size_x=pixel_size_x,
        pixel_size_y=pixel_size_y,
        z_spacing=z_spacing,
        channel_colors=colors,
        channel_names=names,
    )

    # Build status message
    msg = "OK"
    if warnings:
        msg += " (warnings: %s)" % "; ".join(warnings)

    voxel_str = "%.4f x %.4f x %.4f um" % (
        pixel_size_x, pixel_size_y, z_spacing)
    detail = ("%dx%d, %dch, %d-bit, voxel=%s"
              % (width, height, n_channels, bits, voxel_str))

    return output_path, msg, detail


# ==========================================================================
#  GUI
# ==========================================================================

class ConverterApp(object):

    BG           = "#1e1e2e"
    BG_FRAME     = "#262637"
    FG           = "#d4d4e8"
    ACCENT       = "#7aa2f7"
    ACCENT_HOVER = "#89b4fa"
    PASS_FG      = "#9ece6a"
    WARN_FG      = "#f7b955"
    FAIL_FG      = "#f7768e"

    def __init__(self, root):
        self.root = root
        self.root.title("TIFF to IMS Converter")
        self.root.configure(bg=self.BG)
        self.root.minsize(780, 520)

        if sys.platform == "win32":
            self.FONT_FAMILY = "Consolas"
        else:
            self.FONT_FAMILY = "Menlo"
        self.FONT_MONO = (self.FONT_FAMILY, 10)

        # Centre on screen
        sw = self.root.winfo_screenwidth()
        sh = self.root.winfo_screenheight()
        w, h = min(int(sw * 0.55), 960), min(int(sh * 0.6), 640)
        x, y = (sw - w) // 2, (sh - h) // 2
        self.root.geometry("%dx%d+%d+%d" % (w, h, x, y))

        self.input_path  = tk.StringVar(value="")
        self.output_path = tk.StringVar(value="")
        self._build_ui()

    def _build_ui(self):
        # Title
        tk.Label(
            self.root, text="TIFF  \u2192  IMS Converter",
            font=(self.FONT_FAMILY, 16, "bold"),
            bg=self.BG, fg=self.ACCENT,
        ).pack(fill="x", padx=16, pady=(14, 2))

        tk.Label(
            self.root,
            text="Converts ImageJ composite TIFFs to Imaris .ims "
                 "(HDF5), bypassing the TIFF reader.",
            font=(self.FONT_FAMILY, 10), bg=self.BG, fg="#8888aa",
            anchor="w",
        ).pack(fill="x", padx=18, pady=(0, 10))

        # Input folder
        self._folder_row("Input folder (.tiff files):",
                         self.input_path, self._browse_input)

        # Output folder
        self._folder_row("Output folder (.ims files):",
                         self.output_path, self._browse_output)

        # Buttons
        btn_frame = tk.Frame(self.root, bg=self.BG)
        btn_frame.pack(fill="x", padx=16, pady=(6, 6))

        self.run_btn = tk.Button(
            btn_frame, text="  Convert  ",
            command=self._run,
            font=(self.FONT_FAMILY, 11, "bold"),
            bg=self.ACCENT, fg="#1e1e2e",
            activebackground=self.ACCENT_HOVER,
            activeforeground="#1e1e2e",
            relief="flat", cursor="hand2", bd=0, padx=14, pady=4,
        )
        self.run_btn.pack(side="left")

        self.clear_btn = tk.Button(
            btn_frame, text="  Clear  ",
            command=self._clear,
            font=(self.FONT_FAMILY, 10),
            bg="#44446a", fg=self.FG,
            activebackground="#565676", activeforeground=self.FG,
            relief="flat", cursor="hand2", bd=0, padx=10, pady=4,
        )
        self.clear_btn.pack(side="left", padx=(10, 0))

        # Progress
        style = ttk.Style()
        style.theme_use("default")
        style.configure(
            "Conv.Horizontal.TProgressbar",
            troughcolor=self.BG_FRAME, background=self.ACCENT,
            darkcolor=self.ACCENT, lightcolor=self.ACCENT,
            borderwidth=0,
        )
        self.progress = ttk.Progressbar(
            self.root, orient="horizontal", mode="determinate",
            style="Conv.Horizontal.TProgressbar",
        )
        self.progress.pack(fill="x", padx=16, pady=(0, 4))

        # Log
        self.text = ScrolledText(
            self.root, wrap="word", font=self.FONT_MONO,
            bg=self.BG_FRAME, fg=self.FG, insertbackground=self.FG,
            relief="flat", highlightthickness=1,
            highlightbackground="#44446a", highlightcolor=self.ACCENT,
            state="disabled", height=14,
        )
        self.text.pack(fill="both", expand=True, padx=16, pady=(0, 6))
        self.text.tag_configure("ok",   foreground=self.PASS_FG)
        self.text.tag_configure("warn", foreground=self.WARN_FG)
        self.text.tag_configure("fail", foreground=self.FAIL_FG)
        self.text.tag_configure("info", foreground=self.ACCENT)

        # Status bar
        self.status_var = tk.StringVar(value="Ready.")
        tk.Label(
            self.root, textvariable=self.status_var,
            font=(self.FONT_FAMILY, 9),
            bg=self.BG, fg="#666680", anchor="w",
        ).pack(fill="x", padx=18, pady=(0, 8))

    def _folder_row(self, label, var, browse_fn):
        frame = tk.Frame(self.root, bg=self.BG)
        frame.pack(fill="x", padx=16, pady=(0, 4))
        tk.Label(
            frame, text=label,
            font=(self.FONT_FAMILY, 10, "bold"),
            bg=self.BG, fg=self.FG,
        ).pack(side="left")
        entry = tk.Entry(
            frame, textvariable=var,
            font=self.FONT_MONO, bg=self.BG_FRAME, fg=self.FG,
            insertbackground=self.FG, relief="flat",
            highlightthickness=1,
            highlightbackground="#44446a",
            highlightcolor=self.ACCENT,
        )
        entry.pack(side="left", fill="x", expand=True,
                   padx=(8, 8), ipady=3)
        tk.Button(
            frame, text=" Browse... ", command=browse_fn,
            font=(self.FONT_FAMILY, 10),
            bg=self.ACCENT, fg="#1e1e2e",
            activebackground=self.ACCENT_HOVER,
            activeforeground="#1e1e2e",
            relief="flat", cursor="hand2", bd=0, padx=8, pady=2,
        ).pack(side="left")

    # ---- Actions ----

    def _browse_input(self):
        d = filedialog.askdirectory(
            title="Select folder containing .tiff files")
        if d:
            self.input_path.set(d)

    def _browse_output(self):
        d = filedialog.askdirectory(
            title="Select output folder for .ims files")
        if d:
            self.output_path.set(d)

    def _clear(self):
        self.text.configure(state="normal")
        self.text.delete("1.0", "end")
        self.text.configure(state="disabled")
        self.progress["value"] = 0
        self.status_var.set("Cleared.")

    def _log(self, msg, tag=""):
        self.text.configure(state="normal")
        if tag:
            self.text.insert("end", msg + "\n", tag)
        else:
            self.text.insert("end", msg + "\n")
        self.text.see("end")
        self.text.configure(state="disabled")
        self.root.update_idletasks()

    def _run(self):
        in_dir  = self.input_path.get().strip()
        out_dir = self.output_path.get().strip()

        if not in_dir or not os.path.isdir(in_dir):
            messagebox.showwarning(
                "No input folder",
                "Please select a valid input folder.")
            return
        if not out_dir:
            messagebox.showwarning(
                "No output folder",
                "Please select an output folder.")
            return
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)

        files = sorted([f for f in os.listdir(in_dir)
                        if f.lower().endswith((".tiff", ".tif"))])
        if not files:
            messagebox.showinfo(
                "No files",
                "No .tiff/.tif files found in the input folder.")
            return

        self._clear()
        self.run_btn.configure(state="disabled")
        total = len(files)
        ok_count = 0
        fail_count = 0

        self._log("Converting %d file(s)..." % total, "info")
        self._log("Input:  %s" % in_dir)
        self._log("Output: %s" % out_dir)
        self._log("")

        for idx, fname in enumerate(files, 1):
            self.progress.configure(maximum=total, value=idx)
            self.status_var.set("Converting %d/%d: %s" %
                                (idx, total, fname))
            self.root.update_idletasks()

            tiff_path = os.path.join(in_dir, fname)
            try:
                out_path, status, detail = convert_one_file(
                    tiff_path, out_dir)
                out_name = os.path.basename(out_path)

                if "warning" in status.lower():
                    self._log("%s -> %s" % (fname, out_name), "warn")
                    self._log("  %s  (%s)" % (status, detail), "warn")
                else:
                    self._log("%s -> %s" % (fname, out_name), "ok")
                    self._log("  %s" % detail)
                ok_count += 1

            except Exception as e:
                self._log("%s -> FAILED" % fname, "fail")
                self._log("  %s" % str(e), "fail")
                fail_count += 1

        self._log("")
        self._log("Done: %d converted, %d failed." %
                  (ok_count, fail_count), "info")
        self.status_var.set(
            "Done. %d converted, %d failed." %
            (ok_count, fail_count))
        self.run_btn.configure(state="normal")


# ==========================================================================
#  ENTRY POINT
# ==========================================================================

def main():
    root = tk.Tk()
    ConverterApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
