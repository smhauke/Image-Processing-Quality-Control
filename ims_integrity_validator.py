#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
IMS Integrity Validator
=======================
Reads Imaris .ims (HDF5) files directly and checks for:
  - Cross-contamination (same pixel data in different files)
  - Phantom dimensions (zero-filled Z-planes or timepoints)
  - Missing or default spatial calibration
  - Metadata inconsistencies across a batch

Designed to audit files produced by Imaris 9.0.1's TIFF reader,
which can silently merge data from separate images into one file.

Requirements:
    h5py   (Python 2.7: pip install h5py==2.10.0)
    numpy  (Python 2.7: pip install numpy==1.16.6)

Usage:
    python ims_integrity_validator.py
"""

from __future__ import print_function, division

import os
import sys
import csv
import hashlib
from datetime import datetime
from collections import OrderedDict

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

# ======================================================================
#  IMS READER
# ======================================================================

def get_display_name(filepath):
    """
    Derive a human-readable name for an IMS file.

    In Imaris Arena storage, every file is called 'dataset.ims'
    (or similar) inside a folder whose name is the image title.
    This function returns the parent folder name when the filename
    is generic, or the filename itself otherwise.
    """
    fname = os.path.basename(filepath)
    stem = os.path.splitext(fname)[0].lower()
    if stem in ("dataset", "data", "image"):
        return os.path.basename(os.path.dirname(filepath))
    return os.path.splitext(fname)[0]
# ======================================================================

def _ims_attr_str(attr_val):
    """Decode an Imaris HDF5 string attribute to a Python string."""
    if attr_val is None:
        return ""
    if isinstance(attr_val, bytes):
        return attr_val.decode("latin-1", "replace")
    if isinstance(attr_val, np.ndarray):
        if attr_val.dtype.kind == "S":
            return b"".join(attr_val).decode("latin-1", "replace")
        return str(attr_val)
    return str(attr_val)


def read_ims_metadata(filepath):
    """
    Read structural metadata from an IMS file without loading
    all pixel data.  Returns a dict with image properties.
    """
    info = OrderedDict()
    info["filepath"] = filepath
    info["filename"] = os.path.basename(filepath)
    info["file_size_mb"] = os.path.getsize(filepath) / 1.0e6

    try:
        with h5py.File(filepath, "r") as f:

            # --- Image-level metadata ---
            if "DataSetInfo/Image" in f:
                img = f["DataSetInfo/Image"]
                for key in ["X", "Y", "Z", "Unit",
                            "ExtMax0", "ExtMax1", "ExtMax2",
                            "ExtMin0", "ExtMin1", "ExtMin2",
                            "RecordingDate"]:
                    info["Image." + key] = _ims_attr_str(
                        img.attrs.get(key, None))
            else:
                info["_error"] = "Missing DataSetInfo/Image group"
                return info

            # --- Count channels ---
            ch_idx = 0
            info["channel_names"] = []
            info["channel_colors"] = []
            while ("DataSetInfo/Channel %d" % ch_idx) in f:
                ci = f["DataSetInfo/Channel %d" % ch_idx]
                info["channel_names"].append(
                    _ims_attr_str(ci.attrs.get("Name", None)))
                info["channel_colors"].append(
                    _ims_attr_str(ci.attrs.get("Color", None)))
                ch_idx += 1
            info["n_channels_meta"] = ch_idx

            # --- Count timepoints ---
            tp_idx = 0
            while ("DataSet/ResolutionLevel 0/"
                   "TimePoint %d" % tp_idx) in f:
                tp_idx += 1
            info["n_timepoints"] = tp_idx

            # --- Count data channels and get shapes ---
            if tp_idx > 0 and info["n_channels_meta"] > 0:
                tp0 = f["DataSet/ResolutionLevel 0/TimePoint 0"]
                ch0_path = "Channel 0/Data"
                if ch0_path in tp0:
                    shape = tp0[ch0_path].shape
                    info["data_shape"] = "x".join(str(s) for s in shape)
                    info["data_z"] = shape[0] if len(shape) >= 3 else 1
                else:
                    info["data_shape"] = "?"
                    info["data_z"] = 0

            # --- TimeInfo ---
            if "DataSetInfo/TimeInfo" in f:
                ti = f["DataSetInfo/TimeInfo"]
                info["TimeInfo.DatasetTimePoints"] = _ims_attr_str(
                    ti.attrs.get("DatasetTimePoints", None))

            # --- Count resolution levels ---
            rl = 0
            while ("DataSet/ResolutionLevel %d" % rl) in f:
                rl += 1
            info["n_resolution_levels"] = rl

    except Exception as e:
        info["_error"] = str(e)

    return info


def compute_plane_hashes(filepath):
    """
    Compute SHA-256 hashes for every pixel data plane in the file.

    Returns a list of tuples:
        (timepoint, channel, z_index, hash_hex, is_zero, min, max, mean)
    """
    planes = []
    try:
        with h5py.File(filepath, "r") as f:
            tp = 0
            while True:
                tp_path = ("DataSet/ResolutionLevel 0/"
                           "TimePoint %d" % tp)
                if tp_path not in f:
                    break
                ch = 0
                while True:
                    ch_path = tp_path + "/Channel %d/Data" % ch
                    if ch_path not in f:
                        break
                    data = f[ch_path][:]
                    # data shape is (Z, Y, X)
                    if data.ndim == 2:
                        data = data.reshape(1, data.shape[0],
                                            data.shape[1])
                    for z in range(data.shape[0]):
                        plane = data[z]
                        raw = plane.tobytes()
                        h = hashlib.sha256(raw).hexdigest()[:16]
                        is_zero = bool(np.all(plane == 0))
                        planes.append((
                            tp, ch, z, h, is_zero,
                            int(plane.min()), int(plane.max()),
                            float(plane.mean()),
                        ))
                    ch += 1
                tp += 1
    except Exception as e:
        planes.append((-1, -1, -1, "ERROR", False, 0, 0, 0.0))
    return planes


# ======================================================================
#  ANALYSIS
# ======================================================================

def analyze_files(file_list, log_fn=None):
    """
    Run the full validation suite on a list of IMS file paths.

    Parameters
    ----------
    file_list : list of str
        Paths to .ims files.
    log_fn : callable(msg, tag) or None
        Callback for progress logging.

    Returns
    -------
    results : list of dict
        One entry per file with all findings.
    cross_hits : list of tuple
        (file_a, plane_a_desc, file_b, plane_b_desc) for each
        cross-contamination match.
    """
    if log_fn is None:
        def log_fn(msg, tag=""):
            print(msg)

    n = len(file_list)
    log_fn("Phase 1: Reading metadata from %d file(s)..." % n, "info")

    # Build display name map
    name_of = {}  # filepath -> display name
    for fp in file_list:
        name_of[fp] = get_display_name(fp)

    all_meta = []
    for i, fp in enumerate(file_list):
        log_fn("  [%d/%d] %s" % (i + 1, n, name_of[fp]))
        meta = read_ims_metadata(fp)
        all_meta.append(meta)

    log_fn("")
    log_fn("Phase 2: Computing pixel fingerprints...", "info")

    all_planes = {}  # filepath -> list of plane tuples
    hash_index = {}  # hash -> list of (filepath, tp, ch, z)

    for i, fp in enumerate(file_list):
        log_fn("  [%d/%d] %s" % (i + 1, n, name_of[fp]))
        planes = compute_plane_hashes(fp)
        all_planes[fp] = planes

        for (tp, ch, z, h, is_zero, pmin, pmax, pmean) in planes:
            if is_zero or h == "ERROR":
                continue
            key = h
            entry = (fp, tp, ch, z)
            if key not in hash_index:
                hash_index[key] = []
            hash_index[key].append(entry)

    log_fn("")
    log_fn("Phase 3: Analyzing results...", "info")

    # --- Build per-file results ---
    results = []
    cross_hits = []

    for i, fp in enumerate(file_list):
        meta = all_meta[i]
        planes = all_planes.get(fp, [])
        dname = name_of[fp]

        r = OrderedDict()
        r["name"] = dname
        r["file_size_mb"] = "%.2f" % meta.get("file_size_mb", 0)

        # Dimensions
        r["x"] = meta.get("Image.X", "?")
        r["y"] = meta.get("Image.Y", "?")
        r["z_meta"] = meta.get("Image.Z", "?")
        r["z_data"] = str(meta.get("data_z", "?"))
        r["timepoints"] = str(meta.get("n_timepoints", "?"))
        r["channels"] = str(meta.get("n_channels_meta", "?"))
        r["data_shape"] = meta.get("data_shape", "?")

        # Calibration
        ext_x = meta.get("Image.ExtMax0", "0")
        ext_y = meta.get("Image.ExtMax1", "0")
        ext_z = meta.get("Image.ExtMax2", "0")
        r["extent_um"] = "%s x %s x %s" % (ext_x, ext_y, ext_z)
        r["unit"] = meta.get("Image.Unit", "?")

        # --- Check: calibration ---
        cal_issues = []
        try:
            ex = float(ext_x)
            ey = float(ext_y)
            if ex <= 1.0 and ey <= 1.0:
                cal_issues.append("ExtMax near default (1 um)")
            if ex == 0 and ey == 0:
                cal_issues.append("ExtMax is zero")
        except ValueError:
            cal_issues.append("ExtMax not numeric")
        r["calibration_status"] = (
            "; ".join(cal_issues) if cal_issues else "OK")

        # --- Check: phantom dimensions ---
        dim_issues = []
        n_zero = sum(1 for p in planes if p[4])  # is_zero
        n_total = len(planes)
        n_real = n_total - n_zero
        tp_count = meta.get("n_timepoints", 1)
        z_data = meta.get("data_z", 1)

        if n_zero > 0:
            dim_issues.append(
                "%d/%d planes are all zeros" % (n_zero, n_total))
        if tp_count > 1:
            dim_issues.append(
                "%d timepoints (expected 1)" % tp_count)
        if z_data > 1:
            dim_issues.append(
                "Z=%d in data (expected 1 for 2D)" % z_data)

        r["dimension_status"] = (
            "; ".join(dim_issues) if dim_issues else "OK")

        # --- Check: cross-contamination ---
        contam_details = []
        for (tp, ch, z, h, is_zero, pmin, pmax, pmean) in planes:
            if is_zero or h == "ERROR":
                continue
            matches = hash_index.get(h, [])
            for (other_fp, other_tp, other_ch, other_z) in matches:
                if other_fp == fp:
                    continue
                other_name = name_of[other_fp]
                desc_here = "TP%d Ch%d Z%d" % (tp, ch, z)
                desc_there = "TP%d Ch%d Z%d" % (
                    other_tp, other_ch, other_z)
                detail = ("%s -> matches %s %s"
                          % (desc_here, other_name, desc_there))
                contam_details.append(detail)
                cross_hits.append(
                    (dname, desc_here, other_name, desc_there))

        # Deduplicate (each pair appears twice)
        seen = set()
        unique_contam = []
        for d in contam_details:
            if d not in seen:
                seen.add(d)
                unique_contam.append(d)

        if unique_contam:
            r["contamination_status"] = "CONTAMINATED"
            r["contamination_details"] = "; ".join(unique_contam)
        else:
            r["contamination_status"] = "OK"
            r["contamination_details"] = ""

        # --- Overall verdict ---
        problems = []
        if r["calibration_status"] != "OK":
            problems.append("CALIBRATION")
        if r["dimension_status"] != "OK":
            problems.append("DIMENSIONS")
        if r["contamination_status"] != "OK":
            problems.append("CONTAMINATION")
        r["verdict"] = (
            ", ".join(problems) if problems else "PASS")

        results.append(r)

    # --- Log summary ---
    log_fn("")
    n_pass = sum(1 for r in results if r["verdict"] == "PASS")
    n_fail = len(results) - n_pass
    log_fn("=" * 60, "info")
    log_fn("SUMMARY: %d files checked, %d passed, %d flagged"
           % (len(results), n_pass, n_fail), "info")
    log_fn("=" * 60, "info")
    log_fn("")

    for r in results:
        v = r["verdict"]
        dname = r["name"]
        if v == "PASS":
            log_fn("  PASS  %s" % dname, "ok")
            log_fn("         %s, %sch, Z=%s, %sTP, extent=%s"
                   % (r["data_shape"], r["channels"],
                      r["z_data"], r["timepoints"],
                      r["extent_um"]))
        else:
            log_fn("  FAIL  %s  [%s]" % (dname, v), "fail")
            if r["calibration_status"] != "OK":
                log_fn("         Calibration: %s"
                       % r["calibration_status"], "warn")
            if r["dimension_status"] != "OK":
                log_fn("         Dimensions: %s"
                       % r["dimension_status"], "warn")
            if r["contamination_details"]:
                for line in r["contamination_details"].split("; "):
                    log_fn("         Cross-contam: %s" % line,
                           "warn")

    # Deduplicate cross_hits (A->B and B->A)
    seen_pairs = set()
    unique_cross = []
    for (a, ad, b, bd) in cross_hits:
        pair_key = tuple(sorted([(a, ad), (b, bd)]))
        if pair_key not in seen_pairs:
            seen_pairs.add(pair_key)
            unique_cross.append((a, ad, b, bd))

    return results, unique_cross


def write_csv_report(results, cross_hits, output_path):
    """Write the validation results to a CSV file."""
    fieldnames = list(results[0].keys()) if results else []

    with open(output_path, "w") as f:
        # Use newline="" on Py3, not available on Py2
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            writer.writerow(r)

        # Append cross-contamination section
        f.write("\n")
        f.write("Cross-Contamination Matches\n")
        f.write("File A,Plane A,File B,Plane B\n")
        for (a, ad, b, bd) in cross_hits:
            f.write("%s,%s,%s,%s\n" % (a, ad, b, bd))


# ======================================================================
#  GUI  (Tkinter loaded only when main() runs)
# ======================================================================


def _import_tk():
    """Lazy-import Tkinter so the analysis functions can be used
    headlessly (e.g. in a script or notebook)."""
    global tk, filedialog, messagebox, ScrolledText, ttk
    try:
        import Tkinter as tk
        import tkFileDialog as filedialog
        import tkMessageBox as messagebox
        import ScrolledText as _st
        ScrolledText = _st.ScrolledText
        import ttk
    except ImportError:
        import tkinter as tk
        import tkinter.filedialog as filedialog
        import tkinter.messagebox as messagebox
        import tkinter.scrolledtext as _st
        ScrolledText = _st.ScrolledText
        import tkinter.ttk as ttk


class ValidatorApp(object):

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
        self.root.title("IMS Integrity Validator")
        self.root.configure(bg=self.BG)
        self.root.minsize(860, 580)

        if sys.platform == "win32":
            self.FONT_FAMILY = "Consolas"
        else:
            self.FONT_FAMILY = "Menlo"

        self.FONT_MONO = (self.FONT_FAMILY, 10)

        sw = self.root.winfo_screenwidth()
        sh = self.root.winfo_screenheight()
        w, h = min(int(sw * 0.6), 1020), min(int(sh * 0.65), 700)
        x, y = (sw - w) // 2, (sh - h) // 2
        self.root.geometry("%dx%d+%d+%d" % (w, h, x, y))

        self.input_path = tk.StringVar(value="")
        self.filter_var = tk.StringVar(value="")
        self.recurse_var = tk.BooleanVar(value=True)
        self._build_ui()

    def _build_ui(self):
        tk.Label(
            self.root,
            text="IMS Integrity Validator",
            font=(self.FONT_FAMILY, 16, "bold"),
            bg=self.BG, fg=self.ACCENT,
        ).pack(fill="x", padx=16, pady=(14, 2))

        tk.Label(
            self.root,
            text=("Detects cross-contamination, phantom dimensions, "
                  "and calibration errors in Imaris .ims files."),
            font=(self.FONT_FAMILY, 10),
            bg=self.BG, fg="#8888aa", anchor="w",
        ).pack(fill="x", padx=18, pady=(0, 10))

        # Input folder row
        self._entry_row(
            "Arena data folder:",
            self.input_path, self._browse)

        # Experiment filter row
        self._entry_row(
            "Experiment filter (keyword):",
            self.filter_var, None)

        # Help text for filter
        tk.Label(
            self.root,
            text=("  Leave blank to include all .ims files.  "
                  "Otherwise only folders whose name contains "
                  "the keyword are checked."),
            font=(self.FONT_FAMILY, 9),
            bg=self.BG, fg="#666680", anchor="w",
        ).pack(fill="x", padx=18, pady=(0, 6))

        # Options row
        opt_fr = tk.Frame(self.root, bg=self.BG)
        opt_fr.pack(fill="x", padx=16, pady=(0, 6))
        tk.Checkbutton(
            opt_fr, text="Include subfolders",
            variable=self.recurse_var,
            font=(self.FONT_FAMILY, 10),
            bg=self.BG, fg=self.FG, selectcolor=self.BG_FRAME,
            activebackground=self.BG, activeforeground=self.FG,
        ).pack(side="left")

        # Buttons
        btn_fr = tk.Frame(self.root, bg=self.BG)
        btn_fr.pack(fill="x", padx=16, pady=(0, 6))

        self.run_btn = tk.Button(
            btn_fr, text="  Validate  ",
            command=self._run,
            font=(self.FONT_FAMILY, 11, "bold"),
            bg=self.ACCENT, fg="#1e1e2e",
            activebackground=self.ACCENT_HOVER,
            activeforeground="#1e1e2e",
            relief="flat", cursor="hand2", bd=0, padx=14, pady=4,
        )
        self.run_btn.pack(side="left")

        self.clear_btn = tk.Button(
            btn_fr, text="  Clear  ",
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
            "Val.Horizontal.TProgressbar",
            troughcolor=self.BG_FRAME, background=self.ACCENT,
            darkcolor=self.ACCENT, lightcolor=self.ACCENT,
            borderwidth=0,
        )
        self.progress = ttk.Progressbar(
            self.root, orient="horizontal", mode="determinate",
            style="Val.Horizontal.TProgressbar",
        )
        self.progress.pack(fill="x", padx=16, pady=(0, 4))

        # Log
        self.text = ScrolledText(
            self.root, wrap="word", font=self.FONT_MONO,
            bg=self.BG_FRAME, fg=self.FG, insertbackground=self.FG,
            relief="flat", highlightthickness=1,
            highlightbackground="#44446a",
            highlightcolor=self.ACCENT,
            state="disabled", height=18,
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

    def _entry_row(self, label, var, browse_fn):
        fr = tk.Frame(self.root, bg=self.BG)
        fr.pack(fill="x", padx=16, pady=(0, 4))
        tk.Label(
            fr, text=label,
            font=(self.FONT_FAMILY, 10, "bold"),
            bg=self.BG, fg=self.FG,
        ).pack(side="left")
        entry = tk.Entry(
            fr, textvariable=var,
            font=self.FONT_MONO, bg=self.BG_FRAME, fg=self.FG,
            insertbackground=self.FG, relief="flat",
            highlightthickness=1,
            highlightbackground="#44446a",
            highlightcolor=self.ACCENT,
        )
        entry.pack(side="left", fill="x", expand=True,
                   padx=(8, 8), ipady=3)
        if browse_fn:
            tk.Button(
                fr, text=" Browse... ", command=browse_fn,
                font=(self.FONT_FAMILY, 10),
                bg=self.ACCENT, fg="#1e1e2e",
                activebackground=self.ACCENT_HOVER,
                activeforeground="#1e1e2e",
                relief="flat", cursor="hand2", bd=0, padx=8, pady=2,
            ).pack(side="left")

    def _browse(self):
        d = filedialog.askdirectory(
            title="Select folder containing .ims files")
        if d:
            self.input_path.set(d)

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

    def _find_ims_files(self, folder, recurse, keyword):
        found = []
        if recurse:
            for root_dir, dirs, files in os.walk(folder):
                for f in sorted(files):
                    if f.lower().endswith(".ims"):
                        full = os.path.join(root_dir, f)
                        found.append(full)
        else:
            for f in sorted(os.listdir(folder)):
                full = os.path.join(folder, f)
                if f.lower().endswith(".ims"):
                    found.append(full)

        # Apply experiment keyword filter
        if keyword:
            kw = keyword.lower()
            filtered = []
            for fp in found:
                display = get_display_name(fp).lower()
                if kw in display:
                    filtered.append(fp)
            return filtered
        return found

    def _run(self):
        in_dir = self.input_path.get().strip()
        if not in_dir or not os.path.isdir(in_dir):
            messagebox.showwarning(
                "No folder",
                "Please select a valid folder.")
            return

        recurse = self.recurse_var.get()
        keyword = self.filter_var.get().strip()
        files = self._find_ims_files(in_dir, recurse, keyword)

        if not files:
            messagebox.showinfo(
                "No files",
                "No .ims files found in the selected folder.")
            return

        self._clear()
        self.run_btn.configure(state="disabled")

        self._log("IMS Integrity Validator", "info")
        self._log("Folder: %s" % in_dir)
        if keyword:
            self._log("Filter: '%s'" % keyword)
        self._log("Found %d .ims file(s)" % len(files))
        self._log("")

        # Progress: 2 passes (metadata + hashes), so max = 2*n
        self.progress.configure(maximum=len(files) * 2)

        step = [0]
        orig_log = self._log

        def progress_log(msg, tag=""):
            if msg.startswith("  ["):
                step[0] += 1
                self.progress["value"] = step[0]
                self.status_var.set(msg.strip())
                self.root.update_idletasks()
            orig_log(msg, tag)

        results, cross_hits = analyze_files(files, progress_log)

        # Write CSV report next to the input folder
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        csv_name = "ims_validation_%s.csv" % timestamp
        csv_path = os.path.join(in_dir, csv_name)
        try:
            write_csv_report(results, cross_hits, csv_path)
            self._log("")
            self._log("Report saved: %s" % csv_path, "info")
        except Exception as e:
            self._log("Could not save CSV: %s" % str(e), "warn")

        self.progress["value"] = self.progress["maximum"]
        self.status_var.set("Done. %d files checked." % len(files))
        self.run_btn.configure(state="normal")


# ======================================================================
#  ENTRY POINT
# ======================================================================

def main():
    _import_tk()
    root = tk.Tk()
    ValidatorApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
