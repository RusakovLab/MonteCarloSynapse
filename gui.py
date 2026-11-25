"""Desktop GUI for the glutamate diffusion + receptor-kinetics package.

Built on Tkinter (ships with Python -- no extra install) with an embedded
matplotlib canvas. Intended for students and graduate students: set parameters,
run the diffusion Monte Carlo, then run the AMPA/NMDA receptor models on the
result, and view the figures in the same window.

Run with:   python gui.py
"""
from __future__ import annotations

import os
import queue
import sys
import threading
import traceback

# Make "python gui.py" work from any working directory, and from a folder
# unzipped straight off GitHub (no installation required).
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def _preflight():
    """Report anything missing in plain language, instead of a traceback."""
    problems = []
    missing = [name for name in ("numpy", "scipy", "matplotlib")
               if not _importable(name)]
    if missing:
        problems.append(
            "Missing Python package(s): %s\nInstall them with:\n"
            "    %s -m pip install -r requirements.txt"
            % (", ".join(missing), _python_name()))
    if not _importable("tkinter"):
        problems.append(
            "Tkinter (the graphical toolkit) is not available in this Python.\n"
            "On Linux install it with:  sudo apt install python3-tk\n"
            "On macOS/Windows reinstall Python from https://www.python.org/ "
            "(Tkinter is included there).")
    if not problems:
        return
    msg = "\n\n".join(problems)
    print(msg, file=sys.stderr)
    try:
        import tkinter
        from tkinter import messagebox as _mb
        _root = tkinter.Tk()
        _root.withdraw()
        _mb.showerror("Cannot start the GUI", msg)
        _root.destroy()
    except Exception:
        pass
    raise SystemExit(1)


def _importable(name):
    try:
        __import__(name)
    except ImportError:
        return False
    return True


def _python_name():
    return os.path.basename(sys.executable) or "python"


_preflight()

import numpy as np

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, simpledialog

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

from glutamate import default_config, quick_config
from glutamate import diffusion, receptors_ampa, receptors_nmda, io_files, project, video


# --- colour palette ---
BG = "#f2f4f7"
CARD = "#ffffff"
ACCENT = "#2d6cdf"
ACCENT_DARK = "#1e51b0"
TEXT = "#243447"


class ToolTip:
    """Small hover popup that describes a widget."""
    def __init__(self, widget, text, wrap=340):
        self.widget = widget
        self.text = text
        self.wrap = wrap
        self.tip = None
        widget.bind("<Enter>", self._show, add="+")
        widget.bind("<Leave>", self._hide, add="+")
        widget.bind("<ButtonPress>", self._hide, add="+")

    def _show(self, _e=None):
        if self.tip or not self.text:
            return
        x = self.widget.winfo_rootx() + 18
        y = self.widget.winfo_rooty() + self.widget.winfo_height() + 6
        self.tip = tw = tk.Toplevel(self.widget)
        tw.wm_overrideredirect(True)
        tw.wm_geometry("+%d+%d" % (x, y))
        tk.Label(tw, text=self.text, justify="left", background="#2b3a4a",
                 foreground="white", relief="solid", borderwidth=0,
                 font=("Segoe UI", 9), wraplength=self.wrap, padx=8, pady=6).pack()

    def _hide(self, _e=None):
        if self.tip:
            self.tip.destroy()
            self.tip = None


# (label, group, attribute, kind)  kind in {int,float,str,bool,list}
FORM = [
    ("Domain & core", [
        ("Particles (np)", "basic", "NumberOfParticles", "int"),
        ("Total time (ms)", "basic", "TotalTime", "float"),
        ("Trials", "basic", "Trials", "int"),
        ("Release source", "basic", "ReleaseSource", "choice:synapse,astrocyte"),
        ("Diffusion coef (nm^2/ms)", "computation", "DifCoefINIT", "float"),
        ("Brownian step (ms)", "computation", "GTimeStep", "float"),
        ("Distance bins", "computation", "DistanceBins", "int"),
        ("Time bins", "computation", "TimeBins", "int"),
        ("Cube size (nm)", "cube", "Size", "float"),
        ("Snapshot times (ms, csv)", "computation", "SnapshotTimes", "list"),
        ("Diff. window begin (frac)", "basic", "BeginTime", "float"),
        ("Diff. window end (frac)", "basic", "EndTime", "float"),
    ]),
    ("Release schedule", [
        ("Release times (ms, csv)", "basic", "ReleaseTimes", "list"),
        ("Release fractions (csv)", "basic", "ReleaseFractions", "list"),
        ("Astrocyte radius (nm)", "basic", "AstrocyteRadius", "float"),
        ("Astrocyte centre dist (nm)", "basic", "AstrocyteCenterDistance", "float"),
        ("Release patch radius (nm)", "basic", "ReleasePatchRadius", "float"),
    ]),
    ("Obstacles", [
        ("Number of obstacles", "basic", "Numbersphere", "int"),
        ("Min ball radius (nm)", "balls", "MinBallRadius", "float"),
        ("Max ball radius (nm)", "balls", "MaxBallRadius", "float"),
        ("Placement mode", "balls", "PlacementMode", "choice:overlapping,nonoverlapping"),
        ("Adhesive fraction", "balls", "ProbabilityofAstrocytes", "float"),
        ("Occupancy samples", "balls", "OccupancySampleCount", "int"),
        ("Adhesion decay", "balls", "UseAdhesionDecay", "bool"),
        ("Adhesion decay time (ms)", "balls", "AdhesionDecayTime", "float"),
    ]),
    ("Binding & removal", [
        ("Max adhesion prob", "adhesion", "MaxProbAdhesive", "float"),
        ("Adhesion time const (ms)", "adhesion", "TimeINsideAdhesiveZone", "float"),
        ("Unbind prob", "basic", "UnboundProb", "float"),
        ("Exponential removal", "basic", "UseExponentialDecay", "bool"),
        ("Removal time (ms)", "basic", "TimeDecay", "float"),
    ]),
    ("Cleft geometry", [
        ("Cleft radius (nm)", "cleft", "RadiusofCleft", "float"),
        ("Cleft width (nm)", "cleft", "CleftWidth", "float"),
        ("Synaptic envelope (nm)", "cleft", "SynapticEnvelopeRadius", "float"),
        ("Spillover radius (nm)", "cleft", "SpilloverRadius", "float"),
    ]),
    ("Receptors", [
        ("AMPA duration (ms)", "ampa", "TotalTime_ms", "float"),
        ("AMPA output step (ms)", "ampa", "OutputTimeStep_ms", "float"),
        ("Plot distance (um)", "ampa", "MaxDistance_um", "float"),
        ("NMDA duration (ms)", "nmda", "TotalTime_ms", "float"),
    ]),
]


# Hover descriptions, keyed by (group, attribute).
DESCRIPTIONS = {
    ("basic", "NumberOfParticles"): "Total number of glutamate molecules released (summed over all release events).",
    ("basic", "TotalTime"): "Total simulated diffusion time (ms). Number of steps = TotalTime / Brownian step.",
    ("basic", "Trials"): "Number of independent Monte Carlo trials. Receptor models run on the average over all trials.",
    ("basic", "ReleaseSource"): "Where glutamate is released: 'synapse' (point source in the cleft) or 'astrocyte' (surface patch).",
    ("computation", "DifCoefINIT"): "Free diffusion coefficient of glutamate (nm^2/ms). For example = 0.4e6 = 0.4 um^2/ms.",
    ("computation", "GTimeStep"): "Integration time step (ms). Smaller = more accurate but slower. RMS step = sqrt(6*D*dt).",
    ("computation", "DistanceBins"): "Number of radial shells (cleft radius -> domain half-size) for the concentration and receptor maps.",
    ("computation", "TimeBins"): "Number of time samples of the radial bound/free histograms over the total time.",
    ("cube", "Size"): "Side length of the cubic simulation domain (nm). Half of it is the radial analysis limit.",
    ("computation", "SnapshotTimes"): "Times (ms, comma-separated) at which particle-position snapshots (PD files) are saved.",
    ("basic", "BeginTime"): "Start of the window (fraction of Total time) used to estimate the effective diffusion coefficient.",
    ("basic", "EndTime"): "End of the window (fraction of Total time) for the effective-diffusion estimate.",
    ("basic", "ReleaseTimes"): "Times (ms) of each glutamate release, comma-separated. Must start at 0 and strictly increase.",
    ("basic", "ReleaseFractions"): "Fraction of particles released at each release time. Must be positive and sum to 1.",
    ("basic", "AstrocyteRadius"): "Radius of the primary astrocyte sphere (astrocyte source only).",
    ("basic", "AstrocyteCenterDistance"): "Distance from the origin to the astrocyte centre (astrocyte source only).",
    ("basic", "ReleasePatchRadius"): "Radius of the release patch on the astrocyte surface (astrocyte source only).",
    ("basic", "Numbersphere"): "Number of spherical extracellular obstacles. More obstacles -> higher occupancy (less free space).",
    ("balls", "MinBallRadius"): "Minimum obstacle radius (nm). Radii are drawn uniformly between min and max.",
    ("balls", "MaxBallRadius"): "Maximum obstacle radius (nm).",
    ("balls", "PlacementMode"): "'overlapping' (obstacles may intersect, standard) or 'nonoverlapping' (enforce a minimum gap).",
    ("balls", "ProbabilityofAstrocytes"): "Fraction of obstacles that are adhesive (can bind glutamate), 0..1.",
    ("balls", "OccupancySampleCount"): "Number of Monte Carlo points used to estimate the occupied/free volume fraction.",
    ("balls", "UseAdhesionDecay"): "If on, adhesive obstacles gradually lose their adhesiveness over time.",
    ("balls", "AdhesionDecayTime"): "Time constant (ms) for adhesive obstacles losing adhesiveness.",
    ("adhesion", "MaxProbAdhesive"): "Asymptotic binding probability inside an adhesive obstacle (1 = certain, 0 = none).",
    ("adhesion", "TimeINsideAdhesiveZone"): "Time constant (ms) to bind while inside an adhesive obstacle. Shorter = binds faster.",
    ("basic", "UnboundProb"): "Probability that an adhered molecule can later unbind (0 = permanent until removal).",
    ("basic", "UseExponentialDecay"): "If on, bound glutamate is irreversibly removed (uptake) with an exponential hazard.",
    ("basic", "TimeDecay"): "Time constant (ms) of the exponential removal of bound glutamate.",
    ("cleft", "RadiusofCleft"): "Radius of the cylindrical synaptic cleft (nm). The receptor shells start here.",
    ("cleft", "CleftWidth"): "Distance between the two cleft membranes (nm). Molecules reflect off them in z.",
    ("cleft", "SynapticEnvelopeRadius"): "Radius of the solid pre/postsynaptic exclusion sphere; obstacles and molecules stay outside it (except the cleft gap).",
    ("cleft", "SpilloverRadius"): "Radial distance defining spillover: molecules beyond it are counted as escaped.",
    ("ampa", "TotalTime_ms"): "Integration horizon (ms) for the AMPA receptor ODEs.",
    ("ampa", "OutputTimeStep_ms"): "Temporal resolution (ms) of the saved AMPA response grid.",
    ("ampa", "MaxDistance_um"): "Maximum radial distance (um) shown on the receptor maps.",
    ("nmda", "TotalTime_ms"): "Integration horizon (ms) for the NMDA receptor ODEs.",
}


class App:
    def __init__(self, root):
        self.root = root
        root.title("Glutamate Diffusion & Receptor Kinetics")
        root.geometry("1400x800")
        self.cfg = default_config()
        self.vars = {}
        self.q = queue.Queue()
        self.diff_result = None
        self.trials_out = None
        self.results_folder = None
        self.last_result = None
        self._pending_frame = None
        self._field_balls = None
        self.receptor_results = {}
        self._build()
        self._load_from_cfg()
        self.root.after(120, self._poll)

    # ---------------- theme ----------------
    def _style(self):
        style = ttk.Style()
        try:
            style.theme_use("clam")
        except Exception:
            pass
        self.root.configure(bg=BG)
        base = ("Segoe UI", 10)
        style.configure(".", background=BG, foreground=TEXT, font=base)
        style.configure("TFrame", background=BG)
        style.configure("TLabel", background=BG, foreground=TEXT)
        style.configure("TCheckbutton", background=BG)
        style.configure("TLabelframe", background=BG, relief="solid", borderwidth=1)
        style.configure("TLabelframe.Label", background=BG, foreground=ACCENT,
                        font=("Segoe UI", 10, "bold"))
        style.configure("TButton", padding=(10, 5))
        style.map("TButton", background=[("active", "#e6edf9")])
        style.configure("Accent.TButton", padding=(12, 6),
                        font=("Segoe UI", 10, "bold"), foreground="white",
                        background=ACCENT)
        style.map("Accent.TButton",
                  background=[("pressed", ACCENT_DARK), ("active", ACCENT_DARK),
                              ("disabled", "#a9bce0")],
                  foreground=[("disabled", "#eef2fb")])
        style.configure("Header.TLabel", background=BG, foreground="#16233a",
                        font=("Segoe UI", 15, "bold"))
        style.configure("Sub.TLabel", background=BG, foreground="#6b7787",
                        font=("Segoe UI", 9))
        style.configure("TProgressbar", background=ACCENT, troughcolor="#dfe4ea")
        style.configure("TCombobox", padding=3)
        style.configure("TEntry", padding=3)

    # ---------------- layout ----------------
    def _build(self):
        self._style()
        header = ttk.Frame(self.root)
        header.pack(fill="x", padx=12, pady=(8, 2))
        ttk.Label(header, text="Glutamate Diffusion & Receptor Kinetics",
                  style="Header.TLabel").pack(side="left")
        ttk.Label(header, text="Monte Carlo field  +  AMPA / NMDA kinetics",
                  style="Sub.TLabel").pack(side="left", padx=14)
        ttk.Separator(self.root, orient="horizontal").pack(fill="x", padx=10, pady=(2, 4))

        paned = ttk.PanedWindow(self.root, orient="horizontal")
        paned.pack(fill="both", expand=True)

        left = ttk.Frame(paned, width=660)
        paned.add(left, weight=0)
        right = ttk.Frame(paned)
        paned.add(right, weight=1)

        # ---- scrollable two-column parameter form ----
        canvas = tk.Canvas(left, borderwidth=0, width=650, background=BG, highlightthickness=0)
        vs = ttk.Scrollbar(left, orient="vertical", command=canvas.yview)
        form = ttk.Frame(canvas)
        form.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=form, anchor="nw")
        canvas.configure(yscrollcommand=vs.set)
        canvas.pack(side="left", fill="both", expand=True)
        vs.pack(side="right", fill="y")
        canvas.bind_all("<MouseWheel>", lambda e: canvas.yview_scroll(int(-e.delta / 120), "units"))

        columns = [ttk.Frame(form), ttk.Frame(form)]
        columns[0].grid(row=0, column=0, sticky="n", padx=(0, 4))
        columns[1].grid(row=0, column=1, sticky="n")

        total = sum(len(f) for _, f in FORM)
        running = 0
        col = 0
        for section, fields in FORM:
            if col == 0 and running >= total / 2:
                col = 1
            running += len(fields)
            lf = ttk.LabelFrame(columns[col], text=section)
            lf.pack(fill="x", padx=6, pady=4)
            for i, (label, grp, attr, kind) in enumerate(fields):
                lab = ttk.Label(lf, text=label)
                lab.grid(row=i, column=0, sticky="w", padx=4, pady=2)
                key = (grp, attr)
                if kind == "bool":
                    var = tk.BooleanVar()
                    w = ttk.Checkbutton(lf, variable=var)
                elif kind.startswith("choice:"):
                    var = tk.StringVar()
                    opts = kind.split(":", 1)[1].split(",")
                    w = ttk.Combobox(lf, textvariable=var, values=opts, width=13, state="readonly")
                else:
                    var = tk.StringVar()
                    w = ttk.Entry(lf, textvariable=var, width=15)
                w.grid(row=i, column=1, sticky="w")
                self.vars[key] = (var, kind)
                desc = DESCRIPTIONS.get(key, "")
                if desc:
                    ToolTip(lab, desc)
                    ToolTip(w, desc)

        # ---- right: project row + controls + plot + log ----
        self.folder_var = tk.StringVar(value="(no results folder - runs won't be saved)")
        proj = ttk.Frame(right)
        proj.pack(fill="x", padx=6, pady=(6, 0))
        ttk.Button(proj, text="Results folder...", command=self._choose_folder).pack(side="left")
        ttk.Button(proj, text="New folder...", command=self._new_folder).pack(side="left", padx=4)
        ttk.Button(proj, text="Load params", command=self._load_params).pack(side="left", padx=4)
        ttk.Button(proj, text="Save params", command=self._save_params).pack(side="left")
        ttk.Button(proj, text="Load results", command=self._load_results).pack(side="left", padx=(10, 0))
        ttk.Label(proj, textvariable=self.folder_var).pack(side="left", padx=8)

        ctl = ttk.Frame(right)
        ctl.pack(fill="x", padx=6, pady=6)
        ttk.Button(ctl, text="Quick preset", command=self._preset_quick).pack(side="left")
        ttk.Button(ctl, text="Full preset", command=self._preset_full).pack(side="left", padx=4)
        self.seed_var = tk.StringVar(value="")
        ttk.Label(ctl, text="seed:").pack(side="left", padx=(12, 2))
        ttk.Entry(ctl, textvariable=self.seed_var, width=8).pack(side="left")
        self.live_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(ctl, text="Live 3D field (slow)",
                        variable=self.live_var).pack(side="left", padx=(12, 0))
        self.ballshow_var = tk.StringVar(value="both")
        ttk.Label(ctl, text="balls:").pack(side="left", padx=(8, 2))
        ttk.Combobox(ctl, textvariable=self.ballshow_var, width=12, state="readonly",
                     values=["both", "adhesive", "non-adhesive"]).pack(side="left")

        run = ttk.Frame(right)
        run.pack(fill="x", padx=6)
        self.btn_diff = ttk.Button(run, text="Run diffusion", command=self._run_diffusion,
                                   style="Accent.TButton")
        self.btn_diff.pack(side="left")
        self.model_var = tk.StringVar(value="AMPA")
        ttk.Combobox(run, textvariable=self.model_var, values=["AMPA", "AMPA1", "AMPA2"],
                     width=7, state="readonly").pack(side="left", padx=(12, 2))
        self.btn_ampa = ttk.Button(run, text="Run AMPA", command=self._run_ampa,
                                   state="disabled", style="Accent.TButton")
        self.btn_ampa.pack(side="left")
        self.btn_nmda = ttk.Button(run, text="Run NMDA", command=self._run_nmda,
                                   state="disabled", style="Accent.TButton")
        self.btn_nmda.pack(side="left", padx=4)
        ttk.Button(run, text="Save results...", command=self._save).pack(side="left", padx=4)
        self.btn_video = ttk.Button(run, text="Export video...", command=self._export_video)
        self.btn_video.pack(side="left", padx=4)

        prog = ttk.Frame(right)
        prog.pack(fill="x", padx=6, pady=4)
        self.progress = ttk.Progressbar(prog, mode="determinate", maximum=100)
        self.progress.pack(side="left", fill="x", expand=True)
        self.view_var = tk.StringVar(value="Concentration")
        self.view_box = ttk.Combobox(prog, textvariable=self.view_var, width=16, state="readonly",
                                     values=["Concentration", "Open state", "Desensitization",
                                             "Spillover", "Particles"])
        self.view_box.pack(side="left", padx=6)
        self.view_box.bind("<<ComboboxSelected>>", lambda e: self._on_view_change())
        ttk.Label(prog, text="receptor:").pack(side="left", padx=(10, 2))
        self.recv_var = tk.StringVar(value="")
        self.recv_box = ttk.Combobox(prog, textvariable=self.recv_var, width=8,
                                     state="readonly", values=[])
        self.recv_box.pack(side="left")
        self.recv_box.bind("<<ComboboxSelected>>", lambda e: self._show_receptor())

        # ---- colour-scale (colorbar min/max) controls ----
        self._scale_manual = False
        scale = ttk.Frame(right)
        scale.pack(fill="x", padx=6)
        ttk.Label(scale, text="colour scale").pack(side="left")
        ttk.Label(scale, text="min:").pack(side="left", padx=(10, 2))
        self.cmin_var = tk.StringVar(value="")
        emin = ttk.Entry(scale, textvariable=self.cmin_var, width=12)
        emin.pack(side="left")
        ttk.Label(scale, text="max:").pack(side="left", padx=(8, 2))
        self.cmax_var = tk.StringVar(value="")
        emax = ttk.Entry(scale, textvariable=self.cmax_var, width=12)
        emax.pack(side="left")
        emin.bind("<Return>", lambda e: self._apply_scale())
        emax.bind("<Return>", lambda e: self._apply_scale())
        ttk.Button(scale, text="Apply", command=self._apply_scale).pack(side="left", padx=6)
        ttk.Button(scale, text="Auto", command=self._auto_scale).pack(side="left")
        self.logscale_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(scale, text="log scale", variable=self.logscale_var,
                        command=self._apply_scale).pack(side="left", padx=(12, 0))
        ToolTip(emin, "Lower limit of the colour scale for the current map. "
                      "Edit and press Enter (or Apply) to fix it; press Auto to recompute.")
        ToolTip(emax, "Upper limit of the colour scale for the current map. "
                      "Edit and press Enter (or Apply) to fix it; press Auto to recompute.")

        self.fig = Figure(figsize=(6.5, 4.6))
        self.canvas = FigureCanvasTkAgg(self.fig, master=right)
        self.canvas.get_tk_widget().pack(fill="both", expand=True, padx=6, pady=4)
        NavigationToolbar2Tk(self.canvas, right)

        self.log = tk.Text(right, height=6)
        self.log.pack(fill="x", padx=6, pady=4)

    # ---------------- config <-> form ----------------
    def _grp(self, name):
        return getattr(self.cfg, name)

    def _load_from_cfg(self):
        for (grp, attr), (var, kind) in self.vars.items():
            val = getattr(self._grp(grp), attr)
            if kind == "bool":
                var.set(bool(val))
            elif kind == "list":
                var.set(", ".join(str(x) for x in val))
            else:
                var.set(str(val))

    def _apply_to_cfg(self):
        for (grp, attr), (var, kind) in self.vars.items():
            raw = var.get()
            if kind == "int":
                val = int(float(raw))
            elif kind == "float":
                val = float(raw)
            elif kind == "bool":
                val = bool(var.get())
            elif kind == "list":
                val = [float(x) for x in raw.replace(";", ",").split(",") if x.strip()]
            else:
                val = raw
            setattr(self._grp(grp), attr, val)
        # recompute derived + linked values
        self.cfg.computation.__post_init__()
        self.cfg.cube.__post_init__()
        self.cfg.ampa.TimeGlu_ms = self.cfg.basic.TotalTime
        self.cfg.ampa.BinGlu = self.cfg.computation.DistanceBins

    # ---------------- presets ----------------
    def _preset_quick(self):
        self.cfg = quick_config()
        self._load_from_cfg()
        self._say("Quick preset loaded (fast, for trying the workflow).")

    def _preset_full(self):
        self.cfg = default_config()
        self._load_from_cfg()
        self._say("Full preset loaded (paper-scale; slow).")

    # ---------------- run (threaded) ----------------
    def _busy(self, on):
        state = "disabled" if on else "normal"
        self.btn_diff.config(state=state)
        self.btn_ampa.config(state="disabled" if on or self.diff_result is None else "normal")
        self.btn_nmda.config(state="disabled" if on or self.diff_result is None else "normal")

    def _seed(self):
        s = self.seed_var.get().strip()
        return int(s) if s else None

    def _run_diffusion(self):
        try:
            self._apply_to_cfg()
        except Exception as e:
            messagebox.showerror("Invalid parameters", str(e)); return
        self._busy(True)
        self.progress["value"] = 0
        n = int(self.cfg.basic.Trials)
        folder = self.results_folder
        live = self.live_var.get()
        self._field_balls = None
        self._say("Running %d trial(s)%s%s..." % (
            n, (" -> " + folder) if folder else " (choose a Results folder to save them)",
            " [live 3D field]" if live else ""))

        def work():
            try:
                out = project.run_trials(
                    self.cfg, folder=folder, seed=self._seed(), save=bool(folder),
                    progress_cb=lambda p: self.q.put(("progress", p)),
                    frame_cb=(lambda f: self.q.put(("frame", f))) if live else None)
                self.q.put(("diffusion_done", out))
            except Exception:
                self.q.put(("error", traceback.format_exc()))
        threading.Thread(target=work, daemon=True).start()

    def _run_ampa(self):
        if self.trials_out is None:
            return
        self._busy(True); self._say("Running %s (mean of %d trial(s))..." % (
            self.model_var.get(), self.trials_out["n_trials"]))
        meta = self.trials_out["metadata"]
        free = self.trials_out["average_free"]
        model = self.model_var.get()

        def work():
            try:
                res = receptors_ampa.run_ampa_spatial(
                    self.cfg, model_type=model, initial_distribution=free,
                    shell_volume_um3=meta["ShellVolume_um3"],
                    distance_um=meta["OuterRadius_um"])
                self.q.put(("ampa_done", res))
            except Exception:
                self.q.put(("error", traceback.format_exc()))
        threading.Thread(target=work, daemon=True).start()

    def _run_nmda(self):
        if self.trials_out is None:
            return
        self._busy(True); self._say("Running NMDA (mean of %d trial(s))..." % self.trials_out["n_trials"])
        meta = self.trials_out["metadata"]
        free = self.trials_out["average_free"]

        def work():
            try:
                res = receptors_nmda.run_nmda_spatial(
                    self.cfg, initial_distribution=free,
                    shell_volume_um3=meta["ShellVolume_um3"],
                    glu_time_ms=meta["TimePoints_ms"],
                    distance_um=meta["OuterRadius_um"])
                self.q.put(("nmda_done", res))
            except Exception:
                self.q.put(("error", traceback.format_exc()))
        threading.Thread(target=work, daemon=True).start()

    # ---------------- video export ----------------
    def _video_dialog(self):
        top = tk.Toplevel(self.root)
        top.title("Export video")
        top.transient(self.root)
        top.grab_set()
        top.configure(bg=BG)
        fields = [("Frames", "frames", "80"),
                  ("Size (px, square)", "size", "540"),
                  ("FPS", "fps", "20"),
                  ("Target size (MB)", "mb", "1.0"),
                  ("Camera pan (deg)", "span", "90"),
                  ("Scale bar (um, 0=off)", "scale", "1")]
        ents = {}
        for i, (lab, key, default) in enumerate(fields):
            ttk.Label(top, text=lab).grid(row=i, column=0, sticky="w", padx=10, pady=4)
            v = tk.StringVar(value=default)
            ttk.Entry(top, textvariable=v, width=12).grid(row=i, column=1, padx=10, pady=4)
            ents[key] = v
        rot = tk.BooleanVar(value=True)
        ttk.Checkbutton(top, text="Rotate camera", variable=rot).grid(
            row=len(fields), column=0, columnspan=2, sticky="w", padx=10, pady=(2, 0))
        tim = tk.BooleanVar(value=True)
        ttk.Checkbutton(top, text="Time counter", variable=tim).grid(
            row=len(fields) + 1, column=0, columnspan=2, sticky="w", padx=10, pady=(0, 6))
        result = {}

        def ok():
            try:
                result.update(frames=int(float(ents["frames"].get())),
                              size=int(float(ents["size"].get())),
                              fps=int(float(ents["fps"].get())),
                              mb=float(ents["mb"].get()),
                              span=float(ents["span"].get()),
                              scale=float(ents["scale"].get()),
                              rotate=bool(rot.get()),
                              show_time=bool(tim.get()),
                              ball_mode=self.ballshow_var.get())
            except ValueError:
                messagebox.showerror("Invalid", "Please enter numbers.")
                return
            top.destroy()

        btnf = ttk.Frame(top)
        btnf.grid(row=len(fields) + 2, column=0, columnspan=2, pady=10)
        ttk.Button(btnf, text="OK", command=ok, style="Accent.TButton").pack(side="left", padx=4)
        ttk.Button(btnf, text="Cancel", command=top.destroy).pack(side="left", padx=4)
        top.wait_window()
        return result or None

    def _export_video(self):
        try:
            self._apply_to_cfg()
        except Exception as e:
            messagebox.showerror("Invalid parameters", str(e)); return
        opts = self._video_dialog()
        if not opts:
            return
        init = self.results_folder or ""
        path = filedialog.asksaveasfilename(
            title="Save video as", initialdir=init,
            initialfile="StochasticSynapse.mp4", defaultextension=".mp4",
            filetypes=[("MP4 video", "*.mp4"), ("WebM", "*.webm"), ("GIF", "*.gif")])
        if not path:
            return
        self.btn_video.config(state="disabled")
        self.progress["value"] = 0
        self._say("Rendering video (%d frames, %d px, ~%.1f MB target)... this runs a fresh diffusion."
                  % (opts["frames"], opts["size"], opts["mb"]))
        cfg = self.cfg
        seed = self._seed()

        def work():
            try:
                size = video.export_diffusion_video(
                    cfg, path, n_frames=opts["frames"], size_px=opts["size"],
                    fps=opts["fps"], target_mb=opts["mb"], seed=seed,
                    rotate=opts["rotate"], azim_span=opts["span"],
                    ball_mode=opts["ball_mode"], scale_um=opts["scale"],
                    show_time=opts["show_time"],
                    progress_cb=lambda p: self.q.put(("progress", p)))
                self.q.put(("video_done", (path, size)))
            except Exception:
                self.q.put(("error", traceback.format_exc()))
        threading.Thread(target=work, daemon=True).start()

    # ---------------- queue polling ----------------
    def _poll(self):
        try:
            while True:
                kind, payload = self.q.get_nowait()
                if kind == "progress":
                    self.progress["value"] = payload
                elif kind == "diffusion_done":
                    self.trials_out = payload
                    self.diff_result = payload["results"][0]
                    self.progress["value"] = 100
                    self._busy(False)
                    self._say("Diffusion done: %d trial(s). Occupancy=%.2f%% D_mean=%.3g. "
                              "Now run a receptor model." % (payload["n_trials"],
                              100 * payload["occupied_fraction"], payload["mean_diffusion"]))
                    self.view_var.set("Particles"); self._redraw()
                elif kind == "ampa_done":
                    self._busy(False)
                    self._register_receptor(payload.model_type, payload)
                    saved = self._save_receptor(payload, payload.model_type)
                    self._say("%s done. open max=%.4f%s" % (
                        payload.model_type, payload.open_state.max(), " (saved)" if saved else ""))
                    self.view_var.set("Open state"); self._redraw()
                elif kind == "nmda_done":
                    self._busy(False)
                    self._register_receptor("NMDA", payload)
                    saved = self._save_receptor(payload, "NMDA")
                    self._say("NMDA done. open max=%.4f%s" % (
                        payload.open_state.max(), " (saved)" if saved else ""))
                    self.view_var.set("Open state"); self._redraw()
                elif kind == "video_done":
                    path, size = payload
                    self.progress["value"] = 100
                    if hasattr(self, "btn_video"):
                        self.btn_video.config(state="normal")
                    self._say("Video saved: %s (%.2f MB)" % (path, size / 1e6))
                    messagebox.showinfo("Video exported",
                                        "Saved:\n%s\n\nSize: %.2f MB" % (path, size / 1e6))
                elif kind == "error":
                    self._busy(False)
                    if hasattr(self, "btn_video"):
                        self.btn_video.config(state="normal")
                    self._say("ERROR:\n" + payload)
                    messagebox.showerror("Run failed", payload.splitlines()[-1])
                elif kind == "frame":
                    self._pending_frame = payload
        except queue.Empty:
            pass
        if self._pending_frame is not None:
            fr = self._pending_frame
            self._pending_frame = None
            self._draw_field(fr)
        self.root.after(120, self._poll)

    def _draw_field(self, frame):
        try:
            if "balls" in frame:
                self._field_balls = frame["balls"]
            self.fig.clear()
            ax = self.fig.add_subplot(111, projection="3d")
            b = self._field_balls
            if b is not None and b["xyz"].shape[0] > 0:
                xyz = b["xyz"]; r = b["radius"]; adh = b["adhesive"]
                lim = self.cfg.cube.AxisSize
                # Marker size ~ physical diameter (a 400 nm ball looks large,
                # a 50 nm one small) relative to the +/-lim cube.
                ppn = (self.fig.get_size_inches()[0] * 72.0 * 0.55) / (2.0 * lim)
                s = (2.0 * r * ppn) ** 2
                mode = self.ballshow_var.get()
                if mode in ("both", "adhesive"):
                    ax.scatter(xyz[adh, 0], xyz[adh, 1], xyz[adh, 2], s=s[adh],
                               c="red", alpha=0.18, edgecolors="none")
                if mode in ("both", "non-adhesive"):
                    ax.scatter(xyz[~adh, 0], xyz[~adh, 1], xyz[~adh, 2], s=s[~adh],
                               c="royalblue", alpha=0.12, edgecolors="none")
            x = frame["x"]; y = frame["y"]; z = frame["z"]
            act = frame["active"]; bnd = frame["bound"]
            free = act & ~bnd
            ax.scatter(x[free], y[free], z[free], s=8, c="blue")
            ax.scatter(x[act & bnd], y[act & bnd], z[act & bnd], s=16, c="red")
            lim = self.cfg.cube.AxisSize
            ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
            ax.set_xlabel("X (nm)"); ax.set_ylabel("Y (nm)"); ax.set_zlabel("Z (nm)")
            ax.set_title("t = %.3f ms | active = %d | bound = %d"
                         % (frame["t"], int(act.sum()), int((act & bnd).sum())))
            self.canvas.draw()
        except Exception:
            self._say("Field draw error:\n" + traceback.format_exc())

    # ---------------- colour-scale helpers ----------------
    def _on_view_change(self):
        # switching the view returns to automatic colour limits
        self._scale_manual = False
        self._redraw()

    def _apply_scale(self):
        # user typed explicit limits -> honour them
        self._scale_manual = True
        self._redraw()

    def _auto_scale(self):
        self._scale_manual = False
        self._redraw()

    def _color_limits(self, auto_vmin, auto_vmax):
        """Return (vmin, vmax) for the current map.

        In auto mode the entry boxes are refreshed with the computed values so
        the user can see and then tweak them. In manual mode the typed values
        are used (falling back to auto on blank/invalid input).
        """
        if not self._scale_manual:
            self.cmin_var.set("%.4g" % auto_vmin)
            self.cmax_var.set("%.4g" % auto_vmax)
            return auto_vmin, auto_vmax
        try:
            vmin = float(self.cmin_var.get())
        except (ValueError, TypeError):
            vmin = auto_vmin
        try:
            vmax = float(self.cmax_var.get())
        except (ValueError, TypeError):
            vmax = auto_vmax
        if vmax <= vmin:
            vmin, vmax = auto_vmin, auto_vmax
        return vmin, vmax

    # ---------------- plotting ----------------
    def _redraw(self):
        self.fig.clear()
        view = self.view_var.get()
        r = self.last_result
        try:
            if view == "Particles" and self.diff_result is not None:
                snaps = self.diff_result.snapshots
                if snaps:
                    arr = list(snaps.values())[-1]
                    ax = self.fig.add_subplot(111, projection="3d")
                    X, Y, Z, S = arr
                    free = S == 1
                    ax.scatter(X[~free], Y[~free], Z[~free], s=5, alpha=0.15, c="steelblue")
                    ax.scatter(X[free], Y[free], Z[free], s=16, alpha=0.9, c="crimson")
                    lim = self.cfg.cube.AxisSize
                    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
                    ax.set_title("Particles (red=free, blue=bound)")
            elif view == "Spillover" and self.diff_result is not None:
                ax = self.fig.add_subplot(111)
                d = self.diff_result
                ax.plot(d.time_ms, d.fraction_spillover, label="all")
                ax.plot(d.time_ms, d.fraction_free_spillover, label="free")
                ax.set_xlabel("Time (ms)"); ax.set_ylabel("Spillover fraction")
                ax.legend(); ax.grid(True)
            elif view == "Concentration" and self.trials_out is not None:
                from matplotlib.colors import LogNorm, Normalize
                from glutamate.concentration import counts_to_uM
                meta = self.trials_out["metadata"]
                free = self.trials_out["average_free"]
                conc = counts_to_uM(free[:, 1:], np.asarray(meta["ShellVolume_um3"], float))
                dist = np.asarray(meta["OuterRadius_um"], float)
                tp = np.asarray(meta["TimePoints_ms"], float)
                ax = self.fig.add_subplot(111)
                pos = conc[conc > 0]
                auto_vmax = float(pos.max()) if pos.size else 1.0
                auto_vmin = auto_vmax * 1e-4
                vmin, vmax = self._color_limits(auto_vmin, auto_vmax)
                logscale = self.logscale_var.get()
                if logscale:
                    vmin = max(vmin, vmax * 1e-12)  # LogNorm needs vmin > 0
                    norm = LogNorm(vmin=vmin, vmax=vmax)
                    plot = np.where(conc > vmin, conc, np.nan)
                else:
                    norm = Normalize(vmin=vmin, vmax=vmax)
                    plot = conc
                cmap = matplotlib.colormaps["viridis"].copy()
                cmap.set_bad(cmap(0.0))
                im = ax.imshow(plot, aspect="auto", origin="lower", cmap=cmap,
                               extent=[0, tp[-1], dist[0], dist[-1]], norm=norm)
                cb = self.fig.colorbar(im, ax=ax)
                cb.set_label("Glutamate concentration (uM)")
                ax.set_title("Glutamate Concentration Distribution")
                ax.set_xlabel("Time (ms)"); ax.set_ylabel("Distance (um)")
            elif view in ("Open state", "Desensitization") and r is not None:
                from matplotlib.colors import Normalize
                ax = self.fig.add_subplot(111)
                if view == "Open state":
                    data, title = r.open_state, "Open state"
                else:
                    data, title = r.desensitization, "Desensitization"
                auto_vmax = float(np.nanmax(data)) if np.isfinite(data).any() else 1.0
                auto_vmin = float(np.nanmin(data)) if np.isfinite(data).any() else 0.0
                vmin, vmax = self._color_limits(auto_vmin, auto_vmax)
                ext = [0, r.time_ms[-1], r.distance_um[0], r.distance_um[-1]]
                im = ax.imshow(data, aspect="auto", origin="lower", cmap="viridis",
                               extent=ext, norm=Normalize(vmin=vmin, vmax=vmax))
                self.fig.colorbar(im, ax=ax)
                ax.set_title(title)
                ax.set_xlabel("Time (ms)"); ax.set_ylabel("Distance (um)")
            else:
                ax = self.fig.add_subplot(111)
                ax.text(0.5, 0.5, "Run diffusion, then a receptor model.",
                        ha="center", va="center")
                ax.axis("off")
        except Exception:
            self._say("Plot error:\n" + traceback.format_exc())
        self.canvas.draw()

    # ---------------- receptor results ----------------
    def _register_receptor(self, key, res):
        self.receptor_results[key] = res
        self.recv_box["values"] = list(self.receptor_results.keys())
        self.recv_var.set(key)
        self.last_result = res

    def _show_receptor(self):
        key = self.recv_var.get()
        if key in self.receptor_results:
            self.last_result = self.receptor_results[key]
            if self.view_var.get() not in ("Open state", "Desensitization"):
                self.view_var.set("Open state")
            self._scale_manual = False
            self._redraw()

    def _save_receptor(self, res, model):
        if not self.results_folder:
            return False
        try:
            project.save_receptor_result(res, self.results_folder, model)
            return True
        except Exception:
            self._say("Save receptor failed:\n" + traceback.format_exc())
            return False

    def _save(self):
        if self.diff_result is None:
            messagebox.showinfo("Nothing to save", "Run a diffusion simulation first.")
            return
        folder = filedialog.askdirectory(title="Choose output folder")
        if not folder:
            return
        paths = io_files.save_receptor_input(self.diff_result, folder)
        io_files.save_snapshots(self.diff_result, folder)
        self._say("Saved DistanceFree/Bound, metadata, and snapshots to:\n" + folder)

    def _choose_folder(self):
        init = os.path.dirname(self.results_folder) if self.results_folder else ""
        folder = filedialog.askdirectory(title="Choose or create results folder",
                                         mustexist=False, initialdir=init)
        if not folder:
            return
        os.makedirs(folder, exist_ok=True)
        self.results_folder = folder
        self.folder_var.set(folder)
        self._say("Results folder: " + folder)

    def _new_folder(self):
        init = os.path.dirname(self.results_folder) if self.results_folder else ""
        parent = filedialog.askdirectory(title="Choose the parent folder for the new run",
                                         initialdir=init)
        if not parent:
            return
        name = simpledialog.askstring("New results folder",
                                      "Name for the new results folder:", parent=self.root)
        if not name:
            return
        folder = os.path.join(parent, name.strip())
        try:
            os.makedirs(folder, exist_ok=True)
        except Exception as e:
            messagebox.showerror("Cannot create folder", str(e)); return
        self.results_folder = folder
        self.folder_var.set(folder)
        self._say("New results folder: " + folder)

    def _save_params(self):
        try:
            self._apply_to_cfg()
        except Exception as e:
            messagebox.showerror("Invalid parameters", str(e)); return
        folder = self.results_folder or filedialog.askdirectory(title="Folder for Parameters.json")
        if not folder:
            return
        os.makedirs(folder, exist_ok=True)
        path = project.save_config(self.cfg, folder)
        if not self.results_folder:
            self.results_folder = folder
            self.folder_var.set(folder)
        self._say("Saved parameters -> " + path)

    def _load_params(self):
        path = filedialog.askopenfilename(title="Load Parameters.json",
                                          filetypes=[("JSON", "*.json"), ("All files", "*.*")])
        if not path:
            return
        try:
            self.cfg = project.load_config(path)
        except Exception as e:
            messagebox.showerror("Load failed", str(e)); return
        self._load_from_cfg()
        self._say("Loaded parameters <- " + path)

    def _load_results(self):
        folder = self.results_folder or filedialog.askdirectory(title="Results folder to load")
        if not folder:
            return
        try:
            data = project.load_results(folder)
        except Exception as e:
            messagebox.showerror("Load failed", str(e)); return
        if data.get("average_free") is None:
            messagebox.showinfo("No results", "No DistanceFree files found in this folder.")
            return
        if data.get("cfg") is not None:
            self.cfg = data["cfg"]
            self._load_from_cfg()
        self.results_folder = folder
        self.folder_var.set(folder)
        self.trials_out = dict(average_free=data["average_free"],
                               average_bound=data["average_bound"],
                               metadata=data["metadata"], n_trials=data["n_trials"],
                               occupied_fraction=data["occupied_fraction"],
                               mean_diffusion=data["mean_diffusion"], results=[])
        from types import SimpleNamespace
        sp = data["spillover"]
        self.diff_result = SimpleNamespace(
            snapshots=data["snapshots"],
            time_ms=(sp["time_ms"] if sp else np.array([0.0])),
            fraction_spillover=(sp["fraction"] if sp else np.array([0.0])),
            fraction_free_spillover=(sp["fraction_free"] if sp else np.array([0.0])),
            metadata=data["metadata"])
        self.receptor_results = {}
        for model, rd in data.get("receptors", {}).items():
            self.receptor_results[model] = SimpleNamespace(
                model_type=model, open_state=rd["open_state"],
                desensitization=rd["desensitization"], time_ms=rd["time_ms"],
                distance_um=rd["distance_um"], glu_time_ms=rd["glu_time_ms"])
        if self.receptor_results:
            self.recv_box["values"] = list(self.receptor_results.keys())
            first = next(iter(self.receptor_results))
            self.recv_var.set(first)
            self.last_result = self.receptor_results[first]
        else:
            self.recv_box["values"] = []
            self.recv_var.set("")
            self.last_result = None
        self._busy(False)
        rec = ", ".join(self.receptor_results.keys()) if self.receptor_results else "none saved (Run AMPA/NMDA)"
        self._say("Loaded results from %s (%d-trial mean, occupancy %.1f%%). "
                  "Receptors available: %s." % (folder, data["n_trials"],
                  100 * data["occupied_fraction"], rec))
        self.view_var.set("Concentration")
        self._redraw()

    def _say(self, msg):
        self.log.insert("end", msg + "\n")
        self.log.see("end")


def main():
    root = tk.Tk()
    App(root)
    root.mainloop()


if __name__ == "__main__":
    main()
