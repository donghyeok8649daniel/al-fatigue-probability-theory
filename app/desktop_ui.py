"""Lightweight Pre/Solve/Post desktop UI for the probability-PDE solver."""
from __future__ import annotations

import os
import queue
import sys
import threading
import time
import tkinter as tk
from pathlib import Path
from tkinter import messagebox, ttk

import matplotlib

matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib.ticker import ScalarFormatter
import numpy as np

from .solver_adapter import UIAnalysisConfig, physical_load_conversion, run_ui_analysis


APP_BG = "#eef1f4"
PANEL_BG = "#ffffff"
HEADER_BG = "#253746"
SIDEBAR_BG = "#e4e9ed"
ACCENT = "#1677c8"
TEXT = "#202830"
MUTED = "#65717c"
WINDOW_TITLE = "Al Fatigue — Probability PDE"
_INSTANCE_MUTEX = None


FIELD_LABELS = {
    "stress": "Applied stress",
    "strain_components": "Strain components",
    "strain": "Total axial strain",
    "normal_strain": "Normal opening",
    "intrawell_strain": "Intrawell registry",
    "plastic_strain": "Well-index plastic",
    "survival": "Survival",
    "initiation_probability": "Initiation probability",
    "first_passage_flux": "First-passage flux",
    "mass_diagnostics": "Mass diagnostics",
}


def acquire_single_instance() -> bool:
    """Prevent accidental duplicate desktop windows on Windows."""

    global _INSTANCE_MUTEX
    if os.name != "nt":
        return True
    import ctypes

    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    user32 = ctypes.WinDLL("user32", use_last_error=True)
    kernel32.CreateMutexW.argtypes = (ctypes.c_void_p, ctypes.c_bool, ctypes.c_wchar_p)
    kernel32.CreateMutexW.restype = ctypes.c_void_p
    kernel32.CloseHandle.argtypes = (ctypes.c_void_p,)
    user32.FindWindowW.argtypes = (ctypes.c_wchar_p, ctypes.c_wchar_p)
    user32.FindWindowW.restype = ctypes.c_void_p
    handle = kernel32.CreateMutexW(None, False, "Local\\AlFatigueProbabilityPDE")
    if not handle:
        return True
    if ctypes.get_last_error() == 183:
        window = user32.FindWindowW(None, WINDOW_TITLE)
        if window:
            user32.ShowWindow(window, 9)
            user32.SetForegroundWindow(window)
        kernel32.CloseHandle(handle)
        return False
    _INSTANCE_MUTEX = handle
    return True


class DesktopApp:
    """Responsive engineering workspace with explicit Pre/Solve/Post stages."""

    PARAMS = (
        ("young_gpa", "Young's modulus", "69", "GPa"),
        ("stress_mean_mpa", "Mean axial stress", "50", "MPa"),
        ("stress_amplitude_mpa", "Stress amplitude", "100", "MPa"),
        ("model_frequency", "Model frequency", "25", "cycles / model time"),
        ("cycles", "Load cycles", "1", "cycles"),
        ("steps_per_cycle", "Output resolution", "40", "steps/cycle"),
        ("tensile_direction", "Axial direction", "1 0 0", "[x y z]"),
    )

    def __init__(self) -> None:
        self.root = tk.Tk()
        self.root.title(WINDOW_TITLE)
        self.root.configure(bg=APP_BG)
        self.root.minsize(940, 620)
        self._center(1180, 760)

        self.entries: dict[str, ttk.Entry] = {}
        self.spatial_backend = tk.StringVar(value="FVM")
        self.analysis_quality = tk.StringVar(value="Preview (21 x 31, explicit)")
        self.field = tk.StringVar(value="strain_components")
        self.status = tk.StringVar(value="Ready · probability PDE always active")
        self.busy = False
        self.stop_event = threading.Event()
        self.result: dict[str, object] | None = None
        self.live_records: list[dict[str, float]] = []
        self._queue: queue.Queue[tuple[str, object]] = queue.Queue(maxsize=512)
        self._poll_job = None
        self._last_draw = 0.0
        self._last_field: str | None = None
        self._view_limits: dict[str, tuple[tuple[float, float], tuple[float, float]]] = {}
        self._pan_origin = None

        self._styles()
        self._header()
        self._workspace()
        self._statusbar()
        self._connect_plot_events()
        self.root.protocol("WM_DELETE_WINDOW", self._close)
        self._plot()

    def _center(self, width: int, height: int) -> None:
        sw, sh = self.root.winfo_screenwidth(), self.root.winfo_screenheight()
        width = min(width, max(900, sw - 60))
        height = min(height, max(600, sh - 90))
        x, y = max(0, (sw - width) // 2), max(0, (sh - height) // 2)
        self.root.geometry(f"{width}x{height}+{x}+{y}")

    def _styles(self) -> None:
        style = ttk.Style(self.root)
        if "clam" in style.theme_names():
            style.theme_use("clam")
        style.configure("App.TFrame", background=APP_BG)
        style.configure("Panel.TFrame", background=PANEL_BG)
        style.configure("Side.TFrame", background=SIDEBAR_BG)
        style.configure(
            "Header.TLabel",
            background=HEADER_BG,
            foreground="white",
            font=("Segoe UI", 13, "bold"),
        )
        style.configure(
            "SubHeader.TLabel",
            background=HEADER_BG,
            foreground="#cbd5dd",
            font=("Segoe UI", 9),
        )
        style.configure(
            "Section.TLabel", background=PANEL_BG, foreground=TEXT,
            font=("Segoe UI", 10, "bold")
        )
        style.configure("Property.TLabel", background=PANEL_BG, foreground=TEXT)
        style.configure("Unit.TLabel", background=PANEL_BG, foreground=MUTED)
        style.configure(
            "Accent.TButton", background=ACCENT, foreground="white",
            font=("Segoe UI", 9, "bold"), padding=(14, 7)
        )
        style.configure(
            "Field.Toolbutton", background="#e7edf2", foreground=TEXT,
            padding=(8, 5), relief="flat"
        )
        style.map(
            "Field.Toolbutton",
            background=[("selected", ACCENT), ("active", "#cfe4f5")],
            foreground=[("selected", "white")],
        )
        style.configure("TNotebook.Tab", font=("Segoe UI", 10, "bold"), padding=(18, 7))

    def _header(self) -> None:
        header = tk.Frame(self.root, bg=HEADER_BG, height=56)
        header.pack(fill="x")
        header.pack_propagate(False)
        ttk.Label(header, text="AL FATIGUE", style="Header.TLabel").pack(
            side="left", padx=(18, 10), pady=15
        )
        ttk.Label(
            header,
            text="Direct probability PDE · axial loading · Theory always active",
            style="SubHeader.TLabel",
        ).pack(side="left", pady=18)

    def _workspace(self) -> None:
        body = ttk.Frame(self.root, style="App.TFrame")
        body.pack(fill="both", expand=True)
        sidebar = ttk.Frame(body, style="Side.TFrame", width=205)
        sidebar.pack(side="left", fill="y")
        sidebar.pack_propagate(False)
        tk.Label(
            sidebar, text="STUDY TREE", bg=SIDEBAR_BG, fg=MUTED,
            font=("Segoe UI", 8, "bold")
        ).pack(anchor="w", padx=14, pady=(14, 5))
        tree = ttk.Treeview(sidebar, show="tree", height=20)
        project = tree.insert("", "end", text="  Axial fatigue study", open=True)
        pre = tree.insert(project, "end", text="  Pre", open=True)
        tree.insert(pre, "end", text="  Material / loading")
        tree.insert(pre, "end", text="  Axial direction")
        solve = tree.insert(project, "end", text="  Solve", open=True)
        tree.insert(solve, "end", text="  Probability PDE (always on)")
        tree.insert(solve, "end", text="  Spatial display: FVM / FEM")
        post = tree.insert(project, "end", text="  Post", open=True)
        for label in FIELD_LABELS.values():
            tree.insert(post, "end", text=f"  {label}")
        tree.pack(fill="both", expand=True, padx=8, pady=(0, 8))

        main = ttk.Frame(body, style="App.TFrame")
        main.pack(side="left", fill="both", expand=True)
        self.notebook = ttk.Notebook(main)
        self.notebook.pack(fill="both", expand=True, padx=8, pady=8)
        self.pre_tab = ttk.Frame(self.notebook, style="Panel.TFrame")
        self.solve_tab = ttk.Frame(self.notebook, style="Panel.TFrame")
        self.post_tab = ttk.Frame(self.notebook, style="Panel.TFrame")
        self.notebook.add(self.pre_tab, text="1  PRE")
        self.notebook.add(self.solve_tab, text="2  SOLVE")
        self.notebook.add(self.post_tab, text="3  POST")
        self._pre_tab()
        self._solve_tab()
        self._post_tab()

    def _pre_tab(self) -> None:
        ttk.Label(self.pre_tab, text="Physical axial-load inputs", style="Section.TLabel").grid(
            row=0, column=0, columnspan=3, sticky="w", padx=20, pady=(18, 10)
        )
        for row, (key, label, default, unit) in enumerate(self.PARAMS, start=1):
            ttk.Label(self.pre_tab, text=label, style="Property.TLabel").grid(
                row=row, column=0, sticky="w", padx=(20, 8), pady=7
            )
            entry = ttk.Entry(self.pre_tab, width=20)
            entry.insert(0, default)
            entry.grid(row=row, column=1, sticky="ew", padx=4, pady=7)
            ttk.Label(self.pre_tab, text=unit, style="Unit.TLabel").grid(
                row=row, column=2, sticky="w", padx=(6, 20), pady=7
            )
            self.entries[key] = entry
        self.pre_tab.columnconfigure(1, weight=1)
        note = ttk.LabelFrame(self.pre_tab, text="Scope", padding=12)
        note.grid(row=len(self.PARAMS) + 2, column=0, columnspan=3, sticky="ew", padx=20, pady=14)
        ttk.Label(
            note,
            text=(
                "Physical stress is reduced by E and then mapped with the verified relaxed axial κ.\n"
                "The current mechanics applies stress only along the entered axis. Transverse stress is zero."
            ),
            justify="left",
        ).pack(anchor="w")

    def _solve_tab(self) -> None:
        left = ttk.Frame(self.solve_tab, style="Panel.TFrame")
        left.pack(side="left", fill="y", padx=20, pady=18)
        ttk.Label(left, text="Probability PDE", style="Section.TLabel").pack(anchor="w")
        ttk.Label(
            left, text="Always active", foreground=ACCENT,
            background=PANEL_BG, font=("Segoe UI", 9, "bold")
        ).pack(anchor="w", pady=(3, 15))
        ttk.Label(left, text="Spatial display backend", style="Section.TLabel").pack(anchor="w")
        ttk.Combobox(
            left, textvariable=self.spatial_backend, values=("FVM", "FEM"),
            state="readonly", width=22
        ).pack(anchor="w", pady=(6, 15))
        ttk.Label(left, text="Probability grid quality", style="Section.TLabel").pack(
            anchor="w"
        )
        ttk.Combobox(
            left,
            textvariable=self.analysis_quality,
            values=(
                "Preview (21 x 31, explicit)",
                "Resolved (81 x 91, implicit)",
            ),
            state="readonly",
            width=29,
        ).pack(anchor="w", pady=(6, 15))
        ttk.Label(
            left,
            text=(
                "P(a,s,t) is evolved directly.\n"
                "Initial state: Gibbs at σ(t=0).\n"
                "No Monte Carlo counting.\n"
                "Model time is not calibrated seconds."
            ),
            style="Property.TLabel",
            justify="left",
        ).pack(anchor="w", pady=(0, 16))
        self.run_button = ttk.Button(
            left, text="RUN ANALYSIS", style="Accent.TButton", command=self._start_solve
        )
        self.run_button.pack(fill="x")
        self.progress = ttk.Progressbar(left, mode="indeterminate", length=235)
        self.progress.pack(fill="x", pady=12)
        self.summary = tk.Text(
            self.solve_tab, wrap="word", relief="flat", bg="#f7f9fa",
            fg=TEXT, font=("Consolas", 9), padx=14, pady=12
        )
        self.summary.pack(side="left", fill="both", expand=True, padx=(0, 20), pady=18)
        self._set_summary("Ready. Probability and all strain components will come from solver_v1.\n")

    def _post_tab(self) -> None:
        toolbar = ttk.Frame(self.post_tab, style="Panel.TFrame")
        toolbar.pack(fill="x", padx=12, pady=(10, 3))
        for index, (key, label) in enumerate(FIELD_LABELS.items()):
            row, column = divmod(index, 5)
            ttk.Radiobutton(
                toolbar, text=label, value=key, variable=self.field,
                style="Field.Toolbutton", command=self._plot
            ).grid(row=row, column=column, sticky="ew", padx=2, pady=2)
        for column in range(5):
            toolbar.columnconfigure(column, weight=1, uniform="field")
        self.figure = Figure(figsize=(8, 5), dpi=96, facecolor=PANEL_BG)
        self.ax = self.figure.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.figure, master=self.post_tab)
        self.canvas.get_tk_widget().pack(fill="both", expand=True, padx=10, pady=(4, 0))
        nav = ttk.Frame(self.post_tab, style="Panel.TFrame")
        nav.pack(fill="x", padx=10, pady=(0, 7))
        NavigationToolbar2Tk(self.canvas, nav, pack_toolbar=False).pack(side="left")
        ttk.Label(
            nav,
            text="Wheel: zoom · Shift: X only · Ctrl: Y only · drag: pan · double-click/Home: reset",
            foreground=MUTED,
        ).pack(side="right", padx=8)

    def _statusbar(self) -> None:
        tk.Label(
            self.root, textvariable=self.status, anchor="w", bg="#dbe2e7",
            fg=TEXT, font=("Segoe UI", 9), padx=10, pady=4
        ).pack(fill="x", side="bottom")

    def _config(self) -> UIAnalysisConfig:
        direction = self.entries["tensile_direction"].get().replace(",", " ").split()
        if len(direction) != 3 or np.linalg.norm([float(v) for v in direction]) == 0.0:
            raise ValueError("Axial direction requires three values and must be nonzero")
        resolved = self.analysis_quality.get().startswith("Resolved")
        return UIAnalysisConfig(
            young_gpa=float(self.entries["young_gpa"].get()),
            stress_mean_mpa=float(self.entries["stress_mean_mpa"].get()),
            stress_amplitude_mpa=float(self.entries["stress_amplitude_mpa"].get()),
            model_frequency=float(self.entries["model_frequency"].get()),
            cycles=float(self.entries["cycles"].get()),
            steps_per_cycle=int(self.entries["steps_per_cycle"].get()),
            grid_n_a=81 if resolved else 21,
            grid_n_s=91 if resolved else 31,
            integration_method="implicit" if resolved else "explicit",
            analysis_quality="resolved" if resolved else "preview",
        )

    def _start_solve(self) -> None:
        if self.busy:
            self.stop_event.set()
            self.status.set("Stopping after the current PDE step…")
            self.run_button.configure(text="STOPPING…", state="disabled")
            return
        try:
            config = self._config()
            config.validate()
            conversion = physical_load_conversion(config)
        except Exception as exc:
            messagebox.showerror("Invalid input", str(exc), parent=self.root)
            return
        self.busy = True
        self.stop_event.clear()
        self.result = None
        self.live_records.clear()
        self._view_limits.clear()
        self.run_button.configure(text="STOP ANALYSIS", state="normal")
        self.progress.start(12)
        self.status.set("Solving direct probability PDE in background…")
        self._set_summary(
            f"Spatial display: {self.spatial_backend.get()}\n"
            f"Probability solve quality: {config.analysis_quality} "
            f"({config.grid_n_a} x {config.grid_n_s}, {config.integration_method})\n"
            f"Stress: {config.stress_min_mpa:.6g} to {config.stress_max_mpa:.6g} MPa\n"
            f"Reduced stress: {conversion['reduced_stress_min']:.6g} to "
            f"{conversion['reduced_stress_max']:.6g}\n"
            f"Model period: {conversion['model_period']:.8g}\n"
            f"Local pristine De_fast / De_slow: {conversion['de_fast']:.6g} / "
            f"{conversion['de_slow']:.6g}\n"
            f"Relaxed axial κ: {conversion['relaxed_axial_kappa']:.12g}\n"
            f"Initial Gibbs preload force: {conversion['preload_force']:.12g}\n"
            "Computing…\n"
        )
        threading.Thread(target=self._solve_worker, args=(config,), daemon=True).start()
        if self._poll_job is None:
            self._poll_job = self.root.after(30, self._drain_queue)

    def _solve_worker(self, config: UIAnalysisConfig) -> None:
        try:
            def emit(record: dict[str, float]) -> None:
                while not self.stop_event.is_set():
                    try:
                        self._queue.put(("record", record), timeout=0.05)
                        return
                    except queue.Full:
                        continue

            result = run_ui_analysis(
                config, record_callback=emit, stop_requested=self.stop_event.is_set
            )
            self._queue.put(("done", result))
        except Exception as exc:
            self._queue.put(("error", str(exc)))

    def _drain_queue(self) -> None:
        self._poll_job = None
        for _ in range(128):
            try:
                kind, payload = self._queue.get_nowait()
            except queue.Empty:
                break
            if kind == "record":
                self.live_records.append(payload)
                now = time.monotonic()
                if now - self._last_draw > 0.20:
                    self._last_draw = now
                    self._plot()
                self.status.set(
                    f"PDE model time {payload['model_time']:.5g} · "
                    f"S={payload['survival']:.12g} · ε={payload['strain']:.6g}"
                )
            elif kind == "done":
                self._solve_done(payload)
                return
            elif kind == "error":
                self._solve_failed(str(payload))
                return
        if self.busy or not self._queue.empty():
            self._poll_job = self.root.after(30, self._drain_queue)

    def _solve_done(self, result: dict[str, object]) -> None:
        self.result = result
        self.busy = False
        self.progress.stop()
        self.run_button.configure(text="RUN ANALYSIS", state="normal")
        self.notebook.select(self.post_tab)
        self.status.set(
            f"Complete · {len(result['model_time'])} records · "
            f"S={float(np.asarray(result['survival'])[-1]):.12g}"
        )
        self._set_summary(
            f"Probability source: {result['probability_source']}\n"
            f"Solve quality: {result['analysis_quality']} "
            f"{result['grid_shape']}, {result['integration_method']}\n"
            f"Initial condition: {result['initial_condition']}\n"
            f"Initial stress: {result['initial_stress_mpa']:.8g} MPa\n"
            f"Model frequency: {result['model_frequency']:.8g} cycles/model time\n"
            f"Model period: {result['model_period']:.8g}\n"
            f"Local pristine tau_fast / tau_slow: {result['tau_fast']:.8g} / "
            f"{result['tau_slow']:.8g}\n"
            f"Local pristine omega*tau_fast / omega*tau_slow: "
            f"{result['de_fast']:.8g} / {result['de_slow']:.8g}\n"
            f"Local linear |G| / phase: {result['linear_transfer_magnitude']:.8g} / "
            f"{result['linear_transfer_phase_degrees']:.6g} deg\n"
            f"Relaxed axial κ: {result['relaxed_axial_kappa']:.12g}\n"
            f"Final survival: {float(np.asarray(result['survival'])[-1]):.12g}\n"
            f"Final initiation: {float(np.asarray(result['initiation_probability'])[-1]):.12g}\n"
            f"Max |mass residual|: {float(np.max(np.abs(result['mass_balance_residual']))):.3e}\n"
            f"Time status: {result['solver_time_status']}\n"
        )
        self._plot(force_auto=True)

    def _solve_failed(self, detail: str) -> None:
        self.busy = False
        self.progress.stop()
        self.run_button.configure(text="RUN ANALYSIS", state="normal")
        self.status.set("Solve failed")
        messagebox.showerror("PDE solver error", detail, parent=self.root)

    def _set_summary(self, text: str) -> None:
        self.summary.configure(state="normal")
        self.summary.delete("1.0", "end")
        self.summary.insert("end", text)
        self.summary.configure(state="disabled")

    def _plot_data(self) -> dict[str, np.ndarray] | None:
        if self.result is not None:
            return {key: np.asarray(value) for key, value in self.result.items() if isinstance(value, np.ndarray)}
        if not self.live_records:
            return None
        keys = self.live_records[0].keys()
        return {key: np.asarray([row[key] for row in self.live_records], dtype=float) for key in keys}

    def _plot(self, force_auto: bool = False) -> None:
        field = self.field.get()
        if self._last_field is not None and self._last_field != field and self.ax.has_data():
            self._view_limits[self._last_field] = (self.ax.get_xlim(), self.ax.get_ylim())
        self._last_field = field
        data = self._plot_data()
        self.ax.clear()
        if data is None or "model_time" not in data:
            self.ax.text(
                0.5, 0.5, "Run an analysis to view PDE results", ha="center", va="center",
                transform=self.ax.transAxes, color=MUTED
            )
            self.ax.set_axis_off()
            self.canvas.draw_idle()
            return
        self.ax.set_axis_on()
        x = np.asarray(data["model_time"], dtype=float)
        if field == "strain_components":
            for key, label, color in (
                ("normal_strain", "Normal opening", "#4e79a7"),
                ("intrawell_strain", "Intrawell registry", "#59a14f"),
                ("plastic_strain", "Well-index plastic", "#e15759"),
                ("strain", "Total axial strain", "#111827"),
            ):
                self.ax.plot(x, data[key], label=label, color=color, linewidth=1.8)
            ylabel = "Axial strain [-]"
            self.ax.legend(loc="best", fontsize=8, frameon=False)
        elif field == "mass_diagnostics":
            for key, label in (
                ("intact_probability_mass", "Intact mass"),
                ("cumulative_absorbed_mass", "Absorbed mass"),
                ("mass_balance_residual", "Mass-balance residual"),
                ("negative_mass_correction", "Max negative correction"),
            ):
                self.ax.plot(x, data[key], label=label, linewidth=1.5)
            ylabel = "Probability mass [-]"
            self.ax.legend(loc="best", fontsize=8, frameon=False)
        else:
            key = "applied_stress_mpa" if field == "stress" else field
            self.ax.plot(x, data[key], color=ACCENT, linewidth=1.8)
            ylabel = {
                "stress": "Axial stress [MPa]",
                "strain": "Axial strain [-]",
                "normal_strain": "Normal strain [-]",
                "intrawell_strain": "Intrawell strain [-]",
                "plastic_strain": "Plastic strain [-]",
                "survival": "Survival probability [-]",
                "initiation_probability": "Initiation probability [-]",
                "first_passage_flux": "First-passage flux [1/model time]",
            }[field]
        self.ax.set_xlabel("Dimensionless solver time")
        self.ax.set_ylabel(ylabel)
        self.ax.set_title(FIELD_LABELS[field], loc="left", fontsize=12, fontweight="bold")
        self.ax.grid(True, color="#d5dde3", linewidth=0.7, alpha=0.8)
        formatter = ScalarFormatter(useMathText=True, useOffset=True)
        formatter.set_powerlimits((-3, 3))
        self.ax.yaxis.set_major_formatter(formatter)
        if force_auto:
            self._view_limits.pop(field, None)
        if field in self._view_limits:
            self.ax.set_xlim(*self._view_limits[field][0])
            self.ax.set_ylim(*self._view_limits[field][1])
        else:
            self.ax.relim()
            self.ax.autoscale_view()
        self.figure.tight_layout(pad=1.1)
        self.canvas.draw_idle()

    def _connect_plot_events(self) -> None:
        self.canvas.mpl_connect("scroll_event", self._on_scroll)
        self.canvas.mpl_connect("button_press_event", self._on_press)
        self.canvas.mpl_connect("button_release_event", self._on_release)
        self.canvas.mpl_connect("motion_notify_event", self._on_motion)
        self.canvas.mpl_connect("key_press_event", self._on_key)

    def _remember_view(self) -> None:
        self._view_limits[self.field.get()] = (self.ax.get_xlim(), self.ax.get_ylim())

    def _on_scroll(self, event) -> None:
        if event.inaxes is not self.ax or event.xdata is None or event.ydata is None:
            return
        factor = 0.80 if event.button == "up" else 1.25
        key = (event.key or "").lower()
        zoom_x = "control" not in key and "ctrl" not in key
        zoom_y = "shift" not in key
        if zoom_x:
            left, right = self.ax.get_xlim()
            self.ax.set_xlim(
                event.xdata - (event.xdata - left) * factor,
                event.xdata + (right - event.xdata) * factor,
            )
        if zoom_y:
            low, high = self.ax.get_ylim()
            self.ax.set_ylim(
                event.ydata - (event.ydata - low) * factor,
                event.ydata + (high - event.ydata) * factor,
            )
        self._remember_view()
        self.canvas.draw_idle()

    def _on_press(self, event) -> None:
        if event.dblclick and event.inaxes is self.ax:
            self._reset_view()
            return
        if event.button == 1 and event.inaxes is self.ax and event.xdata is not None:
            self._pan_origin = (
                event.xdata, event.ydata, self.ax.get_xlim(), self.ax.get_ylim()
            )

    def _on_motion(self, event) -> None:
        if self._pan_origin is None or event.inaxes is not self.ax:
            return
        if event.xdata is None or event.ydata is None:
            return
        x0, y0, xlim, ylim = self._pan_origin
        dx, dy = event.xdata - x0, event.ydata - y0
        self.ax.set_xlim(xlim[0] - dx, xlim[1] - dx)
        self.ax.set_ylim(ylim[0] - dy, ylim[1] - dy)
        self._remember_view()
        self.canvas.draw_idle()

    def _on_release(self, _event) -> None:
        self._pan_origin = None

    def _on_key(self, event) -> None:
        if (event.key or "").lower() == "home":
            self._reset_view()

    def _reset_view(self) -> None:
        self._view_limits.pop(self.field.get(), None)
        self._plot(force_auto=True)

    def _close(self) -> None:
        self.stop_event.set()
        self.root.destroy()

    def run(self) -> None:
        self.root.mainloop()


def startup_smoke() -> dict[str, float]:
    """Exercise imports and canonical conversion without opening a window."""

    conversion = physical_load_conversion(UIAnalysisConfig())
    return {
        "a0": conversion["a0"],
        "relaxed_axial_kappa": conversion["relaxed_axial_kappa"],
        "model_frequency": conversion["model_frequency"],
    }


def main() -> None:
    if "--smoke" in sys.argv:
        values = startup_smoke()
        print(
            f"desktop-ui smoke OK: a0={values['a0']:.16g}, "
            f"kappa={values['relaxed_axial_kappa']:.16g}"
        )
        return
    if not acquire_single_instance():
        return
    DesktopApp().run()


if __name__ == "__main__":
    main()
