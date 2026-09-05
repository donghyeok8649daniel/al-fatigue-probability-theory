"""Lightweight Pre/Solve/Post UI for uniaxial theory and 3D geometry."""
from __future__ import annotations

import os
import queue
import sys
import threading
import time
import tkinter as tk
from collections import deque
from pathlib import Path
from tkinter import filedialog, messagebox, ttk

import matplotlib

matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import numpy as np


APP_BG = "#f3f5f7"
PANEL_BG = "#ffffff"
HEADER_BG = "#263746"
SIDEBAR_BG = "#e9edf1"
ACCENT = "#1677c8"
TEXT = "#202830"
MUTED = "#66727d"
WINDOW_TITLE = "Al Fatigue — Theory Core v1"
_INSTANCE_MUTEX = None


def acquire_single_instance() -> bool:
    """Keep one Windows UI instance and focus it when launched again."""
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
    user32.ShowWindow.argtypes = (ctypes.c_void_p, ctypes.c_int)
    user32.SetForegroundWindow.argtypes = (ctypes.c_void_p,)
    handle = kernel32.CreateMutexW(None, False, "Local\\AlFatigueDesktopApp")
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
        ("length_mm", "Specimen length", "50", "mm"),
        ("diameter_mm", "Cylinder diameter", "6", "mm"),
        ("young_gpa", "Young's modulus", "69", "GPa"),
        ("poisson_ratio", "Poisson ratio", "0.33", "-"),
        ("loading_direction", "Crystal direction", "1 0 0", "[h k l]"),
        ("tensile_direction", "Tensile stress direction", "1 0 0", "[x y z]"),
        ("elements", "Control volumes", "40", "cells"),
        ("stress_mean_mpa", "Mean normal stress", "50", "MPa"),
        ("stress_amplitude_mpa", "Normal stress amplitude", "100", "MPa"),
        ("frequency_hz", "Loading frequency", "20", "Hz"),
        ("cycles", "Spatial preview", "2", "cycles"),
        ("steps_per_cycle", "Resolution", "80", "steps/cycle"),
        ("deformation_scale", "Display deformation", "1", "x"),
    )

    FIELD_LABELS = {
        "stress": "Normal stress",
        "strain": "Axial strain",
        "diameter": "Diameter change",
        "initiation": "Specimen first passage",
        "survival": "Specimen survival",
        "hazard": "First-passage flux",
        "life": "First-passage mass",
    }

    def __init__(self, project_path: Path | None = None) -> None:
        self.root = tk.Tk()
        self.root.title(WINDOW_TITLE)
        self.root.configure(bg=APP_BG)
        self.root.minsize(980, 640)
        self._center(1180, 760)

        self.output_dir = self._default_output_dir()
        self.geometry_path = self._default_cylinder_path()
        self.spatial_backend = tk.StringVar(value="FVM")
        self.field = tk.StringVar(value="stress")
        self.entries: dict[str, ttk.Entry] = {}
        self.nodes: np.ndarray | None = None
        self.elements: np.ndarray | None = None
        self.current_config = None
        self.busy = False
        self.stop_event = threading.Event()
        self.live_cycles = 0.0
        self.live_records = deque(maxlen=5000)
        self.live_events: list[dict] = []
        self._live_queue: queue.Queue[tuple[str, object]] = queue.Queue(maxsize=512)
        self._live_poll_job = None
        self._last_live_draw = 0.0
        self.playing = False
        self.play_position = tk.DoubleVar(value=0.0)
        self._play_job = None
        self._cursor_line = None
        self._cursor_domain = None
        self._cursor_x = None
        self._cursor_y = None
        self._cursor_value_name = None
        self._last_plot_field = None

        self._styles()
        self._header()
        self._workspace()
        self._statusbar()
        self.root.protocol("WM_DELETE_WINDOW", self._close)
        self._show_empty_plot()
        if project_path is not None:
            self.root.after(0, lambda: self._load_project(Path(project_path)))

    @staticmethod
    def _default_output_dir() -> Path:
        if getattr(sys, "frozen", False):
            local_data = os.environ.get("LOCALAPPDATA")
            base = Path(local_data) if local_data else Path.home() / "AppData" / "Local"
            return base / "AlFatigue" / "results" / "desktop_session"
        return Path("results/data/desktop_session")

    @staticmethod
    def _default_cylinder_path() -> Path:
        root = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parents[1]))
        return root / "examples" / "meshes" / "default_tensile_cylinder.stl"

    def _center(self, width: int, height: int) -> None:
        sw = self.root.winfo_screenwidth()
        sh = self.root.winfo_screenheight()
        width = min(width, max(900, sw - 80))
        height = min(height, max(600, sh - 100))
        self.root.geometry(f"{width}x{height}+{max(0, (sw-width)//2)}+{max(0, (sh-height)//2)}")

    def _styles(self) -> None:
        style = ttk.Style(self.root)
        if "clam" in style.theme_names():
            style.theme_use("clam")
        style.configure("App.TFrame", background=APP_BG)
        style.configure("Panel.TFrame", background=PANEL_BG)
        style.configure("Side.TFrame", background=SIDEBAR_BG)
        style.configure("Header.TLabel", background=HEADER_BG, foreground="white", font=("Segoe UI", 13, "bold"))
        style.configure("SubHeader.TLabel", background=HEADER_BG, foreground="#cbd5dd", font=("Segoe UI", 9))
        style.configure("Section.TLabel", background=PANEL_BG, foreground=TEXT, font=("Segoe UI", 10, "bold"))
        style.configure("Property.TLabel", background=PANEL_BG, foreground=TEXT, font=("Segoe UI", 9))
        style.configure("Unit.TLabel", background=PANEL_BG, foreground=MUTED, font=("Segoe UI", 8))
        style.configure("Accent.TButton", background=ACCENT, foreground="white", font=("Segoe UI", 9, "bold"), padding=(14, 7))
        style.map("Accent.TButton", background=[("active", "#0f65ac")])
        style.configure("Field.Toolbutton", background="#e7edf2", foreground=TEXT, padding=(10, 6), relief="flat")
        style.map("Field.Toolbutton", background=[("selected", ACCENT), ("active", "#cfe4f5")], foreground=[("selected", "white")])
        style.configure("TNotebook", background=APP_BG, borderwidth=0)
        style.configure("TNotebook.Tab", font=("Segoe UI", 10, "bold"), padding=(20, 8))

    def _header(self) -> None:
        header = tk.Frame(self.root, bg=HEADER_BG, height=58)
        header.pack(fill="x")
        header.pack_propagate(False)
        brand = ttk.Label(header, text="AL FATIGUE", style="Header.TLabel")
        brand.pack(side="left", padx=(18, 8), pady=(10, 0), anchor="n")
        subtitle = ttk.Label(header, text="Theory Core v1 always active · uniaxial stress / 3D geometry", style="SubHeader.TLabel")
        subtitle.pack(side="left", padx=4, pady=(18, 0), anchor="n")
        ttk.Button(header, text="Open project", command=self._open_project).pack(side="right", padx=(4, 16), pady=12)
        ttk.Button(header, text="Save project", command=self._save_project).pack(side="right", padx=4, pady=12)

    def _workspace(self) -> None:
        body = ttk.Frame(self.root, style="App.TFrame")
        body.pack(fill="both", expand=True)

        sidebar = ttk.Frame(body, style="Side.TFrame", width=210)
        sidebar.pack(side="left", fill="y")
        sidebar.pack_propagate(False)
        tk.Label(sidebar, text="MODEL", bg=SIDEBAR_BG, fg=MUTED, font=("Segoe UI", 8, "bold")).pack(anchor="w", padx=14, pady=(14, 5))
        self.tree = ttk.Treeview(sidebar, show="tree", selectmode="browse", height=18)
        project = self.tree.insert("", "end", text="  Fatigue Study", open=True)
        pre = self.tree.insert(project, "end", text="  Pre-processing", open=True)
        self.tree.insert(pre, "end", text="  3D cylinder / FVM mesh")
        self.tree.insert(pre, "end", text="  Material / Loading")
        solve = self.tree.insert(project, "end", text="  Solvers", open=True)
        self.tree.insert(solve, "end", text="  Theory Core v1 (always on)")
        self.tree.insert(solve, "end", text="  Spatial backend: FVM / FEM")
        post = self.tree.insert(project, "end", text="  Results", open=True)
        for label in self.FIELD_LABELS.values():
            self.tree.insert(post, "end", text=f"  {label}")
        self.tree.pack(fill="both", expand=True, padx=8, pady=(0, 8))

        main = ttk.Frame(body, style="App.TFrame")
        main.pack(side="left", fill="both", expand=True)
        self.notebook = ttk.Notebook(main)
        self.notebook.pack(fill="both", expand=True, padx=8, pady=8)
        self.pre_tab = ttk.Frame(self.notebook, style="Panel.TFrame")
        self.solve_tab = ttk.Frame(self.notebook, style="Panel.TFrame")
        self.post_tab = ttk.Frame(self.notebook, style="Panel.TFrame")
        self.notebook.add(self.pre_tab, text="1  PRE / MESH")
        self.notebook.add(self.solve_tab, text="2  SOLVE")
        self.notebook.add(self.post_tab, text="3  POST")
        self._pre_tab()
        self._solve_tab()
        self._post_tab()

    def _pre_tab(self) -> None:
        ttk.Label(self.pre_tab, text="Model properties", style="Section.TLabel").grid(row=0, column=0, columnspan=6, sticky="w", padx=18, pady=(16, 10))
        split = (len(self.PARAMS) + 1) // 2
        for i, (key, label, default, unit) in enumerate(self.PARAMS):
            column_group = 0 if i < split else 3
            row = i % split + 1
            ttk.Label(self.pre_tab, text=label, style="Property.TLabel").grid(row=row, column=column_group, sticky="w", padx=(18, 8), pady=7)
            entry = ttk.Entry(self.pre_tab, width=14)
            entry.insert(0, default)
            entry.grid(row=row, column=column_group + 1, sticky="ew", padx=4, pady=7)
            ttk.Label(self.pre_tab, text=unit, style="Unit.TLabel").grid(row=row, column=column_group + 2, sticky="w", padx=(4, 22), pady=7)
            self.entries[key] = entry
        self.pre_tab.columnconfigure(1, weight=1)
        self.pre_tab.columnconfigure(4, weight=1)
        mesh = ttk.LabelFrame(self.pre_tab, text="Geometry and mesh", padding=10)
        mesh.grid(row=split + 2, column=0, columnspan=6, sticky="ew", padx=18, pady=10)
        ttk.Button(mesh, text="VIEW DEFAULT CYLINDER", command=self._view_default_geometry).pack(side="left")
        ttk.Button(mesh, text="OPEN MESH", command=self._open_geometry).pack(side="left", padx=(6, 0))
        self.mesh_label = ttk.Label(
            mesh,
            text="default_tensile_cylinder.stl · 3D triangular surface · axial-only FVM",
            foreground=MUTED,
        )
        self.mesh_label.pack(side="left", padx=14)

    def _solve_tab(self) -> None:
        left = ttk.Frame(self.solve_tab, style="Panel.TFrame")
        left.pack(side="left", fill="y", padx=20, pady=18)
        ttk.Label(left, text="Theory Core v1", style="Section.TLabel").pack(anchor="w", pady=(0, 4))
        ttk.Label(left, text="Always active", foreground=ACCENT, background=PANEL_BG, font=("Segoe UI", 9, "bold")).pack(anchor="w", pady=(0, 16))
        ttk.Label(left, text="Spatial discretization", style="Section.TLabel").pack(anchor="w", pady=(0, 8))
        combo = ttk.Combobox(left, textvariable=self.spatial_backend, values=("FVM", "FEM"), state="readonly", width=24)
        combo.pack(anchor="w", pady=(0, 14))
        ttk.Label(left, text="Theory uses deterministic finite-chain mechanics.\nThe baseline initial measure is mu0 = delta(lambda=1,c=0).\nNormal loading enters exactly as q(tau) = sigma_n(t)/E.\nFVM/FEM applies stress only along the entered axis.\nNo noise, mobility closure, or fitted life distribution is used.", style="Property.TLabel", justify="left").pack(anchor="w", pady=(0, 18))
        self.run_button = ttk.Button(left, text="RUN ANALYSIS", style="Accent.TButton", command=self._start_solve)
        self.run_button.pack(anchor="w", fill="x")
        self.progress = ttk.Progressbar(left, mode="indeterminate", length=250)
        self.progress.pack(anchor="w", fill="x", pady=12)
        self.solve_summary = tk.Text(self.solve_tab, wrap="word", relief="flat", bg="#f8fafb", fg=TEXT, font=("Consolas", 9), padx=14, pady=12)
        self.solve_summary.pack(side="left", fill="both", expand=True, padx=(0, 20), pady=18)
        self.solve_summary.insert("end", "Ready. Theory Core v1 will compute hysteresis, survival, hazard and first passage.\n")
        self.solve_summary.configure(state="disabled")

    def _post_tab(self) -> None:
        toolbar = ttk.Frame(self.post_tab, style="Panel.TFrame")
        toolbar.pack(fill="x", padx=14, pady=(10, 4))
        ttk.Label(toolbar, text="Result field", style="Section.TLabel").grid(
            row=0, column=0, columnspan=4, sticky="w", padx=2, pady=(0, 5)
        )
        for index, (key, label) in enumerate(self.FIELD_LABELS.items()):
            row, column = divmod(index, 4)
            ttk.Radiobutton(
                toolbar,
                text=label,
                value=key,
                variable=self.field,
                style="Field.Toolbutton",
                command=self._plot,
            ).grid(row=row + 1, column=column, sticky="ew", padx=2, pady=2)
        for column in range(4):
            toolbar.columnconfigure(column, weight=1, uniform="post_fields")
        timeline = ttk.Frame(self.post_tab, style="Panel.TFrame")
        timeline.pack(fill="x", padx=14, pady=(2, 2))
        self.play_button = ttk.Button(timeline, text="PLAY", command=self._toggle_playback)
        self.play_button.pack(side="left", padx=(2, 8))
        ttk.Button(timeline, text="RESET", command=self._reset_playback).pack(side="left", padx=(0, 8))
        self.timeline = ttk.Scale(
            timeline,
            from_=0.0,
            to=1.0,
            variable=self.play_position,
            command=self._on_timeline,
        )
        self.timeline.pack(side="left", fill="x", expand=True, padx=4)
        self.timeline_label = ttk.Label(
            timeline, text="Time: 0", style="Property.TLabel", width=38, anchor="e"
        )
        self.timeline_label.pack(side="right", padx=(8, 2))
        self.figure = Figure(figsize=(8, 5), dpi=96, facecolor=PANEL_BG)
        self.ax = self.figure.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.figure, master=self.post_tab)
        self.canvas.get_tk_widget().pack(fill="both", expand=True, padx=10, pady=(4, 0))
        nav = ttk.Frame(self.post_tab, style="Panel.TFrame")
        nav.pack(fill="x", padx=10, pady=(0, 8))
        NavigationToolbar2Tk(self.canvas, nav, pack_toolbar=False).pack(side="left")

    def _statusbar(self) -> None:
        self.status = tk.StringVar(value="Ready · Theory Core v1 + FVM")
        bar = tk.Label(self.root, textvariable=self.status, anchor="w", bg="#dfe5ea", fg=TEXT, font=("Segoe UI", 9), padx=10, pady=4)
        bar.pack(fill="x", side="bottom")

    def _config(self) -> TensionRunConfig:
        from simulations.fem_tension_app import TensionRunConfig

        crystal_direction = self.entries["loading_direction"].get().replace(",", " ").split()
        if len(crystal_direction) != 3:
            raise ValueError("Crystal direction requires three integers: h k l")
        h, k, l = (int(x) for x in crystal_direction)
        tensile_direction = self.entries["tensile_direction"].get().replace(",", " ").split()
        if len(tensile_direction) != 3:
            raise ValueError("Tensile stress direction requires three values: x y z")
        axis_x, axis_y, axis_z = (float(x) for x in tensile_direction)
        diameter = float(self.entries["diameter_mm"].get())
        return TensionRunConfig(
            length_mm=float(self.entries["length_mm"].get()),
            width_mm=diameter,
            thickness_mm=diameter,
            section_shape="circular",
            diameter_mm=diameter,
            young_gpa=float(self.entries["young_gpa"].get()),
            poisson_ratio=float(self.entries["poisson_ratio"].get()),
            loading_h=h, loading_k=k, loading_l=l, elements=int(self.entries["elements"].get()),
            tensile_axis_x=axis_x, tensile_axis_y=axis_y, tensile_axis_z=axis_z,
            stress_mean_mpa=float(self.entries["stress_mean_mpa"].get()),
            stress_amplitude_mpa=float(self.entries["stress_amplitude_mpa"].get()),
            theory_stress_scale_mpa=None,
            frequency_hz=float(self.entries["frequency_hz"].get()), cycles=int(self.entries["cycles"].get()),
            steps_per_cycle=int(self.entries["steps_per_cycle"].get()),
            deformation_scale=float(self.entries["deformation_scale"].get()),
        )

    def _start_solve(self) -> None:
        if self.busy:
            self.stop_event.set()
            self.run_button.configure(text="STOPPING...", state="disabled")
            self.status.set("Stopping after the current solver operation...")
            return
        try:
            config = self._config()
        except Exception as exc:
            messagebox.showerror("Invalid input", str(exc), parent=self.root)
            return
        self.busy = True
        self.stop_event.clear()
        self.live_cycles = 0.0
        self.live_records.clear()
        self.live_events.clear()
        try:
            while True:
                self._live_queue.get_nowait()
        except queue.Empty:
            pass
        self._last_live_draw = 0.0
        self._stop_playback()
        self.current_config = config
        self.run_button.configure(text="STOP ANALYSIS", state="normal")
        self.progress.start(12)
        spatial = self.spatial_backend.get()
        self.status.set(f"Solving Theory Core v1 + {spatial}…")
        if self._live_poll_job is None:
            self._live_poll_job = self.root.after(20, self._drain_live_queue)
        threading.Thread(target=self._solve_worker, args=(config, spatial), daemon=True).start()

    def _solve_worker(self, config: TensionRunConfig, spatial_backend: str) -> None:
        try:
            from simulations.fem_tension_app import (
                run_live_theory_solver,
                run_theory_spatial_solver,
            )
            from simulations.fem_tension_ui import load_fem_history
            spatial_dir = self.output_dir / "spatial"
            run_theory_spatial_solver(
                config,
                self.output_dir,
                spatial_backend,
                stop_requested=self.stop_event.is_set,
            )
            nodes, elements = load_fem_history(spatial_dir)
            self.root.after(0, lambda: self._solve_done(
                nodes, elements, spatial_backend
            ))
            def stream_record(record: dict) -> None:
                while not self.stop_event.is_set():
                    try:
                        self._live_queue.put(("record", dict(record)), timeout=0.05)
                        break
                    except queue.Full:
                        continue
                # Keep every phase sample, while pacing the producer so the
                # plot does not alias a sinusoid into a triangle or beat.
                time.sleep(0.003)

            result = run_live_theory_solver(
                config,
                stream_record,
                self.stop_event.is_set,
            )
            self._live_queue.put(("done", result))
        except InterruptedError:
            self.root.after(0, lambda: self._finish_live_analysis("Stopped by user"))
        except Exception as exc:
            detail = str(exc)
            self.root.after(0, lambda detail=detail: self._solve_failed(detail))

    def _solve_done(
        self,
        nodes: np.ndarray,
        elements: np.ndarray,
        spatial_backend: str,
    ) -> None:
        self.nodes, self.elements = nodes, elements
        self.progress.stop()
        steps = len(np.unique(elements["step"]))
        axis = self.current_config.tensile_unit_vector if self.current_config is not None else np.array([1.0, 0.0, 0.0])
        self._summary(
            f"Theory core: Theory Core v1\n"
            f"Spatial backend: {spatial_backend}\n"
            f"Tensile axis: [{axis[0]:.4g}, {axis[1]:.4g}, {axis[2]:.4g}]\n"
            f"Stress ratio R: {self.current_config.stress_ratio:.4g}\n"
            f"Phase-space P: exact spatial counting over {self.current_config.elements} chain spacings\n"
            f"Specimen survival: deterministic push-forward of the declared mu0\n"
            f"Initial measure: mu0 = delta(lambda=1,c=0)\n"
            f"Load map: q(tau) = full signed sigma_n(t) / E\n"
            f"First passage: lambda_i reaches phi''(lambda_c)=0 boundary\n"
            f"Noise / mobility / external fitted distribution: none\n"
            f"Time status: dimensionless solver time and applied-load cycles\n"
            f"Laboratory fatigue-life calibration: not yet established\n"
            f"Applied transverse stress: 0 Pa\n"
            f"Time records: {steps}\n"
            f"Control volumes: {len(np.unique(elements['element']))}\n"
            f"\nThe deterministic finite-chain solve continues until STOP."
        )
        self.notebook.select(self.post_tab)
        if self.stop_event.is_set():
            self._finish_live_analysis("Stopped")
            return
        self.status.set(f"LIVE | deterministic Theory Core + {spatial_backend}")
        self._plot()

    def _solve_failed(self, detail: str) -> None:
        self.busy = False
        self.progress.stop()
        self.run_button.configure(text="RUN ANALYSIS", state="normal")
        self.status.set("Solve failed")
        messagebox.showerror("Solver error", detail, parent=self.root)

    def _summary(self, text: str) -> None:
        self.solve_summary.configure(state="normal")
        self.solve_summary.delete("1.0", "end")
        self.solve_summary.insert("end", text)
        self.solve_summary.configure(state="disabled")

    def _show_empty_plot(self) -> None:
        self.ax.clear()
        self.ax.text(0.5, 0.5, "Run an analysis to view results", ha="center", va="center", transform=self.ax.transAxes, color=MUTED)
        self.ax.set_axis_off()
        self.canvas.draw_idle()

    def _series_by_step(self, data: np.ndarray, column: str) -> tuple[np.ndarray, np.ndarray]:
        steps = np.unique(data["step"]).astype(int)
        t = np.array([np.mean(data[data["step"] == step]["time_s"]) for step in steps])
        y = np.array([np.mean(data[data["step"] == step][column]) for step in steps])
        return t, y

    def _set_cursor_domain(
        self,
        start: float,
        end: float,
        *,
        logarithmic: bool,
        unit: str,
        x_values=None,
        y_values=None,
        value_name: str | None = None,
    ) -> None:
        start = max(float(start), np.finfo(float).tiny) if logarithmic else float(start)
        end = max(float(end), start * (1.0 + 1.0e-12))
        self._cursor_domain = (start, end, logarithmic, unit)
        self._cursor_x = None if x_values is None else np.asarray(x_values, dtype=float)
        self._cursor_y = None if y_values is None else np.asarray(y_values, dtype=float)
        self._cursor_value_name = value_name
        position = float(self.play_position.get())
        if logarithmic:
            value = 10.0 ** (np.log10(start) + position * (np.log10(end) - np.log10(start)))
        else:
            value = start + position * (end - start)
        self._cursor_line = self.ax.axvline(value, color="#d64b3c", linewidth=1.2, alpha=0.9)
        self._update_timeline_label(value, unit)

    def _update_timeline_label(self, value: float, unit: str) -> None:
        suffix = ""
        if self._cursor_x is not None and self._cursor_y is not None and self._cursor_x.size:
            order = np.argsort(self._cursor_x)
            x_values = self._cursor_x[order]
            y_values = self._cursor_y[order]
            index = int(np.clip(np.searchsorted(x_values, value, side="right") - 1, 0, len(x_values) - 1))
            suffix = f"  |  {self._cursor_value_name} = {y_values[index]:.4g}"
        if unit == "cycles":
            self.timeline_label.configure(text=f"N = {value:.5g}{suffix}")
        elif unit == "model_time":
            self.timeline_label.configure(text=f"model time = {value:.5g}{suffix}")
        else:
            self.timeline_label.configure(text=f"t = {value:.5g} s{suffix}")

    def _on_timeline(self, _value=None) -> None:
        if self._cursor_domain is None or self._cursor_line is None:
            return
        start, end, logarithmic, unit = self._cursor_domain
        position = float(self.play_position.get())
        if logarithmic:
            value = 10.0 ** (np.log10(start) + position * (np.log10(end) - np.log10(start)))
        else:
            value = start + position * (end - start)
        self._cursor_line.set_xdata([value, value])
        self._update_timeline_label(value, unit)
        self.canvas.draw_idle()

    def _toggle_playback(self) -> None:
        if self.busy:
            return
        if self.playing:
            self._stop_playback()
        else:
            self._start_playback()

    def _start_playback(self) -> None:
        if self._cursor_domain is None:
            return
        if float(self.play_position.get()) >= 1.0:
            self.play_position.set(0.0)
        self.playing = True
        self.play_button.configure(text="PAUSE")
        if self._play_job is None:
            self._play_job = self.root.after(80, self._play_tick)

    def _stop_playback(self) -> None:
        self.playing = False
        self.play_button.configure(text="PLAY")
        if self._play_job is not None:
            self.root.after_cancel(self._play_job)
            self._play_job = None

    def _reset_playback(self) -> None:
        self._stop_playback()
        self.play_position.set(0.0)
        self._on_timeline()

    def _play_tick(self) -> None:
        self._play_job = None
        if not self.playing:
            return
        position = min(1.0, float(self.play_position.get()) + 0.008)
        self.play_position.set(position)
        self._on_timeline()
        if position >= 1.0:
            self._stop_playback()
        else:
            self._play_job = self.root.after(80, self._play_tick)

    def _accept_live_record(self, record: dict) -> None:
        """Accept one aggregate emitted by the deterministic chain solver."""
        self.live_records.append(record)
        survival = float(record["survival"])
        if not self.live_events or survival != float(self.live_events[-1]["survival"]):
            self.live_events.append(record)
        self.live_cycles = float(record["cycle"])
        self.status.set(
            f"LIVE | t*={float(record['model_time']):.6g} | "
            f"N={self.live_cycles:.7g} | S_spec={survival:.5f} | "
            f"eps={float(record['strain']):.4g} | "
            f"min phi''={float(record['min_opening_eigenvalue']):.4g}"
        )
        now = time.monotonic()
        if now - self._last_live_draw >= 0.20:
            self._last_live_draw = now
            self.play_position.set(1.0)
            self._plot()

    def _drain_live_queue(self) -> None:
        """Move an ordered batch of complete solver records onto the Tk thread."""
        self._live_poll_job = None
        for _ in range(128):
            try:
                kind, payload = self._live_queue.get_nowait()
            except queue.Empty:
                break
            if kind == "record":
                self._accept_live_record(payload)
            elif kind == "done":
                self._native_solver_stopped(payload)
                return
        if self.busy or not self._live_queue.empty():
            self._live_poll_job = self.root.after(20, self._drain_live_queue)

    def _native_solver_stopped(self, _result: dict) -> None:
        self._finish_live_analysis("Native solver stopped")

    def _finish_live_analysis(self, message: str) -> None:
        self.busy = False
        self.progress.stop()
        self.run_button.configure(text="RUN ANALYSIS", state="normal")
        model_time = float(self.live_records[-1]["model_time"]) if self.live_records else 0.0
        self.status.set(f"{message} | model time={model_time:.7g} | N={self.live_cycles:.7g}")
        self._plot()

    def _close(self) -> None:
        self.stop_event.set()
        self.root.destroy()

    def _plot(self) -> None:
        if self.elements is None:
            self._show_empty_plot()
            return
        field = self.field.get()
        if field != self._last_plot_field:
            self._stop_playback()
            self.play_position.set(1.0 if self.busy else 0.0)
            self._last_plot_field = field
        self.ax.clear()
        self.ax.set_axis_on()
        self._cursor_line = None
        self._cursor_domain = None
        self._cursor_x = None
        self._cursor_y = None
        self._cursor_value_name = None
        if field in {"initiation", "survival", "hazard", "life"}:
            if not self.live_events:
                self.ax.text(0.5, 0.5, "Waiting for deterministic Theory Core records", ha="center", va="center", transform=self.ax.transAxes, color=MUTED)
                self.ax.set_axis_off(); self.canvas.draw_idle(); return
            records = list(self.live_events)
            if self.live_records and self.live_records[-1] is not records[-1]:
                records.append(self.live_records[-1])
            model_time = np.asarray([row["model_time"] for row in records], dtype=float)
            cycles = np.asarray([row["cycle"] for row in records], dtype=float)
            survival = np.asarray([row["survival"] for row in records], dtype=float)
            initiation = 1.0 - survival
            previous_survival = np.concatenate(([1.0], survival[:-1]))
            event_mass = np.maximum(0.0, previous_survival - survival)
            event_resolution = float(records[-1].get("probability_resolution", 1.0))
            probability_ylim = None
            if field == "life":
                event = event_mass > 0.0
                if not np.any(event):
                    x = np.asarray([cycles[0], cycles[-1]], dtype=float)
                    y = np.zeros(2, dtype=float)
                    self.ax.plot(x, y, color=ACCENT, linewidth=1.4)
                    self.ax.text(
                        0.5, 0.72,
                        "Exact zero first-passage mass for the declared initial measure\n"
                        "A nontrivial specimen probability requires a physical mu0",
                        ha="center", va="center", transform=self.ax.transAxes, color=MUTED,
                    )
                    probability_ylim = (0.0, min(1.0, 1.2*event_resolution))
                else:
                    x = cycles[event]
                    y = event_mass[event]
                    self.ax.vlines(x, 0.0, y, color="#79a9cf", linewidth=1.0)
                    self.ax.scatter(x, y, color=ACCENT, s=24, zorder=3)
                ylabel = "Specimen first-passage mass [-]"
                title = "Deterministic push-forward first-passage mass"
                xlabel = "Applied-load cycles N [-]"
                unit = "cycles"
            elif field == "initiation":
                x = model_time
                y = initiation
                ylabel = "Cumulative specimen first-passage mass [-]"
                title = "Deterministic first passage under entered sigma(t)"
                xlabel = "Solver model time"
                unit = "model_time"
                self.ax.step(x, y, where="post", color=ACCENT, linewidth=2.0)
                if float(np.max(y)) == 0.0:
                    self.ax.text(
                        0.02, 0.90,
                        "Pinit = 0 exactly for the declared discrete mu0 so far",
                        transform=self.ax.transAxes, color=MUTED, fontsize=9,
                    )
                    probability_ylim = (0.0, min(1.0, 1.2*event_resolution))
            elif field == "survival":
                x = model_time
                y = survival
                ylabel = "Specimen survival mass [-]"
                title = "Deterministic specimen survival under entered sigma(t)"
                xlabel = "Solver model time"
                unit = "model_time"
                self.ax.step(x, y, where="post", color=ACCENT, linewidth=2.0)
                if float(np.min(y)) == 1.0:
                    self.ax.text(
                        0.02, 0.10,
                        "Sspec = 1 exactly for the declared discrete mu0 so far",
                        transform=self.ax.transAxes, color=MUTED, fontsize=9,
                    )
                    probability_ylim = (max(0.0, 1.0 - 1.2*event_resolution), 1.001)
            else:
                event = event_mass > 0.0
                x = model_time[event]
                y = np.divide(
                    event_mass[event],
                    previous_survival[event],
                    out=np.zeros(np.count_nonzero(event), dtype=float),
                    where=previous_survival[event] > 0.0,
                )
                ylabel = "Discrete conditional first-passage flux ratio [-]"
                title = "Deterministic first-passage flux ratio under entered sigma(t)"
                xlabel = "Solver model time"
                unit = "model_time"
                self.ax.vlines(x, 0.0, y, color="#79a9cf", linewidth=1.0)
                self.ax.scatter(x, y, color=ACCENT, s=24, zorder=3)
                if not np.any(event):
                    self.ax.text(
                        0.5, 0.72,
                        "Exact zero crossing flux for the declared initial measure",
                        ha="center", va="center", transform=self.ax.transAxes, color=MUTED,
                    )
                    probability_ylim = (0.0, min(1.0, 1.2*event_resolution))
            if probability_ylim is None:
                finite_y = y[np.isfinite(y)]
                if field == "survival":
                    low = float(np.min(finite_y)); high = float(np.max(finite_y))
                    pad = max(0.15*(high-low), 0.5*event_resolution)
                    probability_ylim = (max(0.0, low-pad), min(1.001, high+pad))
                else:
                    high = float(np.max(finite_y)) if finite_y.size else 0.0
                    probability_ylim = (0.0, min(1.0, high + max(0.2*high, 0.5*event_resolution)))
            self.ax.set_ylim(*probability_ylim)
            self.ax.set_xlabel(xlabel)
            self.ax.set_ylabel(ylabel)
            self.ax.set_title(title, loc="left", fontsize=12, fontweight="bold")
            self.ax.grid(True, color="#d9e0e5", linewidth=0.7, alpha=0.8)
            cursor_name = "PMF" if field == "life" else {"initiation": "Pinit", "survival": "S", "hazard": "h"}[field]
            axis_start = float(np.min(x)) if x.size else 0.0
            axis_end = float(np.max(model_time if field != "life" else cycles))
            self._set_cursor_domain(
                axis_start,
                axis_end,
                logarithmic=False,
                unit=unit,
                x_values=x,
                y_values=y,
                value_name=cursor_name,
            )
            self.figure.tight_layout(pad=1.2)
            self.canvas.draw_idle()
            return
        live_history = bool(self.live_records) and self.current_config is not None
        if live_history:
            config = self.current_config
            records = list(self.live_records)
            t = np.asarray([row["model_time"] for row in records], dtype=float)
            stress_mpa = np.asarray([row["applied_stress_mpa"] for row in records], dtype=float)
            axial_strain = np.asarray([row["strain"] for row in records], dtype=float)
            normal_strain = np.asarray([row["normal_strain"] for row in records], dtype=float)
            intrawell_strain = np.asarray([row["intrawell_strain"] for row in records], dtype=float)
            plastic_strain = np.asarray([row["plastic_strain"] for row in records], dtype=float)
        if field == "stress":
            if live_history:
                y = stress_mpa
            else:
                t, y = self._series_by_step(self.elements, "stress_pa")
                y = y / 1.0e6 if np.nanmax(np.abs(y)) > 1.0e5 else y
            ylabel = "Normal stress [MPa]"
        elif field == "strain":
            if live_history:
                y = axial_strain
                self.ax.plot(t, normal_strain, color="#4e79a7", linewidth=1.0, label="normal opening")
                self.ax.plot(t, intrawell_strain, color="#59a14f", linewidth=1.0, label="intrawell registry")
                self.ax.plot(t, plastic_strain, color="#e15759", linewidth=1.2, label="well-index plastic")
            else:
                t, y = self._series_by_step(self.elements, "strain")
            ylabel = "Axial strain [-]"
        elif field == "diameter":
            config = self.current_config or self._config()
            if live_history:
                diameter = config.diameter_m * (1.0 - config.poisson_ratio * axial_strain)
            elif "diameter_m" in (self.elements.dtype.names or ()):
                t, diameter = self._series_by_step(self.elements, "diameter_m")
            else:
                t, axial_strain = self._series_by_step(self.elements, "strain")
                diameter = config.diameter_m * (1.0 - config.poisson_ratio * axial_strain)
            y = (diameter - config.diameter_m) * 1.0e6
            ylabel = "Diameter change [µm]"
        self.ax.plot(t, y, color=ACCENT, linewidth=1.8)
        self.ax.fill_between(t, y, alpha=0.12, color=ACCENT)
        if field == "strain" and live_history:
            self.ax.lines[-1].set_label("total")
            self.ax.legend(loc="best", fontsize=8, frameon=False)
        time_unit = "model_time" if live_history else "seconds"
        self.ax.set_xlabel("Solver model time" if live_history else "Time [s]")
        self.ax.set_ylabel(ylabel)
        self.ax.set_title(self.FIELD_LABELS[field], loc="left", fontsize=12, fontweight="bold")
        self.ax.grid(True, color="#d9e0e5", linewidth=0.7, alpha=0.8)
        self._set_cursor_domain(
            float(np.min(t)),
            float(np.max(t)),
            logarithmic=False,
            unit=time_unit,
            x_values=t,
            y_values=y,
            value_name=self.FIELD_LABELS[field],
        )
        self.figure.tight_layout(pad=1.2)
        self.canvas.draw_idle()

    def _open_geometry(self) -> None:
        path = filedialog.askopenfilename(parent=self.root, filetypes=[("Mesh files", "*.obj *.stl *.ply *.vtk")])
        if not path:
            return
        self.geometry_path = Path(path)
        self._view_geometry(self.geometry_path)

    def _view_default_geometry(self) -> None:
        self.geometry_path = self._default_cylinder_path()
        self._view_geometry(self.geometry_path)

    def _view_geometry(self, path: Path) -> None:
        from simulations.mesh_viewer import MeshViewport, load_mesh, orient_local_x

        try:
            mesh = load_mesh(path)
            loading_axis = self._config().tensile_unit_vector
            if path.resolve() == self._default_cylinder_path().resolve():
                mesh = orient_local_x(mesh, loading_axis)
            self.mesh_label.configure(text=f"{path.name} · {len(mesh.vertices)} vertices · {mesh.dimension}D")
            MeshViewport(mesh, display_dimension=3, loading_axis=loading_axis).figure.show()
        except Exception as exc:
            messagebox.showerror("Mesh error", str(exc), parent=self.root)

    def _open_project(self) -> None:
        path = filedialog.askopenfilename(parent=self.root, filetypes=[("Al Fatigue project", "*.ftgsim")])
        if not path:
            return
        self._load_project(Path(path))

    def _load_project(self, path: Path) -> None:
        try:
            from simulations.fem_tension_app import config_from_ftgsim
            from simulations.fem_tension_ui import load_fem_history
            from simulations.ftgsim_format import extract_geometry, extract_results, open_ftgsim

            config, _geometry, _display = config_from_ftgsim(path)
            self.current_config = config
            values = {
                "length_mm": config.length_mm, "diameter_mm": config.diameter_mm,
                "young_gpa": config.young_gpa, "poisson_ratio": config.poisson_ratio,
                "loading_direction": f"{config.loading_h} {config.loading_k} {config.loading_l}",
                "tensile_direction": " ".join(f"{value:.8g}" for value in config.tensile_unit_vector),
                "elements": config.elements, "stress_mean_mpa": config.stress_mean_mpa,
                "stress_amplitude_mpa": config.stress_amplitude_mpa,
                "frequency_hz": config.frequency_hz, "cycles": config.cycles,
                "steps_per_cycle": config.steps_per_cycle, "deformation_scale": config.deformation_scale,
            }
            for key, value in values.items():
                if key in self.entries:
                    self.entries[key].delete(0, "end"); self.entries[key].insert(0, str(value))
            bundle = open_ftgsim(path)
            bundle_id = str(bundle.manifest.get("bundle_id", path.stem))
            self.output_dir = self._default_output_dir() / "opened" / bundle_id
            extract_results(bundle, self.output_dir)
            geometry_files = extract_geometry(bundle, self.output_dir / "geometry")
            if geometry_files:
                self.geometry_path = geometry_files[0]
                self.mesh_label.configure(text=f"{self.geometry_path.name} · saved 3D geometry")
            self.nodes, self.elements = load_fem_history(self.output_dir)
            self.status.set(f"Opened {path.name}")
            self.notebook.select(self.post_tab); self._plot()
        except Exception as exc:
            messagebox.showerror("Open project failed", str(exc), parent=self.root)

    def _save_project(self) -> None:
        if self.elements is None:
            messagebox.showinfo("Save project", "Run an analysis first.", parent=self.root)
            return
        path = filedialog.asksaveasfilename(parent=self.root, defaultextension=".ftgsim", filetypes=[("Al Fatigue project", "*.ftgsim")])
        if not path:
            return
        try:
            from simulations.fem_tension_app import save_tension_ftgsim

            save_tension_ftgsim(
                Path(path),
                self._config(),
                self.output_dir,
                view="3D",
                field=self.field.get(),
                geometry_source=self.geometry_path if self.geometry_path.is_file() else None,
            )
            self.status.set(f"Saved {Path(path).name}")
        except Exception as exc:
            messagebox.showerror("Save project failed", str(exc), parent=self.root)

    def run(self) -> None:
        self.root.mainloop()


def main() -> None:
    if not acquire_single_instance():
        return
    project = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    DesktopApp(project_path=project).run()


if __name__ == "__main__":
    main()
