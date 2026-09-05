"""Lightweight Pre/Solve/Post UI for uniaxial theory and 3D geometry."""
from __future__ import annotations

import os
import sys
import threading
import tkinter as tk
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
        ("theory_stress_scale_mpa", "Theory stress scale", "40", "MPa / model force"),
        ("frequency_hz", "Loading frequency", "20", "Hz"),
        ("cycles", "Loading cycles", "2", "cycles"),
        ("steps_per_cycle", "Resolution", "80", "steps/cycle"),
        ("deformation_scale", "Display deformation", "1", "x"),
    )

    FIELD_LABELS = {
        "stress": "Normal stress",
        "strain": "Axial strain",
        "diameter": "Diameter change",
        "initiation": "First passage",
        "survival": "Survival",
        "hazard": "Hazard",
        "life": "Life distribution",
        "sn": "S-N initiation",
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
        self.initiation: np.ndarray | None = None
        self.sn_curve: np.ndarray | None = None
        self.life_distribution: np.ndarray | None = None
        self.current_config = None
        self.busy = False

        self._styles()
        self._header()
        self._workspace()
        self._statusbar()
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
        ttk.Label(left, text="Theory computes probability and first passage.\nFVM/FEM applies stress only along the entered axis.\nFCC slip projection feeds a separate dislocation-based S-N bridge.\nTheory stress scale remains an uncalibrated mechanism-screening input.\nDiameter change is Poisson kinematics; transverse stress is zero.", style="Property.TLabel", justify="left").pack(anchor="w", pady=(0, 18))
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
            theory_stress_scale_mpa=float(self.entries["theory_stress_scale_mpa"].get()),
            frequency_hz=float(self.entries["frequency_hz"].get()), cycles=int(self.entries["cycles"].get()),
            steps_per_cycle=int(self.entries["steps_per_cycle"].get()),
            deformation_scale=float(self.entries["deformation_scale"].get()),
        )

    def _start_solve(self) -> None:
        if self.busy:
            return
        try:
            config = self._config()
        except Exception as exc:
            messagebox.showerror("Invalid input", str(exc), parent=self.root)
            return
        self.busy = True
        self.current_config = config
        self.run_button.configure(state="disabled")
        self.progress.start(12)
        spatial = self.spatial_backend.get()
        self.status.set(f"Solving Theory Core v1 + {spatial}…")
        threading.Thread(target=self._solve_worker, args=(config, spatial), daemon=True).start()

    def _solve_worker(self, config: TensionRunConfig, spatial_backend: str) -> None:
        try:
            from simulations.fem_tension_app import run_theory_spatial_solver
            from simulations.fem_tension_ui import load_fem_history
            from simulations.visualize_fem1d import load_numeric_csv

            theory_dir = self.output_dir / "theory"
            spatial_dir = self.output_dir / "spatial"
            run_theory_spatial_solver(config, self.output_dir, spatial_backend)
            nodes, elements = load_fem_history(spatial_dir)
            init_path = theory_dir / "initiation_elements.csv"
            initiation = load_numeric_csv(init_path) if init_path.is_file() else None
            sn_path = self.output_dir / "sn_curve.csv"
            sn_curve = load_numeric_csv(sn_path) if sn_path.is_file() else None
            life_path = self.output_dir / "life_distribution.csv"
            life_distribution = load_numeric_csv(life_path) if life_path.is_file() else None
            self.root.after(0, lambda: self._solve_done(
                nodes, elements, initiation, life_distribution, sn_curve, spatial_backend
            ))
        except Exception as exc:
            detail = str(exc)
            self.root.after(0, lambda detail=detail: self._solve_failed(detail))

    def _solve_done(
        self,
        nodes: np.ndarray,
        elements: np.ndarray,
        initiation: np.ndarray | None,
        life_distribution: np.ndarray | None,
        sn_curve: np.ndarray | None,
        spatial_backend: str,
    ) -> None:
        self.nodes, self.elements, self.initiation, self.life_distribution, self.sn_curve = (
            nodes, elements, initiation, life_distribution, sn_curve
        )
        self.busy = False
        self.progress.stop()
        self.run_button.configure(state="normal")
        steps = len(np.unique(elements["step"]))
        fp = "n/a"
        first_passage = "none in simulated interval"
        if initiation is not None:
            fp = f"{float(np.nanmax(initiation['initiation_probability'])):.3f}"
            crossed = initiation[initiation["initiation_probability"] > 0.0]
            if crossed.size:
                first_time_s = float(crossed[0]["time_s"])
                first_passage = (
                    f"{first_time_s:.6g} s "
                    f"({first_time_s * self.current_config.frequency_hz:.4g} cycles)"
                )
        axis = self.current_config.tensile_unit_vector if self.current_config is not None else np.array([1.0, 0.0, 0.0])
        active_slip = self.current_config.fcc_slip_systems[0]
        max_tau_a = active_slip.schmid_factor * self.current_config.stress_amplitude_mpa
        tmw = self.current_config.tmw_initiation
        tmw_life = (
            f"{tmw.cycles_to_initiation:.4g} cycles"
            if np.isfinite(tmw.cycles_to_initiation)
            else "not activated in this mechanism"
        )
        self._summary(
            f"Theory core: Theory Core v1\n"
            f"Spatial backend: {spatial_backend}\n"
            f"Tensile axis: [{axis[0]:.4g}, {axis[1]:.4g}, {axis[2]:.4g}]\n"
            f"Stress ratio R: {self.current_config.stress_ratio:.4g}\n"
            f"Maximum FCC Schmid factor: {active_slip.schmid_factor:.4g}\n"
            f"Maximum resolved-shear amplitude: {max_tau_a:.4g} MPa\n"
            f"FCC dislocation-initiation N50 scale: {tmw_life}\n"
            f"Life law: Theory Core v1 empirical first passage (right-censored)\n"
            f"S-N bridge status: conditional physical scale; no life-distribution fit\n"
            f"Applied transverse stress: 0 Pa\n"
            f"Time records: {steps}\n"
            f"Control volumes: {len(np.unique(elements['element']))}\n"
            f"First detected passage: {first_passage}\n"
            f"Final/maximum first passage: {fp}\n\nSolve completed successfully."
        )
        self.status.set(f"Solved · Theory Core v1 + {spatial_backend} · {steps} records")
        self.notebook.select(self.post_tab)
        self._plot()

    def _solve_failed(self, detail: str) -> None:
        self.busy = False
        self.progress.stop()
        self.run_button.configure(state="normal")
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

    def _plot(self) -> None:
        if self.elements is None:
            self._show_empty_plot()
            return
        field = self.field.get()
        self.ax.clear()
        self.ax.set_axis_on()
        if field == "life":
            if self.life_distribution is None:
                self.ax.text(0.5, 0.5, "Life distribution requires an analysis", ha="center", va="center", transform=self.ax.transAxes, color=MUTED)
                self.ax.set_axis_off(); self.canvas.draw_idle(); return
            cycles = np.asarray(self.life_distribution["initiation_cycles"], dtype=float)
            cdf = np.asarray(self.life_distribution["cumulative_probability"], dtype=float)
            finite = np.isfinite(cycles) & np.isfinite(cdf)
            if not np.any(finite):
                self.ax.text(0.5, 0.5, "No finite first-passage life for this mechanism", ha="center", va="center", transform=self.ax.transAxes, color=MUTED)
                self.ax.set_axis_off(); self.canvas.draw_idle(); return
            order = np.argsort(cycles[finite])
            self.ax.step(cycles[finite][order], cdf[finite][order], where="post", color=ACCENT, linewidth=2.0, label="Initiation CDF")
            self.ax.step(cycles[finite][order], 1.0 - cdf[finite][order], where="post", color="#66727d", linewidth=1.6, label="Survival")
            self.ax.set_xscale("log")
            self.ax.set_ylim(0.0, 1.02)
            self.ax.set_xlabel("Cycles to crack initiation")
            self.ax.set_ylabel("Probability [-]")
            self.ax.set_title("Life first-passage distribution", loc="left", fontsize=12, fontweight="bold")
            self.ax.legend(frameon=False)
            self.ax.grid(True, which="both", color="#d9e0e5", linewidth=0.7, alpha=0.8)
            self.figure.tight_layout(pad=1.2)
            self.canvas.draw_idle()
            return
        if field == "sn":
            if self.sn_curve is None:
                self.ax.text(0.5, 0.5, "S-N result requires an analysis", ha="center", va="center", transform=self.ax.transAxes, color=MUTED)
                self.ax.set_axis_off(); self.canvas.draw_idle(); return
            finite = self.sn_curve[np.isfinite(self.sn_curve["n50_cycles"])]
            finite = finite[finite["n50_cycles"] > 0.0]
            if not finite.size:
                self.ax.text(0.5, 0.5, "This mechanism is below its effective shear threshold", ha="center", va="center", transform=self.ax.transAxes, color=MUTED)
                self.ax.set_axis_off(); self.canvas.draw_idle(); return
            cycles = np.asarray(finite["n50_cycles"], dtype=float)
            n10 = np.asarray(finite["n10_cycles"], dtype=float)
            n80 = np.asarray(finite["n80_cycles"], dtype=float)
            amplitude = np.asarray(finite["axial_stress_amplitude_mpa"], dtype=float)
            order = np.argsort(cycles)
            self.ax.semilogx(n10[order], amplitude[order], color="#79a9cf", linewidth=1.3, label="N10")
            self.ax.semilogx(cycles[order], amplitude[order], color=ACCENT, linewidth=2.2, label="N50")
            self.ax.semilogx(n80[order], amplitude[order], color="#263746", linewidth=1.3, label="N80")
            estimate = self.current_config.tmw_initiation
            if np.isfinite(estimate.cycles_to_initiation):
                self.ax.scatter(
                    [estimate.cycles_to_initiation],
                    [self.current_config.stress_amplitude_mpa],
                    s=52, color="#d64b3c", zorder=3, label="Current load",
                )
            self.ax.legend(frameon=False)
            self.ax.set_xlabel("Cycles to crack initiation, Nc")
            self.ax.set_ylabel("Axial stress amplitude [MPa]")
            self.ax.set_title("FCC slip-band initiation probability quantiles", loc="left", fontsize=12, fontweight="bold")
            self.ax.grid(True, which="both", color="#d9e0e5", linewidth=0.7, alpha=0.8)
            self.figure.tight_layout(pad=1.2)
            self.canvas.draw_idle()
            return
        if field == "stress":
            t, y = self._series_by_step(self.elements, "stress_pa")
            y = y / 1.0e6 if np.nanmax(np.abs(y)) > 1.0e5 else y
            ylabel = "Normal stress [MPa]" if np.nanmax(np.abs(y)) > 1.0e-3 else "Normalized tensile force"
        elif field == "strain":
            t, y = self._series_by_step(self.elements, "strain")
            ylabel = "Axial strain [-]"
        elif field == "diameter":
            config = self.current_config or self._config()
            if "diameter_m" in (self.elements.dtype.names or ()):
                t, diameter = self._series_by_step(self.elements, "diameter_m")
            else:
                t, axial_strain = self._series_by_step(self.elements, "strain")
                diameter = config.diameter_m * (1.0 - config.poisson_ratio * axial_strain)
            y = (diameter - config.diameter_m) * 1.0e6
            ylabel = "Diameter change [µm]"
        else:
            if self.initiation is None:
                self.ax.text(0.5, 0.5, "Probability fields require a Theory Core v1 solve", ha="center", va="center", transform=self.ax.transAxes, color=MUTED)
                self.ax.set_axis_off(); self.canvas.draw_idle(); return
            column = {"initiation": "initiation_probability", "survival": "survival", "hazard": "hazard_per_s"}[field]
            t, y = self._series_by_step(self.initiation, column)
            ylabel = {"initiation": "First-passage probability [-]", "survival": "Survival probability [-]", "hazard": "Hazard [1/s]"}[field]
        self.ax.plot(t, y, color=ACCENT, linewidth=1.8)
        self.ax.fill_between(t, y, alpha=0.12, color=ACCENT)
        self.ax.set_xlabel("Time [s]")
        self.ax.set_ylabel(ylabel)
        self.ax.set_title(self.FIELD_LABELS[field], loc="left", fontsize=12, fontweight="bold")
        self.ax.grid(True, color="#d9e0e5", linewidth=0.7, alpha=0.8)
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
            from simulations.visualize_fem1d import load_numeric_csv

            config, _geometry, _display = config_from_ftgsim(path)
            self.current_config = config
            values = {
                "length_mm": config.length_mm, "diameter_mm": config.diameter_mm,
                "young_gpa": config.young_gpa, "poisson_ratio": config.poisson_ratio,
                "loading_direction": f"{config.loading_h} {config.loading_k} {config.loading_l}",
                "tensile_direction": " ".join(f"{value:.8g}" for value in config.tensile_unit_vector),
                "elements": config.elements, "stress_mean_mpa": config.stress_mean_mpa,
                "stress_amplitude_mpa": config.stress_amplitude_mpa,
                "theory_stress_scale_mpa": config.theory_stress_scale_mpa,
                "frequency_hz": config.frequency_hz, "cycles": config.cycles,
                "steps_per_cycle": config.steps_per_cycle, "deformation_scale": config.deformation_scale,
            }
            for key, value in values.items():
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
            init_path = self.output_dir / "initiation_elements.csv"
            self.initiation = load_numeric_csv(init_path) if init_path.is_file() else None
            sn_path = self.output_dir / "sn_curve.csv"
            self.sn_curve = load_numeric_csv(sn_path) if sn_path.is_file() else None
            life_path = self.output_dir / "life_distribution.csv"
            self.life_distribution = load_numeric_csv(life_path) if life_path.is_file() else None
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
