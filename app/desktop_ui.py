"""Lightweight Pre/Solve/Post desktop UI for Theory Core v1 and 1D solvers."""
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


class DesktopApp:
    """Responsive engineering workspace with explicit Pre/Solve/Post stages."""

    PARAMS = (
        ("length_mm", "Specimen length", "50", "mm"),
        ("width_mm", "Specimen width", "10", "mm"),
        ("thickness_mm", "Specimen thickness", "1", "mm"),
        ("young_gpa", "Young's modulus", "69", "GPa"),
        ("loading_direction", "Crystal direction", "1 0 0", "[h k l]"),
        ("elements", "Control volumes", "40", "cells"),
        ("stress_mean_mpa", "Mean normal stress", "50", "MPa"),
        ("stress_amplitude_mpa", "Normal stress amplitude", "100", "MPa"),
        ("frequency_hz", "Loading frequency", "20", "Hz"),
        ("cycles", "Loading cycles", "2", "cycles"),
        ("steps_per_cycle", "Resolution", "80", "steps/cycle"),
        ("deformation_scale", "Display deformation", "1", "x"),
    )

    FIELD_LABELS = {
        "stress": "Normal stress",
        "strain": "Axial strain",
        "initiation": "First passage",
        "survival": "Survival",
        "hazard": "Hazard",
    }

    def __init__(self, project_path: Path | None = None) -> None:
        self.root = tk.Tk()
        self.root.title("Al Fatigue — Theory Core v1")
        self.root.configure(bg=APP_BG)
        self.root.minsize(980, 640)
        self._center(1180, 760)

        self.output_dir = self._default_output_dir()
        self.spatial_backend = tk.StringVar(value="FVM")
        self.field = tk.StringVar(value="stress")
        self.entries: dict[str, ttk.Entry] = {}
        self.nodes: np.ndarray | None = None
        self.elements: np.ndarray | None = None
        self.initiation: np.ndarray | None = None
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
        subtitle = ttk.Label(header, text="Theory Core v1 always active · axial normal loading", style="SubHeader.TLabel")
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
        self.tree.insert(pre, "end", text="  Geometry / Mesh")
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
        for i, (key, label, default, unit) in enumerate(self.PARAMS):
            column_group = 0 if i < 6 else 3
            row = i % 6 + 1
            ttk.Label(self.pre_tab, text=label, style="Property.TLabel").grid(row=row, column=column_group, sticky="w", padx=(18, 8), pady=7)
            entry = ttk.Entry(self.pre_tab, width=14)
            entry.insert(0, default)
            entry.grid(row=row, column=column_group + 1, sticky="ew", padx=4, pady=7)
            ttk.Label(self.pre_tab, text=unit, style="Unit.TLabel").grid(row=row, column=column_group + 2, sticky="w", padx=(4, 22), pady=7)
            self.entries[key] = entry
        self.pre_tab.columnconfigure(1, weight=1)
        self.pre_tab.columnconfigure(4, weight=1)
        mesh = ttk.LabelFrame(self.pre_tab, text="Geometry and mesh", padding=12)
        mesh.grid(row=8, column=0, columnspan=6, sticky="ew", padx=18, pady=16)
        ttk.Button(mesh, text="Open OBJ / STL / PLY / VTK", command=self._open_geometry).pack(side="left")
        self.mesh_label = ttk.Label(mesh, text="Uniform 1D control-volume mesh", foreground=MUTED)
        self.mesh_label.pack(side="left", padx=14)

    def _solve_tab(self) -> None:
        left = ttk.Frame(self.solve_tab, style="Panel.TFrame")
        left.pack(side="left", fill="y", padx=20, pady=18)
        ttk.Label(left, text="Theory Core v1", style="Section.TLabel").pack(anchor="w", pady=(0, 4))
        ttk.Label(left, text="Always active", foreground=ACCENT, background=PANEL_BG, font=("Segoe UI", 9, "bold")).pack(anchor="w", pady=(0, 16))
        ttk.Label(left, text="Spatial discretization", style="Section.TLabel").pack(anchor="w", pady=(0, 8))
        combo = ttk.Combobox(left, textvariable=self.spatial_backend, values=("FVM", "FEM"), state="readonly", width=24)
        combo.pack(anchor="w", pady=(0, 14))
        ttk.Label(left, text="Theory computes probability, plastic memory,\nsurvival, hazard and first passage.\nFVM/FEM computes the spatial reference field.", style="Property.TLabel", justify="left").pack(anchor="w", pady=(0, 18))
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
        toolbar.pack(fill="x", padx=12, pady=(10, 4))
        ttk.Label(toolbar, text="Result field", style="Section.TLabel").pack(side="left", padx=(4, 12))
        for key, label in self.FIELD_LABELS.items():
            ttk.Radiobutton(toolbar, text=label, value=key, variable=self.field, style="Field.Toolbutton", command=self._plot).pack(side="left", padx=2)
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

        direction = self.entries["loading_direction"].get().replace(",", " ").split()
        if len(direction) != 3:
            raise ValueError("Crystal direction requires three integers: h k l")
        h, k, l = (int(x) for x in direction)
        return TensionRunConfig(
            length_mm=float(self.entries["length_mm"].get()), width_mm=float(self.entries["width_mm"].get()),
            thickness_mm=float(self.entries["thickness_mm"].get()), young_gpa=float(self.entries["young_gpa"].get()),
            loading_h=h, loading_k=k, loading_l=l, elements=int(self.entries["elements"].get()),
            stress_mean_mpa=float(self.entries["stress_mean_mpa"].get()),
            stress_amplitude_mpa=float(self.entries["stress_amplitude_mpa"].get()),
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
            self.root.after(0, lambda: self._solve_done(nodes, elements, initiation, spatial_backend))
        except Exception as exc:
            detail = str(exc)
            self.root.after(0, lambda detail=detail: self._solve_failed(detail))

    def _solve_done(self, nodes: np.ndarray, elements: np.ndarray, initiation: np.ndarray | None, spatial_backend: str) -> None:
        self.nodes, self.elements, self.initiation = nodes, elements, initiation
        self.busy = False
        self.progress.stop()
        self.run_button.configure(state="normal")
        steps = len(np.unique(elements["step"]))
        fp = "n/a"
        if initiation is not None:
            fp = f"{float(np.nanmax(initiation['initiation_probability'])):.3f}"
        self._summary(f"Theory core: Theory Core v1\nSpatial backend: {spatial_backend}\nTime records: {steps}\nControl volumes: {len(np.unique(elements['element']))}\nFinal/maximum first passage: {fp}\n\nSolve completed successfully.")
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
        if field == "stress":
            t, y = self._series_by_step(self.elements, "stress_pa")
            y = y / 1.0e6 if np.nanmax(np.abs(y)) > 1.0e5 else y
            ylabel = "Normal stress [MPa]" if np.nanmax(np.abs(y)) > 1.0e-3 else "Normalized tensile force"
        elif field == "strain":
            t, y = self._series_by_step(self.elements, "strain")
            ylabel = "Axial strain [-]"
        else:
            if self.initiation is None:
                self.ax.text(0.5, 0.5, "Probability fields require a Theory Core v1 solve", ha="center", va="center", transform=self.ax.transAxes, color=MUTED)
                self.ax.set_axis_off(); self.canvas.draw_idle(); return
            column = {"initiation": "initiation_probability", "survival": "survival", "hazard": "hazard_per_s"}[field]
            t, y = self._series_by_step(self.initiation, column)
            ylabel = {"initiation": "First-passage probability [-]", "survival": "Survival probability [-]", "hazard": "Hazard [1 / model time]"}[field]
        self.ax.plot(t, y, color=ACCENT, linewidth=1.8)
        self.ax.fill_between(t, y, alpha=0.12, color=ACCENT)
        self.ax.set_xlabel("Model time")
        self.ax.set_ylabel(ylabel)
        self.ax.set_title(self.FIELD_LABELS[field], loc="left", fontsize=12, fontweight="bold")
        self.ax.grid(True, color="#d9e0e5", linewidth=0.7, alpha=0.8)
        self.figure.tight_layout(pad=1.2)
        self.canvas.draw_idle()

    def _open_geometry(self) -> None:
        from simulations.mesh_viewer import MeshViewport, load_mesh

        path = filedialog.askopenfilename(parent=self.root, filetypes=[("Mesh files", "*.obj *.stl *.ply *.vtk")])
        if not path:
            return
        try:
            mesh = load_mesh(Path(path))
            self.mesh_label.configure(text=f"{Path(path).name} · {len(mesh.vertices)} vertices · {mesh.dimension}D")
            MeshViewport(mesh).figure.show()
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
            from simulations.ftgsim_format import extract_results, open_ftgsim
            from simulations.visualize_fem1d import load_numeric_csv

            config, _geometry, _display = config_from_ftgsim(path)
            values = {
                "length_mm": config.length_mm, "width_mm": config.width_mm,
                "thickness_mm": config.thickness_mm, "young_gpa": config.young_gpa,
                "loading_direction": f"{config.loading_h} {config.loading_k} {config.loading_l}",
                "elements": config.elements, "stress_mean_mpa": config.stress_mean_mpa,
                "stress_amplitude_mpa": config.stress_amplitude_mpa,
                "frequency_hz": config.frequency_hz, "cycles": config.cycles,
                "steps_per_cycle": config.steps_per_cycle, "deformation_scale": config.deformation_scale,
            }
            for key, value in values.items():
                self.entries[key].delete(0, "end"); self.entries[key].insert(0, str(value))
            bundle = open_ftgsim(path)
            bundle_id = str(bundle.manifest.get("bundle_id", path.stem))
            self.output_dir = self._default_output_dir() / "opened" / bundle_id
            extract_results(bundle, self.output_dir)
            self.nodes, self.elements = load_fem_history(self.output_dir)
            init_path = self.output_dir / "initiation_elements.csv"
            self.initiation = load_numeric_csv(init_path) if init_path.is_file() else None
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

            save_tension_ftgsim(Path(path), self._config(), self.output_dir, view="2D", field=self.field.get())
            self.status.set(f"Saved {Path(path).name}")
        except Exception as exc:
            messagebox.showerror("Save project failed", str(exc), parent=self.root)

    def run(self) -> None:
        self.root.mainloop()


def main() -> None:
    project = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    DesktopApp(project_path=project).run()


if __name__ == "__main__":
    main()
