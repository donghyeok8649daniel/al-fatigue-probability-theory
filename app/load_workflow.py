"""Mesh-face boundary-condition preparation, intentionally separate from PDE."""
from dataclasses import dataclass
import tkinter as tk
from tkinter import ttk, messagebox

import numpy as np
from .tensor_load import default_tensor_expressions, compile_tensor_matrix, evaluate_tensor


@dataclass(frozen=True)
class FaceLoad:
    region: str
    normal_mean_mpa: float
    normal_amplitude_mpa: float
    shear_mean_mpa: float
    shear_amplitude_mpa: float
    face_indices: tuple[int, ...]
    area_mm2: float
    tensor_expression: str


class LoadWorkflow:
    REGIONS = ("all", "top", "bottom", "lateral")

    def __init__(self, app, tab):
        self.app = app
        self.region_code = "top"
        self.region = tk.StringVar(master=app.root, value=app._tr("load.region_top"))
        self.normal_mean = tk.StringVar(master=app.root, value="0")
        self.normal_amplitude = tk.StringVar(master=app.root, value="0")
        self.shear_mean = tk.StringVar(master=app.root, value="0")
        self.shear_amplitude = tk.StringVar(master=app.root, value="0")
        self.summary = tk.StringVar(master=app.root)
        self.status = tk.StringVar(master=app.root)
        self.tensor_text = tk.StringVar(master=app.root, value=default_tensor_expressions())
        self.tensor_preview = tk.StringVar(master=app.root)
        app._bind_text(ttk.Label(tab, style="Section.TLabel"), "section.face_load").pack(anchor="w", padx=20, pady=(20, 8))
        app._bind_text(ttk.Label(tab, wraplength=700), "load.scope").pack(anchor="w", padx=20, pady=(0, 12))
        row = ttk.Frame(tab); row.pack(fill="x", padx=20, pady=5)
        app._bind_text(ttk.Label(row), "load.face_region").pack(side="left")
        self.selector = ttk.Combobox(row, textvariable=self.region, state="readonly", width=20)
        self.selector.pack(side="left", padx=8)
        self.selector.bind("<<ComboboxSelected>>", self._selected)
        fields = (("load.normal_mean", self.normal_mean), ("load.normal_amplitude", self.normal_amplitude),
                  ("load.shear_mean", self.shear_mean), ("load.shear_amplitude", self.shear_amplitude))
        for key, variable in fields:
            row = ttk.Frame(tab); row.pack(fill="x", padx=20, pady=5)
            app._bind_text(ttk.Label(row), key).pack(side="left")
            ttk.Entry(row, textvariable=variable, width=16).pack(side="right")
            app._bind_text(ttk.Label(row), "unit.mpa").pack(side="right", padx=8)
        app._bind_text(ttk.Button(tab, command=self.apply), "load.apply").pack(anchor="w", padx=20, pady=12)
        app._bind_text(ttk.Label(tab, style="Section.TLabel"), "load.tensor_section").pack(anchor="w", padx=20, pady=(14, 5))
        app._bind_text(ttk.Label(tab, wraplength=700), "load.tensor_help").pack(anchor="w", padx=20, pady=(0, 5))
        ttk.Entry(tab, textvariable=self.tensor_text, width=105).pack(fill="x", padx=20, pady=4)
        app._bind_text(ttk.Button(tab, command=self.preview_tensor), "load.tensor_preview_button").pack(anchor="w", padx=20, pady=6)
        ttk.Label(tab, textvariable=self.tensor_preview, wraplength=700).pack(anchor="w", padx=20, pady=4)
        ttk.Label(tab, textvariable=self.summary, wraplength=700).pack(anchor="w", padx=20, pady=4)
        ttk.Label(tab, textvariable=self.status, wraplength=700).pack(anchor="w", padx=20, pady=4)
        app._bind_text(ttk.Button(tab, command=lambda: app.notebook.select(app.pre_tab)), "geometry.next_pre").pack(anchor="e", padx=20, pady=12)
        self.refresh()

    def preview_tensor(self):
        try:
            compiled = compile_tensor_matrix(self.tensor_text.get())
            frequency = float(self.app.entries.get("model_frequency", tk.StringVar(value="25")).get())
            values = evaluate_tensor(compiled, t=0.0, frequency=frequency,
                normal_mean=float(self.normal_mean.get()), normal_amplitude=float(self.normal_amplitude.get()),
                shear_mean=float(self.shear_mean.get()), shear_amplitude=float(self.shear_amplitude.get()))
        except (ValueError, SyntaxError, TypeError):
            self.tensor_preview.set(self.app._tr("load.tensor_invalid")); return
        self.tensor_preview.set(self.app._tr("load.tensor_result", matrix=np.array2string(values, precision=5)))

    def _selected(self, _event=None):
        labels = {self.app._tr("load.region_all"): "all", self.app._tr("load.region_top"): "top",
                  self.app._tr("load.region_bottom"): "bottom", self.app._tr("load.region_lateral"): "lateral"}
        self.region_code = labels.get(self.region.get(), "top")
        self.refresh()

    def _indices(self):
        mesh = self.app.geometry_workflow.mesh
        if mesh is None:
            return np.empty(0, dtype=int)
        tri = mesh.vertices[mesh.faces]
        z = tri[:, :, 2].mean(axis=1)
        lo, hi = mesh.vertices[:, 2].min(), mesh.vertices[:, 2].max()
        tol = max((hi-lo)*1e-8, 1e-12)
        top, bottom = z >= hi-tol, z <= lo+tol
        if self.region_code == "top": mask = top
        elif self.region_code == "bottom": mask = bottom
        elif self.region_code == "lateral": mask = ~(top | bottom)
        else: mask = np.ones(len(tri), dtype=bool)
        return np.flatnonzero(mask)

    def refresh(self):
        values = tuple(self.app._tr(f"load.region_{x}") for x in self.REGIONS)
        self.selector.configure(values=values)
        self.region.set(self.app._tr(f"load.region_{self.region_code}"))
        ids = self._indices()
        mesh = self.app.geometry_workflow.mesh
        area = 0.0 if mesh is None else float(mesh.areas[ids].sum())
        self.summary.set(self.app._tr("load.selection_summary", faces=len(ids), area=f"{area:.6g}"))
        if hasattr(self, "applied") and self.applied is not None:
            self.status.set(self.app._tr("load.prepared", faces=len(self.applied.face_indices)))
        else:
            self.status.set(self.app._tr("load.not_prepared"))

    def apply(self):
        try:
            values = [float(x.get()) for x in (self.normal_mean, self.normal_amplitude, self.shear_mean, self.shear_amplitude)]
            if not all(np.isfinite(x) for x in values): raise ValueError
        except ValueError:
            messagebox.showerror(self.app._tr("tab.load"), self.app._tr("load.invalid"), parent=self.app.root)
            return
        ids = self._indices()
        mesh = self.app.geometry_workflow.mesh
        area = 0.0 if mesh is None else float(mesh.areas[ids].sum())
        try:
            compile_tensor_matrix(self.tensor_text.get())
        except (ValueError, SyntaxError):
            self.tensor_preview.set(self.app._tr("load.tensor_invalid"))
            return
        self.applied = FaceLoad(self.region_code, *values, tuple(map(int, ids)), area, self.tensor_text.get())
        self.status.set(self.app._tr("load.prepared", faces=len(ids)))
