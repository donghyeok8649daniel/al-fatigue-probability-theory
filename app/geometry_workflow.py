"""Model/mesh stages, isolated from probability-state and solver inputs."""
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from .specimen_mesh import cylinder, load_surface, refine_surface


class GeometryWorkflow:
    def __init__(self, app, model_tab, mesh_tab):
        self.app = app
        self.geometry = cylinder()
        self.cylinder_dimensions = (5., 30.)
        self.mesh = None
        self.radius = tk.StringVar(master=app.root, value='5')
        self.length = tk.StringVar(master=app.root, value='30')
        self.unit_scale = tk.StringVar(master=app.root, value='1')
        self.target = tk.StringVar(master=app.root, value='3')
        self.description = tk.StringVar(master=app.root)
        self.mesh_description = tk.StringVar(master=app.root)
        for key, variable in [('geometry.radius', self.radius), ('geometry.length', self.length),
                              ('geometry.scale', self.unit_scale)]:
            row = ttk.Frame(model_tab); row.pack(fill='x', padx=20, pady=8)
            app._bind_text(ttk.Label(row), key).pack(side='left')
            ttk.Entry(row, textvariable=variable, width=15).pack(side='right')
        row = ttk.Frame(model_tab); row.pack(fill='x', padx=20, pady=10)
        app._bind_text(ttk.Button(row, command=self.use_cylinder), 'geometry.cylinder').pack(side='left', padx=4)
        app._bind_text(ttk.Button(row, command=self.import_model), 'geometry.import').pack(side='left', padx=4)
        ttk.Label(model_tab, textvariable=self.description, wraplength=650).pack(anchor='w', padx=20, pady=12)
        app._bind_text(ttk.Label(model_tab, wraplength=650), 'geometry.formats').pack(anchor='w', padx=20, pady=8)
        app._bind_text(ttk.Button(model_tab, command=lambda: app.notebook.select(mesh_tab)),
                       'geometry.next_mesh').pack(anchor='e', padx=20, pady=12)
        row = ttk.Frame(mesh_tab); row.pack(fill='x', padx=12, pady=8)
        app._bind_text(ttk.Label(row), 'geometry.target').pack(side='left')
        ttk.Entry(row, textvariable=self.target, width=9).pack(side='left', padx=8)
        app._bind_text(ttk.Button(row, command=self.generate), 'geometry.generate').pack(side='left')
        app._bind_text(ttk.Button(row, command=lambda: app.notebook.select(app.pre_tab)),
                       'geometry.next_pre').pack(side='right')
        ttk.Label(mesh_tab, textvariable=self.mesh_description, wraplength=680).pack(fill='x', padx=12)
        app._bind_text(ttk.Label(mesh_tab, wraplength=680), 'geometry.scope').pack(fill='x', padx=12, pady=5)
        self.figure = Figure(figsize=(6, 4), dpi=90)
        self.ax = self.figure.add_subplot(111, projection='3d')
        self.canvas = FigureCanvasTkAgg(self.figure, master=mesh_tab)
        self.canvas.get_tk_widget().pack(fill='both', expand=True)
        self.refresh()
        self.draw()

    def refresh(self):
        tr = self.app._tr
        self.description.set(tr('geometry.summary', source=(tr('geometry.default') if
            self.geometry.source == 'cylinder' else self.geometry.source),
            vertices=len(self.geometry.vertices), faces=len(self.geometry.faces)))
        self.mesh_description.set(tr('geometry.unmeshed') if self.mesh is None else tr(
            'geometry.mesh_summary', faces=len(self.mesh.faces), area=f'{sum(self.mesh.areas):.6g}',
            closed=tr('geometry.closed' if self.mesh.closed else 'geometry.open')))

    def error(self, exc):
        key = str(exc)
        known = {'invalid_mesh', 'positive_dimensions', 'mesh_limit', 'triangles_only', 'supported_formats'}
        messagebox.showerror(self.app._tr('tab.model'), self.app._tr(
            'geometry.'+key if key in known else 'geometry.invalid_mesh'), parent=self.app.root)

    def use_cylinder(self):
        try:
            dimensions = (float(self.radius.get()), float(self.length.get()))
            geometry = cylinder(*dimensions)
        except ValueError as exc:
            self.error(exc); return
        self.geometry = geometry; self.mesh = None
        self.cylinder_dimensions = dimensions
        self.refresh(); self.draw()

    def import_model(self):
        path = filedialog.askopenfilename(parent=self.app.root,
            title=self.app._tr('geometry.import'), filetypes=[('STL / OBJ', '*.stl *.obj')])
        if not path:
            return
        try:
            geometry = load_surface(path, float(self.unit_scale.get()))
        except (ValueError, OSError, UnicodeError) as exc:
            self.error(exc); return
        self.geometry = geometry; self.mesh = None
        self.cylinder_dimensions = None
        self.refresh(); self.draw()

    def generate(self):
        try:
            target = float(self.target.get())
            geometry = (cylinder(*self.cylinder_dimensions, target_mm=target)
                        if self.cylinder_dimensions is not None else self.geometry)
            mesh = refine_surface(geometry, target)
        except ValueError as exc:
            self.error(exc); return
        self.mesh = mesh
        self.refresh(); self.draw()

    def draw(self):
        mesh = self.mesh if self.mesh is not None else self.geometry
        self.ax.clear()
        triangles = mesh.vertices[mesh.faces]
        self.ax.add_collection3d(Poly3DCollection(triangles, facecolor='#b5d3e8',
            edgecolor='#38566b', linewidth=.25, alpha=.85))
        lo, hi = mesh.vertices.min(axis=0), mesh.vertices.max(axis=0)
        center = (lo+hi)/2; radius = max(float(np.max(hi-lo))/2, 1e-6)
        self.ax.set_xlim(center[0]-radius, center[0]+radius)
        self.ax.set_ylim(center[1]-radius, center[1]+radius)
        self.ax.set_zlim(center[2]-radius, center[2]+radius)
        self.ax.set_box_aspect((1, 1, 1))
        self.ax.set_xlabel('X [mm]'); self.ax.set_ylabel('Y [mm]'); self.ax.set_zlabel('Z [mm]')
        self.canvas.draw_idle()
