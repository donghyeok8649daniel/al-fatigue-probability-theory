"""Actual section stress / local-PDE fields on their source specimen mesh."""
import tkinter as tk
from tkinter import ttk
import numpy as np
from matplotlib import colormaps
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from .face_picker import FacePicker
from .axial_specimen import validate_spatial_result
from .materials import ALUMINUM


def extrapolate_patches(local_probability, areas, correlation_area):
    p, ac = float(local_probability), float(correlation_area)
    areas = np.asarray(areas, dtype=float)
    if not np.isfinite(p) or not 0 <= p <= 1 or not np.isfinite(ac) or ac <= 0:
        raise ValueError('invalid_extrapolation')
    if not areas.size or np.any(~np.isfinite(areas)) or np.any(areas <= 0): raise ValueError('invalid_area')
    if p == 0: log_s = np.zeros_like(areas)
    elif p == 1: log_s = np.full_like(areas, -np.inf)
    else:
        # Log the positive hazard to avoid an overflowing area/count ratio.
        with np.errstate(over='ignore', under='ignore'):
            log_s = -np.exp(np.log(areas)-np.log(ac)+np.log(-np.log1p(-p)))
    return np.exp(log_s), -np.expm1(log_s), log_s


class MeshProbabilityView:
    def __init__(self, app, mesh):
        self.app, self.mesh = app, mesh
        self.window = tk.Toplevel(app.root)
        self.window.title(app._tr('map.title'))
        self.window.geometry('920x760')
        self.note = ttk.Label(self.window, text=app._tr('map.assumption'), wraplength=850)
        self.note.pack(fill='x', padx=12, pady=6)
        self.info = tk.StringVar(master=self.window)
        ttk.Label(self.window, textvariable=self.info, wraplength=850).pack(fill='x', padx=12)
        controls = ttk.Frame(self.window); controls.pack(fill='x', padx=12)
        self.field = tk.StringVar(master=self.window, value='initiation')
        self.text_widgets = [(self.note, 'map.assumption')]
        for code in ('stress', 'survival', 'initiation'):
            button = ttk.Radiobutton(controls, text=app._tr('map.'+code), value=code,
                            variable=self.field, command=self.update)
            button.pack(side='left')
            self.text_widgets.append((button,'map.'+code))
        button = ttk.Button(controls, text=app._tr('map.refresh'), command=self.update)
        button.pack(side='right')
        self.text_widgets.append((button, 'map.refresh'))
        self.solid_controls = ttk.Frame(self.window)
        for code in ('mises', 'xx', 'yy', 'zz', 'xy', 'xz', 'yz'):
            button = ttk.Radiobutton(self.solid_controls, text=app._tr('map.'+code), value=code,
                                    variable=self.field, command=self.update)
            button.pack(side='left')
            self.text_widgets.append((button, 'map.'+code))
        self.volume_controls = ttk.Frame(self.window)
        label = ttk.Label(self.volume_controls, text=app._tr('risk.volume'))
        label.pack(side='left')
        self.text_widgets.append((label, 'risk.volume'))
        ttk.Entry(self.volume_controls, textvariable=app.correlation_volume, width=12).pack(side='left', padx=4)
        self.risk_controls = ttk.Frame(self.window)
        for code in ('hazard_density', 'cell_risk', 'unit_risk'):
            button = ttk.Radiobutton(self.risk_controls, text=app._tr('map.'+code), value=code,
                                    variable=self.field, command=self.update)
            button.pack(side='left')
            self.text_widgets.append((button, 'map.'+code))
        self.slice_controls = ttk.Frame(self.window)
        label = ttk.Label(self.slice_controls, text=app._tr('solid.slice'))
        label.pack(side='left')
        self.text_widgets.append((label, 'solid.slice'))
        self.slice_fraction = tk.DoubleVar(master=self.window, value=1.)
        tk.Scale(self.slice_controls, from_=0, to=1, resolution=.05, orient='horizontal',
                 variable=self.slice_fraction, command=lambda _: self.update(), length=250).pack(side='left')
        self._solid_surface = None
        self._surface_key = None
        self._volume_trace = app.correlation_volume.trace_add('write', lambda *_: self.update())
        self.window.bind('<Destroy>', self._destroy, add='+')
        actions = ttk.Frame(self.window)
        actions.pack(fill='x', padx=12, pady=3)
        button = ttk.Button(actions, text=app._tr('solid.run'),
                            command=lambda: app._start_solve(solid=True))
        button.pack(side='left')
        self.text_widgets.append((button, 'solid.run'))
        button = ttk.Button(actions, text=app._tr('axial.run'),
                            command=lambda: app._start_solve(spatial=True))
        button.pack(side='left', padx=8)
        self.text_widgets.append((button, 'axial.run'))
        self.index = tk.IntVar(master=self.window, value=0)
        self.slider = tk.Scale(self.window, from_=0, to=0, orient='horizontal',
                               variable=self.index, command=lambda _: self.update(), showvalue=True)
        self.slider.pack(fill='x', padx=12)
        self.picker = FacePicker(self.window, lambda *_: None)
        # Reserve a fixed colorbar lane; never steal space on each update.
        self.picker.ax.set_position([.08, .10, .70, .80])
        self.picker.draw(mesh, [])
        self.colorbar = None
        self._result_id = None
        self.update()

    def _destroy(self, event):
        if event.widget is self.window:
            self.app.correlation_volume.trace_remove('write', self._volume_trace)

    def localize(self):
        self.window.title(self.app._tr('map.title'))
        for widget, key in self.text_widgets: widget.configure(text=self.app._tr(key))
        self.update()

    def unavailable(self, message):
        self.info.set(message)
        if self.colorbar is not None:
            self.colorbar.remove()
            self.colorbar = None
        self.picker.draw(self.mesh, [])

    def update(self):
        app, mesh = self.app, self.mesh
        if mesh is None:
            self.unavailable(app._tr('load.empty_selection')); return
        try: volume = f'{mesh.enclosed_volume_mm3:.8g}'
        except ValueError: volume = app._tr('map.volume_unknown')
        geometry = app._tr('map.geometry', area=f'{mesh.areas.sum():.8g}', volume=volume)
        if app.result is not None:
            geometry = app._result_material_text(app.result.get('material_id', ALUMINUM)) + '\n' + geometry
        solid = (app.result or {}).get('solid_specimen')
        if solid is not None:
            self.update_solid(solid, geometry)
            return
        self.solid_controls.pack_forget()
        self.volume_controls.pack_forget()
        self.risk_controls.pack_forget()
        self.slice_controls.pack_forget()
        self.window.title(app._tr('map.title'))
        self.note.configure(text=app._tr('map.assumption'))
        if self.field.get() not in ('stress', 'survival', 'initiation'):
            self.field.set('initiation')
        data = (app.result or {}).get('axial_specimen')
        if data is None:
            self.unavailable(geometry+'\n'+app._tr('map.need_result')); return
        try:
            validate_spatial_result(data, mesh)
            self.slider.configure(to=len(data['model_time'])-1)
            if app.result is not self._result_id:
                self._result_id = app.result
                self.index.set(len(data['model_time'])-1)
            index = max(0, min(self.index.get(), len(data['model_time'])-1))
        except (ValueError, KeyError) as exc:
            key = str(exc) if str(exc).startswith('axial.') else 'axial.incomplete'
            self.unavailable(geometry+'\n'+app._tr(key)); return
        key = {'stress': 'stress_mpa', 'survival': 'local_survival_probability',
               'initiation': 'local_initiation_probability'}[self.field.get()]
        history = np.asarray(data[key])
        values = history[index, np.asarray(data['face_sections'])]
        lower, upper = float(values.min()), float(values.max())
        # Fixed across the stored history: changing time does not stretch colors.
        vmin = float(history.min()) if self.field.get() == 'survival' else min(0., float(history.min()))
        vmax = 1. if self.field.get() == 'survival' else max(0., float(history.max()))
        norm = Normalize(vmin=vmin if vmin != vmax else 0., vmax=vmax if vmin != vmax else 1.)
        cmap = colormaps['viridis']
        self.picker.draw(mesh, [])
        self.picker.ax.collections[0].set_facecolor(cmap(norm(values)))
        if self.colorbar is None:
            color_axes = self.picker.figure.add_axes([.84, .20, .025, .60])
            self.colorbar = self.picker.figure.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=color_axes)
        else:
            self.colorbar.mappable.set_norm(norm)
            self.colorbar.mappable.set_cmap(cmap)
        self.colorbar.set_label(app._tr('map.'+self.field.get()))
        self.info.set(geometry+'\n'+app._tr('map.values',
            time=f"{data['model_time'][index]:.8g}",
            minimum=f'{lower:.8e}', maximum=f'{upper:.8e}',
            force=f"{data['axial_force_n'][index]:.8g}",
            area=f"{data['loaded_area_mm2']:.8g}",
            axis=np.array2string(np.asarray(data['axis']), precision=4),
            count=len(data['areas_mm2']),
            floor=f"{np.max(data['local_rare_event_floor']):.3e}"))
        self.picker.canvas.draw_idle()

    def update_solid(self, data, geometry):
        from .solid_mechanics import validate_solid_result, solid_field, visible_cell_surface
        from .volume_risk import volume_risk
        from types import SimpleNamespace
        app = self.app
        self.window.title(app._tr('solid.map_title'))
        self.note.configure(text=app._tr('solid.map_scope'))
        self.solid_controls.pack(fill='x', padx=12, before=self.slider)
        self.volume_controls.pack(fill='x', padx=12, before=self.slider)
        self.risk_controls.pack(fill='x', padx=12, before=self.slider)
        self.slice_controls.pack(fill='x', padx=12, before=self.slider)
        try:
            validate_solid_result(data, self.mesh)
        except (ValueError, KeyError) as exc:
            key = str(exc) if str(exc).startswith('axial.') else 'axial.incomplete'
            self.unavailable(geometry+'\n'+app._tr(key)); return
        self.slider.configure(to=len(data['model_time'])-1)
        if app.result is not self._result_id:
            self._result_id = app.result
            self.index.set(len(data['model_time'])-1)
        index = max(0, min(self.index.get(), len(data['model_time'])-1))
        risk_field = self.field.get() in ('hazard_density', 'cell_risk', 'unit_risk')
        if risk_field:
            try:
                vc = float(app.correlation_volume.get())
                density, cell_risk, unit_risk = volume_risk(solid_field(data, 'initiation', index),
                                                          data['cell_volumes_mm3'], vc)
                values = dict(hazard_density=density, cell_risk=cell_risk, unit_risk=unit_risk)[self.field.get()]
            except ValueError:
                self.unavailable(geometry+'\n'+app._tr('risk.need_volume')); return
        else:
            values = solid_field(data, self.field.get(), index)
        if self.field.get() in ('survival', 'initiation'):
            source = np.asarray(data['local_survival_probability' if self.field.get() == 'survival'
                                     else 'local_initiation_probability'])
            low = float(source.min()) if self.field.get() == 'survival' else 0.
            high = 1. if self.field.get() == 'survival' else float(source.max())
        elif risk_field:
            # Absorbed local P is nondecreasing; the final frame gives fixed
            # limits for the whole history at the current statistical V_c.
            final_risk = volume_risk(solid_field(data, 'initiation', -1), data['cell_volumes_mm3'], vc)
            final_values = final_risk[('hazard_density', 'cell_risk', 'unit_risk').index(self.field.get())]
            finite = final_values[np.isfinite(final_values)]
            low, high = 0., float(finite.max()) if finite.size else 1.
        else:
            finite = values[np.isfinite(values)]
            low = min(0., float(finite.min())) if finite.size else 0.
            high = max(0., float(finite.max())) if finite.size else 1.
        norm = Normalize(vmin=low if low != high else 0., vmax=high if low != high else 1.)
        cmap = colormaps['viridis'].with_extremes(over='#d01c8b')
        surface_key = (id(data), self.slice_fraction.get())
        if self._surface_key != surface_key:
            faces, owners = visible_cell_surface(data, self.slice_fraction.get())
            self._solid_surface = SimpleNamespace(vertices=data['nodes_mm'], faces=faces)
            self._surface_cells = owners
            self._surface_key = surface_key
        # A display cut changes only visible cells. Preserve the user's camera.
        limits = self.picker.ax.get_xlim3d(), self.picker.ax.get_ylim3d(), self.picker.ax.get_zlim3d()
        self.picker.draw(self._solid_surface, [])
        self.picker.ax.set(xlim=limits[0], ylim=limits[1], zlim=limits[2])
        self.picker.ax.collections[0].set_facecolor(cmap(norm(values[self._surface_cells])))
        from .face_picker import EDGE_DISPLAY_FACE_LIMIT
        edges = len(self._surface_cells) <= EDGE_DISPLAY_FACE_LIMIT
        self.picker.ax.collections[0].set_edgecolor((.15, .20, .25, .20) if edges else 'none')
        self.picker.ax.collections[0].set_linewidth(.15 if edges else 0.)
        if self.colorbar is None:
            color_axes = self.picker.figure.add_axes([.84, .20, .025, .60])
            self.colorbar = self.picker.figure.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=color_axes)
        else:
            self.colorbar.mappable.set_norm(norm)
            self.colorbar.mappable.set_cmap(cmap)
        self.colorbar.set_label(app._tr('map.'+self.field.get()))
        self.info.set(geometry+'\n'+app._tr('solid.map_values', time=f"{data['model_time'][index]:.8g}",
            low=f'{values.min():.6e}', high=f'{values.max():.6e}', cells=len(data['tetrahedra']),
            count=len(data['sample_cells']), axis=np.array2string(data['probability_axis'], precision=4),
            balance=np.array2string(data['force_torque_balance'][index], precision=3),
            residual=f"{np.max(data['relative_linear_residual']):.3e}",
            mismatch=f"{np.max(data['sampling_stress_mismatch_mpa']):.6g}",
            floor=f"{np.max(data['local_rare_event_floor']):.3e}"))
        if risk_field:
            self.info.set(self.info.get()+'\n'+app._tr('risk.assumption', volume=f'{vc:.8g}',
                                                       infinite=int(np.isinf(density).sum())))
        self.picker.canvas.draw_idle()
