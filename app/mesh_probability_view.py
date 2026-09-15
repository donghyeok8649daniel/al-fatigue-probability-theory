"""Explicitly hypothetical surface-patch extrapolation of one local result."""
import tkinter as tk
from tkinter import ttk
import numpy as np
from matplotlib import colormaps
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from .face_picker import FacePicker


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
        self.window.geometry('920x720')
        self.note = ttk.Label(self.window, text=app._tr('map.assumption'), wraplength=850)
        self.note.pack(fill='x', padx=12, pady=6)
        self.info = tk.StringVar(master=self.window)
        ttk.Label(self.window, textvariable=self.info, wraplength=850).pack(fill='x', padx=12)
        controls = ttk.Frame(self.window); controls.pack(fill='x', padx=12)
        self.field = tk.StringVar(master=self.window, value='survival')
        self.text_widgets = [(self.note, 'map.assumption')]
        for code in ('survival', 'initiation'):
            button = ttk.Radiobutton(controls, text=app._tr('map.'+code), value=code,
                            variable=self.field, command=self.update)
            button.pack(side='left')
            self.text_widgets.append((button,'map.'+code))
        button = ttk.Button(controls, text=app._tr('map.refresh'), command=self.update)
        button.pack(side='right')
        self.text_widgets.append((button, 'map.refresh'))
        self.index = tk.IntVar(master=self.window, value=0)
        self.slider = tk.Scale(self.window, from_=0, to=0, orient='horizontal',
                               variable=self.index, command=lambda _: self.update(), showvalue=True)
        self.slider.pack(fill='x', padx=12)
        self.picker = FacePicker(self.window, lambda *_: None)
        self.picker.draw(mesh, [])
        self.colorbar = None
        self._result_id = None
        self.update()

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
        data = app._plot_data()
        if data is None or 'local_initiation_probability' not in data:
            self.unavailable(geometry+'\n'+app._tr('map.need_result')); return
        try:
            ac = float(app.entries['correlation_area_mm2'].get())
            probabilities = np.asarray(data['local_initiation_probability'])
            self.slider.configure(to=len(probabilities)-1)
            if app.result is not self._result_id:
                self._result_id = app.result
                self.index.set(len(probabilities)-1)
            index = min(self.index.get(), len(probabilities)-1)
            survival, initiation, log_s = extrapolate_patches(probabilities[index], mesh.areas, ac)
        except (ValueError, KeyError):
            self.unavailable(geometry+'\n'+app._tr('map.need_area')); return
        values = survival if self.field.get() == 'survival' else initiation
        lower, upper = float(values.min()), float(values.max())
        norm = Normalize(vmin=lower if lower != upper else 0., vmax=upper if lower != upper else 1.)
        cmap = colormaps['viridis']
        self.picker.draw(mesh, [])
        self.picker.ax.collections[0].set_facecolor(cmap(norm(values)))
        if self.colorbar is not None: self.colorbar.remove()
        self.colorbar = self.picker.figure.colorbar(ScalarMappable(norm=norm,cmap=cmap), ax=self.picker.ax, shrink=.65)
        self.colorbar.set_label(app._tr('map.'+self.field.get()))
        total = -np.expm1(log_s.sum())
        self.info.set(geometry+'\n'+app._tr('map.values', p=f'{probabilities[index]:.8e}',
            minimum=f'{lower:.8e}', maximum=f'{upper:.8e}', total=f'{total:.8e}'))
        self.picker.canvas.draw_idle()
