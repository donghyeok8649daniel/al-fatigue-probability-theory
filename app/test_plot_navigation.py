"""Headless checks of actual UI handlers; never rerun the PDE."""
from types import SimpleNamespace

import numpy as np
import pytest
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.backend_bases import MouseEvent

from app.desktop_ui import DesktopApp
from app.i18n import Localizer, tr
from app.solver_adapter import UIAnalysisConfig


def app_stub():
    app = DesktopApp.__new__(DesktopApp)
    app.localizer = Localizer()
    app.last_config = UIAnalysisConfig()
    app.figure = Figure()
    app.canvas = FigureCanvasAgg(app.figure)
    app.ax = app.figure.add_subplot()
    app.result = {"local_initiation_probability": np.array([0., 1.e-10])}
    app.ax.plot([0., 1.], app.result["local_initiation_probability"])
    app.canvas.draw()
    app._pan_origin = None
    app._view_limits = {}
    app.field = SimpleNamespace(get=lambda: "local_initiation_probability")
    return app


def test_pan_uses_fixed_pixel_transform_without_feedback_or_data_mutation():
    app = app_stub()
    original = app.result["local_initiation_probability"].copy()
    xlim, ylim = app.ax.get_xlim(), app.ax.get_ylim()
    x, y = app.ax.transAxes.transform((.5, .5))
    app._on_press(MouseEvent("button_press_event", app.canvas, x, y, button=1))
    for displacement in (15, 30, 45):
        app._on_motion(MouseEvent("motion_notify_event", app.canvas, x+displacement, y))
        expected = np.array(xlim) - displacement/app.ax.bbox.width*(xlim[1]-xlim[0])
        np.testing.assert_allclose(app.ax.get_xlim(), expected)
        np.testing.assert_allclose(app.ax.get_ylim(), ylim)
    np.testing.assert_array_equal(app.result["local_initiation_probability"], original)
    app._on_release(None)
    assert app._pan_origin is None


@pytest.mark.parametrize("mode", ["pan/zoom", "zoom rect"])
def test_toolbar_gesture_is_not_handled_twice(mode):
    app = app_stub()
    app.canvas.toolbar = SimpleNamespace(mode=mode)
    x, y = app.ax.transAxes.transform((.5, .5))
    app._on_press(MouseEvent("button_press_event", app.canvas, x, y, button=1))
    assert app._pan_origin is None


def test_reset_clears_only_current_view_and_requests_autoscale():
    app = app_stub()
    app._view_limits = {"local_initiation_probability": ((2,3),(-1,-.5)), "strain": ((0,1),(0,.1))}
    calls = []
    app._plot = lambda **kwargs: calls.append(kwargs)
    app._reset_view()
    assert set(app._view_limits) == {"strain"}
    assert calls == [{"force_auto": True}]


def test_actual_probability_plot_preserves_absorbed_array_and_resets_bad_view():
    app = app_stub()
    probability = np.array([0., 1.e-20, 1.e-20, 2.e-20])
    app.result = {
        "model_time": np.arange(4.),
        "local_initiation_probability": probability,
        "time_basis": "model",
    }
    app.localizer = Localizer('en')
    app._last_field = None
    app._view_limits = {"local_initiation_probability": ((1.,2.),(-1.e-8,-.5e-8))}
    app._plot()
    np.testing.assert_array_equal(app.ax.lines[0].get_ydata(), probability)
    np.testing.assert_allclose(app.ax.get_ylim(), (-1.e-8,-.5e-8))
    app._reset_view()
    low, high = app.ax.get_ylim()
    assert low <= probability.min() <= probability.max() <= high
    np.testing.assert_array_equal(app.ax.lines[0].get_ydata(), probability)
    assert app.result["local_initiation_probability"] is probability
    assert app.figure._suptitle.get_text() == 'Result material: Aluminum (Al)'


@pytest.mark.parametrize('language', ['ko', 'en'])
def test_live_specimen_plot_and_labels_are_updated_together(language):
    app=app_stub(); app.result=None
    app.localizer=Localizer(language)
    app.time_basis_code='model'; app._last_field=None
    app.field=SimpleNamespace(get=lambda:'specimen_probability_extrapolation')
    app.entries={
        'correlation_area_mm2':SimpleNamespace(get=lambda:'1'),
        'stressed_area_mm2':SimpleNamespace(get=lambda:'100'),
    }
    app.live_records=[dict(model_time=0.,local_initiation_probability=0.),
                      dict(model_time=1.,local_initiation_probability=1e-16)]
    labels={}
    for name in ('specimen_N_eff','local_floor_display','plastic_floor_display',
                 'configurational_barrier_display','opening_barrier_display',
                 'plasticity_detail_display','local_probability_display',
                 'specimen_extrapolation_display','specimen_certified_display',
                 'probability_status_display'):
        setattr(app,name,SimpleNamespace(set=lambda value,name=name:labels.__setitem__(name,value)))
    app._update_specimen_probability()
    expected=-np.expm1(100*np.log1p(-np.array([0.,1e-16])))
    np.testing.assert_array_equal(app.ax.lines[0].get_ydata(),expected)
    assert '1e-14' in labels['specimen_extrapolation_display']
    assert any('1e-14' in text.get_text() for text in app.ax.texts)
    assert '100' in labels['specimen_N_eff']
    assert '—' in labels['specimen_certified_display']
    assert tr('status.specimen_live',language) in labels['probability_status_display']
    assert app.result is None
    assert app.figure._suptitle.get_text() == app._result_material_text(app.last_config.material_id)
