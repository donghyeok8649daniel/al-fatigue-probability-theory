from types import SimpleNamespace
import numpy as np
import pytest
from app.desktop_ui import DesktopApp


def stub():
    app = DesktopApp.__new__(DesktopApp)
    app.result = dict(local_initiation_probability=np.array([0.,1e-16]),
                      local_rare_event_floor=1e-15, probability_resolution_certified=False)
    app._refresh_specimen_labels = lambda: None
    app._plot = lambda: None
    app.entries = {}
    return app


def areas(app, correlation, stressed):
    app.entries = {"correlation_area_mm2": SimpleNamespace(get=lambda: correlation),
                   "stressed_area_mm2": SimpleNamespace(get=lambda: stressed)}


@pytest.mark.parametrize('bad', ['', 'abc', '0', '-1', 'nan', 'inf'])
def test_invalid_area_clears_stale_results_then_recovers_without_pde(bad):
    app = stub()
    local = app.result['local_initiation_probability']
    areas(app, '1', '100')
    app._update_specimen_probability()
    assert app.result['N_eff'] == 100
    assert np.isfinite(app.result['specimen_probability_extrapolation']).all()
    assert np.isnan(app.result['specimen_initiation_probability']).all()
    areas(app, bad, '100')
    app._update_specimen_probability()
    assert 'N_eff' not in app.result
    assert not any(k.startswith('specimen_') for k in app.result)
    assert app._specimen_status_key == 'status.area_required'
    areas(app, '2', '100')
    app._update_specimen_probability()
    assert app.result['N_eff'] == 50
    assert app.result['local_initiation_probability'] is local
    expected = -np.expm1(50*np.log1p(-local[-1]))
    assert app.result['specimen_probability_extrapolation'][-1] == expected


@pytest.mark.parametrize('bad', ['-1', 'nan', 'inf', ''])
def test_invalid_stressed_area_is_handled(bad):
    app = stub()
    areas(app,'1',bad)
    app._update_specimen_probability()
    assert app._specimen_status_key == 'status.area_required'


def test_zero_stressed_area_is_valid_mathematical_zero():
    app=stub()
    areas(app,'1','0')
    app._update_specimen_probability()
    np.testing.assert_array_equal(app.result['specimen_probability_extrapolation'],[0,0])
