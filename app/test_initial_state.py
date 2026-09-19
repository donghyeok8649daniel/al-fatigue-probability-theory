"""Initial ensemble choice is explicit, including an absent loaded basin."""
from dataclasses import replace
import numpy as np
import pytest
from app.solver_adapter import UIAnalysisConfig, run_ui_analysis, physical_load_conversion
from solver_v1.kinetic_calibration_workflow import build_time_basis_model
from solver_v1.probability_pde_2d import (
    Grid2D, InitialGibbsBasinError, initial_gibbs_density, cyclic_load_from_sigma_over_E,
)


@pytest.fixture(scope='module')
def model():
    return build_time_basis_model('two_row_lj_reference')


def test_1200_1500_nominal_initializes_but_neck_reports_actual_unbound_stress(model):
    config = UIAnalysisConfig(stress_mean_mpa=1200, stress_amplitude_mpa=1500,
                              model_frequency=1000., cycles=.002)
    nominal = run_ui_analysis(config, _prepared_model=model)
    assert nominal['initialization'] == 'loaded_gibbs'
    assert nominal['preload_stress_mpa'] == 1200
    peak = 5213.923003084407
    with pytest.raises(ValueError) as caught:
        run_ui_analysis(config, axial_stress_function=lambda _: peak, _prepared_model=model)
    payload = caught.value.ui_error_data
    assert payload['key'] == 'initial.unbound'
    assert float(payload['values']['stress']) == pytest.approx(peak)
    assert float(payload['values']['limit']) == pytest.approx(4199.2211409)
    assert config.initialization == 'loaded_gibbs'  # no automatic change of preparation


def test_explicit_zero_load_preparation_records_instant_opening_and_undefined_strain(model):
    config = UIAnalysisConfig(stress_mean_mpa=1200, stress_amplitude_mpa=1500,
        model_frequency=1000., cycles=.002, initialization='zero_load_gibbs')
    conversion = physical_load_conversion(config)
    assert conversion['preload_force'] == 0.
    assert conversion['reduced_stress_initial'] == 1200/config.young_mpa
    result = run_ui_analysis(config, axial_stress_function=lambda _: 5213.923003084407,
                             _prepared_model=model)
    assert result['preload_force'] == 0.
    assert result['initial_force'] > 6.
    assert result['initial_absorbed_mass'][0] == pytest.approx(1.)
    assert result['local_initiation_probability'][0] == pytest.approx(1.)
    assert result['intact_probability_mass'][0] == 0.
    assert abs(result['mass_balance_residual'][0]) < 1e-14
    assert np.isnan(result['strain'][0])  # survivor strain is undefined after extinction


def test_missing_grid_cells_below_opening_limit_is_not_physical_opening(model):
    grid = Grid2D(a=np.linspace(20, 21, 5), s=np.linspace(-.1,.1,5), da=.25, ds=.05)
    with pytest.raises(InitialGibbsBasinError) as caught:
        initial_gibbs_density(model, grid, preload_force=1.)
    assert caught.value.reason == 'unresolved'


def test_reusing_prepared_model_preserves_independent_pde_histories(model):
    config = UIAnalysisConfig(model_frequency=1000., cycles=.002)
    shared = run_ui_analysis(config, _prepared_model=model)
    fresh = run_ui_analysis(config)
    for key in ('model_time','local_initiation_probability','strain','mass_balance_residual','final_density'):
        if key in fresh:
            np.testing.assert_array_equal(shared[key], fresh[key])
    later = run_ui_analysis(replace(config, stress_mean_mpa=150.), _prepared_model=model)
    assert later['initial_stress_mpa'] == 150.
    assert shared['initial_stress_mpa'] == 50.


def test_force_bridge_reuses_the_same_fixed_load_scale(model, monkeypatch):
    scale = model.sigma_over_E_force_scale()
    load = cyclic_load_from_sigma_over_E(model, sigma_over_E_min=-.01, sigma_over_E_max=.02,
        value_function=lambda t: .003+np.sin(t)*.01)
    monkeypatch.setattr(model, 'force_from_sigma_over_E', lambda _: pytest.fail('recomputed fixed tangent'))
    for time in (0., .2, 1.):
        assert load.value(time) == scale*(.003+np.sin(time)*.01)
