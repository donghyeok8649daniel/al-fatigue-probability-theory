"""No tests promote a source line clock to the probability PDE."""
import json

import numpy as np
import pytest

from .dislocation_line_kinetics import PinnedLineKinetics, fit_velocity_response
from .run_line_kinetic_benchmark import REFERENCE, velocity_points
from .fatigue_validation_reference import normalized_fatigue_records, fatigue_comparison_gate
from .physical_time import load_time_calibration, model_frequency_to_hz


def test_velocity_fit_is_deterministic_and_heldout_does_not_enter_calibration():
    stress = np.arange(1, 9)*1e5
    v = 1.3e-5*stress
    train = np.array([True]*5+[False]*3)
    first = fit_velocity_response(stress, v, train)
    v[~train] *= 3
    second = fit_velocity_response(stress, v, train)
    assert first['velocity_per_stress_m_per_Pa_s'] == second['velocity_per_stress_m_per_Pa_s']
    assert first['velocity_per_stress_m_per_Pa_s'] == pytest.approx(1.3e-5)
    assert second['heldout_rmse_m_s'] > 0
    assert not second['production_clock_calibrated']


def test_digitized_source_units_and_temperature_are_not_guessed():
    data = json.loads(REFERENCE.read_bytes())
    values, uncertainty = velocity_points(data)
    assert data['temperature_K'] == 23+273.15
    assert data['stress_unit_to_Pa'] == 1e5
    assert data['velocity_unit_to_m_s'] == .01
    assert .39e6 < values[:, 0].min() < .40e6
    assert 1.38e6 < values[:, 0].max() < 1.39e6
    assert np.all(uncertainty > 0)
    assert not data['full_figure_population_recovered']


def test_line_dimensions_and_source_span_scaling():
    model = PinnedLineKinetics(2e-5, 1e-9, 2.86e-10, 1e-6)
    assert model.slowest_seconds == pytest.approx(2e-8/np.pi**2)
    longer = PinnedLineKinetics(2e-5, 1e-9, 2.86e-10, 2e-6)
    assert longer.slowest_seconds/model.slowest_seconds == pytest.approx(4)
    assert longer.static_mean_bow_m(1e5)/model.static_mean_bow_m(1e5) == pytest.approx(4)
    with pytest.raises(ValueError):
        PinnedLineKinetics(-1., 1., 1., 1.)


def test_infinite_transfer_tail_and_low_high_frequency_limits():
    model = PinnedLineKinetics(1., 1., 1., np.pi)
    omega = np.array([0., .1, 1., 10., 100.])
    data = model.frequency_response(omega/(2*np.pi), tail_tolerance=1e-8)
    reference = model.frequency_response(omega/(2*np.pi), tail_tolerance=1e-12)
    assert abs(data['transfer'][0]-1) < data['absolute_tail_bound']
    assert np.max(abs(data['transfer']-reference['transfer'])) < data['absolute_tail_bound']
    assert np.all(np.diff(abs(data['transfer'])) < 0)
    assert np.all(np.angle(data['transfer'][1:]) < 0)


def test_independent_spatial_refinement_and_energy_dissipation():
    model = PinnedLineKinetics(2e-5, 1e-9, 2.86e-10, 1e-6)
    frequency = 1/(2*np.pi*model.slowest_seconds)
    exact = model.frequency_response(frequency)['transfer']
    errors = []
    for n in (16, 32, 64, 128):
        data = model.discrete_harmonic(frequency, 1e4, segments=n)
        errors.append(abs(data['transfer']-exact))
        assert data['work_per_cycle_J'] > 0
        assert abs(data['energy_balance_residual_J']) < 2e-11*data['work_per_cycle_J']
    assert np.all(np.array(errors[:-1])/errors[1:] > 3.9)


def test_time_refinement_zero_load_hold_and_no_permanent_slip():
    model = PinnedLineKinetics(1., 1., 1., np.pi)
    residues = []
    for step in (1/16, 1/32, 1/64):
        t = np.arange(0, 36+step/2, step)
        stress = np.where(t <= 16, 1., 0.)
        result = model.integrate_backward_euler(t, stress, segments=64)
        index = int(16/step)
        residues.append(result['mean_bow_m'][-1]/result['mean_bow_m'][index])
        assert not result['permanent_slip_generated']
        assert result['opening_probability'] is None
        assert np.all(np.diff(result['mean_bow_m'][index:]) <= 0)
    assert residues[-1] < residues[0] < 1e-8
    assert residues[-1] > 0  # finite relaxation tail is not clipped to zero


def test_fatigue_stress_range_is_not_amplitude_or_initiation_life():
    rows = normalized_fatigue_records()
    stress = [r for r in rows if r['control'] == 'axial_stress']
    assert [r['amplitude'] for r in stress] == [25., 31.]
    assert [r['laboratory_failure_time_min_seconds'] for r in stress] == [172000., 53000.]
    assert all(r['initiation_cycles'] is None for r in rows)
    assert all(r['local_absorbed_probability'] is None for r in rows)
    assert len(rows) == 7


def test_fatigue_endpoint_and_physical_clock_are_separate_requirements():
    gate = fatigue_comparison_gate(prediction_endpoint='local_atomic_opening_absorption',
        physical_clock_calibrated=False, specimen_mapping_validated=True,
        control_protocol_matched=True, microstructure_matched=True,
        probability_resolution_certified=True)
    assert not gate['comparable'] and len(gate['reasons']) == 2


def test_line_benchmark_cannot_enable_production_physical_hz():
    from pathlib import Path
    calibration = load_time_calibration(Path(__file__).with_name('data')/'aluminum_kinetic_calibration.json')
    assert not calibration.calibrated
    with pytest.raises(ValueError):
        model_frequency_to_hz(1., calibration)
