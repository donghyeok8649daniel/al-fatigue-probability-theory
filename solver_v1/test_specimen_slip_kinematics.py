import json
from pathlib import Path

import numpy as np
import pytest

from .finite_source_reference import line_energy_coefficient, solve_pinned_graph
from .nonlocal_interface_elasticity import cubic_elastic_tensor
from .specimen_slip_kinematics import (
    slip_distortion, symmetric_slip_strain, axial_slip_strain, transport_history,
    uniform_stress_work, pinned_swept_area, pinned_population_strain_budget,
)


def geometry():
    return dict(burgers_vectors_m=[[3e-10, 0, 0]], plane_normals=[[0, 1, 0]],
                specimen_volume_m3=1e-12)


def test_signed_tensor_strain_and_engineering_shear():
    g = geometry()
    areas = np.array([[4e-13], [-4e-13]])
    beta = slip_distortion(areas, **g)
    eps = symmetric_slip_strain(areas, **g)
    assert beta[0, 0, 1] == pytest.approx(1.2e-10)
    np.testing.assert_allclose(eps[1], -eps[0], atol=0)
    assert eps[0, 0, 1] == pytest.approx(6e-11)
    assert np.trace(eps[0]) == 0
    axial = axial_slip_strain(eps, [1/np.sqrt(2), 1/np.sqrt(2), 0])
    np.testing.assert_allclose(axial, [6e-11, -6e-11], atol=1e-25)


def test_burgers_and_plane_sign_gauges():
    g = geometry(); area = [2e-13]
    original = symmetric_slip_strain(area, **g)
    reverse_b = {**g, 'burgers_vectors_m': [[-3e-10, 0, 0]]}
    reverse_n = {**g, 'plane_normals': [[0, -1, 0]]}
    np.testing.assert_allclose(symmetric_slip_strain([-area[0]], **reverse_b), original)
    np.testing.assert_allclose(symmetric_slip_strain([-area[0]], **reverse_n), original)


def test_multiple_slip_systems_tensor_superposition_and_rotation():
    g = dict(burgers_vectors_m=[[3e-10, 0, 0], [0, 0, 3e-10]],
             plane_normals=[[0, 1, 0], [0, 1, 0]], specimen_volume_m3=2e-12)
    areas = np.array([4e-13, -7e-13])
    eps = symmetric_slip_strain(areas, **g)
    angle = .37
    Q = np.array([[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
    rotated = {**g, 'burgers_vectors_m': np.asarray(g['burgers_vectors_m'])@Q.T,
               'plane_normals': np.asarray(g['plane_normals'])@Q.T}
    np.testing.assert_allclose(symmetric_slip_strain(areas, **rotated), Q@eps@Q.T, atol=2e-26)


def test_existing_slip_is_not_new_flow_and_recoverable_bow_has_gross_activity():
    initial = 2e-13
    areas = np.array([[initial], [initial+3e-13], [initial], [initial-3e-13], [initial]])
    result = transport_history(areas, **geometry())
    np.testing.assert_array_equal(result['signed_slip_strain_tensor'][[0, -1]], 0)
    assert result['forward_area_m2'][-1, 0] == pytest.approx(6e-13)
    assert result['backward_area_m2'][-1, 0] == pytest.approx(6e-13)
    assert result['gross_area_m2'][-1, 0] == pytest.approx(1.2e-12)
    assert np.max(abs(result['area_balance_residual_m2'])) < 1e-27
    assert not result['residual_plasticity_certified']


def test_kinematic_retained_area_not_automatically_certified_material_plasticity():
    result = transport_history([[0], [4e-13], [3e-13]], **geometry())
    assert result['signed_slip_strain_tensor'][-1, 0, 1] > 0
    assert not result['residual_plasticity_certified']


def test_uniform_mixed_stress_virtual_work_two_independent_forms():
    stress = np.array([[12, 7, -2], [7, -5, 3], [-2, 3, 8]])*1e6
    g = dict(burgers_vectors_m=[[3e-10, 0, 0], [0, 0, 3e-10]],
             plane_normals=[[0, 1, 0], [0, 1, 0]], specimen_volume_m3=2e-12)
    result = uniform_stress_work(stress, [[4e-13, -2e-13], [-1e-13, 7e-13]], **g)
    exact = np.array([7e6*3e-10*4e-13 - 3e6*3e-10*2e-13,
                      -7e6*3e-10*1e-13 + 3e6*3e-10*7e-13])
    np.testing.assert_allclose(result['tensor_work_J'], exact, rtol=3e-15, atol=1e-30)
    assert np.max(abs(result['residual_J'])) < 1e-30


def test_volume_tiling_changes_neither_strain_nor_work_density():
    eps = symmetric_slip_strain([4e-13], **geometry())
    g = {**geometry(), 'specimen_volume_m3': 7e-12}
    np.testing.assert_allclose(symmetric_slip_strain([7*4e-13], **g), eps)


def test_sampled_gross_area_requires_time_resolution():
    coarse = transport_history([[0], [0]], **geometry())
    resolved = transport_history([[0], [1e-13], [0]], **geometry())
    assert coarse['gross_area_m2'][-1, 0] == 0
    assert resolved['gross_area_m2'][-1, 0] == pytest.approx(2e-13)
    np.testing.assert_array_equal(coarse['signed_slip_strain_tensor'][-1],
                                  resolved['signed_slip_strain_tensor'][-1])


def test_nontangent_motion_and_invalid_stress_are_not_slip():
    with pytest.raises(ValueError):
        symmetric_slip_strain([1e-13], **{**geometry(), 'burgers_vectors_m': [[0, 3e-10, 0]]})
    with pytest.raises(ValueError):
        uniform_stress_work([[0, 1, 0], [0, 0, 0], [0, 0, 0]], [1e-13], **geometry())
    with pytest.raises(ValueError):
        axial_slip_strain(np.zeros((3, 3)), [2, 0, 0])


@pytest.fixture(scope='module')
def constant_line():
    # nu=0 gives an exactly orientation-independent line coefficient.
    return line_energy_coefficient(cubic_elastic_tensor(50e9, 0, 25e9), [3e-10, 0, 0])


def test_analytic_circle_area_and_fold_limit(constant_line):
    for angle in (.2, 1., np.pi/2):
        L = 1e-6
        result = pinned_swept_area(constant_line, angle, span_m=L, outer_log_ratio=7.)
        radius = L/(2*np.sin(angle))
        expected = radius**2*(angle-np.sin(angle)*np.cos(angle))
        assert result['area_m2'] == pytest.approx(expected, rel=2e-12)
        assert result['estimated_quadrature_error_m2'] < 1e-9*expected


def test_area_matches_independent_discrete_line_and_converges(constant_line):
    exact = pinned_swept_area(constant_line, 1., span_m=1e-6, outer_log_ratio=7.)
    errors = []
    for n in (32, 64, 128):
        result = solve_pinned_graph(constant_line, span_m=1e-6, outer_log_ratio=7.,
                                   shear_Pa=exact['applied_shear_Pa'], segments=n)
        assert result['converged']
        errors.append(abs(result['swept_area_m2']-exact['area_m2']))
    assert errors[1] < .3*errors[0] and errors[2] < .3*errors[1]


def test_population_budget_is_not_a_calibrated_density_or_yield(constant_line):
    area = pinned_swept_area(constant_line, np.pi/2, span_m=1e-6, outer_log_ratio=7.)['area_m2']
    budget = pinned_population_strain_budget(area_per_source_m2=area, span_m=1e-6,
        burgers_m=3e-10, schmid_factor=.5, axial_criterion=.002)
    rho = budget['required_initial_pinned_line_density_m2']
    assert .5*3e-10*rho*area/1e-6 == pytest.approx(.002)
    assert budget['required_overlap_parameter'] > 1
    assert not budget['population_calibrated'] and budget['experimental_yield_MPa'] is None


@pytest.mark.parametrize('updates', [
    {'plane_normals': [[0, 2, 0]]}, {'plane_normals': [[1, 0, 0]]},
    {'burgers_vectors_m': [[0, 0, 0]]}, {'specimen_volume_m3': 0},
    {'specimen_volume_m3': np.inf}, {'burgers_vectors_m': [[np.nan, 0, 0]]},
])
def test_invalid_geometry_is_not_silently_normalized(updates):
    with pytest.raises(ValueError):
        symmetric_slip_strain([1e-13], **{**geometry(), **updates})


def test_no_statistical_area_or_fitted_stress_scale_in_api():
    with pytest.raises(TypeError):
        symmetric_slip_strain([1e-13], A_c=1e-12, **geometry())


def test_existing_production_kinetic_status_unchanged():
    path = Path(__file__).parent/'data/aluminum_kinetic_calibration.json'
    record = json.loads(path.read_text(encoding='utf-8'))
    assert record['calibrated'] is False
    assert record['t0_seconds'] is None
