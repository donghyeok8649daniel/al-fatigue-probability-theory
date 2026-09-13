"""Research channel extraction must preserve energy/derivative consistency."""
import numpy as np
import pytest

from .coordination_screening import CoordinationScreenedInterface, CoordinationScreenedCache
from .run_radial_channels_v32 import channel_shape, conditional_svd
from .vector_material_calibration import MaterialObservation
from .vector_interface_reference import VectorInterfaceEvaluation


SHAPE = (4.9, 10., 11., 1000., 8.4, 1.)


@pytest.mark.parametrize('kind,column', [('screened_rank1', 6), ('rank3', 5), ('quadratic_rank2', 7),
                                        ('density_linear', 3), ('density_quadratic', 4)])
def test_extracted_column_is_one_analytic_energy(kind, column):
    shape, index = channel_shape(SHAPE, kind, 5.)
    assert index == column
    coefficients = np.ones(10)
    model = CoordinationScreenedInterface(shape, coefficients, law='power')
    component = {6: 'angular_1', 5: 'angular_3', 7: 'angular_2', 3: 'embedding', 4: 'embedding'}[index]
    zero = coefficients.copy(); zero[index] = 0.
    without = CoordinationScreenedInterface(shape, zero, law='power') if index in (3, 4) else None
    def evaluate(state):
        result = model.evaluate(state)
        jet = result.components[component]
        if index in (3, 4):
            jet = jet-without.evaluate(state).components[component]
        return VectorInterfaceEvaluation.from_jet(jet, {component: jet}, result.diagnostics)
    q = np.array([.879, .267, .113])
    value = evaluate(q)
    jet = np.r_[value.energy, value.gradient, value.hessian[0], value.hessian[1, 1:], value.hessian[2, 2]]
    obs = [MaterialObservation(f'jet{i}', 0., 1., 'test', 'heldout', tuple(q), tuple(w))
           for i, w in enumerate(np.eye(10))]
    matrix = CoordinationScreenedCache(obs, law='power').matrix(shape)
    np.testing.assert_allclose(matrix[:, index], jet, atol=2e-9, rtol=2e-10)
    step = 2e-5
    high = [evaluate(q+step*e) for e in np.eye(3)]
    low = [evaluate(q-step*e) for e in np.eye(3)]
    gradient = [(h.energy-l.energy)/(2*step) for h, l in zip(high, low)]
    hessian = np.column_stack([(h.gradient-l.gradient)/(2*step) for h, l in zip(high, low)])
    np.testing.assert_allclose(gradient, value.gradient, atol=1e-7, rtol=2e-6)
    np.testing.assert_allclose(hessian, value.hessian, atol=2e-6, rtol=3e-6)
    np.testing.assert_array_equal(value.hessian, value.hessian.T)
    repeat = evaluate(q+[0., 1., 0.])
    np.testing.assert_allclose(repeat.hessian, value.hessian, atol=1e-9)


def test_shape_selection_does_not_mutate_parent_and_rejects_invalid():
    original = np.array(SHAPE)
    shape, column = channel_shape(original, 'quadratic_rank2', 3.)
    np.testing.assert_array_equal(original, SHAPE)
    assert shape[3] == 0 and shape[2] == 3 and column == 7
    for kind, decay in [('unknown', 3), ('rank3', 0), ('rank3', np.nan)]:
        with pytest.raises(ValueError):
            channel_shape(original, kind, decay)


def test_duplicate_channel_cannot_be_called_independent_information():
    obs = [MaterialObservation('anchor', 1., 1., 'test', 'exact'),
           MaterialObservation('fit1', 1., 1., 'test', 'fit'),
           MaterialObservation('fit2', 1., 1., 'test', 'fit')]
    M = np.array([[1., 0.], [0., 1.], [1., 2.]])
    base = conditional_svd(M, obs, [1, 2])
    duplicate = conditional_svd(np.column_stack([M, M[:, 1]]), obs, [1, 2])
    assert base['conditional_rank'] == duplicate['conditional_rank'] == 1
    assert duplicate['null_dimension'] == 2
    assert duplicate['singular_values'][-1] < 1e-14


def test_study_variants_cannot_be_silently_combined(tmp_path):
    from .run_radial_channels_v32 import run
    with pytest.raises(ValueError, match='separate'):
        run(tmp_path/'missing', tmp_path/'out', powers_only=True, density_only=True)
    for options in (dict(rank_four_decay=4.),dict(rank_four=True,rank_four_decay=-1.),
                    dict(rank_four=True,spatial_gradient=True)):
        with pytest.raises(ValueError):
            run(tmp_path/'missing',tmp_path/'out',**options)
    assert not (tmp_path/'out').exists()


def test_positive_density_vector_shape_is_sitewise_square():
    from .coordination_power_research import positive_power_site_jet
    density = np.zeros((3, 10)); density[:, 0] = [-.7, -.3, .2]
    moments = np.zeros((3, 10, 3)); moments[:, 0] = np.arange(9).reshape(3, 3)/10
    for t in (-.5, .5, 1., 2., 4., 8.):
        powers = np.array([positive_power_site_jet(density, moments, p) for p in (0., 1., 2.)])
        result = np.array([(1-t)**2, 2*t*(1-t), t*t])@powers
        expected = np.sum((1+t*density[:, 0])**2*np.sum(moments[:, 0]**2, axis=1))
        assert result[0] == pytest.approx(expected, abs=1e-13)
        assert result[0] >= 0


@pytest.mark.parametrize('k1', [8.4, 10.])
def test_zero_vector_amplitude_is_valid_exact_limit(k1):
    from .rank_one_range_material import RankOneRangeInterface
    shape = (*SHAPE[:4], k1, SHAPE[5])
    c = np.ones(10); c[6] = 0.
    old = RankOneRangeInterface(shape[:5], c)
    new = CoordinationScreenedInterface(shape, c, law='power')
    q = (.879, .267, .113)
    left, right = old.evaluate(q), new.evaluate(q)
    assert left.energy == right.energy
    np.testing.assert_array_equal(left.gradient, right.gradient)
    np.testing.assert_array_equal(left.hessian, right.hessian)
    assert not any(moment.invariant.rank == 1 for _, moment in new.base.moments)
    assert np.all(np.isfinite(new.screened_jet(q)[0]))


@pytest.mark.parametrize('power', [0., .5, 1., 2., 4.])
def test_positive_power_nested_identity_and_coordinate_derivatives(power):
    from .coordination_power_research import positive_power_site_jet
    from .coordination_screening import screened_site_jet
    model = CoordinationScreenedInterface(SHAPE, np.ones(10), law='power')
    q = np.array([.879, .267, .113])
    def jet(state):
        density, moments, _ = model.site_inputs(tuple(state))
        return 2*positive_power_site_jet(density, moments, power)
    result = jet(q)
    if power <= 1:
        density, moments, _ = model.site_inputs(tuple(q))
        np.testing.assert_allclose(result, 2*screened_site_jet(density, moments, power, law='power'), atol=1e-12)
    step = 2e-5
    gradient = [(jet(q+step*e)[0]-jet(q-step*e)[0])/(2*step) for e in np.eye(3)]
    H = np.column_stack([(jet(q+step*e)[1:4]-jet(q-step*e)[1:4])/(2*step) for e in np.eye(3)])
    from .vector_interface_reference import HESSIAN_INDICES
    np.testing.assert_allclose(gradient, result[1:4], atol=1e-7, rtol=2e-6)
    np.testing.assert_allclose(H, result[HESSIAN_INDICES], atol=2e-6, rtol=3e-6)
    # Independent real-space scalar/vector sums, nonlinear function per atom.
    if power in (0., 1.):
        assert result[0] == pytest.approx(model.direct_screened_energy(q, radius=16, layers=16)
            if power == 1 else model.direct_rank1_energy(q, radius=16, layers=16), abs=1e-11)
