"""Independent identities for the v19 static hypothesis, never Al acceptance."""
import numpy as np
import pytest

from .coordination_screening import (
    CoordinationScreenedInterface, CoordinationScreenedBulk,
    screening_factor, screened_site_jet,
)
from .rank_one_range_material import RankOneRangeInterface, RankOneRangeBulk
from .vector_interface_reference import HESSIAN_INDICES, _square_jet


SHAPE = (2.4, 6.1, 8., 30., 6.3)
C = np.array([.018, .10, 1.4, 1.6, 2.9, 37., 4.7, 2., 1.5, 32000.])


@pytest.fixture(scope='module')
def model():
    return CoordinationScreenedInterface((*SHAPE, .8), C)


def test_factor_nested_limits_derivatives_and_refusal():
    x = np.array([.5, 1., 1.5])
    np.testing.assert_array_equal(screening_factor(x, 0)[0], np.ones(3))
    np.testing.assert_allclose(screening_factor(x, 1)[0], 1/x)
    z = .7; step = 1e-5
    g, first, second = screening_factor(x, z)
    np.testing.assert_allclose((screening_factor(x+step, z)[0]-screening_factor(x-step, z)[0])/(2*step), first, rtol=3e-10)
    np.testing.assert_allclose((screening_factor(x+step, z)[1]-screening_factor(x-step, z)[1])/(2*step), second, rtol=5e-10)
    for bad_x, bad_z in ((0., 1), (-1., .2), (np.nan, 0), (1., -1), (1., 1.1)):
        with pytest.raises(ValueError):
            screening_factor(bad_x, bad_z)
    with pytest.raises(ValueError, match='isolated'):
        CoordinationScreenedInterface((14., 6., 8., 0., 6., 1.), C)


def test_sitewise_not_global_and_padding():
    rho = np.zeros((2, 10)); rho[:, 0] = [-.4, -.1]
    tensor = np.zeros((3, 10, 3)); tensor[:, 0, 0] = [.2, .3, .05]
    out = screened_site_jet(rho, tensor, .9)
    expected = np.sum(np.array([.04, .09, .0025])/(1+.9*np.r_[rho[:, 0], 0]))
    assert out[0] == pytest.approx(expected)
    assert abs(out[0]-(.04+.09+.0025)/(1+.9*np.mean(rho[:, 0]))) > .005
    np.testing.assert_allclose(screened_site_jet(rho, tensor, 0), _square_jet(tensor))


def test_zero_screening_exact_previous_model():
    new = CoordinationScreenedInterface((*SHAPE, 0.), C)
    old = RankOneRangeInterface(SHAPE, C)
    for q in ((new.h, 0., 0.), (1.2*new.h, .23, -.05)):
        a, b = new.evaluate(q), old.evaluate(q)
        assert a.energy == b.energy
        np.testing.assert_array_equal(a.gradient, b.gradient)
        np.testing.assert_array_equal(a.hessian, b.hessian)


def test_perfect_tangent_unchanged_without_coefficient_refit(model):
    q = (model.h, 0., 0.)
    old = RankOneRangeInterface(SHAPE, C).evaluate(q)
    new = model.evaluate(q)
    assert abs(new.energy) < 1e-12
    np.testing.assert_allclose(new.gradient, old.gradient, atol=1e-12)
    np.testing.assert_allclose(new.hessian, old.hessian, atol=1e-11)
    # Arbitrary test coefficients are not force-equilibrated material data.


@pytest.mark.parametrize('q', [(1.12, .21, -.07), (.84, .47, .09)])
def test_coordinate_gradient_hessian_and_periodicity(model, q):
    q = np.array(q); result = model.evaluate(q); errors = []
    for step in (4e-5, 2e-5):
        gradient, H = [], []
        for direction in np.eye(3):
            hi, lo = model.evaluate(q+step*direction), model.evaluate(q-step*direction)
            gradient.append((hi.energy-lo.energy)/(2*step))
            H.append((hi.gradient-lo.gradient)/(2*step))
        errors.append(np.max(abs(np.array(H).T-result.hessian)))
        np.testing.assert_allclose(gradient, result.gradient, rtol=2e-6, atol=1e-7)
    assert errors[-1] < 3e-5
    assert errors[-1] < .4*errors[0]
    np.testing.assert_array_equal(result.hessian, result.hessian.T)
    repeat = model.evaluate(q+[0., 1., 0.])
    np.testing.assert_allclose(repeat.hessian, result.hessian, atol=2e-10, rtol=1e-10)


def test_direct_per_site_and_reciprocal_refinement(model):
    q = (1.5*model.h, .24, -.08)
    exact = model.screened_jet(q)[0][0]
    errors = [abs(model.direct_screened_energy(q, radius=r, layers=r)-exact) for r in (6, 10, 16)]
    assert errors[-1] < 1e-11
    assert errors[1] < errors[0]
    finer = CoordinationScreenedInterface((*SHAPE, .8), C, tolerance=2e-13)
    np.testing.assert_allclose(finer.evaluate(q).hessian, model.evaluate(q).hessian, atol=3e-10)


def test_full_finite_q_reference_and_dilated_factor():
    for stretch in (1., .99, 1.01):
        old = RankOneRangeBulk(SHAPE, stretch=stretch, radius=8.)
        new = CoordinationScreenedBulk((*SHAPE, .8), stretch=stretch, radius=8.)
        H, err = old.evaluate((.3, .2, .1)); screened, bounds = new.evaluate((.3, .2, .1))
        expected = H.copy(); expected[6] *= screening_factor(old.x, .8)[0]
        np.testing.assert_allclose(screened, expected, rtol=1e-14)
        assert np.all(bounds >= 0)
        if stretch == 1:
            np.testing.assert_array_equal(screened, H)
            np.testing.assert_array_equal(bounds, err)


@pytest.mark.parametrize('exponent', [-1., .6, 1.])
def test_power_law_derivatives_and_same_nested_reference(exponent):
    q = np.array([1.04, .27, -.063])
    obj = CoordinationScreenedInterface((*SHAPE, exponent), C, law='power')
    actual = obj.evaluate(q); step = 2e-5
    H = np.column_stack([(obj.evaluate(q+step*e).gradient-obj.evaluate(q-step*e).gradient)/(2*step)
                         for e in np.eye(3)])
    np.testing.assert_allclose(H, actual.hessian, atol=3e-6, rtol=2e-6)
    assert obj.direct_screened_energy(q, radius=16, layers=16) == pytest.approx(
        obj.screened_jet(q)[0][0], rel=1e-10, abs=1e-12)
    if exponent == -1:
        inverse = CoordinationScreenedInterface((*SHAPE, 1.), C)
        np.testing.assert_allclose(actual.hessian, inverse.evaluate(q).hessian, atol=1e-12)


@pytest.mark.parametrize('law,parameter', [('rational', .85), ('power', .7)])
def test_displaced_site_energy_recovers_dilated_finite_q(law, parameter):
    obj = CoordinationScreenedBulk((*SHAPE, parameter), radius=12., stretch=1.007, law=law)
    columns, tail = obj.evaluate((1/3, 1/3, 1/3))
    H = np.einsum('c,cij->ij', C, columns)
    for v in np.eye(3):
        options = dict(planes=6, mode=2, polarization_plane=v, validation_radius=12.,
                       quadrupole_saturation=SHAPE[3])
        zero = obj.direct_sinusoidal_energy(C, amplitude=0., **options)
        second = []
        for h in (1e-4, 5e-5, 2.5e-5):
            hi = obj.direct_sinusoidal_energy(C, amplitude=h, **options)
            lo = obj.direct_sinusoidal_energy(C, amplitude=-h, **options)
            second.append(2*(hi+lo-2*zero)/h**2)
        fourth = [(4*b-a)/3 for a, b in zip(second, second[1:])]
        sixth = (16*fourth[1]-fourth[0])/15
        assert abs(sixth-v@H@v) < 2e-5+tail@abs(C)


def test_new_validation_excluded_and_fixed_coefficient_matrix(model):
    from .run_coordination_screening import declared_new_observations
    from .run_vector_material_calibration import source_and_targets
    from .coordination_screening import CoordinationScreenedCache
    source, _, _ = source_and_targets()
    first = declared_new_observations(source)
    second = declared_new_observations(source, followup=True)
    assert len(first) == len(second) == 40
    assert all(o.role == 'heldout' for o in first+second)
    assert not {o.state for o in first} & {o.state for o in second}
    subset = second[:4]
    cache = CoordinationScreenedCache(subset, law='power')
    obj = CoordinationScreenedInterface((*SHAPE, .6), C, law='power')
    M = cache.matrix((*SHAPE, .6))
    for row, o in zip(M, subset):
        actual = obj.evaluate(o.state)
        jet = np.r_[actual.energy, actual.gradient, actual.hessian[0], actual.hessian[1, 1:], actual.hessian[2, 2]]
        assert row@C == pytest.approx(np.asarray(o.jet_weights)@jet, abs=3e-10, rel=1e-10)


def test_heldout_targets_cannot_change_coefficient_loss():
    from dataclasses import replace
    from .vector_material_calibration import MaterialObservation
    from .tail_constrained_material import convex_profile
    # Two coefficient toy problem exercises the SAME constrained solver, not
    # a statement that every measured source observable must fit perfectly.
    M = np.array([[1., 0.], [0., 1.], [0., 2.], [1., 1.]])
    obs = [MaterialObservation('exact', 1., 1., 'test', 'exact'),
           MaterialObservation('fit1', 2., 1., 'test', 'fit'),
           MaterialObservation('fit2', 4.1, 1., 'test', 'fit'),
           MaterialObservation('held', 9., 1., 'test', 'heldout')]
    first = convex_profile(M, obs, nonnegative=(0, 1))
    altered = obs[:-1]+[replace(obs[-1], target=-1e9)]
    second = convex_profile(M, altered, nonnegative=(0, 1))
    np.testing.assert_array_equal(first['coefficients'], second['coefficients'])
    assert first['squared_loss'] == second['squared_loss']


def test_research_loader_binds_explicit_law_and_refuses_zero_pair(tmp_path):
    import json
    from .validate_coordination_screening import load_candidate, independent_states
    result = tmp_path/'candidate'; result.mkdir()
    (tmp_path/'definition.json').write_text(json.dumps(dict(screening_law='power')))
    record = dict(completed=True, shape=[*SHAPE, -.7], screening_law='power',
                  best=dict(strictly_positive_LJ=True, coefficients=C.tolist()))
    (result/'calibration.json').write_text(json.dumps(record))
    obj, shape, coefficients, law, _ = load_candidate(result)
    assert law == obj.law == 'power' and obj.screening == -.7
    np.testing.assert_array_equal(coefficients, C)
    assert len(independent_states()) == len(set(independent_states())) == 10
    record['best']['strictly_positive_LJ'] = False
    (result/'calibration.json').write_text(json.dumps(record))
    with pytest.raises(ValueError, match='positive-LJ'):
        load_candidate(result)


def test_cubic_export_targets_units_and_residuals_stay_together():
    from .run_vector_material_calibration import source_and_targets
    from .run_coordination_shape_refinement import cubic_observations, residual_table
    _, original, _ = source_and_targets()
    converted = cubic_observations(original)
    predictions = np.array([o.target for o in converted])
    prescribed = np.linspace(-.01, .01, len(original))
    predictions += prescribed*np.array([o.scale for o in converted])
    rows = residual_table(original, dict(predictions=predictions, residuals=prescribed))
    assert [r['observable'] for r in rows[2:5]] == ['C11_GPa', 'C12_GPa', 'C44_GPa']
    assert all(r['units'] == 'GPa' for r in rows[2:5])
    np.testing.assert_allclose([r['target'] for r in rows[2:5]], [114., 62., 32.])
    for row in rows:
        assert (row['prediction']-row['target'])/row['scale'] == pytest.approx(row['normalized_residual'])
    with pytest.raises(AssertionError):
        residual_table(original, dict(predictions=predictions, residuals=prescribed+1.))


def test_saved_v19_fit_tables_have_consistent_units_and_residuals():
    import csv
    from pathlib import Path
    root = Path(__file__).resolve().parents[1]/'results/fcc111_active_interface/coordination_screening_v19'
    tables = list(root.glob('**/residuals.csv'))
    assert len(tables) == 6
    for table in tables:
        with table.open(encoding='utf-8', newline='') as stream:
            rows = list(csv.DictReader(stream))
        assert [r['observable'] for r in rows[2:5]] == ['C11_GPa', 'C12_GPa', 'C44_GPa']
        assert all(r['units'] == 'GPa' for r in rows[2:5])
        for row in rows:
            recomputed = (float(row['prediction'])-float(row['target']))/float(row['scale'])
            assert recomputed == pytest.approx(float(row['normalized_residual']), abs=1e-10)
