import numpy as np
import pytest

from .quadrupole_saturation import (
    SaturatedQuadrupoleInterface, SaturatedQuadrupoleCache, saturated_site_norm_jet,
)
from .quartic_angular_material import QuarticSymmetryInterface, QuarticSymmetryObservationCache


C = np.array([.3, .5, 2., 1., .2, 3., -.1, 2., .03, 1200.])
D = (3., 6., 3.5)


def test_site_counting_rotation_gauge_and_no_pole():
    z = np.zeros((2, 10, 3)); z[:, 0] = [[1., 0., 0.], [0., 2., 0.]]
    v = saturated_site_norm_jet(z, alpha=2., reference_curvature=3.)
    assert v[0] == pytest.approx(1/(1+2/3) + 4/(1+8/3))
    assert v[0] != pytest.approx(5/(1+10/3))
    R = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
    np.testing.assert_allclose(v, saturated_site_norm_jet(z@R, alpha=2., reference_curvature=3.))
    np.testing.assert_allclose(v, saturated_site_norm_jet(7*z, alpha=2., reference_curvature=147.)/49)
    for a in (0., 1e-9, 1., 1e8):
        result = saturated_site_norm_jet(z, alpha=a, reference_curvature=3.)
        assert 0 < result[0] <= 5
    with pytest.raises(ValueError):
        saturated_site_norm_jet(z, alpha=-1., reference_curvature=3.)


def test_zero_shape_exactly_recovers_old_and_reference_hessian_is_preserved():
    old = QuarticSymmetryInterface(D, C)
    zero = SaturatedQuadrupoleInterface((*D, 0.), C)
    changed = SaturatedQuadrupoleInterface((*D, 30.), C)
    q = (1.1*old.h, .16, -.05)
    assert zero.evaluate(q).energy == old.evaluate(q).energy
    np.testing.assert_array_equal(zero.evaluate(q).hessian, old.evaluate(q).hessian)
    q = (old.h, 0., 0.)
    a, b = changed.evaluate(q), old.evaluate(q)
    np.testing.assert_allclose(a.hessian, b.hessian, rtol=1e-13, atol=1e-12)
    np.testing.assert_allclose(a.gradient, b.gradient, rtol=0, atol=1e-13)
    assert abs(a.energy) < 1e-13


@pytest.mark.parametrize('q', [(1.04*np.sqrt(2/3), .23, .07),
                               (1.23*np.sqrt(2/3), .37, -.08)])
def test_analytic_gradient_hessian_and_periodicity(q):
    model = SaturatedQuadrupoleInterface((*D, 12.), C)
    q = np.array(q); v = model.evaluate(q)
    errors = []
    for h in (4e-5, 2e-5):
        gradient, hessian = [], []
        for axis in range(3):
            step = np.eye(3)[axis]*h
            p, m = model.evaluate(q+step), model.evaluate(q-step)
            gradient.append((p.energy-m.energy)/(2*h))
            hessian.append((p.gradient-m.gradient)/(2*h))
        errors.append(max(np.max(abs(np.array(gradient)-v.gradient)),
                          np.max(abs(np.array(hessian).T-v.hessian))))
    assert errors[1] < errors[0]*.35
    assert errors[1] < 1e-5
    np.testing.assert_array_equal(v.hessian, v.hessian.T)
    repeat = model.evaluate(q+[0., 1., 0.])
    assert repeat.energy == pytest.approx(v.energy, abs=2e-12)
    np.testing.assert_allclose(repeat.hessian, v.hessian, rtol=2e-11, atol=2e-11)


def test_exact_coefficient_matrix_bulk_unchanged_and_new_holdouts_excluded():
    from .run_vector_material_calibration import source_and_targets
    from .interface_even_development_targets import even_development_observations
    from .interface_development_targets import development_observations
    source, old_obs, states = source_and_targets()
    previous, _ = development_observations(source, old_obs, states)
    obs, provenance = even_development_observations(source, old_obs, states)
    previous_hold = {o.name for o in previous if o.role == 'heldout'}
    assert all(o.role == 'fit' for o in obs if o.name in previous_hold)
    hold = [o for o in obs if o.role == 'heldout']
    assert len(hold) == 36 and all(o.name.startswith('v16_') for o in hold)
    assert set(provenance) == {o.name for o in obs}
    shape = (*D, 10.)
    mat = SaturatedQuadrupoleCache(obs).matrix(shape)
    base = QuarticSymmetryObservationCache(obs).matrix(D)
    np.testing.assert_array_equal(mat[:5], base[:5])
    model = SaturatedQuadrupoleInterface(shape, C)
    for o, row in zip(obs[::9], mat[::9]):
        if o.bulk_index is not None:
            continue
        v = model.evaluate(o.state)
        jet = np.r_[v.energy, v.gradient, v.hessian[0], v.hessian[1, 1:], v.hessian[2, 2]]
        assert row@C == pytest.approx(np.asarray(o.jet_weights)@jet, abs=2e-9, rel=2e-9)


def test_direct_lattice_and_reciprocal_neighborhood_convergence():
    model = SaturatedQuadrupoleInterface((*D, 8.), C)
    q = (1.08*model.h, .27, .08)
    jets = model.even_jets(q)
    expected = np.array([v[0] for v in jets[:2]])
    coarse = model.direct_even_energy(q, radius=6, layers=6)
    fine = model.direct_even_energy(q, radius=16, layers=16)
    assert np.linalg.norm(fine-expected) < np.linalg.norm(coarse-expected)
    np.testing.assert_allclose(fine, expected, rtol=2e-10, atol=2e-12)
    tight = SaturatedQuadrupoleInterface((*D, 8.), C, tolerance=2e-13)
    np.testing.assert_allclose(model.evaluate(q).hessian, tight.evaluate(q).hessian,
                               rtol=3e-11, atol=3e-11)


def test_unchanged_finite_q_is_verified_from_nonzero_shape_site_energy():
    from .isotropic_bulk_validation import IsotropicBulkBasis
    obj = IsotropicBulkBasis(D, stretch=1.007, radius=10.)
    columns, _ = obj.evaluate((1/3, 1/3, 1/3))
    H = np.einsum('c,cij->ij', C, columns)
    for v in np.eye(3):
        options = dict(planes=6, mode=2, polarization_plane=v,
                       validation_radius=10., quadrupole_saturation=20.)
        zero = obj.direct_sinusoidal_energy(C, amplitude=0., **options)
        approx = []
        # The rational term can have large fourth derivatives despite its
        # unchanged harmonic limit. Resolve that limit before roundoff takes
        # over; do not relax the finite-q agreement tolerance.
        for h in (1e-4, 5e-5, 2.5e-5):
            p = obj.direct_sinusoidal_energy(C, amplitude=h, **options)
            m = obj.direct_sinusoidal_energy(C, amplitude=-h, **options)
            approx.append(2*(p+m-2*zero)/h**2)
        fourth = [(4*b-a)/3 for a, b in zip(approx, approx[1:])]
        if v[2]:
            assert abs(fourth[1]-v@H@v) < .1*abs(fourth[0]-v@H@v)
        sixth = (16*fourth[1]-fourth[0])/15
        assert abs(sixth-v@H@v) < 1e-5


def test_curvature_report_keeps_all_levels_and_nonlinear_parameters():
    from .run_even_environment_validation import curvature_refinement

    class ExactPolynomialReference:
        radius = 12.

        def evaluate(self, q):
            return np.eye(3)[None]*7., np.array([.002])

        def direct_sinusoidal_energy(self, c, *, amplitude, **options):
            assert options['quadrupole_saturation'] == 3.
            assert options['validation_radius'] == self.radius
            return 7.*amplitude**2/4 + 1e6*amplitude**4

    rows = curvature_refinement(ExactPolynomialReference(), np.ones(1),
                                alpha=3., polarization=[0., 0., 1.])
    assert len(rows) == 6 and rows[0]['fourth_order'] is None
    assert rows[0]['central_error'] > rows[1]['central_error'] > 0
    assert rows[1]['sixth_order'] is None
    assert max(abs(r['sixth_error']) for r in rows[2:]) < 1e-13
    assert all(r['harmonic_tail_bound'] == .002 for r in rows)


def test_geometric_opening_brackets_find_narrow_extrema_missed_by_uniform_grids():
    from types import SimpleNamespace
    from .run_even_environment_validation import resolved_opening_extrema

    class NarrowAnalyticReference:
        h = 1.

        def evaluate(self, q):
            d = q[0]-1.
            H = (d-1e-4)*(d-2e-4)*1e8
            g = (d**3/3-3e-4*d**2/2+2e-8*d)*1e8
            return SimpleNamespace(energy=0., gradient=np.array([g, 0., 0.]),
                                   hessian=np.diag([H, 1., 1.]))

    model = NarrowAnalyticReference()
    for samples in (129, 257):
        coarse = [model.evaluate([a]).hessian[0, 0] for a in np.linspace(1., 5., samples)]
        assert min(coarse) > 0  # both old grids miss the negative-curvature pocket
        roots = resolved_opening_extrema(model, samples)
        assert len(roots) == 2
        np.testing.assert_allclose([r['a_over_h'] for r in roots], [1.0001, 1.0002], atol=3e-12, rtol=0)
