"""Independent transform/quadrature/real-space checks of a research primitive."""
import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import j0
from .polynomial_exponential_density import QuadraticEnvelopeDensity, positive_exponential_log_curvature, plane_sum
from .fcc111_geometry import fcc111_geometry_from_b


def test_log_curvature_obstruction_and_positive_kernel():
    assert positive_exponential_log_curvature(2., [1.], [3.]) == 0.
    assert positive_exponential_log_curvature(2., [1., 3.], [2., 5.]) > 0
    model = QuadraticEnvelopeDensity(.7, 2., 1.8, .3)
    r = np.linspace(0., 8., 201)
    f, fp, fpp = [model.radial(r, n) for n in (0, 1, 2)]
    assert np.all(f > 0)
    assert np.min(fpp/f-(fp/f)**2) < -1
    with pytest.raises(ValueError):
        QuadraticEnvelopeDensity(1., 1., 1., 0.)


@pytest.mark.parametrize('d,G', [(.7, 0.), (.7, 3.), (1.2, 7.)])
def test_fourier_hankel_quadrature_and_normal_derivatives(d, G):
    model = QuadraticEnvelopeDensity(.7, 2., 1.8, .3)
    reference, error = quad(lambda r: 2*np.pi*r*model.radial(np.hypot(d, r))*j0(G*r),
                            0, np.inf, epsabs=2e-12, epsrel=2e-12, limit=300)
    assert abs(model.plane_transform(d, G)-reference) < 3e-11+error
    h = 2e-5
    first = (model.plane_transform(d+h, G)-model.plane_transform(d-h, G))/(2*h)
    second = (model.plane_transform(d+h, G, 1)-model.plane_transform(d-h, G, 1))/(2*h)
    assert first == pytest.approx(model.plane_transform(d, G, 1), rel=2e-7, abs=1e-8)
    assert second == pytest.approx(model.plane_transform(d, G, 2), rel=2e-7, abs=1e-8)


def test_length_unit_conversion():
    physical = QuadraticEnvelopeDensity(.7, 2., 1.8, .3)
    L = 2.8637824638
    reduced = QuadraticEnvelopeDensity(physical.C*L*L, physical.k*L, physical.center/L, physical.width/L)
    for order in (0, 1, 2):
        assert reduced.plane_transform(.9/L, 2.1*L, order) == pytest.approx(
            L**(order-2)*physical.plane_transform(.9, 2.1, order), abs=2e-12)


def test_infinite_plane_direct_sum_full_jet_and_refinement():
    geometry = fcc111_geometry_from_b(1.)
    model = QuadraticEnvelopeDensity(.7, 2., 1.8, .3)
    d, delta = .81, np.array([.17, -.23])
    result = plane_sum(model, d, delta, geometry)
    errors = []
    for radius in (8, 16, 32):
        m, n = np.meshgrid(np.arange(-radius, radius+1), np.arange(-radius, radius+1), indexing='ij')
        xy = m.ravel()[:, None]*geometry.a1+n.ravel()[:, None]*geometry.a2+delta
        R = np.column_stack([np.full(len(xy), d), xy]); r = np.linalg.norm(R, axis=1); e = R/r[:, None]
        f, fp, fpp = [model.radial(r, order) for order in (0, 1, 2)]
        gradient = np.sum(fp[:, None]*e, axis=0)
        H = np.einsum('n,ni,nj->ij', fpp-fp/r, e, e)+np.eye(3)*np.sum(fp/r)
        errors.append(abs(f.sum()-result.value))
    assert errors[-1] < 2e-11 and errors[-1] < errors[0]
    np.testing.assert_allclose(gradient, np.r_[result.d_d, result.grad_delta], atol=2e-11)
    assembled = np.block([[np.array([[result.d2_dd]]), result.mixed_d_delta[None]],
                          [result.mixed_d_delta[:, None], result.hess_delta]])
    np.testing.assert_allclose(H, assembled, atol=3e-11)
    tight = plane_sum(model, d, delta, geometry, tolerance=1e-14)
    np.testing.assert_allclose(tight.hess_delta, result.hess_delta, atol=2e-12)
    assert plane_sum(model, d, delta+geometry.a1, geometry).value == pytest.approx(result.value, abs=2e-12)
