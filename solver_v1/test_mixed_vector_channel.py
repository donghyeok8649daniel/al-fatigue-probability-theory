"""Coherent radial fields are summed before their per-site quadratic energy."""
import numpy as np
import pytest
from .coordination_screening import CoordinationScreenedInterface, CoordinationScreenedBulk, screened_site_jet
from .mixed_vector_channel import mix_site_moments, mixed_vector_bulk_column
from .tail_constrained_material import normalized_amplitude
from .vector_interface_reference import HESSIAN_INDICES


SHAPE = (4.9, 10., 11., 1000., 8.4, 1.)
OTHER = (*SHAPE[:4], 5., SHAPE[-1])


def test_site_cancellation_and_coherent_not_incoherent():
    a = np.zeros((2, 10, 3)); a[:, 0] = [[1., 2., 3.], [.1, -.2, .3]]
    density = np.zeros((2, 10))
    np.testing.assert_array_equal(mix_site_moments(a, a, -1.), 0.)
    whole = screened_site_jet(density, mix_site_moments(a, a, 1.), 1., law='power')
    single = screened_site_jet(density, a, 1., law='power')
    np.testing.assert_allclose(whole, 4*single)
    assert whole[0] != 2*single[0]
    with pytest.raises(ValueError):
        mix_site_moments(a, a, np.nan)


def test_mixed_interface_analytic_jet():
    first = CoordinationScreenedInterface(SHAPE, np.ones(10), law='power')
    second = CoordinationScreenedInterface(OTHER, np.ones(10), law='power')
    def jet(q, weight=-1.):
        rho, a, _ = first.site_inputs(tuple(q))
        _, b, _ = second.site_inputs(tuple(q))
        return 2*screened_site_jet(rho, mix_site_moments(a, b, weight), 1., law='power')
    q = np.array([.879, .267, .113]); result = jet(q); h = 2e-5
    np.testing.assert_allclose(jet(q, 0), first.screened_jet(q)[0], atol=1e-13)
    g = [(jet(q+h*e)[0]-jet(q-h*e)[0])/(2*h) for e in np.eye(3)]
    H = np.column_stack([(jet(q+h*e)[1:4]-jet(q-h*e)[1:4])/(2*h) for e in np.eye(3)])
    np.testing.assert_allclose(g, result[1:4], atol=2e-7, rtol=2e-6)
    np.testing.assert_allclose(H, result[HESSIAN_INDICES], atol=3e-6, rtol=3e-6)


def test_bulk_nested_direct_energy_and_tail():
    first = CoordinationScreenedBulk(SHAPE, radius=12., law='power')
    second = CoordinationScreenedBulk(OTHER, radius=12., law='power')
    q = (1/3,)*3
    H0, _ = mixed_vector_bulk_column(first, second, q, 0.)
    np.testing.assert_array_equal(H0, first.evaluate(q)[0][6])
    H, tail = mixed_vector_bulk_column(first, second, q, -1.)
    assert np.linalg.eigvalsh(H).min() >= -1e-12
    R = first.R; layers = np.rint(R[:, 2]/first.geometry.h111)
    def direct(amplitude, v):
        values = []
        for site in range(6):
            shift = amplitude*(np.cos(2*np.pi*(site+layers)/3)-np.cos(2*np.pi*site/3))
            positions = R+shift[:, None]*v; r = np.linalg.norm(positions, axis=1)
            rho = np.sum(normalized_amplitude(SHAPE[0])*np.exp(-SHAPE[0]*r))
            moments = [(normalized_amplitude(k)*np.exp(-k*r))@positions for k in (SHAPE[4], OTHER[4])]
            Q = moments[0]-moments[1]
            values.append(rho*(Q@Q))
        return np.mean(values)
    for v in np.eye(3):
        slopes = [2*(direct(h, v)+direct(-h, v)-2*direct(0., v))/h**2 for h in (1e-4, 5e-5)]
        extrapolated = (4*slopes[1]-slopes[0])/3
        assert abs(extrapolated-v@H@v) < 2e-6+tail
