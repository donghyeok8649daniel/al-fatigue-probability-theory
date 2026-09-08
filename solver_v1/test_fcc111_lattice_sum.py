import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import j0

from solver_v1.fcc111_geometry import fcc111_geometry_from_b
from solver_v1.fcc111_lattice_sum import (
    exponential_kernel_fourier_2d,
    plane_exponential_sum_direct,
    plane_exponential_sum_reciprocal,
    plane_power_sum_direct,
    plane_power_sum_reciprocal,
    power_kernel_fourier_2d,
    same_plane_power_sum_direct,
    triangular_epstein_zeta,
    triangular_reciprocal_shells,
)


def test_reciprocal_shell_degeneracies_are_generated_not_assumed():
    geometry = fcc111_geometry_from_b(1.0)
    shells = triangular_reciprocal_shells(geometry, 8)
    observed = {shell.quadratic_index: shell.degeneracy for shell in shells}
    assert observed[1] == 6
    assert observed[3] == 6
    assert observed[4] == 6
    assert observed[7] == 12
    for shell in shells:
        np.testing.assert_allclose(
            np.linalg.norm(shell.vectors, axis=1), shell.magnitude, rtol=2e-15
        )


def test_power_and_exponential_fourier_transforms_match_hankel_integrals():
    d = 0.83
    wave = 2.4
    exponent = 3.0
    decay = 1.7
    power_numeric = 2.0 * math.pi * quad(
        lambda radius: radius * j0(wave * radius) * (d * d + radius * radius) ** (-exponent),
        0.0,
        80.0,
        epsabs=2e-11,
        limit=1000,
    )[0]
    exponential_numeric = 2.0 * math.pi * quad(
        lambda radius: radius
        * j0(wave * radius)
        * math.exp(-decay * math.sqrt(d * d + radius * radius)),
        0.0,
        80.0,
        epsabs=2e-11,
        limit=1000,
    )[0]
    assert power_kernel_fourier_2d(d, wave, exponent) == pytest.approx(
        power_numeric, rel=2e-9, abs=2e-10
    )
    assert exponential_kernel_fourier_2d(d, wave, decay) == pytest.approx(
        exponential_numeric, rel=2e-10, abs=2e-11
    )
    assert power_kernel_fourier_2d(d, 0.0, exponent) == pytest.approx(
        math.pi / (exponent - 1.0) * d ** (2.0 - 2.0 * exponent)
    )


@pytest.mark.parametrize("exponent,radius", ((3.0, 300), (6.0, 90)))
def test_plane_power_poisson_sum_matches_large_direct_disk(exponent, radius):
    geometry = fcc111_geometry_from_b(1.0)
    d = 0.73
    delta = np.array([0.17, -0.09])
    reciprocal = plane_power_sum_reciprocal(
        d, delta, p=exponent, geometry=geometry
    )
    direct = plane_power_sum_direct(
        d, delta, p=exponent, geometry=geometry, radius_index=radius
    )
    assert abs(reciprocal.value - direct.value) <= max(
        1.1 * direct.continuum_tail_estimate, 2.0e-12
    )
    np.testing.assert_allclose(
        reciprocal.grad_delta, direct.grad_delta, rtol=2e-10, atol=2e-11
    )
    np.testing.assert_allclose(
        reciprocal.mixed_d_delta,
        direct.mixed_d_delta,
        rtol=2e-10,
        atol=2e-11,
    )
    np.testing.assert_allclose(
        reciprocal.hess_delta, direct.hess_delta, rtol=3e-10, atol=4e-11
    )


def test_plane_exponential_poisson_sum_matches_direct_to_roundoff():
    geometry = fcc111_geometry_from_b(1.0)
    reciprocal = plane_exponential_sum_reciprocal(
        0.73,
        np.array([0.17, -0.09]),
        kappa=2.3,
        amplitude=1.7,
        geometry=geometry,
    )
    direct = plane_exponential_sum_direct(
        0.73,
        np.array([0.17, -0.09]),
        kappa=2.3,
        amplitude=1.7,
        geometry=geometry,
        radius_index=24,
    )
    assert reciprocal.value == pytest.approx(direct.value, abs=2e-14)
    assert reciprocal.d_d == pytest.approx(direct.d_d, abs=2e-14)
    assert reciprocal.d2_dd == pytest.approx(direct.d2_dd, abs=2e-13)
    np.testing.assert_allclose(reciprocal.grad_delta, direct.grad_delta, atol=2e-14)
    np.testing.assert_allclose(
        reciprocal.mixed_d_delta, direct.mixed_d_delta, atol=2e-14
    )
    np.testing.assert_allclose(reciprocal.hess_delta, direct.hess_delta, atol=2e-13)


def test_analytic_plane_derivatives_match_centered_finite_differences():
    geometry = fcc111_geometry_from_b(1.0)
    d = 0.82
    delta = np.array([0.11, 0.07])
    step = 2.0e-5

    def value(distance, shift):
        return plane_power_sum_reciprocal(
            distance, shift, p=3.0, geometry=geometry
        ).value

    result = plane_power_sum_reciprocal(d, delta, p=3.0, geometry=geometry)
    fd_d = (value(d + step, delta) - value(d - step, delta)) / (2.0 * step)
    fd_dd = (value(d + step, delta) - 2.0 * value(d, delta) + value(d - step, delta)) / step**2
    assert result.d_d == pytest.approx(fd_d, rel=3e-8)
    assert result.d2_dd == pytest.approx(fd_dd, rel=2e-6)
    for axis in range(2):
        offset = np.zeros(2)
        offset[axis] = step
        fd = (value(d, delta + offset) - value(d, delta - offset)) / (2.0 * step)
        assert result.grad_delta[axis] == pytest.approx(fd, rel=2e-8, abs=2e-10)


def test_same_plane_epstein_value_matches_large_direct_sum_with_tail_bound():
    geometry = fcc111_geometry_from_b(1.0)
    for exponent, radius in ((3.0, 300), (6.0, 80)):
        analytic = triangular_epstein_zeta(exponent, geometry.b)
        direct, _, tail = same_plane_power_sum_direct(
            p=exponent, geometry=geometry, radius_index=radius
        )
        assert abs(analytic - direct) <= max(1.1 * tail, 8.0e-15)
