import numpy as np
from scipy.integrate import quad

from solver_v1.exponential_density import (
    ExponentialDensityParams,
    direct_exponential_cross_density,
    exponential_cross_density,
    exponential_fourier_transform,
    same_row_exponential_density,
)


def test_exponential_fourier_transform_matches_direct_quadrature() -> None:
    for a, kappa, wave in ((0.4, 1.2, 0.0), (0.8, 3.0, 2.3), (1.7, 0.7, 8.0)):
        numeric = 2.0 * quad(
            lambda x: np.exp(-kappa * np.sqrt(a * a + x * x)) * np.cos(wave * x),
            0.0,
            np.inf,
            epsabs=2.0e-13,
            epsrel=2.0e-13,
            limit=500,
        )[0]
        analytic = float(exponential_fourier_transform(a, wave, kappa))
        np.testing.assert_allclose(analytic, numeric, rtol=2.0e-11, atol=2.0e-12)


def test_poisson_bessel_density_matches_large_real_space_sum() -> None:
    cases = [
        (0.45, -0.31, 0.8, 1.2, 0.9),
        (0.77, 0.0, 1.0, 3.0, 1.0),
        (1.3, 0.27, 1.4, 2.2, 1.1),
        (2.0, 0.49, 0.7, 5.0, 0.85),
    ]
    rng = np.random.default_rng(8649)
    for _ in range(12):
        b = rng.uniform(0.8, 1.3)
        cases.append(
            (
                rng.uniform(0.4, 2.0),
                rng.uniform(-0.5 * b, 0.5 * b),
                rng.uniform(0.5, 1.5),
                rng.uniform(0.8, 6.0),
                b,
            )
        )
    for a, s, rho0, beta, b in cases:
        params = ExponentialDensityParams(
            rho0=rho0, beta_rho=beta, r_e=1.15, tol=1.0e-14,
            max_modes=180,
        )
        reciprocal = exponential_cross_density(a, s, b=b, params=params)
        direct = direct_exponential_cross_density(
            a, s, b=b, params=params, images=800
        )
        np.testing.assert_allclose(reciprocal.value, direct, rtol=2.0e-12, atol=2.0e-12)
        assert reciprocal.modes_used < params.max_modes
        assert reciprocal.estimated_tail_absolute < 1.0e-11


def test_adaptive_reciprocal_truncation_matches_high_accuracy_reference() -> None:
    rng = np.random.default_rng(20260908)
    for _ in range(16):
        a = rng.uniform(0.4, 1.8)
        b = rng.uniform(0.8, 1.3)
        s = rng.uniform(-0.5 * b, 0.5 * b)
        shared = dict(
            rho0=rng.uniform(0.5, 1.5),
            beta_rho=rng.uniform(0.8, 6.0),
            r_e=rng.uniform(0.8, 1.3),
            max_modes=180,
        )
        adaptive = exponential_cross_density(
            a, s, b=b, params=ExponentialDensityParams(tol=1.0e-10, **shared)
        )
        reference = exponential_cross_density(
            a, s, b=b, params=ExponentialDensityParams(tol=1.0e-15, **shared)
        )
        np.testing.assert_allclose(
            adaptive.value, reference.value, rtol=2.0e-10, atol=2.0e-11
        )
        assert adaptive.modes_used <= reference.modes_used


def test_same_row_density_is_exact_geometric_series_without_self() -> None:
    params = ExponentialDensityParams(rho0=1.7, beta_rho=2.4, r_e=0.9)
    exact = same_row_exponential_density(b=1.2, params=params)
    n = np.arange(1, 1000, dtype=float)
    direct = float(2.0 * np.sum(params.C_rho * np.exp(-params.kappa * 1.2 * n)))
    np.testing.assert_allclose(exact, direct, rtol=2.0e-15, atol=0.0)


def test_density_periodicity_and_large_a_decay() -> None:
    params = ExponentialDensityParams(beta_rho=3.2, r_e=1.1)
    first = exponential_cross_density(0.9, 0.17, b=1.0, params=params)
    translated = exponential_cross_density(0.9, 1.17, b=1.0, params=params)
    np.testing.assert_allclose(first.value, translated.value, rtol=0.0, atol=2.0e-14)
    np.testing.assert_allclose(first.d_ds, translated.d_ds, rtol=0.0, atol=2.0e-13)
    far = exponential_cross_density(8.0, 0.17, b=1.0, params=params)
    assert float(far.value) < 1.0e-8 * float(first.value)
    symmetric = exponential_cross_density(0.9, 0.0, b=1.0, params=params)
    assert abs(float(symmetric.d_ds)) < 1.0e-13


def test_large_beta_localizes_environment_density() -> None:
    broad = ExponentialDensityParams(beta_rho=2.0, r_e=1.0)
    localized = ExponentialDensityParams(beta_rho=10.0, r_e=1.0)
    broad_ratio = float(
        exponential_cross_density(1.8, 0.0, b=1.0, params=broad).value
        / exponential_cross_density(0.8, 0.0, b=1.0, params=broad).value
    )
    localized_ratio = float(
        exponential_cross_density(1.8, 0.0, b=1.0, params=localized).value
        / exponential_cross_density(0.8, 0.0, b=1.0, params=localized).value
    )
    assert localized_ratio < 1.0e-3 * broad_ratio


def test_analytic_density_gradient_and_hessian_match_finite_differences() -> None:
    params = ExponentialDensityParams(beta_rho=2.7, r_e=1.05, tol=2.0e-15)
    h1 = 2.0e-6
    h2 = 2.0e-4
    for a, s in ((0.62, -0.23), (0.81, 0.13), (1.25, 0.37)):
        center = exponential_cross_density(a, s, b=1.0, params=params)
        ap = exponential_cross_density(a + h1, s, b=1.0, params=params)
        am = exponential_cross_density(a - h1, s, b=1.0, params=params)
        sp = exponential_cross_density(a, s + h1, b=1.0, params=params)
        sm = exponential_cross_density(a, s - h1, b=1.0, params=params)
        rho_a_fd = (ap.value - am.value) / (2.0 * h1)
        rho_s_fd = (sp.value - sm.value) / (2.0 * h1)
        np.testing.assert_allclose(center.d_da, rho_a_fd, rtol=2.0e-9, atol=2.0e-9)
        np.testing.assert_allclose(center.d_ds, rho_s_fd, rtol=2.0e-9, atol=2.0e-9)

        ap2 = exponential_cross_density(a + h2, s, b=1.0, params=params)
        am2 = exponential_cross_density(a - h2, s, b=1.0, params=params)
        sp2 = exponential_cross_density(a, s + h2, b=1.0, params=params)
        sm2 = exponential_cross_density(a, s - h2, b=1.0, params=params)
        aa_fd = (ap2.value - 2.0 * center.value + am2.value) / h2**2
        ss_fd = (sp2.value - 2.0 * center.value + sm2.value) / h2**2
        mixed_h = 5.0e-5
        mixed_sp = exponential_cross_density(
            a, s + mixed_h, b=1.0, params=params
        )
        mixed_sm = exponential_cross_density(
            a, s - mixed_h, b=1.0, params=params
        )
        as_fd = (mixed_sp.d_da - mixed_sm.d_da) / (2.0 * mixed_h)
        np.testing.assert_allclose(center.d2_daa, aa_fd, rtol=2.0e-6, atol=2.0e-6)
        np.testing.assert_allclose(center.d2_dss, ss_fd, rtol=3.0e-6, atol=3.0e-6)
        np.testing.assert_allclose(center.d2_das, as_fd, rtol=2.0e-7, atol=2.0e-7)


def test_zero_density_amplitude_is_exactly_zero() -> None:
    params = ExponentialDensityParams(rho0=0.0)
    result = exponential_cross_density(
        np.array([0.7, 1.0]), np.array([0.1, 0.2]), b=1.0, params=params
    )
    for value in (
        result.value, result.d_da, result.d_ds, result.d2_daa,
        result.d2_das, result.d2_dss,
    ):
        np.testing.assert_array_equal(value, 0.0)
    assert same_row_exponential_density(b=1.0, params=params) == 0.0
