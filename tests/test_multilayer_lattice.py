from __future__ import annotations

import math

import numpy as np
import pytest

from theory.multilayer_probability import (
    cumulative_hysteresis,
    governing_equations_metrics,
    moving_barrier_outflux,
    plastic_strains,
    well_populations,
)
from theory.registry_lattice import (
    MultilayerPotentialParameters,
    b_q_direct,
    bessel_lambert,
    bessel_lambert_polylog,
    dU_da,
    dU_da_direct,
    dU_ds,
    delta_h_q,
    generalized_lj_coefficient,
    h_q_bessel,
    h_q_direct,
    h_q_delta_derivative_direct,
    h_q_polylog,
    normal_stationary_points,
    plastic_well_index,
    shifted_inverse_power_bessel,
    u0,
    u0_direct,
    v_slip,
)


@pytest.mark.parametrize("q", [6, 12])
@pytest.mark.parametrize("delta,eta", [(0.0, 0.7), (0.23, 1.0), (0.5, 1.4)])
def test_direct_h_matches_exact_bessel_lambert(q: int, delta: float, eta: float) -> None:
    direct = h_q_direct(q, delta, eta, kmax=400, pmax=800)
    analytic = h_q_bessel(q, delta, eta, modes=20, layer_modes=48)
    assert direct == pytest.approx(analytic, rel=1.3e-11, abs=3.0e-12)


def test_h_is_unweighted_sum_of_row_kernels_with_same_registry() -> None:
    q, delta, eta, kmax, pmax = 6.0, 0.27, 0.9, 30, 1000
    direct = h_q_direct(q, delta, eta, kmax, pmax)
    row_sum = sum(b_q_direct(q, delta, k * eta, pmax) for k in range(1, kmax + 1))
    weighted_wrong = sum(
        k * b_q_direct(q, delta, k * eta, pmax) for k in range(1, kmax + 1)
    )
    shifted_wrong = sum(
        b_q_direct(q, k * delta, k * eta, pmax) for k in range(1, kmax + 1)
    )
    assert direct == pytest.approx(row_sum, rel=2.0e-15)
    assert abs(direct - weighted_wrong) > 1.0e-3
    assert abs(direct - shifted_wrong) > 1.0e-6


@pytest.mark.parametrize("q", [6, 12])
def test_polylog_bessel_lambert_coefficients_are_independently_exact(q: int) -> None:
    for x in (0.9, 2.5, 7.0):
        direct = bessel_lambert((q - 1) / 2, x, kmax=120)
        polylog = bessel_lambert_polylog(q, x)
        assert polylog == pytest.approx(direct, rel=4.0e-15, abs=2.0e-15)
    for delta, eta in ((0.11, 0.7), (0.37, 1.2)):
        assert h_q_polylog(q, delta, eta, 24) == pytest.approx(
            h_q_bessel(q, delta, eta, 24, 80), rel=4.0e-15, abs=2.0e-15
        )


def test_periodicity_even_symmetry_and_slip_reference() -> None:
    params = MultilayerPotentialParameters(b=1.3, epsilon_lj=2.1, sigma_lj=0.91)
    a, s = 1.1, 0.22
    assert u0(a, s + params.b, params) == pytest.approx(u0(a, s, params), abs=2.0e-13)
    assert u0(a, -s, params) == pytest.approx(u0(a, s, params), abs=2.0e-13)
    assert v_slip(a, 0.0, params, 0.0) == pytest.approx(0.0, abs=2.0e-15)
    assert delta_h_q(6, 0.0, 0.0, a / params.b) == pytest.approx(0.0, abs=1.0e-15)


def test_u0_direct_and_analytic_physical_energy_agree() -> None:
    params = MultilayerPotentialParameters(b=1.0, epsilon_lj=1.7, sigma_lj=0.83)
    direct = u0_direct(0.94, 0.21, params, kmax=400, pmax=800)
    analytic = u0(0.94, 0.21, params)
    assert analytic == pytest.approx(direct, rel=1.5e-11)
    assert generalized_lj_coefficient(12, 6) == pytest.approx(4.0)


def test_analytic_u0_derivatives_match_finite_differences() -> None:
    params = MultilayerPotentialParameters(b=1.0, epsilon_lj=1.0, sigma_lj=0.82)
    a, s, step = 0.97, 0.19, 2.0e-6
    fd_a = (u0(a + step, s, params) - u0(a - step, s, params)) / (2 * step)
    fd_s = (u0(a, s + step, params) - u0(a, s - step, params)) / (2 * step)
    assert dU_da(a, s, params) == pytest.approx(fd_a, rel=2.0e-8, abs=2.0e-8)
    assert dU_ds(a, s, params) == pytest.approx(fd_s, rel=2.0e-8, abs=2.0e-8)
    assert dU_da(a, s, params) == pytest.approx(
        dU_da_direct(a, s, params, 300, 600), rel=2.0e-9, abs=2.0e-8
    )
    scale = generalized_lj_coefficient(params.m, params.n) * params.epsilon_lj / params.b
    direct_s = scale * (
        (params.sigma_lj / params.b) ** params.m
        * h_q_delta_derivative_direct(params.m, s / params.b, a / params.b, 300, 600)
        - (params.sigma_lj / params.b) ** params.n
        * h_q_delta_derivative_direct(params.n, s / params.b, a / params.b, 300, 600)
    )
    assert dU_ds(a, s, params) == pytest.approx(direct_s, rel=2.0e-9, abs=2.0e-8)


def test_unwrapped_wells_and_plastic_strain_use_z_transition() -> None:
    s = np.array([-0.75, -0.25, 0.25, 0.75, 1.25])
    assert np.array_equal(plastic_well_index(s, 1.0), [-1, 0, 0, 1, 1])
    density = np.array([[0.0, 0.0, 0.0, 0.4, 0.6]])
    populations, mean_z = well_populations(s, density, 1.0, 0.5, 1.0, 0.0)
    assert sum(populations.values()) == pytest.approx(0.5)
    gamma, epsilon = plastic_strains(mean_z, 1.0, 2.0, 0.4)
    assert epsilon == pytest.approx(0.4 * gamma)


def test_four_governing_metrics_separate_u0_from_dissipation() -> None:
    a, s = np.linspace(0.9, 1.1, 5), np.linspace(-0.2, 0.2, 5)
    da, ds = a[1] - a[0], s[1] - s[0]
    density = np.ones((5, 5)) / (25 * da * ds)
    energy = (a[:, None] - 1.0) ** 2 + s[None, :] ** 2
    ja, js = np.full((5, 5), 0.02), np.full((5, 5), -0.01)
    metrics = governing_equations_metrics(a, s, density, energy, ja, js, 2.0, 3.0, 0.0)
    assert metrics.normalization_or_survival == pytest.approx(1.0)
    assert metrics.mean_spacing == pytest.approx(1.0)
    assert metrics.dissipation_rate > 0.0
    cumulative = cumulative_hysteresis(np.array([0.0, 1.0, 2.0]), np.array([0.0, 2.0, 2.0]))
    assert np.array_equal(cumulative, [0.0, 1.0, 3.0])
    assert np.all(np.diff(cumulative) >= 0.0)


def test_moving_barrier_flux_uses_relative_current() -> None:
    s = np.linspace(-1.0, 1.0, 101)
    barrier = 2.0 + 0.2 * s
    velocity = np.full_like(s, 0.3)
    density = np.full_like(s, 0.4)
    ja = np.full_like(s, 0.5)
    js = np.full_like(s, 0.1)
    result = moving_barrier_outflux(s, barrier, velocity, density, ja, js)
    expected_density = 0.5 - 0.1 * 0.2 - 0.4 * 0.3
    assert result == pytest.approx(2.0 * expected_density)


def test_normal_escape_boundary_is_outer_negative_curvature_root() -> None:
    params = MultilayerPotentialParameters(sigma_lj=0.82, bessel_modes=12, layer_modes=24)
    roots = normal_stationary_points(
        generalized_force=0.1,
        s=0.0,
        params=params,
        a_min=0.6,
        a_max=2.8,
        samples=140,
        kmax=80,
        pmax=160,
    )
    assert len(roots) == 2
    assert roots[0][1] > 0.0
    assert roots[1][1] < 0.0
    assert roots[1][0] > roots[0][0]
