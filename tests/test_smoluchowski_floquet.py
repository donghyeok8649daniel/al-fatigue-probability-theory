import numpy as np
import pytest

from theory.smoluchowski_escape import (
    TransportConfig,
    conditional_equilibrium,
    finite_volume_generator,
    step,
    transport_grid,
)
from theory.smoluchowski_floquet import (
    adiabatic_mean_escape_rate,
    asymptotic_survival_prefactor,
    cycle_operator_matrix,
    dense_cycle_spectrum,
    direct_cycle_survival_ratios,
    frozen_principal_escape_rate,
    principal_survival_mode,
    propagate_cycle,
)


def _protocol(period=4.0, points=61):
    time = np.linspace(0.0, period, points)
    force = 0.008 + 0.007 * np.sin(2.0 * np.pi * time / period)
    return time, force


def _config(cells=70):
    return TransportConfig(
        cells=cells,
        inverse_temperature=1000.0,
        boundary="absorbing",
        initiation_definition="tangent_instability",
    )


def test_protocol_must_be_periodic_and_boundary_absorbing():
    time, force = _protocol()
    with pytest.raises(ValueError, match="absorbing"):
        principal_survival_mode(time, force, TransportConfig(cells=50))
    force[-1] += 1e-4
    with pytest.raises(ValueError, match="endpoint"):
        principal_survival_mode(time, force, _config(50))


def test_principal_mode_is_positive_periodic_and_mass_balanced():
    time, force = _protocol()
    c = _config()
    result = principal_survival_mode(
        time, force, c, max_dt=0.04, tolerance=1e-10)
    x, dx = transport_grid(c)
    assert 0.0 < result.multiplier < 1.0
    assert np.min(result.phase_density) >= 0.0
    assert np.sum(result.start_density) * dx == pytest.approx(1.0)
    assert result.phase_survival[-1] == pytest.approx(result.multiplier)
    assert 1.0 - result.multiplier == pytest.approx(
        np.sum(result.phase_outflux[1:] * np.diff(time)), abs=3e-12)
    assert result.integrated_hazard == pytest.approx(-np.log(result.multiplier))
    assert np.sum(result.phase_hazard[1:] * np.diff(time)) == pytest.approx(
        result.integrated_hazard, abs=3e-14)
    assert np.all(result.phase_hazard >= 0.0)
    assert result.mean_hazard_rate == pytest.approx(
        result.integrated_hazard / result.period)
    assert result.residual_l1 < 2e-9
    assert np.sum(np.abs(result.phase_conditional_density[-1]
                         - result.start_density)) * dx < 2e-9


def test_eigenstate_survival_is_exact_geometric_and_generic_ratios_converge():
    time, force = _protocol(points=41)
    c = _config(60)
    result = principal_survival_mode(
        time, force, c, max_dt=0.05, tolerance=1e-11)
    survival, ratios = direct_cycle_survival_ratios(
        result.start_density, time, force, c, cycles=6, max_dt=0.05)
    expected = result.multiplier ** np.arange(7)
    assert np.max(np.abs(survival - expected)) < 3e-10
    assert np.max(np.abs(ratios - result.multiplier)) < 3e-10

    x, dx = transport_grid(c)
    initial = conditional_equilibrium(x, dx, float(force[0]), c)
    _, generic_ratios = direct_cycle_survival_ratios(
        initial, time, force, c, cycles=8, max_dt=0.05)
    assert abs(generic_ratios[-1] - result.multiplier) < 2e-8
    generic = initial.copy()
    for _ in range(6):
        generic, _, _ = propagate_cycle(
            generic, time, force, c, max_dt=0.05)
        generic /= np.sum(generic) * dx
    assert np.sum(np.abs(generic - result.start_density)) * dx < 2e-8


def test_dense_operator_spectral_radius_matches_power_iteration():
    time, force = _protocol(period=2.0, points=21)
    c = _config(40)
    operator = cycle_operator_matrix(time, force, c, max_dt=0.05)
    result = principal_survival_mode(
        time, force, c, max_dt=0.05, tolerance=1e-11)
    eigenvalues = np.linalg.eigvals(operator)
    spectral_radius = float(np.max(np.abs(eigenvalues)))
    assert np.min(operator) >= -1e-14
    assert np.max(np.sum(operator, axis=0)) <= 1.0 + 2e-12
    assert spectral_radius == pytest.approx(result.multiplier, rel=2e-10)


def test_left_perron_mode_predicts_generic_initial_survival_prefactor():
    time, force = _protocol(period=2.0, points=21)
    c = _config(40)
    spectrum = dense_cycle_spectrum(time, force, c, max_dt=0.05)
    x, dx = transport_grid(c)
    assert spectrum.multiplier > spectrum.second_eigenvalue_modulus >= 0.0
    assert 0.0 <= spectrum.spectral_ratio < 1.0
    assert np.min(spectrum.right_density) >= 0.0
    assert np.min(spectrum.left_survival_weight) >= 0.0
    assert np.sum(spectrum.right_density) * dx == pytest.approx(1.0)
    assert np.dot(spectrum.left_survival_weight,
                  spectrum.right_density) * dx == pytest.approx(1.0)
    assert np.sum(np.abs(
        spectrum.operator @ spectrum.right_density
        - spectrum.multiplier * spectrum.right_density)) * dx < 2e-11
    assert np.max(np.abs(
        spectrum.operator.T @ spectrum.left_survival_weight
        - spectrum.multiplier * spectrum.left_survival_weight)) < 2e-11

    initial = conditional_equilibrium(x, dx, float(force[0]), c)
    coefficient = asymptotic_survival_prefactor(initial, spectrum, c)
    survival, _ = direct_cycle_survival_ratios(
        initial, time, force, c, cycles=10, max_dt=0.05)
    scaled = survival / spectrum.multiplier ** np.arange(survival.size)
    assert scaled[-1] == pytest.approx(coefficient, rel=2e-8)


def test_floquet_multiplier_grid_and_timestep_refinement():
    time, force = _protocol(points=81)
    coarse = principal_survival_mode(
        time, force, _config(60), max_dt=0.05, tolerance=1e-10)
    medium = principal_survival_mode(
        time, force, _config(100), max_dt=0.025, tolerance=1e-10)
    fine = principal_survival_mode(
        time, force, _config(140), max_dt=0.0125, tolerance=1e-10)
    coarse_error = abs(coarse.multiplier - fine.multiplier)
    medium_error = abs(medium.multiplier - fine.multiplier)
    assert medium_error < coarse_error
    assert medium_error < 6e-5


def test_cycle_flux_equals_mass_loss_for_arbitrary_density():
    time, force = _protocol(points=41)
    c = _config(60)
    x, dx = transport_grid(c)
    density = conditional_equilibrium(x, dx, float(force[0]), c)
    end, escaped, _ = propagate_cycle(
        density, time, force, c, max_dt=0.05)
    assert escaped == pytest.approx(1.0 - np.sum(end) * dx, abs=3e-12)


def test_principal_multiplier_is_independent_of_cycle_phase_origin():
    period = 4.0
    time = np.linspace(0.0, period, 41)
    c = _config(60)
    multipliers = []
    for phase in (0.0, 0.25, 0.5, 0.75):
        force = 0.008 + 0.007 * np.sin(
            2.0 * np.pi * (time / period + phase))
        multipliers.append(principal_survival_mode(
            time, force, c, max_dt=0.05, tolerance=1e-11).multiplier)
    assert np.ptp(multipliers) < 2e-13


def test_fast_cycle_escape_is_linear_in_period_to_leading_order():
    c = TransportConfig(
        cells=50, inverse_temperature=2000.0, boundary="absorbing",
        initiation_definition="tangent_instability")
    rates = []
    escapes = []
    for period in (0.125, 0.25):
        time = np.linspace(0.0, period, 41)
        force = 0.008 + 0.007 * np.sin(2.0 * np.pi * time / period)
        result = principal_survival_mode(
            time, force, c, max_dt=period / 80.0, tolerance=1e-10)
        rates.append(result.integrated_hazard / period)
        escapes.append(result.escape_per_cycle)
    assert abs(rates[1] - rates[0]) / rates[0] < 7e-4
    assert escapes[1] / escapes[0] == pytest.approx(2.0, rel=8e-4)


def test_continuous_generator_matches_backward_euler_step():
    c = _config(40)
    x, dx = transport_grid(c)
    density = conditional_equilibrium(x, dx, 0.008, c)
    dt = 0.017
    generator = finite_volume_generator(0.008, c)
    off_diagonal = generator.copy()
    np.fill_diagonal(off_diagonal, 0.0)
    assert np.min(off_diagonal) >= -1e-14
    assert np.max(np.sum(generator, axis=0)) <= 2e-12
    matrix_solution = np.linalg.solve(
        np.eye(c.cells) - dt * generator, density)
    flux_solution, _ = step(density, x, dx, dt, 0.008, c)
    assert np.max(np.abs(matrix_solution - flux_solution)) < 2e-13


def test_slow_cycle_hazard_approaches_frozen_generator_average():
    c = TransportConfig(
        cells=40, inverse_temperature=2000.0, boundary="absorbing",
        initiation_definition="tangent_instability")
    phase = np.linspace(0.0, 1.0, 41)
    force = 0.008 + 0.007 * np.sin(2.0 * np.pi * phase)
    adiabatic_rate = adiabatic_mean_escape_rate(force, c)
    slow = principal_survival_mode(
        80.0 * phase, force, c, max_dt=0.05, tolerance=1e-9)
    assert slow.mean_hazard_rate == pytest.approx(adiabatic_rate, rel=5e-3)

    frozen_rate = frozen_principal_escape_rate(0.008, c)
    constant = principal_survival_mode(
        2.0 * phase, np.full_like(phase, 0.008), c,
        max_dt=0.01, tolerance=1e-10)
    assert constant.mean_hazard_rate == pytest.approx(frozen_rate, rel=5e-5)
