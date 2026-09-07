import numpy as np
import pytest

from app.solver_adapter import canonical_model_params
from solver_v1.dynamics_diagnostics import (
    conditional_fast_a_density,
    extract_harmonic_response,
    linear_axial_transfer,
    model_frequency_diagnostics,
    pristine_relaxation_spectrum,
    stable_normal_equilibrium,
    stationary_total_equilibrium_curve,
    stationary_total_equilibrium,
)
from solver_v1.model import TwoRowLJ
import solver_v1.probability_pde_2d as pde


@pytest.fixture(scope="module")
def prepared_model() -> TwoRowLJ:
    model = TwoRowLJ(canonical_model_params())
    model._build_opening_table()
    return model


def test_pristine_relaxation_rates_and_times_are_positive() -> None:
    model = TwoRowLJ(canonical_model_params())
    spectrum = pristine_relaxation_spectrum(model)
    np.testing.assert_allclose(
        spectrum.eigenvalues,
        [2.517379007312628, 122.78696473236548],
        rtol=2.0e-12,
    )
    assert spectrum.tau_fast == pytest.approx(0.00814418698417756)
    assert spectrum.tau_slow == pytest.approx(0.39723855529705393)
    assert np.all(spectrum.eigenvalues > 0.0)
    assert np.all(np.isfinite(1.0 / spectrum.eigenvalues))


def test_linear_transfer_has_relaxed_limit_and_finite_rate_lag() -> None:
    model = TwoRowLJ(canonical_model_params())
    assert linear_axial_transfer(model, 0.0) == pytest.approx(1.0 + 0.0j)
    frequencies = np.array([0.01, 0.1, 1.0, 10.0]) / (
        2.0 * np.pi * pristine_relaxation_spectrum(model).tau_slow
    )
    responses = np.asarray(
        [linear_axial_transfer(model, 2.0 * np.pi * f) for f in frequencies]
    )
    assert np.all(np.diff(np.abs(responses)) < 0.0)
    assert np.all(np.angle(responses) < 0.0)


def test_current_ui_model_frequency_diagnostics_are_dimensionless() -> None:
    model = TwoRowLJ(canonical_model_params())
    values = model_frequency_diagnostics(model, 25.0)
    assert values["model_period"] == pytest.approx(0.04)
    assert values["omega"] == pytest.approx(157.07963267948966)
    assert values["de_fast"] == pytest.approx(1.279285899947692)
    assert values["de_slow"] == pytest.approx(62.39808635219237)


def test_stable_fast_a_branch_and_quasistatic_tangent() -> None:
    model = TwoRowLJ(canonical_model_params())
    for reduced_stress in (-1.0e-5, 1.0e-5):
        force = float(model.force_from_sigma_over_E(reduced_stress))
        normal = stable_normal_equilibrium(model, 0.0, force)
        assert normal is not None
        assert normal.normal_stiffness > 0.0
        assert normal.residual < 1.0e-10
        total = stationary_total_equilibrium(model, reduced_stress)
        assert total["strain"] / reduced_stress == pytest.approx(1.0, abs=2.0e-3)
        assert total["plastic_strain"] == 0.0
    curve = stationary_total_equilibrium_curve(
        model, np.array([-1.0e-5, 0.0, 1.0e-5])
    )
    centered_tangent = (curve["strain"][2] - curve["strain"][0]) / 2.0e-5
    assert centered_tangent == pytest.approx(1.0, abs=2.0e-4)
    np.testing.assert_array_equal(curve["well_index"], 0.0)


def test_conditional_fast_a_projection_preserves_s_marginal(
    prepared_model: TwoRowLJ,
) -> None:
    grid = pde.build_grid(
        prepared_model,
        pde.Grid2DParams(
            n_a=31, n_s=31, a_lower_factor=0.82, a_upper=1.10, s_wells=1
        ),
    )
    density = pde.initial_gibbs_density(prepared_model, grid, preload_force=0.0)
    projected, unbound = conditional_fast_a_density(
        density, prepared_model, grid, force=0.1
    )
    np.testing.assert_allclose(
        np.sum(projected, axis=0) * grid.da,
        np.sum(density, axis=0) * grid.da,
        rtol=2.0e-13,
        atol=2.0e-15,
    )
    assert unbound == pytest.approx(0.0)


def test_implicit_generator_is_the_existing_sg_rhs(prepared_model: TwoRowLJ) -> None:
    grid = pde.build_grid(
        prepared_model,
        pde.Grid2DParams(n_a=11, n_s=15, a_upper=1.10, s_wells=1),
    )
    energy = pde.energy_grid(prepared_model, grid, force=0.0)
    density = np.arange(energy.size, dtype=float).reshape(energy.shape) + 1.0
    rhs, _ = pde._sg_rates_and_rhs(density, energy, prepared_model, grid)
    generator = pde._sg_generator_2d(energy, prepared_model, grid)
    np.testing.assert_allclose(
        generator @ density.ravel(), rhs.ravel(), rtol=2.0e-13, atol=5.0e-11
    )
    assert np.max(np.abs(np.asarray(generator.sum(axis=0)))) < 1.0e-11


def _pde_harmonic_response(
    model: TwoRowLJ,
    *,
    omega: float,
    steps_per_cycle: int,
    n_a: int = 41,
    n_s: int = 31,
):
    period = 2.0 * np.pi / omega
    reduced_amplitude = 1.0e-4
    load = pde.cyclic_load_from_sigma_over_E(
        model,
        sigma_over_E_min=-reduced_amplitude,
        sigma_over_E_max=reduced_amplitude,
        period=period,
        cycles=4.0,
    )
    result = pde.run_probability_pde_2d(
        prepared_model=model,
        grid_params=pde.Grid2DParams(
            n_a=n_a,
            n_s=n_s,
            a_lower_factor=0.82,
            a_upper=1.10,
            s_wells=1,
        ),
        time_params=pde.PDETimeParams(
            max_dt=period / steps_per_cycle,
            cfl=0.40,
            record_interval=period / 50.0,
            integrator="implicit",
        ),
        load=load,
        preload_force=0.0,
    )
    measured = extract_harmonic_response(
        result["time"],
        np.asarray(result["force"]) / model.sigma_over_E_force_scale(),
        result["strain"],
        omega,
        start_time=2.0 * period,
    )
    return measured, result


def test_small_signal_pde_matches_full_linear_response_after_transient(
    prepared_model: TwoRowLJ,
) -> None:
    omega = pristine_relaxation_spectrum(prepared_model).lambda_slow
    analytic = linear_axial_transfer(prepared_model, omega)
    measured, result = _pde_harmonic_response(
        prepared_model, omega=omega, steps_per_cycle=100
    )
    assert measured.magnitude == pytest.approx(abs(analytic), rel=0.05)
    assert measured.phase_degrees == pytest.approx(
        np.degrees(np.angle(analytic)), abs=3.0
    )
    assert np.max(np.abs(result["mass_balance_residual"])) < 1.0e-10
    assert np.min(result["minimum_density"]) >= -2.0e-12


def test_implicit_timestep_refinement_converges_at_ui_model_frequency(
    prepared_model: TwoRowLJ,
) -> None:
    omega = 2.0 * np.pi * 25.0
    measured = []
    mass_errors = []
    final_survival = []
    for steps in (40, 80, 160):
        response, result = _pde_harmonic_response(
            prepared_model,
            omega=omega,
            steps_per_cycle=steps,
            n_a=31,
            n_s=21,
        )
        measured.append(response.transfer)
        mass_errors.append(float(np.max(np.abs(result["mass_balance_residual"]))))
        final_survival.append(float(result["survival"][-1]))
    coarse_change = abs(measured[1] - measured[0])
    fine_change = abs(measured[2] - measured[1])
    assert fine_change < coarse_change
    assert fine_change < 0.025
    assert max(mass_errors) < 1.0e-10
    np.testing.assert_allclose(final_survival, 1.0, rtol=0.0, atol=1.0e-10)
