from dataclasses import replace

import numpy as np
import pytest

from app.solver_adapter import canonical_model_params
from solver_v1.model import TwoRowLJ
from solver_v1.probability_pde_2d import CyclicLoad2D
import solver_v1.probability_pde_2d as full_pde
from solver_v1.reduced_fast_a import (
    ReducedFastAParams,
    bound_basin_conditional,
    conditional_fast_a,
    fast_a_qsd,
    laplace_F_eff,
    reduced_equilibrium,
    run_reduced_fast_a,
)


def _model(*, kT: float = 0.02) -> TwoRowLJ:
    return TwoRowLJ(replace(canonical_model_params(), kT=kT))


def test_fixed_domain_F_eff_gradient_is_conditional_mean_force() -> None:
    model = _model()
    s = 0.13
    force = 0.7
    step = 2.0e-5
    center = conditional_fast_a(
        model, s, force, a_lower=0.68, a_upper=1.20, quadrature_order=192
    )
    plus = conditional_fast_a(
        model, s + step, force, a_lower=0.68, a_upper=1.20,
        quadrature_order=192,
    )
    minus = conditional_fast_a(
        model, s - step, force, a_lower=0.68, a_upper=1.20,
        quadrature_order=192,
    )
    finite_difference = (plus.F_eff - minus.F_eff) / (2.0 * step)
    assert center.boundary_correction_s == 0.0
    assert center.dF_eff_ds == pytest.approx(center.mean_dG_ds, abs=1.0e-14)
    assert finite_difference == pytest.approx(center.mean_dG_ds, rel=2.0e-6)


def test_moving_saddle_derivative_contains_Leibniz_boundary_term() -> None:
    model = _model()
    s = 0.17
    force = 3.5
    step = 1.0e-5
    center = bound_basin_conditional(model, s, force, quadrature_order=256)
    plus = bound_basin_conditional(model, s + step, force, quadrature_order=256)
    minus = bound_basin_conditional(model, s - step, force, quadrature_order=256)
    finite_difference = (plus.F_eff - minus.F_eff) / (2.0 * step)
    assert abs(center.boundary_correction_s) > 1.0e-6
    assert finite_difference == pytest.approx(center.dF_eff_ds, rel=2.0e-5, abs=2.0e-8)
    assert abs(finite_difference - center.mean_dG_ds) > 1.0e-10


def test_conditional_normalization_quadrature_and_additive_constant() -> None:
    model = _model()
    values = [
        conditional_fast_a(
            model, 0.11, 0.8, a_lower=0.68, a_upper=1.20,
            quadrature_order=order,
        )
        for order in (64, 128, 256)
    ]
    for value in values:
        assert np.sum(value.probability_weights) == pytest.approx(1.0, abs=2.0e-14)
    assert abs(values[2].F_eff - values[1].F_eff) < abs(
        values[1].F_eff - values[0].F_eff
    )
    assert abs(values[2].mean_a - values[1].mean_a) < abs(
        values[1].mean_a - values[0].mean_a
    )

    shifted = conditional_fast_a(
        model, 0.11, 0.8, a_lower=0.68, a_upper=1.20,
        quadrature_order=256, energy_offset=3.25,
    )
    assert shifted.F_eff - values[2].F_eff == pytest.approx(3.25, abs=2.0e-12)
    np.testing.assert_allclose(
        shifted.probability_weights, values[2].probability_weights,
        rtol=2.0e-13, atol=2.0e-15,
    )
    assert shifted.dF_eff_ds == pytest.approx(values[2].dF_eff_ds, abs=2.0e-12)


def test_low_temperature_conditional_converges_to_a_star_and_laplace_limit() -> None:
    errors_a = []
    errors_F = []
    for temperature in (0.02, 0.01, 0.005):
        model = _model(kT=temperature)
        exact = conditional_fast_a(
            model, 0.08, 0.6, a_lower=0.68, a_upper=1.20,
            quadrature_order=256,
        )
        laplace = laplace_F_eff(model, 0.08, 0.6)
        errors_a.append(abs(exact.mean_a - laplace["a_star"]))
        errors_F.append(abs(exact.F_eff - laplace["F_eff_laplace"]))
    assert errors_a[2] < errors_a[1] < errors_a[0]
    assert errors_F[2] < errors_F[1] < errors_F[0]


def test_s_periodicity_and_generalized_load_tilt() -> None:
    model = _model()
    s = 0.12
    force = 0.9
    first = conditional_fast_a(
        model, s, force, a_lower=0.68, a_upper=1.20, quadrature_order=192
    )
    shifted = conditional_fast_a(
        model, s + model.p.b, force, a_lower=0.68, a_upper=1.20,
        quadrature_order=192,
    )
    expected_tilt = -force * model.p.chi_axial_projection * model.p.b
    assert shifted.F_eff - first.F_eff == pytest.approx(expected_tilt, abs=2.0e-12)
    assert shifted.mean_a == pytest.approx(first.mean_a, abs=2.0e-13)
    assert shifted.dF_eff_ds == pytest.approx(first.dF_eff_ds, abs=2.0e-12)


def test_qsd_is_normalized_and_not_identical_to_truncated_gibbs() -> None:
    model = _model()
    qsd = fast_a_qsd(
        model, s=0.18, force=3.5, a_lower=0.68, a_upper=1.40, n_a=121
    )
    da = float(qsd.a[1] - qsd.a[0])
    assert np.sum(qsd.density) * da == pytest.approx(1.0, abs=2.0e-12)
    assert np.sum(qsd.truncated_gibbs_density) * da == pytest.approx(
        1.0, abs=2.0e-12
    )
    assert qsd.escape_rate > 0.0
    assert qsd.l1_from_truncated_gibbs > 1.0e-3


def test_reduced_equilibrium_strain_tangent_and_decomposition() -> None:
    model = _model()
    params = ReducedFastAParams(
        n_s=81, s_wells=1, a_upper=1.20, quadrature_order=128
    )
    reduced_stress = 1.0e-5
    minus = reduced_equilibrium(
        model, float(model.force_from_sigma_over_E(-reduced_stress)), params=params
    )
    plus = reduced_equilibrium(
        model, float(model.force_from_sigma_over_E(reduced_stress)), params=params
    )
    tangent = (float(plus["strain"]) - float(minus["strain"])) / (
        2.0 * reduced_stress
    )
    assert tangent == pytest.approx(1.0, rel=0.12)
    for result in (minus, plus):
        assert float(result["strain"]) == pytest.approx(
            float(result["normal_strain"])
            + float(result["intrawell_strain"])
            + float(result["plastic_strain"]),
            abs=2.0e-14,
        )
        assert float(result["plastic_strain"]) == pytest.approx(0.0, abs=1.0e-14)


def test_reduced_s_pde_preserves_stationary_mass_and_positivity() -> None:
    model = _model()
    params = ReducedFastAParams(
        n_s=41,
        s_wells=1,
        a_upper=1.20,
        quadrature_order=64,
        max_dt=0.005,
        record_interval=0.01,
    )
    result = run_reduced_fast_a(
        model,
        params=params,
        load=CyclicLoad2D(force_min=0.0, force_max=0.0, period=1.0, cycles=0.05),
        preload_force=0.0,
    )
    equilibrium = reduced_equilibrium(model, 0.0, params=params)
    assert np.max(np.abs(result["mass_balance_residual"])) < 1.0e-11
    assert np.min(result["minimum_density"]) >= -2.0e-12
    assert float(result["probability_mass"][-1]) == pytest.approx(1.0, abs=1.0e-11)
    np.testing.assert_allclose(
        result["density_s"], equilibrium["density_s"], rtol=2.0e-11, atol=2.0e-13
    )
    np.testing.assert_allclose(
        result["strain"],
        result["normal_strain"]
        + result["intrawell_strain"]
        + result["plastic_strain"],
        rtol=0.0,
        atol=2.0e-13,
    )


def test_full_2d_marginal_converges_toward_reduction_as_Ma_increases() -> None:
    model = _model()
    model._build_opening_table()
    load = full_pde.cyclic_load_from_sigma_over_E(
        model,
        sigma_over_E_min=-900.0 / 69000.0,
        sigma_over_E_max=-100.0 / 69000.0,
        period=0.04,
        cycles=1.0,
    )
    preload = float(model.force_from_sigma_over_E(-500.0 / 69000.0))
    reduced_params = ReducedFastAParams(
        n_s=31,
        s_wells=1,
        a_upper=1.10,
        quadrature_order=64,
        max_dt=0.001,
        record_interval=0.001,
    )
    reduced = run_reduced_fast_a(
        model, params=reduced_params, load=load, preload_force=preload
    )
    grid_params = full_pde.Grid2DParams(
        n_a=41, n_s=31, s_wells=1, a_upper=1.10
    )
    time_params = full_pde.PDETimeParams(
        max_dt=0.001,
        cfl=0.40,
        record_interval=0.001,
        integrator="implicit",
    )
    errors = []
    try:
        for mobility_a in (1.0, 10.0):
            model.p.mobility_a = mobility_a
            full = full_pde.run_probability_pde_2d(
                prepared_model=model,
                grid_params=grid_params,
                time_params=time_params,
                load=load,
                preload_force=preload,
            )
            marginal = np.sum(full["density"], axis=0) * full["grid"].da
            marginal /= float(full["survival"][-1])
            marginal_l1 = float(
                np.sum(np.abs(marginal - reduced["density_s"])) * full["grid"].ds
            )
            strain_error = float(
                np.max(np.abs(np.asarray(full["strain"]) - reduced["strain"]))
            )
            errors.append((marginal_l1, strain_error))
    finally:
        model.p.mobility_a = 1.0
    assert errors[1][0] < errors[0][0]
    assert errors[1][1] < errors[0][1]
