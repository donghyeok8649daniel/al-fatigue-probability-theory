import inspect
import math
import csv
from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import root

from solver_v1 import aluminum_calibration as calibration
from solver_v1.aluminum_calibration import (
    EffectiveHybridParameters,
    barrier_ratio,
    calibrate_linear_hybrid,
    calibrate_square_root_hybrid,
    cell_ev_to_energy_density,
    energy_density_to_cell_ev,
    evaluate_extended_observables,
    fcc_111_geometry,
    mishin_lu_wei_targets,
    normal_curvature_target_ev,
    raw_to_effective_embedding,
    sensitivity_matrix,
    model_from_extended_parameters,
)
from solver_v1.analytic_lj_eam import SquareRootLinearEmbedding
from solver_v1.analytic_lj_eam import AnalyticLJEAM
from solver_v1.exponential_density import ExponentialDensityParams
from solver_v1.model import ModelParams, TwoRowLJ


@pytest.fixture(scope="module")
def linear_fit():
    return calibrate_linear_hybrid(max_nfev=500)


def test_atomistic_unit_conversion_and_crystallographic_area() -> None:
    geometry = fcc_111_geometry(4.05)
    assert geometry.burgers_angstrom == pytest.approx(4.05 / np.sqrt(2.0))
    assert geometry.plane_spacing_angstrom == pytest.approx(4.05 / np.sqrt(3.0))
    assert geometry.atomic_cell_area_angstrom2 == pytest.approx(
        np.sqrt(3.0) * 4.05**2 / 4.0
    )
    assert geometry.atomic_cell_area_angstrom2 == pytest.approx(7.102490842787127)
    value = energy_density_to_cell_ev(0.250, geometry.atomic_cell_area_angstrom2)
    assert value == pytest.approx(0.110825653154936)
    assert cell_ev_to_energy_density(
        value, geometry.atomic_cell_area_angstrom2
    ) == pytest.approx(0.250)


def test_elastic_curvature_mapping_uses_fixed_lateral_111_modulus() -> None:
    geometry, targets = mishin_lu_wei_targets()
    c111 = (114.0 + 2.0 * 62.0 + 4.0 * 32.0) / 3.0
    assert normal_curvature_target_ev(geometry, c111) == pytest.approx(
        targets.normal_curvature_ev
    )
    assert targets.normal_curvature_ev == pytest.approx(18.96905843, rel=2e-10)


def test_raw_density_embedding_gauge_collapses_to_same_two_combinations() -> None:
    first = raw_to_effective_embedding(
        amplitude=0.5, rho0=2.0, beta_rho=3.0,
        r_e_over_b=1.5, rho_ref=4.0,
    )
    second = raw_to_effective_embedding(
        amplitude=0.5, rho0=2.0 * math.e, beta_rho=2.0,
        r_e_over_b=1.0, rho_ref=4.0,
    )
    np.testing.assert_allclose(first, second, rtol=2.0e-15, atol=0.0)


def test_square_root_sensitivity_has_four_identifiable_columns() -> None:
    _, targets = mishin_lu_wei_targets()
    initial = EffectiveHybridParameters(0.20, 0.82, 3.0, 0.30)
    matrix = sensitivity_matrix(np.log(initial.as_array()), targets)
    assert matrix.shape == (5, 4)
    assert np.all(np.isfinite(matrix))
    assert np.linalg.matrix_rank(matrix) == 4


def test_baseline_square_root_fit_quantifies_structural_failure() -> None:
    fit = calibrate_square_root_hybrid(max_nfev=300)
    assert fit.rank == 4
    assert fit.cost > 100.0
    assert abs(fit.normalized_residuals[1]) > 4.0
    assert abs(fit.normalized_residuals[4]) > 10.0
    assert fit.observables.min_hessian_eigenvalue_ev > 0.0


def test_minimal_linear_extension_is_reproducible_and_stable(linear_fit) -> None:
    repeated = calibrate_linear_hybrid(max_nfev=500)
    np.testing.assert_allclose(
        linear_fit.parameters.as_array(), repeated.parameters.as_array(),
        rtol=2.0e-7, atol=2.0e-9,
    )
    assert linear_fit.cost < 0.003
    assert np.max(np.abs(linear_fit.normalized_residuals)) < 0.05
    assert linear_fit.observables.min_hessian_eigenvalue_ev > 0.0
    strict = evaluate_extended_observables(linear_fit.parameters)
    assert strict == linear_fit.observables


def test_best_feasible_hybrid_kappa_has_signed_unit_tangent(linear_fit) -> None:
    model = model_from_extended_parameters(linear_fit.parameters, chi=0.31)
    hessian = model.local_hessian(model.a0, 0.0)
    direction = np.array([1.0, model.p.chi_axial_projection])
    for reduced_stress in (-1.0e-6, 1.0e-6):
        force = model.force_from_sigma_over_E(reduced_stress)
        initial = np.array([model.a0, 0.0]) + np.linalg.solve(
            hessian, force * direction
        )
        stationary = root(
            lambda q: np.array(
                [model.local_deda(q[0], q[1]), model.local_deds(q[0], q[1])]
            ) - force * direction,
            initial,
            options={"xtol": 1.0e-11},
        )
        strain = (
            stationary.x[0] - model.a0
            + model.p.chi_axial_projection * stationary.x[1]
        ) / model.a0
        assert stationary.success
        assert strain / reduced_stress == pytest.approx(1.0, rel=3.0e-4)


def test_linear_extension_reports_ill_conditioning(linear_fit) -> None:
    assert linear_fit.rank == 5
    assert linear_fit.condition_number > 1.0e7
    assert linear_fit.singular_values[-1] < 1.0e-4


def test_linear_embedding_derivatives_are_analytic() -> None:
    embedding = SquareRootLinearEmbedding(amplitude=2.3, linear=0.7, rho_ref=4.2)
    rho = np.array([0.7, 2.0, 5.0])
    h = 2.0e-5
    first_fd = (embedding.value(rho + h) - embedding.value(rho - h)) / (2.0 * h)
    second_fd = (
        embedding.value(rho + h) - 2.0 * embedding.value(rho)
        + embedding.value(rho - h)
    ) / h**2
    np.testing.assert_allclose(
        embedding.first_derivative(rho), first_fd, rtol=2.0e-9, atol=2.0e-10
    )
    np.testing.assert_allclose(
        embedding.second_derivative(rho), second_fd, rtol=2.0e-5, atol=2.0e-6
    )


def test_linear_embedding_is_not_discarded_when_sqrt_amplitude_is_zero() -> None:
    base = ModelParams(n_cells=1)
    model = AnalyticLJEAM(
        base,
        density_params=ExponentialDensityParams(),
        embedding=SquareRootLinearEmbedding(
            amplitude=0.0, linear=1.0e-3, rho_ref=1.0
        ),
    )
    assert not model.embedding_is_zero
    assert abs(float(model.embedding_terms(model.a0, 0.0).embedding_energy)) > 0.0


def test_held_out_observables_are_not_silently_in_fit_loss() -> None:
    _, targets = mishin_lu_wei_targets()
    assert "intrinsic_stacking_fault" not in targets.names
    assert "vacancy_formation" not in targets.names
    assert len(targets.values()) == len(targets.scales) == 5


def test_barrier_ratio_requires_two_positive_barriers() -> None:
    assert barrier_ratio(0.2, 0.4) == pytest.approx(0.5)
    assert math.isnan(barrier_ratio(None, 0.4))
    assert math.isnan(barrier_ratio(0.2, 0.0))


def test_original_lj_reference_is_unchanged() -> None:
    model = TwoRowLJ(
        ModelParams(
            n_cells=1, kT=0.02, mobility_a=1.0, mobility_s=0.05,
            chi_axial_projection=0.20,
        )
    )
    assert model.a0 == pytest.approx(0.7713438268704838, abs=2.0e-15)
    assert model.sigma_over_E_force_scale() == pytest.approx(86.29296488740997)


def test_specimen_correlation_area_cannot_enter_calibration_api() -> None:
    public = (
        calibration.fcc_111_geometry,
        calibration.energy_density_to_cell_ev,
        calibration.model_from_effective_parameters,
        calibration.calibrate_square_root_hybrid,
        calibration.calibrate_linear_hybrid,
    )
    for function in public:
        assert "A_c" not in inspect.signature(function).parameters


def test_committed_fit_table_matches_strict_model_evaluation(linear_fit) -> None:
    path = (
        Path(__file__).resolve().parents[1]
        / "results" / "aluminum_calibration" / "fit_residuals.csv"
    )
    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    fitted = {
        row["target"]: row for row in rows
        if row["model"] == "linear_best_feasible" and row["used_for_fit"] == "yes"
    }
    actual = linear_fit.observables
    expected = {
        "a0_over_b": actual.a0_over_b,
        "cohesive_energy_eV_atom": actual.cohesive_energy_ev_per_atom,
        "normal_curvature_eV": actual.normal_curvature_ev,
        "registry_barrier_eV_cell": actual.registry_barrier_ev_per_cell,
        "work_separation_eV_cell": actual.work_separation_ev_per_cell,
    }
    assert set(fitted) == set(expected)
    for name, value in expected.items():
        assert float(fitted[name]["prediction"]) == pytest.approx(value, rel=2.0e-10)
