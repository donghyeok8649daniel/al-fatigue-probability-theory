import numpy as np
import pytest

from app.specimen_probability import (
    BELOW_RESOLUTION,
    RESOLVED,
    assess_convergence_rare_event,
)
from solver_v1.model import ModelParams, TwoRowLJ
from solver_v1.plasticity_diagnostics import (
    PERSISTENT_RESIDUAL,
    RESOLVED_INTERWELL,
    UNRESOLVED,
    assess_plasticity_convergence,
    classify_zero_load_persistence,
)
from solver_v1.probability_pde_2d import (
    CyclicLoad2D,
    Grid2DParams,
    PDETimeParams,
    run_probability_pde_2d,
)


def _result(signal: float, *, error: float) -> dict[str, np.ndarray]:
    coefficient = 0.25
    mean = np.array([0.0, signal / coefficient])
    outside = max(abs(signal) / coefficient, 0.0)
    return {
        "well_indices": np.array([-1, 0, 1]),
        "well_populations": np.array(
            [[0.0, 1.0, 0.0], [0.0, 1.0 - outside, outside]]
        ),
        "plastic_strain": np.array([0.0, signal]),
        "mean_well_index": mean,
        "cumulative_interwell_net_transfer": np.array(
            [[0.0, 0.0], [0.0, outside]]
        ),
        "cumulative_interwell_gross_transfer": np.array(
            [[0.0, 0.0], [0.0, outside]]
        ),
        "well_population_balance_residual": np.array([[0.0], [error]]),
        "cumulative_negative_mass_correction": np.array([0.0, error]),
    }


def test_tiny_nonconverged_plastic_signal_is_not_called_physical() -> None:
    assessment = assess_plasticity_convergence(
        [_result(0.0, error=2.0e-11), _result(7.0e-11, error=2.0e-11)]
    )
    assert assessment.maximum_plastic_strain == 7.0e-11
    assert assessment.plastic_strain_floor >= 7.0e-11
    assert assessment.classification == UNRESOLVED


def test_converged_well_population_and_flux_can_resolve_transfer() -> None:
    assessment = assess_plasticity_convergence(
        [_result(1.01e-3, error=1.0e-8), _result(1.0e-3, error=1.0e-8)]
    )
    assert assessment.classification == RESOLVED_INTERWELL
    assert assessment.maximum_outside_well_mass > assessment.well_mass_floor
    assert assessment.cumulative_gross_transfer > assessment.interwell_transfer_floor
    assert (
        classify_zero_load_persistence(
            assessment,
            end_of_loading_plastic_strain=1.0e-3,
            end_of_hold_plastic_strain=8.0e-4,
        )
        == PERSISTENT_RESIDUAL
    )


@pytest.fixture(scope="module")
def model() -> TwoRowLJ:
    prepared = TwoRowLJ(
        ModelParams(
            n_cells=1,
            kT=0.02,
            mobility_a=1.0,
            mobility_s=0.05,
            chi_axial_projection=0.20,
        )
    )
    prepared._build_opening_table()
    return prepared


def _physical_run(
    model: TwoRowLJ,
    n_a: int,
    n_s: int,
    *,
    force: float = 3.0,
    duration: float = 0.1,
    value_function=None,
    s_wells: int = 3,
) -> dict[str, object]:
    return run_probability_pde_2d(
        prepared_model=model,
        grid_params=Grid2DParams(
            n_a=n_a, n_s=n_s, s_wells=s_wells, a_upper=1.6
        ),
        time_params=PDETimeParams(
            max_dt=2.0e-3,
            cfl=0.4,
            record_interval=min(0.02, duration / 2.0),
            integrator="implicit",
        ),
        load=CyclicLoad2D(
            force_min=force,
            force_max=force,
            period=1.0,
            cycles=duration,
            value_function=value_function,
        ),
        preload_force=0.0,
    )


def test_apparent_interwell_transfer_shrinks_under_aligned_s_grid_refinement(
    model: TwoRowLJ,
) -> None:
    results = [
        _physical_run(model, 11, 18),
        _physical_run(model, 15, 24),
        _physical_run(model, 19, 30),
    ]
    assessment = assess_plasticity_convergence(results)
    maxima = np.array(
        [np.max(np.abs(result["plastic_strain"])) for result in results]
    )
    assert np.all(np.diff(maxima) < 0.0)
    assert assessment.classification == UNRESOLVED
    assert assessment.maximum_plastic_strain <= assessment.plastic_strain_floor
    for result in results:
        assert np.max(np.abs(result["well_population_balance_residual"])) < 1.0e-12


def test_zero_load_hold_does_not_promote_unresolved_transfer_to_residual_plasticity(
    model: TwoRowLJ,
) -> None:
    def load_then_hold(time: float) -> float:
        if time < 0.1:
            return float(1.5 * (1.0 - np.cos(20.0 * np.pi * time)))
        return 0.0

    results = [
        _physical_run(
            model, 11, 18, duration=0.2, value_function=load_then_hold
        ),
        _physical_run(
            model, 15, 24, duration=0.2, value_function=load_then_hold
        ),
    ]
    assessment = assess_plasticity_convergence(results)
    reference = results[-1]
    end_of_loading = int(np.argmin(np.abs(reference["time"] - 0.1)))
    classification = classify_zero_load_persistence(
        assessment,
        end_of_loading_plastic_strain=float(reference["plastic_strain"][end_of_loading]),
        end_of_hold_plastic_strain=float(reference["plastic_strain"][-1]),
    )
    assert assessment.classification == UNRESOLVED
    assert classification == UNRESOLVED


def test_zero_load_hold_is_not_created_by_outer_s_domain_boundaries(
    model: TwoRowLJ,
) -> None:
    def load_then_hold(time: float) -> float:
        if time < 0.02:
            return float(1.5 * (1.0 - np.cos(100.0 * np.pi * time)))
        return 0.0

    results = [
        _physical_run(
            model,
            11,
            6 * wells,
            duration=0.04,
            value_function=load_then_hold,
            s_wells=wells,
        )
        for wells in (3, 5, 7)
    ]
    final_plastic = np.asarray(
        [result["plastic_strain"][-1] for result in results], dtype=float
    )
    assert np.ptp(final_plastic) < 1.0e-12
    assert max(abs(final_plastic)) < 1.0e-9
    assert np.max(results[-1]["s_truncation_boundary_mass"]) < 1.0e-14


def test_compressive_control_absorption_is_below_convergence_floor(
    model: TwoRowLJ,
) -> None:
    results = [
        _physical_run(model, 11, 18, force=-1.0, duration=0.02),
        _physical_run(model, 15, 24, force=-1.0, duration=0.02),
    ]
    assessment = assess_convergence_rare_event(results)
    assert assessment.signal == 0.0
    assert assessment.status == BELOW_RESOLUTION


def test_extreme_control_opening_loss_is_above_observed_numerical_floor(
    model: TwoRowLJ,
) -> None:
    results = [
        _physical_run(model, 11, 18, force=5.1, duration=0.05),
        _physical_run(model, 15, 24, force=5.1, duration=0.05),
        _physical_run(model, 19, 30, force=5.1, duration=0.05),
    ]
    assessment = assess_convergence_rare_event(results)
    assert assessment.signal > 0.1
    assert assessment.signal > assessment.floor
    assert assessment.status == RESOLVED
