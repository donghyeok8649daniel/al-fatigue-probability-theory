"""Deterministic UI convergence workflow for first-passage certification."""
from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Callable

from .solver_adapter import UIAnalysisConfig, run_ui_analysis
from .specimen_probability import attach_convergence_certificate
from solver_v1.plasticity_diagnostics import assess_plasticity_convergence


@dataclass(frozen=True)
class ConvergenceCheckReport:
    certified_result: dict[str, object]
    run_labels: tuple[str, ...]
    skipped: tuple[str, ...]


def _finer_odd(value: int) -> int:
    candidate = max(value + 4, int(round(1.2 * value)))
    return candidate if candidate % 2 == 1 else candidate + 1


def run_convergence_check(
    config: UIAnalysisConfig,
    *,
    reference_result: dict[str, object],
    runner: Callable[..., dict[str, object]] = run_ui_analysis,
    stop_requested: Callable[[], bool] | None = None,
) -> ConvergenceCheckReport:
    """Run timestep, integrator-when-practical, and spatial refinements.

    A Preview result is deliberately rejected: Preview alone, or a study based
    on it, must not unlock a physical specimen probability.  The returned
    result is the finer spatial run, not the original reference.
    """

    config.validate()
    if config.analysis_quality != "resolved":
        raise ValueError("Preview results are not convergence-certified")
    if reference_result.get("analysis_quality") != "resolved":
        raise ValueError("a resolved reference result is required")

    common = {"stop_requested": stop_requested} if stop_requested is not None else {}
    time_config = replace(config, max_dt=0.5 * config.max_dt)
    time_result = runner(time_config, **common)
    comparisons = [reference_result, time_result]
    labels = ["reference", "dt/2"]
    skipped: list[str] = []

    # The old explicit SG path can require hundreds of thousands of steps on
    # resolved repulsive-LJ grids. Compare integrators only on tractable grids.
    if config.grid_n_a * config.grid_n_s <= 2500:
        alternate = "implicit" if config.integration_method == "explicit" else "explicit"
        alternate_result = runner(
            replace(config, integration_method=alternate, max_dt=0.5 * config.max_dt),
            **common,
        )
        comparisons.append(alternate_result)
        labels.append(f"{alternate}-integrator")
    else:
        skipped.append(
            "explicit/implicit comparison skipped: explicit CFL cost is impractical "
            "on the resolved grid"
        )

    spatial_config = replace(
        config,
        grid_n_a=_finer_odd(config.grid_n_a),
        grid_n_s=_finer_odd(config.grid_n_s),
        max_dt=0.5 * config.max_dt,
    )
    spatial_result = runner(spatial_config, **common)
    attach_convergence_certificate(spatial_result, comparisons)
    plasticity = assess_plasticity_convergence([*comparisons, spatial_result])
    spatial_result.update(
        {
            "plastic_signal_floor": plasticity.plastic_strain_floor,
            "plastic_well_mass_floor": plasticity.well_mass_floor,
            "plastic_interwell_transfer_floor": (
                plasticity.interwell_transfer_floor
            ),
            "plastic_resolution_status": plasticity.classification,
            "plastic_resolution_certified": True,
            "plastic_floor_scope": "grid/time/integrator convergence",
        }
    )
    spatial_result["convergence_check_runs"] = tuple(
        [*labels, f"finer-grid-{spatial_config.grid_n_a}x{spatial_config.grid_n_s}"]
    )
    spatial_result["convergence_check_skipped"] = tuple(skipped)
    return ConvergenceCheckReport(
        certified_result=spatial_result,
        run_labels=spatial_result["convergence_check_runs"],
        skipped=tuple(skipped),
    )
