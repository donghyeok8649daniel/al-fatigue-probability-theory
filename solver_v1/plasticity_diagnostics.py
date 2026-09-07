"""Convergence-based diagnostics for configurational well transfer.

These diagnostics never infer plasticity from a merely nonzero floating-point
``plastic_strain``.  Resolved interwell transfer requires converged well
population and boundary-flux evidence; residual plasticity additionally
requires persistence after unloading and a zero-load hold.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np


UNRESOLVED = "numerical_or_unresolved"
RESOLVED_INTERWELL = "resolved_interwell_transfer"
REVERSIBLE_OR_THERMAL = "reversible_or_thermal"
PERSISTENT_RESIDUAL = "persistent_residual_candidate"


def _array(result: Mapping[str, object], key: str) -> np.ndarray:
    return np.asarray(result[key], dtype=float)


def _max_abs(result: Mapping[str, object], key: str) -> float:
    values = _array(result, key)
    return float(np.nanmax(np.abs(values))) if values.size else 0.0


def _outside_well_mass(result: Mapping[str, object]) -> np.ndarray:
    populations = _array(result, "well_populations")
    wells = np.asarray(result["well_indices"], dtype=int)
    return np.sum(populations[:, wells != 0], axis=1)


@dataclass(frozen=True)
class PlasticityAssessment:
    maximum_plastic_strain: float
    final_plastic_strain: float
    maximum_mean_well_index: float
    maximum_outside_well_mass: float
    final_outside_well_mass: float
    cumulative_net_transfer: float
    cumulative_gross_transfer: float
    plastic_strain_floor: float
    well_mass_floor: float
    interwell_transfer_floor: float
    well_balance_residual: float
    classification: str


def single_run_plasticity_lower_bound(
    result: Mapping[str, object],
) -> dict[str, float | str]:
    """Return an honest lower bound when no refinement ensemble is available."""

    mean_index = _array(result, "mean_well_index")
    plastic = _array(result, "plastic_strain")
    nonzero = np.abs(mean_index) > 10.0 * np.finfo(float).tiny
    coefficient = (
        float(np.median(plastic[nonzero] / mean_index[nonzero]))
        if np.any(nonzero)
        else 1.0
    )
    well_balance = _max_abs(result, "well_population_balance_residual")
    numerical_repair = _max_abs(result, "cumulative_negative_mass_correction")
    identity_error = float(np.max(np.abs(plastic - coefficient * mean_index)))
    floor = max(
        abs(coefficient) * well_balance,
        abs(coefficient) * numerical_repair,
        identity_error,
    )
    return {
        "plastic_signal_floor": float(floor),
        "plastic_signal_maximum": _max_abs(result, "plastic_strain"),
        "plastic_resolution_status": "requires_convergence",
        "plastic_floor_scope": "single-run lower bound",
    }


def assess_plasticity_convergence(
    results: Sequence[Mapping[str, object]],
) -> PlasticityAssessment:
    """Classify the finest result using observed discretization differences.

    The last item is the reference/final-resolution run.  No absolute empirical
    cutoff is used: all floors are built from grid/time/integrator variation,
    well-balance residual, numerical probability repair, and the exact
    relation between mean well index and the reported plastic strain.
    """

    if len(results) < 2:
        raise ValueError("plasticity convergence requires at least two results")
    reference = results[-1]
    plastic_metrics = np.asarray(
        [_max_abs(result, "plastic_strain") for result in results], dtype=float
    )
    outside_metrics = np.asarray(
        [float(np.max(_outside_well_mass(result))) for result in results],
        dtype=float,
    )
    net_metrics = np.asarray(
        [
            _max_abs(result, "cumulative_interwell_net_transfer")
            for result in results
        ],
        dtype=float,
    )
    gross_metrics = np.asarray(
        [
            _max_abs(result, "cumulative_interwell_gross_transfer")
            for result in results
        ],
        dtype=float,
    )

    well_balance = _max_abs(reference, "well_population_balance_residual")
    numerical_repair = _max_abs(reference, "cumulative_negative_mass_correction")
    mean_index = _array(reference, "mean_well_index")
    plastic = _array(reference, "plastic_strain")
    nonzero = np.abs(mean_index) > 10.0 * np.finfo(float).tiny
    if np.any(nonzero):
        coefficient = float(np.median(plastic[nonzero] / mean_index[nonzero]))
    else:
        coefficient = 1.0
    identity_error = float(np.max(np.abs(plastic - coefficient * mean_index)))

    plastic_floor = max(
        float(np.max(np.abs(plastic_metrics[:-1] - plastic_metrics[-1]))),
        abs(coefficient) * well_balance,
        abs(coefficient) * numerical_repair,
        identity_error,
    )
    well_mass_floor = max(
        float(np.max(np.abs(outside_metrics[:-1] - outside_metrics[-1]))),
        well_balance,
        numerical_repair,
    )
    transfer_floor = max(
        float(np.max(np.abs(net_metrics[:-1] - net_metrics[-1]))),
        float(np.max(np.abs(gross_metrics[:-1] - gross_metrics[-1]))),
        well_balance,
        numerical_repair,
    )

    max_plastic = float(plastic_metrics[-1])
    outside = _outside_well_mass(reference)
    max_outside = float(np.max(outside))
    gross = float(gross_metrics[-1])
    net = float(net_metrics[-1])
    resolved = (
        max_plastic > plastic_floor
        and max_outside > well_mass_floor
        and gross > transfer_floor
    )
    return PlasticityAssessment(
        maximum_plastic_strain=max_plastic,
        final_plastic_strain=float(plastic[-1]),
        maximum_mean_well_index=_max_abs(reference, "mean_well_index"),
        maximum_outside_well_mass=max_outside,
        final_outside_well_mass=float(outside[-1]),
        cumulative_net_transfer=net,
        cumulative_gross_transfer=gross,
        plastic_strain_floor=float(plastic_floor),
        well_mass_floor=float(well_mass_floor),
        interwell_transfer_floor=float(transfer_floor),
        well_balance_residual=float(well_balance),
        classification=RESOLVED_INTERWELL if resolved else UNRESOLVED,
    )


def classify_zero_load_persistence(
    assessment: PlasticityAssessment,
    *,
    end_of_loading_plastic_strain: float,
    end_of_hold_plastic_strain: float,
) -> str:
    """Classify persistence without calling a transient crossing residual."""

    floor = assessment.plastic_strain_floor
    if assessment.classification != RESOLVED_INTERWELL:
        return UNRESOLVED
    if abs(float(end_of_hold_plastic_strain)) <= floor:
        return REVERSIBLE_OR_THERMAL
    if abs(float(end_of_loading_plastic_strain)) <= floor:
        return UNRESOLVED
    return PERSISTENT_RESIDUAL
