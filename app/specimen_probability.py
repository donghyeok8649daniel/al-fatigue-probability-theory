"""Rare-event diagnostics and specimen-scale probability aggregation.

The areas in this module are statistical/specimen aggregation inputs only.
They never enter the constitutive stress-to-force conversion.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np


BELOW_RESOLUTION = "below_numerical_resolution"
RESOLVED = "resolved"


@dataclass(frozen=True)
class RareEventAssessment:
    signal: float
    floor: float
    mass_residual: float
    accumulated_repair: float
    flux_discrepancy: float
    resolution_change: float
    status: str


def _last(values) -> float:
    array = np.asarray(values, dtype=float).reshape(-1)
    if array.size == 0:
        raise ValueError("probability diagnostic array is empty")
    return float(array[-1])


def assess_local_rare_event(
    result: Mapping[str, object],
    *,
    resolution_change: float = 0.0,
) -> RareEventAssessment:
    """Estimate a data-derived local rare-event floor for one or more runs.

    ``resolution_change`` is supplied by a grid/time/integrator convergence
    comparison.  With zero it is a single-run lower bound, not a convergence
    certificate.
    """

    signal = _last(result["cumulative_absorbed_mass"])
    mass_residual = float(
        np.max(np.abs(np.asarray(result["mass_balance_residual"], dtype=float)))
    )
    repair_key = (
        "cumulative_negative_mass_correction"
        if "cumulative_negative_mass_correction" in result
        else "negative_mass_correction"
    )
    accumulated_repair = float(
        np.max(np.abs(np.asarray(result[repair_key], dtype=float)))
    )
    flux_discrepancy = float(
        np.max(
            np.abs(
                np.asarray(result.get("flux_consistency_residual", 0.0), dtype=float)
            )
        )
    )
    resolution_change = abs(float(resolution_change))
    floor = max(
        mass_residual,
        accumulated_repair,
        flux_discrepancy,
        resolution_change,
    )
    status = RESOLVED if signal > floor else BELOW_RESOLUTION
    return RareEventAssessment(
        signal=signal,
        floor=floor,
        mass_residual=mass_residual,
        accumulated_repair=accumulated_repair,
        flux_discrepancy=flux_discrepancy,
        resolution_change=resolution_change,
        status=status,
    )


def assess_convergence_rare_event(
    results: Sequence[Mapping[str, object]],
) -> RareEventAssessment:
    """Assess the finest/last result including observed resolution variation."""

    if len(results) < 2:
        raise ValueError("rare-event convergence requires at least two results")
    signals = np.asarray(
        [_last(result["cumulative_absorbed_mass"]) for result in results],
        dtype=float,
    )
    reference = float(signals[-1])
    resolution_change = float(np.max(np.abs(signals[:-1] - reference)))
    internal_floors = [assess_local_rare_event(result).floor for result in results]
    assessment = assess_local_rare_event(
        results[-1],
        resolution_change=max(resolution_change, *internal_floors),
    )
    return assessment


def attach_convergence_certificate(
    result: dict[str, object],
    comparison_results: Sequence[Mapping[str, object]],
) -> RareEventAssessment:
    """Attach a grid/time/integrator assessment to an existing UI result.

    ``comparison_results`` must contain at least one independent run; ``result``
    is evaluated as the final/reference resolution.  Certification means the
    convergence study was performed, not that its signal necessarily resolved.
    """

    assessment = assess_convergence_rare_event([*comparison_results, result])
    result.update(
        {
            "local_rare_event_floor": assessment.floor,
            "probability_resolution_status": assessment.status,
            "probability_resolution_certified": True,
            "rare_event_floor_mass_residual": assessment.mass_residual,
            "rare_event_floor_accumulated_repair": assessment.accumulated_repair,
            "rare_event_floor_flux_discrepancy": assessment.flux_discrepancy,
            "rare_event_floor_resolution_change": assessment.resolution_change,
            "rare_event_floor_scope": "grid/time/integrator convergence",
        }
    )
    return assessment


@dataclass(frozen=True)
class SpecimenProbability:
    N_eff: float
    local_initiation_probability: np.ndarray
    local_survival_probability: np.ndarray
    specimen_initiation_probability: np.ndarray
    specimen_survival_probability: np.ndarray
    mathematical_extrapolation: np.ndarray
    resolved_mask: np.ndarray
    status: str


def aggregate_specimen_probability(
    local_initiation_probability,
    *,
    correlation_area_mm2: float,
    stressed_area_mm2: float,
    local_numerical_floor: float,
    resolution_certified: bool = False,
) -> SpecimenProbability:
    """Aggregate independent equivalent regions using stable log arithmetic.

    Values at or below the local numerical floor are returned as NaN in the
    physical specimen arrays.  Their mathematical extrapolation remains
    available separately and must be labeled unresolved by callers.
    """

    correlation_area = float(correlation_area_mm2)
    stressed_area = float(stressed_area_mm2)
    floor = float(local_numerical_floor)
    if not np.isfinite(correlation_area) or correlation_area <= 0.0:
        raise ValueError("correlation_area_mm2 must be finite and positive")
    if not np.isfinite(stressed_area) or stressed_area < 0.0:
        raise ValueError("stressed_area_mm2 must be finite and nonnegative")
    if not np.isfinite(floor) or floor < 0.0:
        raise ValueError("local_numerical_floor must be finite and nonnegative")

    local = np.atleast_1d(np.asarray(local_initiation_probability, dtype=float))
    if np.any(~np.isfinite(local)) or np.any((local < 0.0) | (local > 1.0)):
        raise ValueError("local initiation probability must lie in [0, 1]")
    N_eff = stressed_area / correlation_area
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        log_survival = N_eff * np.log1p(-local)
        extrapolation = -np.expm1(log_survival)
    if N_eff == 0.0:
        extrapolation = np.zeros_like(local)
    extrapolation = np.asarray(extrapolation, dtype=float)
    # A one-run residual floor is only a lower bound.  Area amplification is
    # physically exposed only after a grid/time/integrator convergence study
    # has certified the local signal.
    resolved_mask = (local > floor) & bool(resolution_certified)
    specimen_initiation = np.where(resolved_mask, extrapolation, np.nan)
    specimen_survival = np.where(resolved_mask, 1.0 - extrapolation, np.nan)
    status = RESOLVED if bool(resolved_mask[-1]) else BELOW_RESOLUTION
    return SpecimenProbability(
        N_eff=float(N_eff),
        local_initiation_probability=local,
        local_survival_probability=1.0 - local,
        specimen_initiation_probability=specimen_initiation,
        specimen_survival_probability=specimen_survival,
        mathematical_extrapolation=extrapolation,
        resolved_mask=resolved_mask,
        status=status,
    )
