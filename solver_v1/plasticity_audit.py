"""Deterministic strain/registry audit helpers.

The helpers summarize existing PDE fields.  They do not alter mobilities,
barriers, temperature, strain, or probability and do not assign physical time.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from .plasticity_diagnostics import PlasticityAssessment, assess_plasticity_convergence


@dataclass(frozen=True)
class PlasticityRunSummary:
    label: str
    grid_a: int
    grid_s: int
    dt: float
    integrator: str
    energy_model: str
    chi: float
    max_epsilon_p: float
    final_epsilon_p: float
    final_epsilon_xi: float
    max_decomposition_residual: float
    max_outside_central_well_mass: float
    final_registry_moment: float
    cumulative_forward_transfer: float
    cumulative_backward_transfer: float
    cumulative_net_registry_transfer: float
    cumulative_gross_registry_activity: float
    opening_absorbed_mass: float
    mass_residual: float
    well_balance_residual: float
    negative_mass_correction: float
    flux_absorption_residual: float
    well_boundary_alignment_error: float

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


def _max_abs(values) -> float:
    array = np.asarray(values, dtype=float)
    return float(np.max(np.abs(array))) if array.size else 0.0


def outside_central_well_mass(result: dict[str, object]) -> np.ndarray:
    populations = np.asarray(result["well_populations"], dtype=float)
    wells = np.asarray(result["well_indices"], dtype=int)
    return np.sum(populations[:, wells != 0], axis=1)


def summarize_plasticity_run(
    result: dict[str, object],
    *,
    label: str,
    dt: float,
    integrator: str,
    energy_model: str,
) -> PlasticityRunSummary:
    """Return machine-readable metrics without interpreting nonzero as plasticity."""

    model = result["model"]
    grid = result["grid"]
    forward = np.asarray(result["cumulative_interwell_forward_transfer"], dtype=float)
    backward = np.asarray(
        result["cumulative_interwell_backward_transfer"], dtype=float
    )
    outside = outside_central_well_mass(result)
    return PlasticityRunSummary(
        label=str(label),
        grid_a=int(grid.a.size),
        grid_s=int(grid.s.size),
        dt=float(dt),
        integrator=str(integrator),
        energy_model=str(energy_model),
        chi=float(model.p.chi_axial_projection),
        max_epsilon_p=_max_abs(result["plastic_strain"]),
        final_epsilon_p=float(np.asarray(result["plastic_strain"])[-1]),
        final_epsilon_xi=float(np.asarray(result["intrawell_strain"])[-1]),
        max_decomposition_residual=_max_abs(
            result["strain_decomposition_residual"]
        ),
        max_outside_central_well_mass=float(np.max(outside)),
        final_registry_moment=float(
            np.asarray(result["unnormalized_registry_moment"])[-1]
        ),
        cumulative_forward_transfer=float(np.sum(forward[-1])),
        cumulative_backward_transfer=float(np.sum(backward[-1])),
        cumulative_net_registry_transfer=float(
            np.asarray(result["accumulated_net_registry_transfer"])[-1]
        ),
        cumulative_gross_registry_activity=float(
            np.sum(np.asarray(result["cumulative_interwell_gross_transfer"])[-1])
        ),
        opening_absorbed_mass=float(
            np.asarray(result["cumulative_absorbed_mass"])[-1]
        ),
        mass_residual=_max_abs(result["mass_balance_residual"]),
        well_balance_residual=_max_abs(result["well_population_balance_residual"]),
        negative_mass_correction=_max_abs(
            result["cumulative_negative_mass_correction"]
        ),
        flux_absorption_residual=_max_abs(result["flux_consistency_residual"]),
        well_boundary_alignment_error=_max_abs(
            result["interwell_boundary_alignment_error"]
        ),
    )


def certify_plasticity_refinement(
    results: list[dict[str, object]],
) -> PlasticityAssessment:
    """Use the existing convergence-based classifier for an ordered run set."""

    return assess_plasticity_convergence(results)


def final_well_populations(result: dict[str, object]) -> dict[int, float]:
    wells = np.asarray(result["well_indices"], dtype=int)
    populations = np.asarray(result["well_populations"], dtype=float)[-1]
    return {int(well): float(value) for well, value in zip(wells, populations)}
