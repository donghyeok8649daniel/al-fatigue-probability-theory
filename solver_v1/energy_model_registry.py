"""Explicit energy-surface choices for the N=1 probability solver.

The registry prevents a dimensionless TwoRowLJ calculation from being labeled
as an Al-target hybrid result.  It changes only the selected energy surface;
the probability equation, mobilities, strain bridge, and opening bookkeeping
remain common.
"""
from __future__ import annotations

from dataclasses import dataclass

from .aluminum_calibration import (
    EffectiveHybridParameters,
    ExtendedHybridParameters,
    model_from_effective_parameters,
    model_from_extended_parameters,
)
from .analytic_lj_eam import AnalyticLJEAM
from .model import ModelParams, TwoRowLJ


TWO_ROW_LJ_REFERENCE = "two_row_lj_reference"
ANALYTIC_LJ_EAM_HYPOTHETICAL = "analytic_lj_eam_hypothetical"
AL_TARGET_BEST_FEASIBLE = "al_target_best_feasible_hybrid"
ENERGY_MODEL_IDS = (
    TWO_ROW_LJ_REFERENCE,
    ANALYTIC_LJ_EAM_HYPOTHETICAL,
    AL_TARGET_BEST_FEASIBLE,
)


@dataclass(frozen=True)
class EnergyModelMetadata:
    model_id: str
    display_key: str
    python_class: str
    family: str
    calibration_status: str
    parameter_source: str
    energy_scale: str
    geometry_scale: str


_METADATA = {
    TWO_ROW_LJ_REFERENCE: EnergyModelMetadata(
        model_id=TWO_ROW_LJ_REFERENCE,
        display_key="energy_model.lj_reference",
        python_class="solver_v1.model.TwoRowLJ",
        family="TwoRowLJ",
        calibration_status="dimensionless reference; not calibrated Al",
        parameter_source="app canonical N=1 reference parameters",
        energy_scale="dimensionless epsilon_LJ=1",
        geometry_scale="dimensionless b=1",
    ),
    ANALYTIC_LJ_EAM_HYPOTHETICAL: EnergyModelMetadata(
        model_id=ANALYTIC_LJ_EAM_HYPOTHETICAL,
        display_key="energy_model.hybrid_hypothetical",
        python_class="solver_v1.analytic_lj_eam.AnalyticLJEAM",
        family="analytic LJ-EAM hybrid",
        calibration_status="hypothetical sensitivity set; uncalibrated",
        parameter_source="results/aluminum_calibration/parameter_sets.csv:hypothetical_v1",
        energy_scale="eV per reduced two-row cell",
        geometry_scale="reduced b=1",
    ),
    AL_TARGET_BEST_FEASIBLE: EnergyModelMetadata(
        model_id=AL_TARGET_BEST_FEASIBLE,
        display_key="energy_model.al_best_feasible",
        python_class="solver_v1.analytic_lj_eam.AnalyticLJEAM",
        family="analytic LJ-EAM hybrid",
        calibration_status="Al-target best-feasible; practically non-identifiable",
        parameter_source="results/aluminum_calibration/parameter_sets.csv:linear_best_feasible",
        energy_scale="eV per reduced two-row cell",
        geometry_scale="b=2.863782463805517 angstrom represented as reduced b=1",
    ),
}


HYPOTHETICAL_PARAMETERS = EffectiveHybridParameters(
    epsilon_lj_ev=0.02,
    sigma_lj_over_b=1.0,
    density_decay_b=3.0,
    embedding_scale_ev=0.08,
)


AL_BEST_FEASIBLE_PARAMETERS = ExtendedHybridParameters(
    epsilon_lj_ev=0.1551326727,
    sigma_lj_over_b=0.8524688688,
    density_decay_b=1.6951301716,
    embedding_scale_ev=6.1512658971,
    linear_ev=3.1077175868,
)


def energy_model_metadata(model_id: str) -> EnergyModelMetadata:
    try:
        return _METADATA[str(model_id)]
    except KeyError as exc:
        raise ValueError(f"unknown energy model: {model_id}") from exc


def build_energy_model(
    model_id: str,
    *,
    chi: float = 0.20,
    kT: float = 0.02,
) -> TwoRowLJ | AnalyticLJEAM:
    """Construct one declared energy model with common kinetic parameters."""

    energy_model_metadata(model_id)
    if model_id == TWO_ROW_LJ_REFERENCE:
        return TwoRowLJ(
            ModelParams(
                n_cells=1,
                kT=float(kT),
                mobility_a=1.0,
                mobility_s=0.05,
                chi_axial_projection=float(chi),
            )
        )
    if model_id == ANALYTIC_LJ_EAM_HYPOTHETICAL:
        return model_from_effective_parameters(
            HYPOTHETICAL_PARAMETERS, chi=float(chi), kT_ev=float(kT)
        )
    return model_from_extended_parameters(
        AL_BEST_FEASIBLE_PARAMETERS, chi=float(chi), kT_ev=float(kT)
    )


def energy_model_result_metadata(model_id: str, model: TwoRowLJ) -> dict[str, object]:
    metadata = energy_model_metadata(model_id)
    return {
        "energy_model_id": metadata.model_id,
        "energy_model_display_key": metadata.display_key,
        "energy_model_python_class": metadata.python_class,
        "energy_model_family": metadata.family,
        "energy_model_calibration_status": metadata.calibration_status,
        "energy_model_parameter_source": metadata.parameter_source,
        "energy_model_energy_scale": metadata.energy_scale,
        "energy_model_geometry_scale": metadata.geometry_scale,
        "energy_model_a0": float(model.a0),
        "energy_model_b": float(model.p.b),
        "energy_model_chi": float(model.p.chi_axial_projection),
        "energy_model_kT": float(model.p.kT),
        "energy_model_kappa_axial": float(model.sigma_over_E_force_scale()),
    }
