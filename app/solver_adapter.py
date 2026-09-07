"""Testable desktop adapter for the canonical N=1 probability PDE.

This module is the sole source of probability and axial-strain histories used
by the desktop UI. It performs no trajectory sampling and applies no empirical
fatigue-life distribution.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from solver_v1.dynamics_diagnostics import model_frequency_diagnostics
from solver_v1.model import ModelParams, TwoRowLJ
from solver_v1.probability_pde_2d import (
    Grid2D,
    Grid2DParams,
    PDETimeParams,
    cyclic_load_from_sigma_over_E,
    run_probability_pde_2d,
)


PDE_RESULT_FIELDS = (
    "normal_strain",
    "intrawell_strain",
    "plastic_strain",
    "strain",
    "survival",
    "survival_probability",
    "initiation_probability",
    "first_passage_flux",
    "intact_probability_mass",
    "cumulative_absorbed_mass",
    "mass_balance_residual",
    "negative_mass_correction",
    "minimum_density",
    "raw_intact_mass",
    "raw_one_minus_survival",
)


@dataclass(frozen=True)
class UIAnalysisConfig:
    """Physical UI inputs plus declared numerical mechanism controls."""

    young_gpa: float = 69.0
    stress_mean_mpa: float = 50.0
    stress_amplitude_mpa: float = 100.0
    model_frequency: float = 25.0
    cycles: float = 1.0
    steps_per_cycle: int = 40
    grid_n_a: int = 21
    grid_n_s: int = 31
    s_wells: int = 3
    a_upper: float = 1.60
    max_dt: float = 1.0e-3
    integration_method: str = "explicit"
    analysis_quality: str = "preview"

    @property
    def young_mpa(self) -> float:
        return float(self.young_gpa * 1.0e3)

    @property
    def stress_min_mpa(self) -> float:
        return float(self.stress_mean_mpa - self.stress_amplitude_mpa)

    @property
    def stress_max_mpa(self) -> float:
        return float(self.stress_mean_mpa + self.stress_amplitude_mpa)

    @property
    def model_period(self) -> float:
        return 1.0 / float(self.model_frequency)

    def stress_mpa(self, model_time):
        phase = 2.0 * np.pi * np.asarray(model_time, dtype=float) / self.model_period
        value = self.stress_mean_mpa + self.stress_amplitude_mpa * np.sin(phase)
        return float(value) if np.ndim(value) == 0 else value

    def validate(self) -> None:
        positive = {
            "young_gpa": self.young_gpa,
            "model_frequency": self.model_frequency,
            "cycles": self.cycles,
            "max_dt": self.max_dt,
        }
        for name, value in positive.items():
            if not np.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive")
        if not np.isfinite(self.stress_mean_mpa):
            raise ValueError("stress_mean_mpa must be finite")
        if not np.isfinite(self.stress_amplitude_mpa) or self.stress_amplitude_mpa < 0.0:
            raise ValueError("stress_amplitude_mpa must be finite and nonnegative")
        if self.steps_per_cycle < 4:
            raise ValueError("steps_per_cycle must be at least 4")
        if self.grid_n_a < 5 or self.grid_n_s < 5:
            raise ValueError("grid dimensions must both be at least 5")
        if self.s_wells < 1 or self.s_wells % 2 == 0:
            raise ValueError("s_wells must be a positive odd integer")
        if self.integration_method not in {"explicit", "implicit"}:
            raise ValueError("integration_method must be 'explicit' or 'implicit'")
        if self.analysis_quality not in {"preview", "resolved"}:
            raise ValueError("analysis_quality must be 'preview' or 'resolved'")


def canonical_model_params() -> ModelParams:
    """Return the verified N=1 UI mechanism parameters."""

    return ModelParams(
        n_cells=1,
        kT=0.02,
        mobility_a=1.0,
        mobility_s=0.05,
        chi_axial_projection=0.20,
    )


def physical_load_conversion(config: UIAnalysisConfig) -> dict[str, float]:
    """Map physical stress inputs through sigma/E and relaxed axial kappa."""

    config.validate()
    model = TwoRowLJ(canonical_model_params())
    kappa = model.sigma_over_E_force_scale()
    reduced_min = config.stress_min_mpa / config.young_mpa
    reduced_max = config.stress_max_mpa / config.young_mpa
    reduced_initial = config.stress_mean_mpa / config.young_mpa
    dynamics = model_frequency_diagnostics(model, config.model_frequency)
    return {
        "a0": float(model.a0),
        "relaxed_axial_kappa": float(kappa),
        "frozen_normal_kappa": float(
            model.frozen_normal_sigma_over_E_force_scale()
        ),
        "reduced_stress_min": float(reduced_min),
        "reduced_stress_max": float(reduced_max),
        "reduced_stress_initial": float(reduced_initial),
        "force_min": float(model.force_from_sigma_over_E(reduced_min)),
        "force_max": float(model.force_from_sigma_over_E(reduced_max)),
        "preload_force": float(model.force_from_sigma_over_E(reduced_initial)),
        **dynamics,
    }


def physical_probability_bookkeeping(
    raw_intact_mass,
    cumulative_absorbed_mass,
) -> dict[str, np.ndarray]:
    """Separate physical first-passage probability from numerical mass error.

    Physical initiation is the mass actually removed at the opening boundary.
    ``1 - raw_intact_mass`` also contains any conservative-solve roundoff and is
    retained only as an explicitly named numerical diagnostic.  No clipping or
    rounding is applied to either quantity.
    """

    raw_intact = np.asarray(raw_intact_mass, dtype=float)
    absorbed = np.asarray(cumulative_absorbed_mass, dtype=float)
    survival = 1.0 - absorbed
    return {
        "survival": survival,
        "survival_probability": survival,
        "initiation_probability": absorbed,
        "raw_intact_mass": raw_intact,
        "raw_one_minus_survival": 1.0 - raw_intact,
    }


def _decorate_record(
    config: UIAnalysisConfig,
    record: dict[str, float],
) -> dict[str, float]:
    model_time = float(record["time"])
    cycle = model_time / config.model_period
    probability = physical_probability_bookkeeping(
        record["survival"], record["cumulative_absorbed_mass"]
    )
    return {
        **record,
        **{key: float(value) for key, value in probability.items()},
        "model_time": model_time,
        "load_cycle": float(cycle),
        "applied_stress_mpa": float(config.stress_mpa(model_time)),
    }


def result_field_mapping(result: dict[str, object]) -> dict[str, np.ndarray]:
    """Return UI plot arrays without reconstructing any PDE observable."""

    required = ("model_time", "load_cycle", "applied_stress_mpa", *PDE_RESULT_FIELDS)
    missing = [name for name in required if name not in result]
    if missing:
        raise KeyError(f"PDE result is missing required UI fields: {missing}")
    return {
        "stress": np.asarray(result["applied_stress_mpa"], dtype=float),
        "normal_strain": np.asarray(result["normal_strain"], dtype=float),
        "intrawell_strain": np.asarray(result["intrawell_strain"], dtype=float),
        "plastic_strain": np.asarray(result["plastic_strain"], dtype=float),
        "strain": np.asarray(result["strain"], dtype=float),
        "survival": np.asarray(result["survival"], dtype=float),
        "survival_probability": np.asarray(
            result["survival_probability"], dtype=float
        ),
        "initiation_probability": np.asarray(
            result["initiation_probability"], dtype=float
        ),
        "first_passage_flux": np.asarray(result["first_passage_flux"], dtype=float),
        "intact_probability_mass": np.asarray(
            result["intact_probability_mass"], dtype=float
        ),
        "cumulative_absorbed_mass": np.asarray(
            result["cumulative_absorbed_mass"], dtype=float
        ),
        "mass_balance_residual": np.asarray(
            result["mass_balance_residual"], dtype=float
        ),
        "negative_mass_correction": np.asarray(
            result["negative_mass_correction"], dtype=float
        ),
        "minimum_density": np.asarray(result["minimum_density"], dtype=float),
        "raw_intact_mass": np.asarray(result["raw_intact_mass"], dtype=float),
        "raw_one_minus_survival": np.asarray(
            result["raw_one_minus_survival"], dtype=float
        ),
    }


def run_ui_analysis(
    config: UIAnalysisConfig,
    *,
    record_callback: Callable[[dict[str, float]], None] | None = None,
    density_record_callback: (
        Callable[[float, float, np.ndarray, TwoRowLJ, Grid2D], None] | None
    ) = None,
    stop_requested: Callable[[], bool] | None = None,
) -> dict[str, object]:
    """Run the canonical PDE path used by the desktop application."""

    config.validate()
    model_params = canonical_model_params()
    calibration_model = TwoRowLJ(model_params)
    load = cyclic_load_from_sigma_over_E(
        calibration_model,
        sigma_over_E_min=config.stress_min_mpa / config.young_mpa,
        sigma_over_E_max=config.stress_max_mpa / config.young_mpa,
        period=config.model_period,
        cycles=config.cycles,
    )
    preload_force = float(
        calibration_model.force_from_sigma_over_E(
            config.stress_mean_mpa / config.young_mpa
        )
    )

    def forward(record: dict[str, float]) -> None:
        if record_callback is not None:
            record_callback(_decorate_record(config, record))

    raw = run_probability_pde_2d(
        model_params=model_params,
        grid_params=Grid2DParams(
            n_a=config.grid_n_a,
            n_s=config.grid_n_s,
            s_wells=config.s_wells,
            a_upper=config.a_upper,
        ),
        time_params=PDETimeParams(
            max_dt=min(config.max_dt, config.model_period / config.steps_per_cycle),
            cfl=0.40,
            record_interval=config.model_period / config.steps_per_cycle,
            integrator=config.integration_method,
        ),
        load=load,
        preload_force=preload_force,
        record_callback=forward,
        density_record_callback=density_record_callback,
        stop_requested=stop_requested,
    )
    model_time = np.asarray(raw["time"], dtype=float)
    load_cycle = model_time / config.model_period
    dynamics = model_frequency_diagnostics(
        calibration_model, config.model_frequency
    )
    probability = physical_probability_bookkeeping(
        raw["survival"], raw["cumulative_absorbed_mass"]
    )
    result: dict[str, object] = {
        **raw,
        **probability,
        "model_time": model_time,
        "load_cycle": load_cycle,
        "applied_stress_mpa": np.asarray(config.stress_mpa(model_time), dtype=float),
        "initial_stress_mpa": float(config.stress_mean_mpa),
        "initial_force": preload_force,
        "initial_condition": "conditional Gibbs at sigma(t=0)",
        "relaxed_axial_kappa": float(calibration_model.sigma_over_E_force_scale()),
        "frozen_normal_kappa": float(
            calibration_model.frozen_normal_sigma_over_E_force_scale()
        ),
        "solver_time_status": "dimensionless model time; not calibrated seconds",
        "probability_source": "N=1 direct Smoluchowski/Fokker-Planck PDE",
        "analysis_quality": config.analysis_quality,
        "integration_method": config.integration_method,
        "grid_shape": (config.grid_n_a, config.grid_n_s),
        **dynamics,
    }
    result_field_mapping(result)
    return result
