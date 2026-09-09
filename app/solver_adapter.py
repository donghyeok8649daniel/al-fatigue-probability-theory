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
from solver_v1.configurational_landscape import bound_configurational_barrier
from solver_v1.energy_model_registry import (
    ENERGY_MODEL_IDS,
    TWO_ROW_LJ_REFERENCE,
    build_energy_model,
    energy_model_result_metadata,
)
from solver_v1.model import ModelParams, TwoRowLJ
from solver_v1.probability_pde_2d import (
    Grid2D,
    Grid2DParams,
    PDETimeParams,
    cyclic_load_from_sigma_over_E,
    run_probability_pde_2d,
)
from solver_v1.plasticity_diagnostics import single_run_plasticity_lower_bound
from solver_v1.physical_time import (
    PhysicalTimeCalibration,
    hz_to_model_frequency,
    model_time_to_seconds,
)
from .specimen_probability import assess_local_rare_event
from solver_v1.kinetic_calibration_workflow import (
    build_time_basis_model, validate_for_energy_model,
)


PDE_RESULT_FIELDS = (
    "normal_strain",
    "intrawell_strain",
    "plastic_strain",
    "strain_decomposition_residual",
    "strain",
    "survival",
    "survival_probability",
    "local_survival_probability",
    "initiation_probability",
    "local_initiation_probability",
    "first_passage_flux",
    "intact_probability_mass",
    "cumulative_absorbed_mass",
    "absorbed_mass_increment",
    "initial_absorbed_mass",
    "integrated_first_passage_flux",
    "flux_consistency_residual",
    "mass_balance_residual",
    "negative_mass_correction",
    "cumulative_negative_mass_correction",
    "minimum_density",
    "raw_intact_mass",
    "raw_one_minus_survival",
    "mean_well_index",
    "plastic_well_activity",
    "unnormalized_registry_moment",
    "accumulated_net_registry_transfer",
    "absorbed_registry_moment",
    "registry_moment_balance_residual",
    "net_interwell_registry_rate",
    "gross_interwell_activity_rate",
    "net_plastic_flow_rate",
    "gross_configurational_slip_activity",
    "selective_opening_plastic_rate",
    "cumulative_forward_registry_activity",
    "cumulative_backward_registry_activity",
    "cumulative_gross_registry_activity",
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
    time_basis: str = "model"
    physical_frequency_hz: float | None = None
    time_calibration: PhysicalTimeCalibration | None = None
    energy_model: str = TWO_ROW_LJ_REFERENCE

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
        return 1.0 / self.effective_model_frequency

    @property
    def effective_model_frequency(self) -> float:
        if self.time_basis == "model":
            return float(self.model_frequency)
        if self.time_basis != "physical":
            raise ValueError("time_basis must be 'model' or 'physical'")
        if self.time_calibration is None or self.physical_frequency_hz is None:
            raise ValueError("physical time requires a kinetic calibration and Hz input")
        return float(
            hz_to_model_frequency(self.physical_frequency_hz, self.time_calibration)
        )

    @property
    def frequency_hz(self) -> float | None:
        if self.time_basis != "physical":
            return None
        self.validate_time_basis()
        return float(self.physical_frequency_hz)

    @property
    def physical_period_seconds(self) -> float | None:
        frequency = self.frequency_hz
        return None if frequency is None else 1.0 / frequency

    @property
    def physical_duration_seconds(self) -> float | None:
        period = self.physical_period_seconds
        return None if period is None else self.cycles * period

    def stress_mpa(self, model_time):
        phase = 2.0 * np.pi * np.asarray(model_time, dtype=float) / self.model_period
        value = self.stress_mean_mpa + self.stress_amplitude_mpa * np.sin(phase)
        return float(value) if np.ndim(value) == 0 else value

    def validate_time_basis(self) -> None:
        if self.time_basis not in {"model", "physical"}:
            raise ValueError("time_basis must be 'model' or 'physical'")
        if self.time_basis == "physical":
            if self.time_calibration is None:
                raise ValueError("physical time requires a kinetic calibration")
            self.time_calibration.require_calibrated()
            validate_for_energy_model(self.time_calibration, self.energy_model)
            if (
                self.physical_frequency_hz is None
                or not np.isfinite(self.physical_frequency_hz)
                or self.physical_frequency_hz <= 0.0
            ):
                raise ValueError("physical_frequency_hz must be finite and positive")

    def validate(self) -> None:
        self.validate_time_basis()
        positive = {
            "young_gpa": self.young_gpa,
            "model_frequency": self.effective_model_frequency,
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
        if self.energy_model not in ENERGY_MODEL_IDS:
            raise ValueError(f"unknown energy model: {self.energy_model}")


def load_interpretation(
    *, young_gpa: float, stress_mean_mpa: float, stress_amplitude_mpa: float
) -> dict[str, float | str]:
    """Return a display-only signed stress range and mechanism-regime label.

    The regime never clamps or modifies solver input.  Thresholds are only UI
    context for distinguishing small-signal checks from extreme reduced-stress
    mechanism probes.
    """

    young_mpa = float(young_gpa) * 1000.0
    mean = float(stress_mean_mpa)
    amplitude = float(stress_amplitude_mpa)
    if not np.isfinite(young_mpa) or young_mpa <= 0.0:
        raise ValueError("young_gpa must be finite and positive")
    if not np.isfinite(mean) or not np.isfinite(amplitude) or amplitude < 0.0:
        raise ValueError("stress inputs must be finite and amplitude nonnegative")
    minimum, maximum = mean - amplitude, mean + amplitude
    reduced_min, reduced_max = minimum / young_mpa, maximum / young_mpa
    maximum_reduced_magnitude = max(abs(reduced_min), abs(reduced_max))
    if maximum < 0.0 and maximum_reduced_magnitude <= 2.0e-2:
        regime = "compression"
    elif maximum_reduced_magnitude <= 1.0e-3:
        regime = "small"
    elif maximum_reduced_magnitude <= 1.0e-2:
        regime = "moderate"
    else:
        regime = "extreme"
    return {
        "stress_min_mpa": float(minimum),
        "stress_max_mpa": float(maximum),
        "reduced_stress_min": float(reduced_min),
        "reduced_stress_max": float(reduced_max),
        "maximum_reduced_stress_magnitude": float(maximum_reduced_magnitude),
        "regime": regime,
    }


def canonical_model_params() -> ModelParams:
    """Return the verified N=1 UI mechanism parameters."""

    return ModelParams(
        n_cells=1,
        kT=0.02,
        mobility_a=1.0,
        mobility_s=0.05,
        chi_axial_projection=0.20,
    )


def physical_load_conversion(config: UIAnalysisConfig) -> dict[str, object]:
    """Map physical stress inputs through sigma/E and relaxed axial kappa."""

    config.validate()
    model = build_time_basis_model(config.energy_model, time_basis=config.time_basis,
                                  calibration=config.time_calibration)
    kappa = model.sigma_over_E_force_scale()
    reduced_min = config.stress_min_mpa / config.young_mpa
    reduced_max = config.stress_max_mpa / config.young_mpa
    reduced_initial = config.stress_mean_mpa / config.young_mpa
    dynamics = model_frequency_diagnostics(model, config.effective_model_frequency)
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
        **energy_model_result_metadata(config.energy_model, model),
        **production_capabilities(),
        "time_basis": config.time_basis,
        "time_unit": "s" if config.time_basis == "physical" else "model time",
        "frequency_unit": (
            "Hz" if config.time_basis == "physical" else "cycles / model time"
        ),
        "frequency_hz": config.frequency_hz,
        "physical_period_seconds": config.physical_period_seconds,
        "physical_duration_seconds": config.physical_duration_seconds,
        "t0_seconds": (
            float(config.time_calibration.t0_seconds)
            if config.time_basis == "physical" and config.time_calibration is not None
            else None
        ),
        "physical_M_a": (config.time_calibration.M_a_phys_m2_per_J_s
                          if config.time_basis == "physical" else None),
        "physical_M_s": (config.time_calibration.M_s_phys_m2_per_J_s
                          if config.time_basis == "physical" else None),
        "kinetic_calibration_source": (config.time_calibration.source
                          if config.time_basis == "physical" else None),
        **dynamics,
    }


def production_capabilities() -> dict[str, object]:
    """Actual UI/PDE capabilities, not research scripts or placeholder controls."""
    return dict(state_coordinates=["a", "s"], probability_state_dimension=2,
                loading_mode="scalar axial sigma/E with fixed chi=0.2",
                independent_shear_input=False, orientation_input_active=False,
                vector_registry_pde=False, spatial_specimen_solver=False)


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
        "local_survival_probability": survival,
        "initiation_probability": absorbed,
        "local_initiation_probability": absorbed,
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
        record["intact_probability_mass"], record["cumulative_absorbed_mass"]
    )
    physical_time = (
        float(model_time_to_seconds(model_time, config.time_calibration))
        if config.time_basis == "physical"
        else None
    )
    return {
        **record,
        **{key: float(value) for key, value in probability.items()},
        "model_time": model_time,
        "physical_time_seconds": physical_time,
        "plot_time": physical_time if physical_time is not None else model_time,
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
        "strain_decomposition_residual": np.asarray(
            result["strain_decomposition_residual"], dtype=float
        ),
        "strain": np.asarray(result["strain"], dtype=float),
        "survival": np.asarray(result["survival"], dtype=float),
        "survival_probability": np.asarray(
            result["survival_probability"], dtype=float
        ),
        "local_survival_probability": np.asarray(
            result["local_survival_probability"], dtype=float
        ),
        "initiation_probability": np.asarray(
            result["initiation_probability"], dtype=float
        ),
        "local_initiation_probability": np.asarray(
            result["local_initiation_probability"], dtype=float
        ),
        "first_passage_flux": np.asarray(result["first_passage_flux"], dtype=float),
        "intact_probability_mass": np.asarray(
            result["intact_probability_mass"], dtype=float
        ),
        "cumulative_absorbed_mass": np.asarray(
            result["cumulative_absorbed_mass"], dtype=float
        ),
        "absorbed_mass_increment": np.asarray(
            result["absorbed_mass_increment"], dtype=float
        ),
        "initial_absorbed_mass": np.asarray(
            result["initial_absorbed_mass"], dtype=float
        ),
        "integrated_first_passage_flux": np.asarray(
            result["integrated_first_passage_flux"], dtype=float
        ),
        "flux_consistency_residual": np.asarray(
            result["flux_consistency_residual"], dtype=float
        ),
        "mass_balance_residual": np.asarray(
            result["mass_balance_residual"], dtype=float
        ),
        "negative_mass_correction": np.asarray(
            result["negative_mass_correction"], dtype=float
        ),
        "cumulative_negative_mass_correction": np.asarray(
            result["cumulative_negative_mass_correction"], dtype=float
        ),
        "minimum_density": np.asarray(result["minimum_density"], dtype=float),
        "raw_intact_mass": np.asarray(result["raw_intact_mass"], dtype=float),
        "raw_one_minus_survival": np.asarray(
            result["raw_one_minus_survival"], dtype=float
        ),
        "mean_well_index": np.asarray(result["mean_well_index"], dtype=float),
        "plastic_well_activity": np.asarray(
            result["plastic_well_activity"], dtype=float
        ),
        "unnormalized_registry_moment": np.asarray(
            result["unnormalized_registry_moment"], dtype=float
        ),
        "accumulated_net_registry_transfer": np.asarray(
            result["accumulated_net_registry_transfer"], dtype=float
        ),
        "absorbed_registry_moment": np.asarray(
            result["absorbed_registry_moment"], dtype=float
        ),
        "registry_moment_balance_residual": np.asarray(
            result["registry_moment_balance_residual"], dtype=float
        ),
        "net_interwell_registry_rate": np.asarray(
            result["net_interwell_registry_rate"], dtype=float
        ),
        "gross_interwell_activity_rate": np.asarray(
            result["gross_interwell_activity_rate"], dtype=float
        ),
        "net_plastic_flow_rate": np.asarray(
            result["net_plastic_flow_rate"], dtype=float
        ),
        "gross_configurational_slip_activity": np.asarray(
            result["gross_configurational_slip_activity"], dtype=float
        ),
        "selective_opening_plastic_rate": np.asarray(
            result["selective_opening_plastic_rate"], dtype=float
        ),
        "cumulative_forward_registry_activity": np.asarray(
            result["cumulative_forward_registry_activity"], dtype=float
        ),
        "cumulative_backward_registry_activity": np.asarray(
            result["cumulative_backward_registry_activity"], dtype=float
        ),
        "cumulative_gross_registry_activity": np.asarray(
            result["cumulative_gross_registry_activity"], dtype=float
        ),
    }


def per_cycle_first_passage_diagnostics(
    result: dict[str, object],
) -> list[dict[str, float]]:
    """Summarize mechanical, registry, and opening diagnostics by cycle."""

    cycle = np.asarray(result["load_cycle"], dtype=float)
    absorbed = np.asarray(result["cumulative_absorbed_mass"], dtype=float)
    survival = np.asarray(result["survival"], dtype=float)
    flux = np.asarray(result["first_passage_flux"], dtype=float)
    force = np.asarray(result["force"], dtype=float)
    total_strain = np.asarray(result["strain"], dtype=float)
    normal_strain = np.asarray(result["normal_strain"], dtype=float)
    intrawell_strain = np.asarray(result["intrawell_strain"], dtype=float)
    plastic_strain = np.asarray(result["plastic_strain"], dtype=float)
    net_transfer = np.asarray(
        result.get("cumulative_forward_registry_activity", np.zeros_like(cycle)),
        dtype=float,
    ) - np.asarray(
        result.get("cumulative_backward_registry_activity", np.zeros_like(cycle)),
        dtype=float,
    )
    gross_activity = np.asarray(
        result.get("cumulative_gross_registry_activity", np.zeros_like(cycle)),
        dtype=float,
    )
    model = result["model"]
    if cycle.size == 0:
        return []
    last_cycle = int(np.ceil(float(cycle[-1]) - 1.0e-12))
    s_probe = np.linspace(-0.5 * model.p.b, 0.5 * model.p.b, 33)
    barrier_cache: dict[float, float] = {}

    def minimum_barrier(load_force: float) -> float:
        key = round(float(load_force), 12)
        if key not in barrier_cache:
            barrier_cache[key] = float(
                min(model.opening_barrier(float(s), key) for s in s_probe)
            )
        return barrier_cache[key]

    rows: list[dict[str, float]] = []
    for number in range(1, last_cycle + 1):
        start = float(number - 1)
        end = min(float(number), float(cycle[-1]))
        mask = (cycle >= start - 1.0e-12) & (cycle <= end + 1.0e-12)
        if not np.any(mask):
            continue
        absorbed_start = float(np.interp(start, cycle, absorbed))
        absorbed_end = float(np.interp(end, cycle, absorbed))
        survival_end = float(np.interp(end, cycle, survival))
        plastic_start = float(np.interp(start, cycle, plastic_strain))
        plastic_end = float(np.interp(end, cycle, plastic_strain))
        net_start = float(np.interp(start, cycle, net_transfer))
        net_end = float(np.interp(end, cycle, net_transfer))
        gross_start = float(np.interp(start, cycle, gross_activity))
        gross_end = float(np.interp(end, cycle, gross_activity))
        rows.append(
            {
                "cycle": float(number),
                "cycle_end": end,
                "mean_total_strain": float(np.mean(total_strain[mask])),
                "total_strain_amplitude": float(
                    0.5 * np.ptp(total_strain[mask])
                ),
                "mean_normal_strain": float(np.mean(normal_strain[mask])),
                "mean_intrawell_strain": float(np.mean(intrawell_strain[mask])),
                "plastic_strain_increment": plastic_end - plastic_start,
                "cumulative_plastic_strain": plastic_end,
                "net_registry_transfer": net_end - net_start,
                "gross_registry_activity": gross_end - gross_start,
                "absorbed_mass": absorbed_end - absorbed_start,
                "minimum_configurational_barrier": float(
                    result.get("configurational_barrier", np.nan)
                ),
                "minimum_opening_barrier": min(
                    minimum_barrier(value) for value in force[mask]
                ),
                "peak_first_passage_flux": float(np.max(flux[mask])),
                "survival_at_cycle_end": survival_end,
            }
        )
    return rows


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
    calibration_model = build_time_basis_model(config.energy_model, time_basis=config.time_basis,
                                              calibration=config.time_calibration)
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
        prepared_model=calibration_model,
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
        calibration_model, config.effective_model_frequency
    )
    tau_fast_seconds = (
        dynamics["tau_fast"] * float(config.time_calibration.t0_seconds)
        if config.time_basis == "physical" and config.time_calibration is not None
        else None
    )
    tau_slow_seconds = (
        dynamics["tau_slow"] * float(config.time_calibration.t0_seconds)
        if config.time_basis == "physical" and config.time_calibration is not None
        else None
    )
    probability = physical_probability_bookkeeping(
        raw["intact_probability_mass"], raw["cumulative_absorbed_mass"]
    )
    peak_force = float(max(load.force_min, load.force_max))
    barrier = bound_configurational_barrier(
        calibration_model,
        peak_force,
        n_s=61,
        quadrature_order=24,
    )
    landscape_diagnostics = {
        "barrier_diagnostic_force": peak_force,
        "configurational_barrier": (
            float(barrier.barrier) if barrier is not None else np.nan
        ),
        "opening_barrier_at_configurational_saddle": (
            float(barrier.opening_barrier_at_saddle)
            if barrier is not None
            else np.nan
        ),
        "configurational_barrier_status": (
            "metastable branch resolved"
            if barrier is not None
            else "no distinct metastable minimum/saddle"
        ),
    }
    if config.time_basis == "physical":
        physical_time = np.asarray(
            model_time_to_seconds(model_time, config.time_calibration), dtype=float
        )
        solver_time_status = "physical seconds from calibrated kinetic mobility"
    else:
        physical_time = None
        solver_time_status = "dimensionless model time; not calibrated seconds"
    result: dict[str, object] = {
        **raw,
        **probability,
        "model_time": model_time,
        "physical_time_seconds": physical_time,
        "plot_time": physical_time if physical_time is not None else model_time,
        "time_basis": config.time_basis,
        "time_unit": "s" if config.time_basis == "physical" else "model time",
        "frequency_unit": (
            "Hz" if config.time_basis == "physical" else "cycles / model time"
        ),
        "frequency_hz": config.frequency_hz,
        "physical_period_seconds": config.physical_period_seconds,
        "physical_duration_seconds": config.physical_duration_seconds,
        "tau_fast_seconds": tau_fast_seconds,
        "tau_slow_seconds": tau_slow_seconds,
        "t0_seconds": (
            float(config.time_calibration.t0_seconds)
            if config.time_basis == "physical" and config.time_calibration is not None
            else None
        ),
        "physical_M_a": (
            config.time_calibration.M_a_phys_m2_per_J_s
            if config.time_basis == "physical" and config.time_calibration is not None
            else None
        ),
        "physical_M_s": (
            config.time_calibration.M_s_phys_m2_per_J_s
            if config.time_basis == "physical" and config.time_calibration is not None
            else None
        ),
        "kinetic_calibration_source": (
            config.time_calibration.source
            if config.time_basis == "physical" and config.time_calibration is not None
            else None
        ),
        "load_cycle": load_cycle,
        "applied_stress_mpa": np.asarray(config.stress_mpa(model_time), dtype=float),
        "initial_stress_mpa": float(config.stress_mean_mpa),
        "initial_force": preload_force,
        "initial_condition": "conditional Gibbs at sigma(t=0)",
        "relaxed_axial_kappa": float(calibration_model.sigma_over_E_force_scale()),
        "frozen_normal_kappa": float(
            calibration_model.frozen_normal_sigma_over_E_force_scale()
        ),
        "solver_time_status": solver_time_status,
        "probability_source": "N=1 direct Smoluchowski/Fokker-Planck PDE",
        **energy_model_result_metadata(config.energy_model, calibration_model),
        **production_capabilities(),
        "kinetic_calibration": (
            config.time_calibration.to_dict() if config.time_basis == "physical" else None
        ),
        "analysis_quality": config.analysis_quality,
        "integration_method": config.integration_method,
        "grid_shape": (config.grid_n_a, config.grid_n_s),
        **landscape_diagnostics,
        **dynamics,
    }
    assessment = assess_local_rare_event(result)
    result.update(
        {
            "cumulative_forward_registry_activity": np.sum(
                np.asarray(result["cumulative_interwell_forward_transfer"]), axis=1
            ),
            "cumulative_backward_registry_activity": np.sum(
                np.asarray(result["cumulative_interwell_backward_transfer"]), axis=1
            ),
            "cumulative_gross_registry_activity": np.sum(
                np.asarray(result["cumulative_interwell_gross_transfer"]), axis=1
            ),
        }
    )
    result.update(
        {
            "local_rare_event_floor": assessment.floor,
            "probability_resolution_status": assessment.status,
            "rare_event_floor_mass_residual": assessment.mass_residual,
            "rare_event_floor_accumulated_repair": assessment.accumulated_repair,
            "rare_event_floor_flux_discrepancy": assessment.flux_discrepancy,
            "rare_event_floor_resolution_change": assessment.resolution_change,
            "rare_event_floor_scope": "single-run lower bound",
            "probability_resolution_certified": False,
        }
    )
    result.update(single_run_plasticity_lower_bound(result))
    result["per_cycle_diagnostics"] = per_cycle_first_passage_diagnostics(result)
    result_field_mapping(result)
    return result
