from dataclasses import replace

import numpy as np
import pytest

from app.convergence_check import run_convergence_check
from app.i18n import plot_strings
from app.solver_adapter import (
    UIAnalysisConfig,
    per_cycle_first_passage_diagnostics,
    physical_load_conversion,
    run_ui_analysis,
)
from app.specimen_probability import aggregate_specimen_probability
from solver_v1.physical_time import calibration_from_mobilities
from solver_v1.physical_time import uncalibrated_time_calibration


def _calibration(t0: float = 0.2):
    length = 2.8e-10
    energy = 1.6e-19
    mobility_a = length**2 / (t0 * energy)
    return calibration_from_mobilities(
        length_scale_m=length,
        energy_scale_J=energy,
        M_a_phys_m2_per_J_s=mobility_a,
        M_s_phys_m2_per_J_s=0.05 * mobility_a,
        source="hypothetical UI regression fixture",
        temperature_K=300.0,
    )


def test_physical_frequency_and_cycles_define_one_consistent_duration() -> None:
    config = UIAnalysisConfig(
        time_basis="physical",
        physical_frequency_hz=5.0,
        time_calibration=_calibration(0.2),
        cycles=3.0,
    )
    config.validate()
    assert config.effective_model_frequency == pytest.approx(1.0)
    assert config.model_period == pytest.approx(1.0)
    assert config.physical_period_seconds == pytest.approx(0.2)
    assert config.physical_duration_seconds == pytest.approx(0.6)


def test_uncalibrated_physical_mode_is_rejected_and_model_mode_is_unchanged() -> None:
    uncalibrated = uncalibrated_time_calibration(
        length_scale_m=2.8e-10,
        energy_scale_J=1.6e-19,
    )
    with pytest.raises(ValueError, match="physical seconds/Hz"):
        UIAnalysisConfig(
            time_basis="physical",
            physical_frequency_hz=5.0,
            time_calibration=uncalibrated,
        ).validate()

    default = UIAnalysisConfig(model_frequency=25.0)
    default.validate()
    assert default.effective_model_frequency == pytest.approx(25.0)
    assert default.frequency_hz is None


def test_physical_and_model_time_inputs_map_to_identical_pde_controls() -> None:
    calibration = _calibration(0.2)
    model = UIAnalysisConfig(model_frequency=1.0, cycles=0.02)
    physical = replace(
        model,
        model_frequency=999.0,
        time_basis="physical",
        physical_frequency_hz=5.0,
        time_calibration=calibration,
    )
    model_values = physical_load_conversion(model)
    physical_values = physical_load_conversion(physical)
    for key in ("force_min", "force_max", "preload_force", "model_period"):
        assert physical_values[key] == pytest.approx(model_values[key])

    common = dict(
        young_gpa=69.0,
        stress_mean_mpa=0.0,
        stress_amplitude_mpa=20.0,
        cycles=0.02,
        steps_per_cycle=4,
        grid_n_a=7,
        grid_n_s=9,
        max_dt=2.0e-3,
    )
    model_result = run_ui_analysis(UIAnalysisConfig(model_frequency=1.0, **common))
    physical_result = run_ui_analysis(
        UIAnalysisConfig(
            model_frequency=999.0,
            time_basis="physical",
            physical_frequency_hz=5.0,
            time_calibration=calibration,
            **common,
        )
    )
    for key in ("model_time", "strain", "cumulative_absorbed_mass"):
        np.testing.assert_allclose(physical_result[key], model_result[key], rtol=0, atol=0)
    np.testing.assert_allclose(
        physical_result["physical_time_seconds"],
        np.asarray(model_result["model_time"]) * 0.2,
    )


def test_plot_metadata_switches_seconds_without_changing_result_arrays() -> None:
    model = plot_strings("first_passage_flux", "en", time_basis="model")
    physical = plot_strings("first_passage_flux", "en", time_basis="physical")
    assert model["xlabel"] == "Dimensionless solver time"
    assert model["ylabel"] == "First-passage flux [1/model time]"
    assert physical["xlabel"] == "Time [s]"
    assert physical["ylabel"] == "First-passage flux [1/s]"


def test_physical_frequency_changes_the_existing_model_dynamics_parameter() -> None:
    calibration = _calibration(0.2)
    low = UIAnalysisConfig(
        time_basis="physical",
        physical_frequency_hz=1.0,
        time_calibration=calibration,
    )
    high = replace(low, physical_frequency_hz=10.0)
    low_values = physical_load_conversion(low)
    high_values = physical_load_conversion(high)
    assert high_values["model_frequency"] == pytest.approx(
        10.0 * low_values["model_frequency"]
    )
    assert high_values["model_period"] == pytest.approx(
        0.1 * low_values["model_period"]
    )
    assert high_values["de_fast"] == pytest.approx(10.0 * low_values["de_fast"])
    assert high_values["de_slow"] == pytest.approx(10.0 * low_values["de_slow"])
    assert high_values["linear_transfer_magnitude"] < low_values[
        "linear_transfer_magnitude"
    ]


def _fake_result(signal: float, quality: str = "resolved") -> dict[str, object]:
    return {
        "analysis_quality": quality,
        "cumulative_absorbed_mass": np.array([0.0, signal]),
        "mass_balance_residual": np.array([0.0, 1.0e-12]),
        "cumulative_negative_mass_correction": np.zeros(2),
        "flux_consistency_residual": np.zeros(2),
    }


def test_convergence_workflow_performs_refinements_before_certifying() -> None:
    calls: list[UIAnalysisConfig] = []

    def runner(config, **_kwargs):
        calls.append(config)
        return _fake_result(1.0e-5 + len(calls) * 1.0e-8)

    config = UIAnalysisConfig(
        model_frequency=1.0,
        cycles=0.1,
        grid_n_a=31,
        grid_n_s=41,
        max_dt=1.0e-3,
        integration_method="implicit",
        analysis_quality="resolved",
    )
    report = run_convergence_check(
        config, reference_result=_fake_result(1.0e-5), runner=runner
    )
    assert len(calls) == 3
    assert calls[0].max_dt == pytest.approx(0.5e-3)
    assert calls[1].integration_method == "explicit"
    assert calls[-1].grid_n_a > config.grid_n_a
    assert calls[-1].grid_n_s > config.grid_n_s
    assert report.certified_result["probability_resolution_certified"] is True


def test_preview_cannot_certify_specimen_probability() -> None:
    preview = UIAnalysisConfig(analysis_quality="preview")
    with pytest.raises(ValueError, match="Preview"):
        run_convergence_check(
            preview, reference_result=_fake_result(1.0e-4, "preview")
        )
    aggregation = aggregate_specimen_probability(
        [1.0e-4],
        correlation_area_mm2=1.0,
        stressed_area_mm2=100.0,
        local_numerical_floor=1.0e-12,
        resolution_certified=False,
    )
    assert np.isfinite(aggregation.mathematical_extrapolation[-1])
    assert np.isnan(aggregation.specimen_initiation_probability[-1])


class _BarrierModel:
    class _Parameters:
        b = 1.0

    p = _Parameters()

    @staticmethod
    def opening_barrier(s: float, force: float) -> float:
        return 2.0 - force + 0.1 * s**2


def test_per_cycle_absorption_table_distinguishes_first_cycle_transient() -> None:
    result = {
        "load_cycle": np.array([0.0, 0.5, 1.0, 1.5, 2.0]),
        "cumulative_absorbed_mass": np.array([0.0, 0.10, 0.12, 0.125, 0.13]),
        "survival": np.array([1.0, 0.90, 0.88, 0.875, 0.87]),
        "first_passage_flux": np.array([0.0, 0.2, 0.04, 0.01, 0.01]),
        "force": np.array([0.0, 1.0, 0.0, 1.0, 0.0]),
        "model": _BarrierModel(),
    }
    rows = per_cycle_first_passage_diagnostics(result)
    assert len(rows) == 2
    assert rows[0]["absorbed_mass"] == pytest.approx(0.12)
    assert rows[1]["absorbed_mass"] == pytest.approx(0.01)
    assert rows[0]["peak_first_passage_flux"] > rows[1]["peak_first_passage_flux"]
    assert rows[0]["minimum_opening_barrier"] == pytest.approx(1.0)
