from functools import lru_cache

import numpy as np

from app.solver_adapter import (
    PDE_RESULT_FIELDS,
    UIAnalysisConfig,
    physical_probability_bookkeeping,
    physical_load_conversion,
    load_interpretation,
    result_field_mapping,
    run_ui_analysis,
)


@lru_cache(maxsize=2)
def _physical_case(name: str) -> tuple[dict[str, object], tuple[dict[str, float], ...]]:
    if name == "compression":
        mean, amplitude = -500.0, 400.0
    elif name == "nonlinear":
        mean, amplitude = 900.0, 2000.0
    else:
        raise ValueError(name)
    records: list[dict[str, float]] = []
    result = run_ui_analysis(
        UIAnalysisConfig(
            young_gpa=69.0,
            stress_mean_mpa=mean,
            stress_amplitude_mpa=amplitude,
            model_frequency=25.0,
            cycles=1.0,
            steps_per_cycle=40,
        ),
        record_callback=records.append,
    )
    return result, tuple(records)


def _assert_common_invariants(result: dict[str, object]) -> None:
    total = np.asarray(result["strain"], dtype=float)
    normal = np.asarray(result["normal_strain"], dtype=float)
    intrawell = np.asarray(result["intrawell_strain"], dtype=float)
    plastic = np.asarray(result["plastic_strain"], dtype=float)
    survival = np.asarray(result["survival"], dtype=float)
    initiation = np.asarray(result["initiation_probability"], dtype=float)
    absorbed = np.asarray(result["cumulative_absorbed_mass"], dtype=float)
    raw_intact = np.asarray(result["raw_intact_mass"], dtype=float)
    np.testing.assert_allclose(total, normal + intrawell + plastic, rtol=0.0, atol=2e-12)
    np.testing.assert_array_equal(initiation, absorbed)
    np.testing.assert_array_equal(survival, 1.0 - absorbed)
    np.testing.assert_array_equal(result["survival_probability"], survival)
    np.testing.assert_array_equal(result["local_survival_probability"], survival)
    np.testing.assert_array_equal(result["local_initiation_probability"], absorbed)
    np.testing.assert_array_equal(result["raw_one_minus_survival"], 1.0 - raw_intact)
    np.testing.assert_allclose(
        result["mass_balance_residual"], raw_intact + absorbed - 1.0,
        rtol=0.0, atol=2e-15,
    )
    assert np.all(np.diff(survival) <= 1e-10)
    assert np.all(np.diff(initiation) >= -1e-10)
    assert np.max(np.abs(np.asarray(result["mass_balance_residual"]))) < 1e-8
    assert np.max(np.abs(np.asarray(result["flux_consistency_residual"]))) < 1e-12
    assert result["plastic_floor_scope"] == "single-run lower bound"
    assert result["plastic_resolution_status"] == "requires_convergence"
    for key in PDE_RESULT_FIELDS:
        assert np.all(np.isfinite(np.asarray(result[key], dtype=float)))


def test_physical_stress_uses_verified_relaxed_axial_force_mapping() -> None:
    config = UIAnalysisConfig(
        young_gpa=69.0,
        stress_mean_mpa=-500.0,
        stress_amplitude_mpa=400.0,
    )
    conversion = physical_load_conversion(config)
    assert conversion["a0"] == 0.7713438268704838
    assert np.isclose(conversion["relaxed_axial_kappa"], 86.29296488740997)
    assert np.isclose(conversion["frozen_normal_kappa"], 94.71096726647392)
    assert np.isclose(conversion["reduced_stress_min"], -900.0 / 69000.0)
    assert np.isclose(conversion["reduced_stress_max"], -100.0 / 69000.0)
    assert np.isclose(
        conversion["force_min"],
        conversion["relaxed_axial_kappa"] * (-900.0 / 69000.0),
    )
    assert conversion["model_frequency"] == 25.0
    assert conversion["model_period"] == 0.04
    assert np.isclose(conversion["de_fast"], 1.279285899947692)
    assert np.isclose(conversion["de_slow"], 62.39808635219237)


def test_load_interpretation_warns_without_modifying_signed_inputs() -> None:
    case_a = load_interpretation(
        young_gpa=69.0,
        stress_mean_mpa=-500.0,
        stress_amplitude_mpa=400.0,
    )
    case_b = load_interpretation(
        young_gpa=69.0,
        stress_mean_mpa=900.0,
        stress_amplitude_mpa=2000.0,
    )
    assert case_a["stress_min_mpa"] == -900.0
    assert case_a["stress_max_mpa"] == -100.0
    assert case_a["regime"] == "compression"
    assert case_b["stress_min_mpa"] == -1100.0
    assert case_b["stress_max_mpa"] == 2900.0
    assert case_b["regime"] == "extreme"


def test_physical_probability_uses_absorbed_mass_not_roundoff_residual() -> None:
    raw_intact = np.array([1.0, 0.800001])
    absorbed = np.array([0.0, 0.2])
    mapped = physical_probability_bookkeeping(raw_intact, absorbed)
    np.testing.assert_array_equal(mapped["initiation_probability"], absorbed)
    np.testing.assert_array_equal(mapped["survival"], 1.0 - absorbed)
    np.testing.assert_array_equal(mapped["survival_probability"], 1.0 - absorbed)
    np.testing.assert_array_equal(mapped["local_survival_probability"], 1.0 - absorbed)
    np.testing.assert_array_equal(mapped["local_initiation_probability"], absorbed)
    np.testing.assert_array_equal(mapped["raw_intact_mass"], raw_intact)
    np.testing.assert_allclose(
        mapped["raw_one_minus_survival"], np.array([0.0, 0.199999]),
        rtol=0.0, atol=2e-16,
    )
    assert mapped["initiation_probability"][1] != mapped["raw_one_minus_survival"][1]


def test_ui_result_registry_maps_required_pde_fields_directly() -> None:
    result, records = _physical_case("compression")
    fields = result_field_mapping(result)
    assert set(PDE_RESULT_FIELDS).issubset(result)
    assert records
    assert records[0]["applied_stress_mpa"] == -500.0
    assert result["initial_condition"] == "conditional Gibbs at sigma(t=0)"
    assert result["model_frequency"] == 25.0
    assert result["solver_time_status"] == (
        "dimensionless model time; not calibrated seconds"
    )
    for source, mapped in (
        ("strain", "strain"),
        ("normal_strain", "normal_strain"),
        ("intrawell_strain", "intrawell_strain"),
        ("plastic_strain", "plastic_strain"),
        ("survival", "survival"),
        ("survival_probability", "survival_probability"),
        ("initiation_probability", "initiation_probability"),
        ("first_passage_flux", "first_passage_flux"),
    ):
        assert np.shares_memory(fields[mapped], np.asarray(result[source]))


def test_compression_only_ui_data_path_is_physically_consistent() -> None:
    result, _ = _physical_case("compression")
    stress = np.asarray(result["applied_stress_mpa"], dtype=float)
    strain = np.asarray(result["strain"], dtype=float)
    assert np.max(stress) < 0.0
    assert np.min(stress) <= -899.0
    assert np.max(stress) >= -101.0
    assert np.all(strain < 0.0)
    assert 1.0e-3 < np.max(np.abs(strain)) < 1.0e-1
    assert float(np.max(result["first_passage_flux"])) < 1.0e-8
    assert float(np.max(result["initiation_probability"])) < 1.0e-8
    assert float(np.min(result["survival"])) > 1.0 - 1.0e-8
    _assert_common_invariants(result)


def test_high_load_ui_data_path_exposes_nonlinear_probability_response() -> None:
    result, _ = _physical_case("nonlinear")
    stress = np.asarray(result["applied_stress_mpa"], dtype=float)
    strain = np.asarray(result["strain"], dtype=float)
    flux = np.asarray(result["first_passage_flux"], dtype=float)
    assert np.min(stress) <= -1099.0
    assert np.max(stress) >= 2899.0
    assert np.ptp(strain) > 1.0e-3
    assert np.max(np.abs(strain)) > 1.0e-3
    assert np.max(flux) > 0.0
    assert float(np.asarray(result["initiation_probability"])[-1]) >= 0.0
    assert float(np.asarray(result["survival"])[-1]) <= 1.0
    _assert_common_invariants(result)
