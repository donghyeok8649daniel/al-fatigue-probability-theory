import numpy as np

from app.specimen_probability import (
    BELOW_RESOLUTION,
    RESOLVED,
    aggregate_specimen_probability,
    attach_convergence_certificate,
    assess_convergence_rare_event,
    assess_local_rare_event,
)


def _diagnostic_result(signal: float, residual: float = 0.0) -> dict[str, object]:
    return {
        "cumulative_absorbed_mass": np.array([0.0, signal]),
        "mass_balance_residual": np.array([0.0, residual]),
        "cumulative_negative_mass_correction": np.array([0.0, 0.0]),
        "flux_consistency_residual": np.array([0.0, 0.0]),
    }


def test_rare_event_floor_separates_absorption_from_mass_residual() -> None:
    assessment = assess_local_rare_event(_diagnostic_result(1.0e-16, 3.0e-16))
    assert assessment.signal == 1.0e-16
    assert assessment.floor == 3.0e-16
    assert assessment.status == BELOW_RESOLUTION

    converged = assess_convergence_rare_event(
        [_diagnostic_result(1.0e-5), _diagnostic_result(1.01e-5)]
    )
    assert np.isclose(converged.floor, 1.0e-7)
    assert converged.status == RESOLVED

    reference = _diagnostic_result(1.01e-5)
    attached = attach_convergence_certificate(
        reference, [_diagnostic_result(1.0e-5)]
    )
    assert reference["probability_resolution_certified"] is True
    assert reference["local_rare_event_floor"] == attached.floor


def test_N_eff_one_returns_local_probability() -> None:
    local = np.array([1.0e-8, 0.2, 0.95])
    result = aggregate_specimen_probability(
        local,
        correlation_area_mm2=4.0,
        stressed_area_mm2=4.0,
        local_numerical_floor=0.0,
        resolution_certified=True,
    )
    np.testing.assert_allclose(result.specimen_initiation_probability, local)
    np.testing.assert_allclose(result.specimen_survival_probability, 1.0 - local)


def test_specimen_probability_increases_monotonically_with_N_eff() -> None:
    values = []
    for stressed_area in (1.0, 5.0, 20.0):
        result = aggregate_specimen_probability(
            [0.01],
            correlation_area_mm2=1.0,
            stressed_area_mm2=stressed_area,
            local_numerical_floor=0.0,
            resolution_certified=True,
        )
        values.append(result.specimen_initiation_probability[0])
    assert values[0] < values[1] < values[2]


def test_log_domain_aggregation_is_stable_in_extreme_regimes() -> None:
    tiny = aggregate_specimen_probability(
        [1.0e-14],
        correlation_area_mm2=1.0,
        stressed_area_mm2=1.0e8,
        local_numerical_floor=0.0,
        resolution_certified=True,
    )
    assert np.isclose(tiny.mathematical_extrapolation[0], 1.0e-6, rtol=1.0e-6)
    near_one = aggregate_specimen_probability(
        [1.0 - 1.0e-14],
        correlation_area_mm2=1.0,
        stressed_area_mm2=1.0e6,
        local_numerical_floor=0.0,
        resolution_certified=True,
    )
    assert near_one.specimen_initiation_probability[0] == 1.0
    assert np.isfinite(near_one.specimen_initiation_probability[0])


def test_small_probability_asymptotic_and_unresolved_safeguard() -> None:
    resolved = aggregate_specimen_probability(
        [1.0e-8],
        correlation_area_mm2=1.0,
        stressed_area_mm2=100.0,
        local_numerical_floor=1.0e-12,
        resolution_certified=True,
    )
    assert np.isclose(
        resolved.specimen_initiation_probability[0], 100.0e-8, rtol=5.0e-7
    )

    unresolved = aggregate_specimen_probability(
        [1.0e-16],
        correlation_area_mm2=1.0e-6,
        stressed_area_mm2=1.0e8,
        local_numerical_floor=2.0e-16,
    )
    assert unresolved.status == BELOW_RESOLUTION
    assert np.isnan(unresolved.specimen_initiation_probability[0])
    assert unresolved.mathematical_extrapolation[0] > 0.0

    uncertified = aggregate_specimen_probability(
        [1.0e-4],
        correlation_area_mm2=1.0,
        stressed_area_mm2=100.0,
        local_numerical_floor=1.0e-12,
    )
    assert uncertified.status == BELOW_RESOLUTION
    assert np.isnan(uncertified.specimen_initiation_probability[0])
    assert uncertified.mathematical_extrapolation[0] > 0.0
