from pathlib import Path

import numpy as np
import pytest

from solver_v1.physical_time import (
    BOLTZMANN_J_PER_K,
    UncalibratedTimeError,
    calibration_from_mobilities,
    hz_to_model_frequency,
    infer_diagonal_mobilities_from_relaxation_times,
    load_time_calibration,
    mobility_from_collective_diffusivity,
    mobility_from_decoupled_relaxation,
    model_frequency_to_hz,
    model_time_to_seconds,
    physical_relaxation_spectrum,
    seconds_to_model_time,
)


L0 = 2.8e-10
E0 = 1.6e-19
T0 = 2.5e-12
MA_STAR = 1.0
MS_STAR = 0.05
MA_PHYS = MA_STAR * L0**2 / (T0 * E0)
MS_PHYS = MS_STAR * L0**2 / (T0 * E0)


@pytest.fixture
def calibrated():
    return calibration_from_mobilities(
        length_scale_m=L0,
        energy_scale_J=E0,
        M_a_phys_m2_per_J_s=MA_PHYS,
        M_s_phys_m2_per_J_s=MS_PHYS,
        source="hypothetical regression fixture",
        temperature_K=300.0,
    )


def test_mobility_dimensions_produce_declared_common_time_scale(calibrated) -> None:
    assert calibrated.t0_seconds == pytest.approx(T0)
    assert calibrated.mobility_ratio == pytest.approx(MS_STAR / MA_STAR)
    assert MA_PHYS * E0 * T0 / L0**2 == pytest.approx(MA_STAR)
    assert MS_PHYS * E0 * T0 / L0**2 == pytest.approx(MS_STAR)


def test_time_and_frequency_round_trips(calibrated) -> None:
    model_time = np.array([0.0, 0.25, 3.0])
    seconds = model_time_to_seconds(model_time, calibrated)
    np.testing.assert_allclose(seconds_to_model_time(seconds, calibrated), model_time)
    model_frequency = np.array([0.2, 10.0, 50.0])
    hz = model_frequency_to_hz(model_frequency, calibrated)
    np.testing.assert_allclose(
        hz_to_model_frequency(hz, calibrated), model_frequency
    )
    assert hz[1] == pytest.approx(model_frequency[1] / T0)


def test_uncalibrated_file_refuses_seconds_and_hertz() -> None:
    path = Path(__file__).with_name("data") / "aluminum_kinetic_calibration.json"
    calibration = load_time_calibration(path)
    assert not calibration.calibrated
    with pytest.raises(UncalibratedTimeError):
        model_time_to_seconds(1.0, calibration)
    with pytest.raises(UncalibratedTimeError):
        hz_to_model_frequency(1.0, calibration)


def test_diffusion_and_decoupled_relaxation_routes_have_consistent_units() -> None:
    diffusivity = 4.0e-20
    temperature = 350.0
    mobility = mobility_from_collective_diffusivity(diffusivity, temperature)
    assert mobility * BOLTZMANN_J_PER_K * temperature == pytest.approx(diffusivity)
    curvature = 8.0
    tau = 0.125
    assert mobility_from_decoupled_relaxation(tau, curvature) == pytest.approx(1.0)


def test_coupled_relaxation_inference_recovers_synthetic_mobilities() -> None:
    hessian = np.array([[12.0, 1.5], [1.5, 5.0]])
    expected = np.array([0.8, 0.06])
    sqrt_m = np.diag(np.sqrt(expected))
    rates = np.linalg.eigvalsh(sqrt_m @ hessian @ sqrt_m)
    inferred = infer_diagonal_mobilities_from_relaxation_times(
        hessian,
        tau_fast_seconds=1.0 / rates[-1],
        tau_slow_seconds=1.0 / rates[0],
        initial_mobilities_m2_per_J_s=(0.7, 0.05),
    )
    assert inferred.M_a_phys_m2_per_J_s == pytest.approx(expected[0], rel=1e-10)
    assert inferred.M_s_phys_m2_per_J_s == pytest.approx(expected[1], rel=1e-10)
    assert inferred.relative_residual < 1e-11


def test_physical_coupled_rates_equal_model_rates_divided_by_t0(calibrated) -> None:
    hessian = np.array([[122.0, 2.0], [2.0, 50.0]])
    model_mobility = np.diag([MA_STAR, MS_STAR])
    sqrt_model = np.diag(np.sqrt(np.diag(model_mobility)))
    expected = np.linalg.eigvalsh(sqrt_model @ hessian @ sqrt_model) / T0
    actual = physical_relaxation_spectrum(hessian, calibrated)
    np.testing.assert_allclose(actual["rates_per_s"], expected, rtol=2e-15)
