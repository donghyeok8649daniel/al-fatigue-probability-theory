"""Physical-time dimensionalization for the overdamped probability model.

This module contains conversions and calibration routes only.  It never
supplies a default kinetic scale and never changes the Smoluchowski generator.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares


BOLTZMANN_J_PER_K = 1.380649e-23


class UncalibratedTimeError(ValueError):
    """Raised when seconds or hertz are requested without kinetic data."""


@dataclass(frozen=True)
class PhysicalTimeCalibration:
    """One common time scale and two physical reduced-coordinate mobilities.

    Both physical coordinates use the same declared length scale ``L0``:
    ``a*=a_phys/L0`` and ``s*=s_phys/L0``.  Consequently

    ``M_i* = t0 M_i_phys E0 / L0^2``.
    """

    calibrated: bool
    source: str | None
    temperature_K: float | None
    length_scale_m: float
    energy_scale_J: float
    M_a_phys_m2_per_J_s: float | None
    M_s_phys_m2_per_J_s: float | None
    t0_seconds: float | None
    mobility_ratio: float | None
    notes: str = ""
    model_mobility_a: float = 1.0
    model_mobility_s: float = 0.05
    energy_model_id: str | None = None
    energy_model_fingerprint: str | None = None
    coordinate_definition: str | None = None

    def validate(self, *, relative_tolerance: float = 2.0e-8) -> None:
        if type(self.calibrated) is not bool:
            raise ValueError("calibrated must be a JSON boolean, not a text label")
        if not np.isfinite(self.length_scale_m) or self.length_scale_m <= 0.0:
            raise ValueError("length_scale_m must be finite and positive")
        if not np.isfinite(self.energy_scale_J) or self.energy_scale_J <= 0.0:
            raise ValueError("energy_scale_J must be finite and positive")
        if any(not np.isfinite(x) or x <= 0.0 for x in
               (self.model_mobility_a, self.model_mobility_s)):
            raise ValueError("model mobilities must be positive")
        if self.temperature_K is not None and (
            not np.isfinite(self.temperature_K) or self.temperature_K <= 0.0
        ):
            raise ValueError("temperature_K must be finite and positive")
        if not self.calibrated:
            return
        values = (
            self.M_a_phys_m2_per_J_s,
            self.M_s_phys_m2_per_J_s,
            self.t0_seconds,
            self.mobility_ratio,
        )
        if any(value is None or not np.isfinite(value) or value <= 0.0 for value in values):
            raise ValueError("calibrated time data must be finite and positive")
        if not isinstance(self.source, str) or not self.source.strip():
            raise ValueError("a calibrated time scale requires a source")
        t0_a = self.model_mobility_a * self.length_scale_m**2 / (
            float(self.M_a_phys_m2_per_J_s) * self.energy_scale_J
        )
        t0_s = self.model_mobility_s * self.length_scale_m**2 / (
            float(self.M_s_phys_m2_per_J_s) * self.energy_scale_J
        )
        if not np.isclose(t0_a, self.t0_seconds, rtol=relative_tolerance, atol=0.0):
            raise ValueError("t0 is inconsistent with the physical a mobility")
        if not np.isclose(t0_s, self.t0_seconds, rtol=relative_tolerance, atol=0.0):
            raise ValueError("t0 is inconsistent with the physical s mobility")
        ratio = float(self.M_s_phys_m2_per_J_s) / float(
            self.M_a_phys_m2_per_J_s
        )
        if not np.isclose(ratio, self.mobility_ratio, rtol=relative_tolerance, atol=0.0):
            raise ValueError("mobility_ratio is inconsistent with physical mobilities")

    def require_calibrated(self) -> None:
        self.validate()
        if not self.calibrated:
            raise UncalibratedTimeError(
                "physical seconds/Hz require calibrated reduced-coordinate mobility"
            )

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class MobilityInference:
    M_a_phys_m2_per_J_s: float
    M_s_phys_m2_per_J_s: float
    fitted_rates_per_s: np.ndarray
    target_rates_per_s: np.ndarray
    relative_residual: float
    note: str


def uncalibrated_time_calibration(
    *,
    length_scale_m: float,
    energy_scale_J: float,
    notes: str = "",
    model_mobility_a: float = 1.0,
    model_mobility_s: float = 0.05,
) -> PhysicalTimeCalibration:
    calibration = PhysicalTimeCalibration(
        calibrated=False,
        source=None,
        temperature_K=None,
        length_scale_m=float(length_scale_m),
        energy_scale_J=float(energy_scale_J),
        M_a_phys_m2_per_J_s=None,
        M_s_phys_m2_per_J_s=None,
        t0_seconds=None,
        mobility_ratio=None,
        notes=notes,
        model_mobility_a=float(model_mobility_a),
        model_mobility_s=float(model_mobility_s),
    )
    calibration.validate()
    return calibration


def calibration_from_mobilities(
    *,
    length_scale_m: float,
    energy_scale_J: float,
    M_a_phys_m2_per_J_s: float,
    M_s_phys_m2_per_J_s: float,
    source: str,
    temperature_K: float | None,
    notes: str = "",
    model_mobility_a: float = 1.0,
    model_mobility_s: float = 0.05,
) -> PhysicalTimeCalibration:
    if not isinstance(source, str) or not source.strip():
        raise ValueError("calibration requires a nonempty source provenance string")
    if any(not np.isfinite(v) or v <= 0 for v in (
        length_scale_m, energy_scale_J, M_a_phys_m2_per_J_s,
        M_s_phys_m2_per_J_s, model_mobility_a, model_mobility_s,
    )):
        raise ValueError("finite positive coordinate/energy scales and mobilities required")
    t0 = float(model_mobility_a) * float(length_scale_m) ** 2 / (
        float(M_a_phys_m2_per_J_s) * float(energy_scale_J)
    )
    calibration = PhysicalTimeCalibration(
        calibrated=True,
        source=str(source),
        temperature_K=None if temperature_K is None else float(temperature_K),
        length_scale_m=float(length_scale_m),
        energy_scale_J=float(energy_scale_J),
        M_a_phys_m2_per_J_s=float(M_a_phys_m2_per_J_s),
        M_s_phys_m2_per_J_s=float(M_s_phys_m2_per_J_s),
        t0_seconds=t0,
        mobility_ratio=float(M_s_phys_m2_per_J_s) / float(M_a_phys_m2_per_J_s),
        notes=notes,
        model_mobility_a=float(model_mobility_a),
        model_mobility_s=float(model_mobility_s),
    )
    calibration.validate()
    return calibration


def load_time_calibration(path: str | Path) -> PhysicalTimeCalibration:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    calibration = PhysicalTimeCalibration(**data)
    calibration.validate()
    return calibration


def save_time_calibration(
    calibration: PhysicalTimeCalibration, path: str | Path
) -> None:
    calibration.validate()
    Path(path).write_text(
        json.dumps(calibration.to_dict(), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def _scaled(value, factor: float):
    array = np.asarray(value, dtype=float) * float(factor)
    return float(array) if array.ndim == 0 else array


def model_time_to_seconds(model_time, calibration: PhysicalTimeCalibration):
    calibration.require_calibrated()
    return _scaled(model_time, float(calibration.t0_seconds))


def seconds_to_model_time(seconds, calibration: PhysicalTimeCalibration):
    calibration.require_calibrated()
    return _scaled(seconds, 1.0 / float(calibration.t0_seconds))


def model_frequency_to_hz(model_frequency, calibration: PhysicalTimeCalibration):
    calibration.require_calibrated()
    return _scaled(model_frequency, 1.0 / float(calibration.t0_seconds))


def hz_to_model_frequency(frequency_hz, calibration: PhysicalTimeCalibration):
    calibration.require_calibrated()
    return _scaled(frequency_hz, float(calibration.t0_seconds))


def mobility_from_collective_diffusivity(
    diffusivity_m2_per_s: float, temperature_K: float
) -> float:
    """Einstein route ``M=D/(k_B T)`` for the same collective coordinate."""

    diffusivity = float(diffusivity_m2_per_s)
    temperature = float(temperature_K)
    if diffusivity <= 0.0 or not np.isfinite(diffusivity):
        raise ValueError("collective-coordinate diffusivity must be positive")
    if temperature <= 0.0 or not np.isfinite(temperature):
        raise ValueError("temperature must be positive")
    return diffusivity / (BOLTZMANN_J_PER_K * temperature)


def mobility_from_decoupled_relaxation(
    relaxation_seconds: float, physical_curvature_J_per_m2: float
) -> float:
    """Return ``M=1/(tau H)`` only for a justified decoupled mode."""

    tau = float(relaxation_seconds)
    curvature = float(physical_curvature_J_per_m2)
    if not np.all(np.isfinite([tau, curvature])) or tau <= 0.0 or curvature <= 0.0:
        raise ValueError("relaxation time and curvature must be positive")
    return 1.0 / (tau * curvature)


def infer_diagonal_mobilities_from_relaxation_times(
    physical_hessian_J_per_m2,
    *,
    tau_fast_seconds: float,
    tau_slow_seconds: float,
    initial_mobilities_m2_per_J_s: tuple[float, float],
) -> MobilityInference:
    """Fit diagonal mobilities to two coupled relaxation rates.

    Relaxation times without mode shapes can admit a coordinate-exchange
    branch. ``initial_mobilities`` selects a deterministic branch; MD mode
    eigenvectors or a fitted relaxation matrix are required to remove that
    ambiguity physically.
    """

    hessian = np.asarray(physical_hessian_J_per_m2, dtype=float)
    if hessian.shape != (2, 2) or not np.allclose(hessian, hessian.T):
        raise ValueError("physical_hessian_J_per_m2 must be symmetric 2x2")
    if np.min(np.linalg.eigvalsh(hessian)) <= 0.0:
        raise ValueError("physical Hessian must be positive definite")
    taus = np.array([tau_slow_seconds, tau_fast_seconds], dtype=float)
    if np.any(~np.isfinite(taus)) or np.any(taus <= 0.0):
        raise ValueError("relaxation times must be finite and positive")
    target = np.sort(1.0 / taus)
    initial = np.asarray(initial_mobilities_m2_per_J_s, dtype=float)
    if initial.shape != (2,) or np.any(initial <= 0.0):
        raise ValueError("two positive initial mobilities are required")

    def rates(log_mobility: np.ndarray) -> np.ndarray:
        mobility = np.exp(log_mobility)
        sqrt_m = np.diag(np.sqrt(mobility))
        return np.linalg.eigvalsh(sqrt_m @ hessian @ sqrt_m)

    solution = least_squares(
        lambda log_m: np.log(rates(log_m) / target),
        np.log(initial),
        method="trf",
        ftol=1.0e-13,
        xtol=1.0e-13,
        gtol=1.0e-13,
    )
    mobility = np.exp(solution.x)
    fitted = rates(solution.x)
    residual = float(np.max(np.abs(fitted / target - 1.0)))
    return MobilityInference(
        M_a_phys_m2_per_J_s=float(mobility[0]),
        M_s_phys_m2_per_J_s=float(mobility[1]),
        fitted_rates_per_s=fitted,
        target_rates_per_s=target,
        relative_residual=residual,
        note="times-only branch; require mode shapes to remove exchange ambiguity",
    )


def physical_relaxation_spectrum(
    reduced_hessian,
    calibration: PhysicalTimeCalibration,
) -> dict[str, np.ndarray | float]:
    """Return coupled physical rates using the declared scales and mobilities."""

    calibration.require_calibrated()
    h_star = np.asarray(reduced_hessian, dtype=float)
    if (h_star.shape != (2, 2) or not np.all(np.isfinite(h_star))
            or not np.allclose(h_star, h_star.T)
            or np.min(np.linalg.eigvalsh(h_star)) <= 0.0):
        raise ValueError("reduced_hessian must be finite symmetric positive definite 2x2")
    h_phys = calibration.energy_scale_J / calibration.length_scale_m**2 * h_star
    mobility = np.diag(
        [calibration.M_a_phys_m2_per_J_s, calibration.M_s_phys_m2_per_J_s]
    ).astype(float)
    sqrt_m = np.diag(np.sqrt(np.diag(mobility)))
    rates = np.linalg.eigvalsh(sqrt_m @ h_phys @ sqrt_m)
    return {
        "rates_per_s": rates,
        "tau_seconds": 1.0 / rates,
        "physical_hessian_J_per_m2": h_phys,
    }
