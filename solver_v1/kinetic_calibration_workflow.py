"""Evidence-bound kinetic import; no default Al clock or trajectory simulation.

The correlation route consumes measured collective-coordinate covariance, not
bulk self diffusion or dislocation-centre mobility. Synthetic data belong in
tests only. Validation cannot establish the truth of a user-supplied source.
"""
from __future__ import annotations

from dataclasses import asdict, replace
import hashlib
import json

import numpy as np
from scipy.linalg import expm, logm

from .energy_model_registry import (
    AL_TARGET_BEST_FEASIBLE, TWO_ROW_LJ_REFERENCE, build_energy_model,
)
from .physical_time import (
    BOLTZMANN_J_PER_K, PhysicalTimeCalibration, calibration_from_mobilities,
)

EV_J = 1.602176634e-19
REDUCED_COORDINATES = "a_phys=L0*a; s_phys=L0*s; one reduced two-row cell"


def energy_fingerprint(model_id: str) -> str:
    """Static parameter signature; temperature/mobility are not energy refits."""
    model = build_energy_model(model_id)
    parameters = {key: getattr(model.p, key) for key in
                  ("n_cells", "b", "epsilon", "sigma_lj", "chi_axial_projection")}
    data = dict(model_id=model_id, parameters=parameters)
    if hasattr(model, "embedding"):
        data.update(embedding=asdict(model.embedding), density=asdict(model.density_params))
    return hashlib.sha256(json.dumps(data, sort_keys=True).encode()).hexdigest()


def bind_to_energy_model(calibration: PhysicalTimeCalibration, model_id: str):
    """Explicit declaration by the data supplier, NOT automatic physical proof."""
    bound = replace(calibration, energy_model_id=model_id,
                    energy_model_fingerprint=energy_fingerprint(model_id),
                    coordinate_definition=REDUCED_COORDINATES)
    validate_for_energy_model(bound, model_id)
    return bound


def validate_for_energy_model(calibration: PhysicalTimeCalibration, model_id: str):
    calibration.require_calibrated()
    if calibration.energy_model_id != model_id:
        raise ValueError("kinetic calibration energy_model_id does not match selected model")
    if calibration.energy_model_fingerprint != energy_fingerprint(model_id):
        raise ValueError("kinetic calibration static energy fingerprint is stale or missing")
    if calibration.coordinate_definition != REDUCED_COORDINATES:
        raise ValueError("kinetic data must describe the same reduced a,s cell coordinates")
    if calibration.temperature_K is None:
        raise ValueError("kinetic calibration requires temperature_K")
    if model_id != TWO_ROW_LJ_REFERENCE and not np.isclose(
        calibration.energy_scale_J, EV_J, rtol=1e-10, atol=0.
    ):
        raise ValueError("hybrid energy numbers use one eV; energy_scale_J must match")
    if model_id == AL_TARGET_BEST_FEASIBLE and not np.isclose(
        calibration.length_scale_m, 4.05 / np.sqrt(2) * 1e-10, rtol=1e-10, atol=0.
    ):
        raise ValueError("Al-target reduced geometry length_scale_m does not match")


def build_time_basis_model(model_id, *, time_basis="model", calibration=None):
    """Same PDE/energy; physical mode uses the declared T and mobility ratio.

    Model mode is exactly the historical M=(1,.05), kT=.02 path. A loaded
    calibration never silently changes an existing model-time result.
    """
    if time_basis == "model":
        return build_energy_model(model_id, chi=.20, kT=.02)
    if time_basis != "physical" or calibration is None:
        raise ValueError("physical time requires a kinetic calibration")
    validate_for_energy_model(calibration, model_id)
    kT = BOLTZMANN_J_PER_K * calibration.temperature_K / calibration.energy_scale_J
    model = build_energy_model(model_id, chi=.20, kT=kT)
    model.p.mobility_a = calibration.model_mobility_a
    model.p.mobility_s = calibration.model_mobility_s
    return model


def infer_from_correlation_matrices(times_seconds, covariance_m2, hessian_J_m2, *,
                                    temperature_K, relative_tolerance):
    """C(t)=exp(-B t) C(0), B=M H, for a stationary harmonic Markov process.

    Input rows include t=0 and at least two positive lags. Column ordering is
    always (a,s). No independent fits of coupled raw coordinate traces. Reject
    oscillatory/non-Markov fits, non-diagonal mobility and equilibrium mismatch.
    ``relative_tolerance`` is a declared measurement acceptance tolerance, not
    a fabricated material uncertainty. Block/window uncertainty must be assessed
    externally on the original MD data before adopting the output.
    """
    t = np.asarray(times_seconds, dtype=float)
    C = np.asarray(covariance_m2, dtype=float)
    H = np.asarray(hessian_J_m2, dtype=float)
    if (t.ndim != 1 or len(t) < 3 or t[0] != 0. or np.any(np.diff(t) <= 0)
            or not np.all(np.isfinite(t)) or C.shape != (len(t), 2, 2)
            or not np.all(np.isfinite(C))):
        raise ValueError("finite increasing lag times starting at zero and (lag,2,2) covariance required")
    if (H.shape != (2, 2) or not np.all(np.isfinite(H))
            or not np.allclose(H, H.T) or np.min(np.linalg.eigvalsh(H)) <= 0):
        raise ValueError("physical Hessian must be symmetric positive definite")
    if not (np.isfinite(temperature_K) and temperature_K > 0 and
            np.isfinite(relative_tolerance) and 0 < relative_tolerance < 1):
        raise ValueError("positive temperature and explicit relative tolerance in (0,1) required")
    if (not np.allclose(C[0], C[0].T, rtol=relative_tolerance, atol=0.)
            or np.min(np.linalg.eigvalsh(C[0])) <= 0):
        raise ValueError("C(0) must be symmetric positive definite")
    expected_C0 = BOLTZMANN_J_PER_K * temperature_K * np.linalg.inv(H)
    equilibrium_error = np.linalg.norm(C[0] - expected_C0) / np.linalg.norm(expected_C0)
    matrices = []
    for lag, current in zip(t[1:], C[1:]):
        propagator = np.linalg.solve(C[0].T, current.T).T
        rate = -logm(propagator) / lag
        if np.iscomplexobj(rate) and np.max(abs(rate.imag)) > 1e-10 * np.linalg.norm(rate):
            raise ValueError("correlation has no real overdamped relaxation logarithm")
        matrices.append(rate.real)
    B = np.mean(matrices, axis=0)
    M = np.linalg.solve(H.T, B.T).T
    diag = np.diag(M)
    if np.any(diag <= 0) or not np.all(np.isfinite(M)):
        raise ValueError("inferred diagonal mobilities must be finite and positive")
    diagonal_error = np.linalg.norm(M - np.diag(diag)) / np.linalg.norm(M)
    rates = np.linalg.eigvalsh(np.diag(np.sqrt(diag)) @ H @ np.diag(np.sqrt(diag)))
    # Compare to the ACTUAL proposed diagonal-M model, not an unrestricted fit.
    B_diagonal = np.diag(diag) @ H
    errors = [np.linalg.norm(expm(-B_diagonal * lag) @ C[0] - current) /
              np.linalg.norm(C[0]) for lag, current in zip(t[1:], C[1:])]
    lag_error = max(np.linalg.norm(r - B_diagonal) / np.linalg.norm(B_diagonal)
                    for r in matrices)
    diagnostics = dict(equilibrium_relative_error=float(equilibrium_error),
                       off_diagonal_mobility_relative_error=float(diagonal_error),
                       covariance_relative_error=float(max(errors)),
                       lag_rate_relative_error=float(lag_error))
    if max(diagnostics.values()) > relative_tolerance:
        raise ValueError(f"collective correlation does not support this diagonal Markov model: {diagnostics}")
    return dict(M_a_phys_m2_per_J_s=float(diag[0]), M_s_phys_m2_per_J_s=float(diag[1]),
                tau_fast_seconds=float(1/rates[-1]), tau_slow_seconds=float(1/rates[0]),
                acceptance_tolerance=float(relative_tolerance), **diagnostics)


def calibration_from_collective_correlations(times_seconds, covariance_m2, *, model_id,
        length_scale_m, energy_scale_J, temperature_K, source, relative_tolerance, notes=""):
    model = build_energy_model(model_id)
    H = model.hessian(model.a0, 0.) * energy_scale_J / length_scale_m**2
    fit = infer_from_correlation_matrices(times_seconds, covariance_m2, H,
                temperature_K=temperature_K, relative_tolerance=relative_tolerance)
    ma, ms = fit["M_a_phys_m2_per_J_s"], fit["M_s_phys_m2_per_J_s"]
    # One clock: choose M_a*=1 and retain the MEASURED ratio, not .05 by fiat.
    calibration = calibration_from_mobilities(length_scale_m=length_scale_m,
        energy_scale_J=energy_scale_J, M_a_phys_m2_per_J_s=ma, M_s_phys_m2_per_J_s=ms,
        source=source, temperature_K=temperature_K, model_mobility_a=1.,
        model_mobility_s=ms/ma, notes=notes)
    return bind_to_energy_model(calibration, model_id), fit
