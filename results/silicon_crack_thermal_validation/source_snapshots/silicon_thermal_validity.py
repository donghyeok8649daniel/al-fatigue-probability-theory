"""Constrained harmonic thermal diagnostics; no material or kinetic calibration."""
from __future__ import annotations

import numpy as np
from scipy.constants import Boltzmann, atomic_mass, electron_volt, hbar

KB = Boltzmann / electron_volt
MASS_UNIT = atomic_mass * 1e4 / electron_volt  # amu -> eV ps^2 / Angstrom^2
HBAR = hbar / electron_volt * 1e12  # eV ps


def oscillator_statistics(curvatures, temperature, mass_amu=28.0855):
    """Equal-mass orthonormal Cartesian modes; covariance is displacement^2.

    Quantum results include zero-point variance, not an effective temperature.
    Positive modes are required; no clipping of unstable or constrained modes.
    """
    lam = np.asarray(curvatures, float)
    if (lam.ndim != 1 or not len(lam) or not np.isfinite(lam).all()
            or np.any(lam <= 0) or not np.isfinite(temperature) or temperature <= 0
            or not np.isfinite(mass_amu) or mass_amu <= 0):
        raise ValueError('positive finite modes, temperature and mass required')
    kt = KB * temperature
    omega = np.sqrt(lam / (mass_amu * MASS_UNIT))
    energy = HBAR * omega
    x = energy / (2 * kt)
    ratio = x / np.tanh(x)
    classical = kt / lam
    quantum = classical * ratio
    f_quantum = energy / 2 + kt * np.log(-np.expm1(-energy / kt))
    f_classical = kt * np.log(energy / kt)
    return dict(omega_rad_ps=omega, quantum_energy_eV=energy,
                classical_variance_A2=classical, quantum_variance_A2=quantum,
                variance_ratio=ratio, quantum_free_energy_eV=f_quantum,
                classical_free_energy_eV=f_classical)


def validate_fixed_cartesian(hessian, basis, free):
    """Reject moving grips and nonorthonormal coordinates before assigning mass."""
    h, b, free = np.asarray(hessian, float), np.asarray(basis, float), np.asarray(free)
    if free.ndim != 1 or free.dtype != np.bool_ or not free.any():
        raise ValueError('boolean free-atom mask required')
    n = 3 * int(free.sum())
    if (h.shape != (n, n) or b.shape != (3 * len(free), n)
            or not np.isfinite(h).all() or not np.isfinite(b).all()):
        raise ValueError('complete fixed-grip Cartesian basis required')
    ortho = float(np.max(abs(b.T @ b - np.eye(n))))
    fixed = float(np.max(abs(b.reshape(len(free), 3, n)[~free]), initial=0))
    asymmetry = float(np.max(abs(h - h.T)))
    if ortho > 1e-10 or fixed > 1e-12 or asymmetry > 1e-8:
        raise ValueError('invalid fixed-grip basis or asymmetric Hessian')
    return (h + h.T) / 2, dict(basis_orthogonality_error=ortho,
                              fixed_basis_maximum=fixed, raw_asymmetry=asymmetry)
