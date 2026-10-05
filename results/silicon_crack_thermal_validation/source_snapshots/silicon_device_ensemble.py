"""Research-only finite loading-device and constrained canonical ensembles.

No device stiffness, atomistic sample, or Monte Carlo index calibrates the
production mobility/clock. Energies are total eV, lengths Angstrom.
"""
from __future__ import annotations
import numpy as np
from scipy.linalg import solve
from scipy.constants import electron_volt
from scipy.spatial.distance import pdist
from .silicon_conditional_research import KB_EV_K
from .silicon_thermal_research import GaussianReference

EV_A2_TO_N_MM = electron_volt / 1e-20 / 1000.


def device_response(hessian, extension_direction, stiffness_eV_A2, temperature_K):
    """Local harmonic response at the SAME equilibrium; Delta = v.T x.

    Add K/2 (X-Delta)^2. At X=Delta_equilibrium the curvature is H+K vv.T.
    K=0 gives the dead-load local ensemble, not a globally normalizable
    tensile-force ensemble. Infinite K is the constrained Gaussian limit.
    """
    h, v = np.asarray(hessian, float), np.asarray(extension_direction, float)
    k, t = float(stiffness_eV_A2), float(temperature_K)
    if (h.ndim != 2 or h.shape != (v.size, v.size) or v.ndim != 1 or not v.size
            or not np.isfinite(h).all() or not np.isfinite(v).all() or not np.any(v)
            or not np.allclose(h, h.T, rtol=0, atol=1e-10)
            or np.isnan(k) or k < 0 or not np.isfinite(t) or t <= 0):
        raise ValueError('symmetric finite positive Hessian, nonzero coordinate, K>=0 and T>0 required')
    np.linalg.cholesky(h)
    compliance_vector = solve(h, v, assume_a='pos')
    compliance = float(v @ compliance_vector)
    if np.isinf(k):
        transmission, variance = 1., 0.
    else:
        transmission = k*compliance/(1+k*compliance)
        variance = KB_EV_K*t*compliance/(1+k*compliance)
    return dict(specimen_stiffness_eV_A2=1/compliance,
                specimen_stiffness_N_mm=EV_A2_TO_N_MM/compliance,
                extension_variance_A2=variance,
                d_extension_d_actuator=transmission,
                d_force_d_actuator_eV_A2=transmission/compliance)


class FixedGripTarget:
    """Exact nonlinear energy correction for a constant-Jacobian Gaussian map.

    A shared absolute Cartesian box and minimum pair distance define the
    computational domain. They are neither an intact basin nor a crack test.
    The evaluator supplies actual U/forces; no quadratic energy substitution.
    """
    def __init__(self, positions, basis, hessian, gradient, temperature_K,
                 common_center, halfwidth_A=3., minimum_pair_A=1.):
        self.positions = np.asarray(positions, float).copy()
        self.basis = np.asarray(basis, float).copy()
        self.center = np.asarray(common_center, float).copy()
        self.halfwidth = float(halfwidth_A)
        self.minimum_pair = float(minimum_pair_A)
        h = np.asarray(hessian, float)
        g = np.asarray(gradient, float)
        if (self.positions.ndim != 2 or self.positions.shape[1] != 3
                or self.center.shape != self.positions.shape
                or self.basis.shape != (self.positions.size, len(g))
                or not np.isfinite(self.positions).all() or not np.isfinite(self.center).all()
                or not np.isfinite(self.basis).all()
                or not np.allclose(self.basis.T@self.basis, np.eye(len(g)), atol=1e-12, rtol=0)
                or not np.isfinite([self.halfwidth, self.minimum_pair]).all()
                or self.halfwidth <= 0 or self.minimum_pair <= 0):
            raise ValueError('finite orthonormal fixed Cartesian geometry/domain required')
        self.free = np.any(self.basis.reshape(len(self.positions), 3, -1) != 0, axis=(1, 2))
        if not np.any(self.free) or not np.allclose(self.center[~self.free], self.positions[~self.free], atol=1e-12, rtol=0):
            raise ValueError('common domain must preserve the same fixed grips')
        self.reference = GaussianReference.from_hessian(h, temperature_K, -solve(h, g, assume_a='pos'))

    def coordinates(self, u):
        d = self.reference.transform(u)
        return self.positions+(self.basis@d).reshape(self.positions.shape)

    def inside(self, r):
        r = np.asarray(r, float)
        return (r.shape == self.positions.shape and np.isfinite(r).all()
                and np.array_equal(r[~self.free], self.positions[~self.free])
                and np.max(abs(r[self.free]-self.center[self.free])) < self.halfwidth
                and pdist(r).min() >= self.minimum_pair)

    def correction(self, u, actual_energy_eV, reference_energy_eV):
        u = np.asarray(u, float)
        return float((actual_energy_eV-reference_energy_eV)/self.reference.thermal_energy-.5*u@u)
