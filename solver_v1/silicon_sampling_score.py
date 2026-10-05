"""Exact force mapping for a constant-Jacobian Gaussian proposal reference.

Research integration coordinates and this dimensionless score are not physical
atomic momentum, a Langevin mobility, or a fracture counter.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import solve_triangular


def nonlinear_correction_gradient(target, coordinates, total_forces_eV_A):
    """d[(U-Uref)/kBT-u.u/2]/du for r=rref+B(mu+sqrt(kBT)L^-T u).

    U is TOTAL energy; F=-grad_r(U). No per-atom rescaling, clipping or force
    normalization. Fixed Cartesian rows of B are zero and exert no contribution.
    """
    u, forces = np.asarray(coordinates, float), np.asarray(total_forces_eV_A, float)
    if (u.shape != (target.basis.shape[1],) or forces.shape != target.positions.shape
            or not np.isfinite(u).all() or not np.isfinite(forces).all()):
        raise ValueError('finite reference coordinates and total-energy Cartesian forces required')
    projected = -target.basis.T@forces.ravel()
    return (solve_triangular(target.reference.lower, projected, lower=True,
                             check_finite=False)/np.sqrt(target.reference.thermal_energy)-u)
