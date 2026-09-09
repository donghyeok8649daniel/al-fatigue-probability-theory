"""Relax the normal coordinate along a declared scalar registry path.

The coupled Schur curvature locates a shear-controlled local instability on
that branch. A fixed-gap traction maximum is NOT a macroscopic yield stress.
No energy, mobility, area, temperature, or physical stress is rescaled to fit.
"""
from functools import lru_cache

import numpy as np
from scipy.optimize import least_squares, brentq


def normal_relaxed_registry(evaluate, *, h, period, units, samples=65,
                            normal_traction_mpa=0., force_tolerance=2e-9):
    if samples < 9 or h <= 0 or period <= 0:
        raise ValueError("positive geometry and at least nine samples required")
    fn = float(units.traction_mpa_to_force(normal_traction_mpa))
    rows = []; previous_a = h
    @lru_cache(maxsize=8192)
    def value(a, s):
        return np.asarray(evaluate(float(a), float(s)), dtype=float)
    def solve(s, seed):
        sol = least_squares(lambda q: [value(q[0],s)[1]-fn], [seed],
            jac=lambda q: [[value(q[0],s)[3]]], bounds=([.60*h],[4*h]),
            ftol=2e-13, xtol=2e-13, gtol=2e-13, max_nfev=80)
        v = value(sol.x[0],s)
        if not sol.success or abs(v[1]-fn)>force_tolerance or v[3]<=0:
            raise ValueError(f"normal branch unresolved at s={s}; not a proven spinodal")
        return float(sol.x[0]), v
    for s in np.linspace(0., .5*period, samples):
        a, v = solve(float(s), previous_a)
        fixed = value(h,float(s))
        rows.append(dict(s_reduced=float(s), a_reduced=a,
            relaxed_energy_ev=float(v[0]-fn*(a-h)), fixed_energy_ev=float(fixed[0]),
            relaxed_shear_mpa=float(units.force_to_traction_mpa(v[2])),
            fixed_shear_mpa=float(units.force_to_traction_mpa(fixed[2])),
            schur_curvature=float(v[5]-v[4]**2/v[3]),
            normal_curvature=float(v[3]), coupling=float(v[4]),
            normal_force_residual=float(abs(v[1]-fn))))
        previous_a = a
    # Locate FIRST local shear spinodal reached from pristine, not a global
    # traction maximum on another branch. Repeat with sample/tolerance refinement.
    spinodal = None
    for left, right in zip(rows[:-1], rows[1:]):
        if left["schur_curvature"]>0 and right["schur_curvature"]<0:
            seed = .5*(left["a_reduced"]+right["a_reduced"])
            def schur(s):
                a, v = solve(float(s),seed)
                return float(v[5]-v[4]**2/v[3])
            s = brentq(schur,left["s_reduced"],right["s_reduced"],xtol=2e-12)
            a,v = solve(s,seed)
            spinodal = dict(s_reduced=float(s),a_reduced=a,
                shear_mpa=float(units.force_to_traction_mpa(v[2])),
                energy_ev=float(v[0]-fn*(a-h)),
                normal_force_residual=float(abs(v[1]-fn)),
                schur_curvature=float(v[5]-v[4]**2/v[3]),
                status="first scalar-path shear spinodal with relaxed normal coordinate")
            break
    return rows, spinodal
