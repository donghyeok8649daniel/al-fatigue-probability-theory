"""Bounded static branch diagnostics for the three-coordinate reference only."""
from functools import lru_cache

import numpy as np
from scipy.optimize import least_squares, minimize


def stationary_state(model, initial, *, force=None, expected_index=None, tolerance=2e-9):
    """Find a LOCAL stationary point and check forces AND the full Morse index.

    Bounded numerical search windows prevent a finite-cutoff SOURCE's flat
    separated region from masquerading as a stationary intact state. Bounds
    are never a physical wall; a solution touching one is rejected.
    """
    initial = np.asarray(initial, dtype=float)
    force = np.zeros(3) if force is None else np.asarray(force, dtype=float)
    if initial.shape != (3,) or force.shape != (3,) or not np.all(np.isfinite([initial, force])):
        raise ValueError('finite three-component initial state and conjugate force required')
    lower = np.r_[.7*model.h, initial[1:]-0.45]
    upper = np.r_[1.7*model.h, initial[1:]+0.45]
    @lru_cache(maxsize=1)
    def value(q):
        return model.evaluate(q)
    solved = least_squares(lambda q: value(tuple(q)).gradient-force, initial,
        jac=lambda q: value(tuple(q)).hessian, bounds=(lower, upper),
        xtol=2e-13, ftol=2e-13, gtol=2e-13, max_nfev=100)
    out = value(tuple(solved.x)); eigenvalues = np.linalg.eigvalsh(out.hessian)
    residual = float(np.max(np.abs(out.gradient-force)))
    index = int(np.count_nonzero(eigenvalues < -tolerance))
    margin = float(np.min(np.r_[solved.x-lower, upper-solved.x]))
    valid = (residual < tolerance and margin > 1e-6
             and min(abs(eigenvalues)) > tolerance
             and (expected_index is None or index == expected_index))
    return dict(q=solved.x, evaluation=out, eigenvalues=eigenvalues, morse_index=index,
                force_residual=residual, search_margin=margin, valid=bool(valid),
                optimizer_success=bool(solved.success), nfev=solved.nfev)


def relax_at_x(model, x, initial=None, *, normal_force=0., transverse_force=0., tolerance=2e-9):
    """Continue local equilibrium of z=(a,y) at fixed registry x.

    This is a local constrained branch, NOT a global minimum over all wells.
    Schur curvature is d^2[W-fa*a-fy*y]/dx^2 along that stationary branch.
    """
    initial = np.array([model.h, 0.]) if initial is None else np.asarray(initial, dtype=float)
    target = np.array([normal_force, transverse_force])
    axes = np.array([0, 2]); lower = np.array([.7*model.h, initial[1]-.4])
    upper = np.array([1.7*model.h, initial[1]+.4])
    @lru_cache(maxsize=1)
    def value(z):
        return model.evaluate([z[0], x, z[1]])
    solved = least_squares(lambda z: value(tuple(z)).gradient[axes]-target, initial,
        jac=lambda z: value(tuple(z)).hessian[np.ix_(axes, axes)], bounds=(lower, upper),
        xtol=3e-12, ftol=3e-12, gtol=3e-12, max_nfev=65)
    out = value(tuple(solved.x)); H = out.hessian
    sub = H[np.ix_(axes, axes)]; residual = float(max(abs(out.gradient[axes]-target)))
    eigen = np.linalg.eigvalsh(sub)
    margin = float(np.min(np.r_[solved.x-lower, upper-solved.x]))
    if residual > tolerance or eigen[0] <= tolerance or margin < 1e-6:
        raise RuntimeError(f'transverse branch unresolved: residual={residual}, eigen={eigen}, margin={margin}')
    derivative = -np.linalg.solve(sub, H[axes, 1])
    schur = float(H[1, 1]+H[1, axes]@derivative)
    return dict(q=np.array([solved.x[0], x, solved.x[1]]), evaluation=out,
                schur_curvature=schur, eliminated_eigenvalues=eigen,
                transverse_force_residual=residual, dq_eliminated_dx=derivative,
                normal_force=normal_force, transverse_force=transverse_force)


def saddle_downhill_endpoints(model, stationary, *, perturbation=.01):
    """Two force-minimized endpoints of the unique negative-curvature mode.

    Establishes basin connectivity for this local saddle, not a global minimum
    energy path or finite dislocation nucleation activation energy.
    """
    if not stationary['valid'] or stationary['morse_index'] != 1:
        raise ValueError('validated index-one saddle required')
    q = stationary['q']; out = stationary['evaluation']
    eigen, vectors = np.linalg.eigh(out.hessian); direction = vectors[:, 0]
    endpoints = []
    for sign in (-1, 1):
        @lru_cache(maxsize=1)
        def value(x):
            return model.evaluate(x)
        bounds = [(.7*model.h, 1.7*model.h), (q[1]-.65, q[1]+.65), (q[2]-.65, q[2]+.65)]
        solved = minimize(lambda x: value(tuple(x)).energy, q+sign*perturbation*direction,
            jac=lambda x: value(tuple(x)).gradient, method='L-BFGS-B', bounds=bounds,
            options=dict(ftol=2e-15, gtol=1e-10, maxiter=180, maxls=35))
        result = stationary_state(model, solved.x, expected_index=0)
        result['energy_decrease'] = out.energy-result['evaluation'].energy
        if not result['valid'] or result['energy_decrease'] <= 0:
            raise RuntimeError('saddle downhill endpoint is not a validated lower minimum')
        endpoints.append(result)
    return endpoints
