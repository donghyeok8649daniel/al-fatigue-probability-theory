"""Explicit calibration controls/status, not additional atomistic observations."""
import numpy as np
from scipy.optimize import minimize

from .vector_material_calibration import MaterialObservation


def append_fixed_pair(matrix, observations, control):
    """Retain an already specified LJ pair for a declared embedding-only test.

    These are exact CONTROL constraints, not two extra measured Al targets.
    No lower bound is invented to rescue an optimizer's zero-pair closure.
    """
    if control is None:
        return matrix, observations
    pair = np.asarray(control['coefficients'], float)
    if pair.shape != (2,) or np.any(~np.isfinite(pair)) or np.any(pair <= 0):
        raise ValueError('two finite positive inherited LJ coefficients required')
    if not control.get('calibration_sha256'):
        raise ValueError('inherited pair provenance binding required')
    rows = [MaterialObservation('control_pair_'+name, value, value, 'eV/basis', 'exact')
            for name, value in zip(('u', 'v'), pair)]
    return np.vstack([matrix, np.eye(matrix.shape[1])[:2]]), list(observations)+rows


def select_profile_result(profiles):
    """A zero-pair closure is diagnostic only, NEVER an accepted LJ model."""
    if not profiles:
        raise ValueError('no successfully certified coefficient profiles')
    admissible = [r for r in profiles if r['strictly_positive_LJ']]
    return (min(admissible or profiles, key=lambda r: r['squared_loss']), bool(admissible))


def feasible_simplex_search(residual_function, initial, bounds, *, max_evaluations):
    """Deterministic bounded search on the VERIFIED profile's domain.

    Undefined/infeasible coefficient profiles return the extended-real +inf
    indicator, never a finite penalty or a repaired unstable parameter vector.
    Caller retains every exception/certificate. This is an optional numerical
    continuation after a trust-region trial left the exact-constraint domain.
    Stopping success is not proof of a global optimum or physical acceptance.
    """
    initial = np.asarray(initial, float)
    limits = np.asarray(bounds, float)
    if (limits.shape != (2, len(initial)) or not np.isfinite(limits).all()
            or np.any(limits[0] >= limits[1]) or not np.isfinite(initial).all()
            or np.any(initial < limits[0]) or np.any(initial > limits[1])
            or int(max_evaluations) != max_evaluations or max_evaluations < len(initial)+1):
        raise ValueError('finite bounds, interior starting point and sufficient evaluation budget required')
    # Explicit valid initial point prevents an all-infinite starting simplex.
    baseline = np.asarray(residual_function(initial), float)
    if not np.isfinite(baseline).all():
        raise ValueError('verified finite starting residual required')
    rejected = []

    def objective(x):
        try:
            residual = np.asarray(residual_function(x), float)
            if not np.isfinite(residual).all():
                raise ArithmeticError('nonfinite coefficient-profile residual')
            return float(residual@residual)
        except (ValueError, ArithmeticError) as error:
            rejected.append(dict(log_shape=np.array(x).tolist(), message=str(error)))
            return np.inf

    simplex = np.tile(initial, (len(initial)+1, 1))
    for i in range(len(initial)):
        direction = 1 if limits[1, i]-initial[i] >= initial[i]-limits[0, i] else -1
        available = (limits[1, i]-initial[i]) if direction > 0 else (initial[i]-limits[0, i])
        simplex[i+1, i] += direction*min(.025, available/2)
    result = minimize(objective, initial, method='Nelder-Mead', bounds=limits.T,
        options=dict(initial_simplex=simplex, adaptive=True, maxfev=max_evaluations,
                     xatol=1e-6, fatol=1e-6))
    result.rejected_profiles = rejected
    return result
