"""Dual envelope derivative of the EXISTING coefficient minimax problem.

This differentiates a calibration objective, not a new energy/force law.
At active-set changes a single dual solution need not define a unique
derivative; independently reprofile finite perturbations before using it.
"""
import numpy as np


def envelope_gradient(matrix_derivatives, observations, profile):
    if not profile.get('completed') or not profile.get('certified_lower_bound_available'):
        raise ValueError('certified coefficient profile required')
    derivative = np.asarray(matrix_derivatives, float)
    c = np.asarray(profile['coefficients'], float)
    scales = np.array([o.scale for o in observations])
    if (derivative.ndim != 3 or derivative.shape[1:] != (len(observations), len(c))
            or np.any(~np.isfinite(derivative)) or np.any(scales <= 0)):
        raise ValueError('finite (shape coordinate, observable, coefficient) derivative required')
    exact = np.asarray(profile['exact_rows'], int)
    rows = np.asarray(profile['tested_rows'], int)
    de = np.einsum('kij,j->ki', derivative[:, exact], c)/scales[exact]
    dr = np.einsum('kij,j->ki', derivative[:, rows], c)/scales[rows]
    eq = np.asarray(profile['equality_duals'])
    net = np.asarray(profile['upper_error_duals'])-np.asarray(profile['lower_error_duals'])
    return -de@eq-dr@net


def difference_stencil(x, index, lower, upper, step):
    """Second-order central or one-sided derivative within declared bounds."""
    x, lower, upper = (np.asarray(v, float) for v in (x, lower, upper))
    if (x.ndim != 1 or lower.shape != x.shape or upper.shape != x.shape
            or not 0 <= index < len(x) or not np.isfinite(step) or step <= 0
            or np.any(~np.isfinite(np.r_[x, lower, upper]))
            or np.any(x < lower) or np.any(x > upper) or np.any(lower >= upper)):
        raise ValueError('finite feasible point, ordered bounds and positive step required')
    h = min(step, (upper[index]-lower[index])/4)
    if x[index]-h >= lower[index] and x[index]+h <= upper[index]:
        offsets, weights = (-h,h), (-1/(2*h),1/(2*h))
    elif x[index]+2*h <= upper[index]:
        offsets, weights = (0.,h,2*h), (-3/(2*h),2/h,-1/(2*h))
    else:
        offsets, weights = (0.,-h,-2*h), (3/(2*h),-2/h,1/(2*h))
    result=[]
    for offset, weight in zip(offsets, weights):
        point=x.copy();point[index]+=offset
        result.append((point,weight))
    return result
