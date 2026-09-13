"""Smooth joint epigraph constraints, without differentiating a profiled LP."""
import numpy as np


def constraints_and_jacobian(matrix, derivatives, targets, scales, exact, tested,
                             coefficients, coefficient_scales, eta):
    """Variables are (shape, c/D, eta); inequalities use >=0 convention."""
    M = np.asarray(matrix, float)
    dM = np.asarray(derivatives, float)
    c = np.asarray(coefficients, float)
    D = np.asarray(coefficient_scales, float)
    scales = np.asarray(scales, float)
    r = (M @ c - targets) / scales
    jr = np.column_stack([(dM @ c).T / scales[:, None],
                          M * D / scales[:, None], np.zeros(len(M))])
    eq = r[exact]
    inequality = np.r_[eta-r[tested], eta+r[tested]]
    ji = np.vstack([-jr[tested], jr[tested]])
    ji[:, -1] = 1.
    return eq, inequality, jr[exact], ji
