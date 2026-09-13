"""Unrestricted nonnegative density exponent, separate static research probe.

The previous [-1,1] interpolation was a search convention, not a physical
upper bound. This module tests x**p |Q1|**2 for p>=0 without changing v19 or
production. It retains per-site infinite density/moment jets and their gauge.
"""
import numpy as np
from .vector_interface_reference import HESSIAN_INDICES


def positive_power_site_jet(density_changes, moment_jets, exponent):
    density, moments = np.asarray(density_changes, float), np.asarray(moment_jets, float)
    p = float(exponent)
    if (not np.isfinite(p) or p < 0 or density.ndim != 2 or density.shape[1] != 10
            or moments.ndim != 3 or moments.shape[1:] != (10, 3)
            or not len(density) or not len(moments)
            or np.any(~np.isfinite(density)) or np.any(~np.isfinite(moments))):
        raise ValueError('finite site jets and nonnegative exponent required')
    count = max(len(density), len(moments))
    x = np.zeros((count, 10)); x[:len(density)] = density; x[:, 0] += 1
    if np.any(x[:, 0] <= 0):
        raise ValueError('strictly positive environmental density required')
    q = np.zeros((count, 10, 3)); q[:len(moments)] = moments
    Q, Qp, Qpq = q[:, 0], q[:, 1:4], q[:, HESSIAN_INDICES]
    I = np.einsum('nc,nc->n', Q, Q)
    Ip = 2*np.einsum('nc,nic->ni', Q, Qp)
    Ipp = 2*(np.einsum('nic,njc->nij', Qp, Qp)+np.einsum('nc,nijc->nij', Q, Qpq))
    g = np.exp(p*np.log(x[:, 0]))
    gx, gxx = p*g/x[:, 0], p*(p-1)*g/x[:, 0]**2
    xp, xpp = x[:, 1:4], x[:, HESSIAN_INDICES]
    gradient = np.einsum('n,ni->i', g, Ip)+np.einsum('n,n,ni->i', gx, I, xp)
    H = (np.einsum('n,nij->ij', g, Ipp)+np.einsum('n,ni,nj->ij', gx, xp, Ip)
         +np.einsum('n,ni,nj->ij', gx, Ip, xp)+np.einsum('n,n,nij->ij', gx, I, xpp)
         +np.einsum('n,n,ni,nj->ij', gxx, I, xp, xp))
    result = np.r_[g@I, gradient, H[0], H[1, 1:], H[2, 2]]
    if np.any(~np.isfinite(result)):
        raise ArithmeticError('unresolved/overflowed analytic site jet')
    return result
