"""Positive spectral quadrature of finite harmonic memory; no damping fit.

Ranks measure representational accuracy, not the number of calibrated material
parameters. Coefficients are derived from the same full bath Hessian.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import eigh

from .silicon_conditional_research import SI_MASS_AMU,AMU_EV_PS2_A2


def positive_memory_quadrature(eigenvalues,weights,rank):
    """Lanczos/Gauss quadrature of the positive scalar memory spectral measure.

    Full reorthogonalization controls roundoff. No pole or weight is fitted or
    clipped. Nodes must stay positive; finite-rank truncation is independently
    compared on the full observation interval.
    """
    values,weights=np.asarray(eigenvalues,float),np.asarray(weights,float)
    if (values.ndim!=1 or values.shape!=weights.shape or not np.isfinite(values).all()
            or not np.isfinite(weights).all() or np.any(values<=0) or np.any(weights<0)
            or weights.sum()<=0 or int(rank)!=rank or not 1<=rank<=len(values)):
        raise ValueError('positive bath spectrum and nonnegative nonzero weights required')
    total=float(weights.sum())
    vectors=[]
    q=np.sqrt(weights/total)
    alphas,betas=[],[]
    previous=np.zeros_like(q);beta=0.
    for index in range(rank):
        vectors.append(q.copy())
        z=values*q-beta*previous
        alpha=float(q@z);z-=alpha*q
        # Two passes keep the Krylov basis orthogonal near spectral clusters.
        basis=np.asarray(vectors)
        for _ in range(2):
            z-=basis.T@(basis@z)
        alphas.append(alpha)
        next_beta=float(np.linalg.norm(z))
        if index==rank-1 or next_beta<=1e-13*max(values):
            break
        betas.append(next_beta)
        previous,q=q,z/next_beta
        beta=next_beta
    matrix=np.diag(alphas)+np.diag(betas,1)+np.diag(betas,-1)
    nodes,rotation=eigh(matrix)
    if nodes[0]<=0:
        raise np.linalg.LinAlgError('nonpositive reduced bath pole')
    quadrature_weights=total*rotation[0]**2
    return dict(nodes_eV_A2=nodes,weights_eV_A2=quadrature_weights,rank=len(nodes),
        requested_rank=int(rank),total_weight_eV_A2=total,
        orthogonality_error=float(np.max(abs(basis@basis.T-np.eye(len(basis))))))


def quadrature_memory(nodes,weights,time_ps,*,mass_amu=SI_MASS_AMU):
    nodes,weights,times=map(lambda x:np.asarray(x,float),(nodes,weights,time_ps))
    if (nodes.shape!=weights.shape or nodes.ndim!=1 or np.any(nodes<=0) or np.any(weights<0)
            or times.ndim!=1 or not all(np.isfinite(x).all() for x in (nodes,weights,times))
            or mass_amu<=0):
        raise ValueError('positive finite spectral representation required')
    return np.cos(np.outer(times,np.sqrt(nodes/(mass_amu*AMU_EV_PS2_A2))))@weights


def quadrature_release(nodes,weights,relaxed_curvature,q_mass,time_ps,*,mass_amu=SI_MASS_AMU):
    """Stable retained-coordinate release from a conditional bath minimum.

    The counterterm Hqq=K+sum(weights) preserves EXACT static marginal K while
    the positive oscillator bath supplies memory and canonical FDT. This does
    not assert agreement with a nonlinear finite-temperature atom model.
    """
    nodes,weights=np.asarray(nodes,float),np.asarray(weights,float)
    if relaxed_curvature<=0 or q_mass<=0 or np.any(nodes<=0) or np.any(weights<0):
        raise ValueError('stable quadratic release and positive masses required')
    coupling=np.sqrt(nodes*weights)
    h=np.diag(np.r_[relaxed_curvature+weights.sum(),nodes])
    h[0,1:]=coupling;h[1:,0]=coupling
    masses=np.r_[q_mass,np.full(len(nodes),mass_amu*AMU_EV_PS2_A2)]
    inverse_sqrt=1/np.sqrt(masses)
    dynamical=h*np.outer(inverse_sqrt,inverse_sqrt)
    values,vectors=eigh(dynamical)
    if values[0]<=0:
        raise np.linalg.LinAlgError('unstable reduced release')
    initial=np.r_[1.,-coupling/nodes]
    mode_weights=(inverse_sqrt[0]*vectors[0])*(vectors.T@(initial*np.sqrt(masses)))
    return np.cos(np.outer(np.asarray(time_ps,float),np.sqrt(values)))@mode_weights
