"""Research validation of fast electronic equilibration in the common SG law.

Finite charge-sector/coordinate generator, reflecting spatial boundaries.
No Si energy, physical rate, mobility calibration or production registration.
All transitions evolve deterministic probability masses; no Monte Carlo.
"""
from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.special import logsumexp

from .probability_pde_2d import _bernoulli


def sg_chain(potential, *, thermal_energy, spacing, face_mobility):
    """Existing SG flux convention on a uniform reflecting 1D mass grid."""
    f = np.asarray(potential, float)
    kt, dx = float(thermal_energy), float(spacing)
    if (f.ndim != 1 or len(f) < 2 or not np.all(np.isfinite(f))
            or not np.isfinite(kt) or not np.isfinite(dx) or kt <= 0 or dx <= 0):
        raise ValueError('finite 1D energies and positive temperature/spacing required')
    m = np.broadcast_to(np.asarray(face_mobility, float), (len(f)-1,))
    if not np.all(np.isfinite(m)) or np.any(m <= 0):
        raise ValueError('positive finite face mobilities required')
    difference = np.diff(f)/kt
    forward = kt*m/dx**2*_bernoulli(difference)
    reverse = kt*m/dx**2*_bernoulli(-difference)
    left, right = np.arange(len(f)-1), np.arange(1, len(f))
    return sparse.coo_matrix((np.concatenate([forward, -forward, reverse, -reverse]),
                             (np.concatenate([right, left, left, right]),
                              np.concatenate([left, left, right, right]))),
                            shape=(len(f), len(f))).tocsr()


def coupled_charge_generator(free_energies, *, excess_counts, chemical_potential,
                             thermal_energy, spacing, mobilities, charge_rate):
    """q-grid SG plus local charge relaxation rate*(w_N*sum_M P_M-P_N).

    Each row is a distinct integer electron sector of the SAME fixed chemical
    configuration. The chosen reset kinetics obey detailed balance but are
    synthetic. Actual electron kinetics need independent data.
    """
    f = np.asarray(free_energies, float)
    n = np.asarray(excess_counts, float)
    mu, kt, rate = map(float, [chemical_potential, thermal_energy, charge_rate])
    if (f.ndim != 2 or min(f.shape) < 1 or f.shape[1] < 2
            or n.shape != (f.shape[0],) or not np.all(np.isfinite(f))
            or not np.all(np.isfinite(n)) or np.any(n != np.rint(n))
            or len(np.unique(n)) != len(n) or not np.all(np.isfinite([mu, kt, rate]))
            or kt <= 0 or rate < 0):
        raise ValueError('finite sector data, distinct integer counts, positive kT and nonnegative charge rate required')
    mobility = np.asarray(mobilities, float)
    if mobility.shape != (f.shape[0],) or np.any(mobility <= 0) or not np.all(np.isfinite(mobility)):
        raise ValueError('one positive mobility per sector required')
    sectors, cells = f.shape
    phi = f-mu*n[:, None]
    log_local = logsumexp(-phi/kt, axis=0)
    weights = np.exp(-phi/kt-log_local)
    omega = -kt*log_local
    equilibrium = np.exp((-phi/kt).ravel()-logsumexp(-phi/kt))
    slow = sparse.block_diag([sg_chain(row, thermal_energy=kt, spacing=spacing, face_mobility=m)
                              for row, m in zip(phi, mobility)], format='csr')
    projection = sparse.hstack([sparse.eye(cells)]*sectors, format='csr')
    lift = sparse.vstack([sparse.diags(w) for w in weights], format='csr')
    fast = rate*(lift@projection-sparse.eye(cells*sectors, format='csr'))
    # Exact slow generator on the local-equilibrium subspace of THIS grid.
    coarse = (projection@slow@lift).tocsr()
    mean_mobility = mobility@weights
    continuum = sg_chain(omega, thermal_energy=kt, spacing=spacing,
                         face_mobility=.5*(mean_mobility[:-1]+mean_mobility[1:]))
    return dict(generator=slow+fast, slow=slow, charge=fast, projection=projection, lift=lift,
                coarse=coarse, continuum=continuum, weights=weights, omega=omega,
                mean_mobility=mean_mobility, equilibrium=equilibrium,
                coarse_equilibrium=projection@equilibrium,
                status='synthetic charge kinetics; no physical Si clock')


def reversible_propagator(generator, equilibrium, time):
    """Independent dense spectral exponential for modest validation grids.

    Detailed balance symmetrizes L. No density clipping or renormalization.
    A failed stationary/detailed-balance check is an error, not repaired.
    """
    from scipy.linalg import eigh
    l = generator.toarray() if sparse.issparse(generator) else np.asarray(generator, float)
    pi = np.asarray(equilibrium, float)
    time = float(time)
    if (l.ndim != 2 or l.shape[0] != l.shape[1] or pi.shape != (len(l),)
            or not np.all(np.isfinite(l)) or not np.all(np.isfinite(pi)) or np.any(pi <= 0)
            or not np.isfinite(time) or time < 0 or not np.isclose(sum(pi), 1, rtol=0, atol=1e-13)):
        raise ValueError('finite generator and normalized positive equilibrium required')
    from .silicon_initiation_probability import check_generator
    check_generator(l)
    # Uniform changes of time units must not hide probability loss or a false
    # equilibrium. A global/unit-sized absolute tolerance also masks a slow
    # disconnected sector next to a fast one; check local flux scales instead.
    stationary=l@pi
    local_scale=abs(l)@pi
    conductance=l*pi[None,:]
    if (not np.all(np.isfinite(conductance)) or not np.all(np.isfinite(local_scale))
            or np.any(abs(stationary)>1e-11*local_scale)
            or np.any(abs(conductance-conductance.T)>1e-11*(abs(conductance)+abs(conductance.T)))):
        raise ValueError('conservative reversible stationary generator required')
    off=l.copy();np.fill_diagonal(off,0)
    if np.any((off>0)&(conductance==0)):
        raise ValueError('equilibrium conductance numerically underflowed')
    root = np.sqrt(pi)
    h = l*root[None, :]/root[:, None]
    if not np.all(np.isfinite(h)):
        raise ValueError('reversible similarity transform numerically unresolved')
    scale = float(np.max(abs(l)))
    if scale==0:return np.eye(len(l))
    values, vectors = eigh(h)
    if values[-1] > 1e-10*scale:
        raise ValueError('generator has a positive mode')
    # Do not clip even the floating-point zero eigenvalue; expose its residual.
    result=root[:, None]*((vectors*np.exp(values*time))@vectors.T)/root[None, :]
    if (not np.all(np.isfinite(result)) or np.min(result)<-1e-11
            or np.max(abs(result.sum(axis=0)-1))>1e-9):
        raise ValueError('spectral propagator probability balance numerically unresolved')
    return result


def charge_memory_blocks(research_generator):
    """Exact fixed-parameter coordinates (marginal p, hidden charge y).

    For sectors j>0, y_j=P_j-w_j*p. With P=T*p+U*y, the full law is
    p'=A*p+B*y; y'=C*p+D*y. Eliminating y gives the memory kernel
    B exp(D*t) C, plus B exp(D*t)y(0). It is not a new constant mobility.
    This identity assumes time-independent weights/chemical potential.
    """
    r = research_generator
    weights = np.asarray(r['weights'], float)
    if weights.ndim != 2 or len(weights) < 2:
        raise ValueError('at least two charge sectors required for hidden modes')
    sectors, cells = weights.shape
    eye = sparse.eye(cells, format='csr')
    zero = sparse.csr_matrix((cells,cells))
    hidden_lift = sparse.bmat([[-eye]*(sectors-1)]+
        [[eye if i==j else zero for j in range(sectors-1)] for i in range(sectors-1)],format='csr')
    selector = sparse.hstack([sparse.csr_matrix((cells*(sectors-1),cells)),
                              sparse.eye(cells*(sectors-1))],format='csr')
    hidden_projection = selector-sparse.vstack([sparse.diags(w) for w in weights[1:]],format='csr')@r['projection']
    s,t,l = r['projection'],r['lift'],r['generator']
    return dict(A=(s@l@t).tocsr(),B=(s@l@hidden_lift).tocsr(),
                C=(hidden_projection@l@t).tocsr(),D=(hidden_projection@l@hidden_lift).tocsr(),
                hidden_lift=hidden_lift,hidden_projection=hidden_projection)
