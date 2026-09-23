"""First passage to specified initiation basins in a deterministic generator.

The basins and physical rates must come from independent material research.
No bond-length threshold, Si rate, initiation law or physical clock is supplied.
Column convention: dP/dt = L P. Existing SG generators obey this convention.
"""
from __future__ import annotations

import numpy as np
import warnings
from scipy import sparse
from scipy.linalg import solve, LinAlgWarning, LinAlgError
from scipy.sparse.linalg import expm_multiply


def check_generator(generator):
    g=sparse.csc_matrix(generator,dtype=float)
    if (g.shape[0]!=g.shape[1] or not g.shape[0] or not np.all(np.isfinite(g.data))):
        raise ValueError('finite nonempty square generator required')
    off=g-sparse.diags(g.diagonal())
    if np.any(off.data<0) or np.any(g.diagonal()>0):
        raise ValueError('generator has negative transition or positive diagonal')
    # A unit-sized/global absolute tolerance can hide missing probability in a
    # slow column, or accept an invalid generator merely by changing time units.
    column_scale=np.asarray(abs(g).sum(axis=0)).ravel()
    column_residual=np.asarray(g.sum(axis=0)).ravel()
    if (not np.all(np.isfinite(column_scale))
            or np.any(abs(column_residual)>1e-12*column_scale)):
        raise ValueError('generator must conserve probability in each column')
    return g


def _indices(values,n):
    raw=np.asarray(values)
    if (raw.ndim!=1 or not len(raw) or not np.issubdtype(raw.dtype,np.integer)
            or np.any(raw<0) or np.any(raw>=n) or len(np.unique(raw))!=len(raw)):
        raise ValueError('distinct in-range integer state indices required')
    return raw.astype(int)


def _all_can_reach(g,targets):
    """Walk predecessor states; g[row destination, column source]."""
    rows=g.tocsr(); reached=set(map(int,targets)); stack=list(reached)
    while stack:
        state=stack.pop(); row=rows.getrow(state)
        for predecessor,rate in zip(row.indices,row.data):
            if rate>0 and int(predecessor) not in reached:
                reached.add(int(predecessor));stack.append(int(predecessor))
    return len(reached)==g.shape[0]


def _reliable_solve(a,b):
    with warnings.catch_warnings():
        warnings.simplefilter('error',LinAlgWarning)
        try:value=solve(a,b,assume_a='gen')
        except (LinAlgWarning,LinAlgError) as exc:
            raise ValueError('first-passage linear system is numerically unresolved') from exc
    if not np.all(np.isfinite(value)):
        raise ValueError('first-passage solve returned nonfinite values')
    return value


def initiation_collectors(generator,basins):
    """Replace each disjoint specified B basin by a first-entry collector.

    basins maps cause names to original state indices. Rates for entering B are
    preserved exactly, while exits from B are removed only in this first-passage
    representation. Overlapping cause labels would double-count and are rejected.
    """
    g=check_generator(generator);n=g.shape[0]
    if not basins or any(not isinstance(k,str) or not k for k in basins):
        raise ValueError('nonempty named initiation basins required')
    names=list(basins);sets=[_indices(basins[k],n) for k in names]
    occupied=np.concatenate(sets)
    if len(np.unique(occupied))!=len(occupied):
        raise ValueError('initiation cause basins must be disjoint')
    transient=np.setdiff1d(np.arange(n),occupied)
    if not len(transient):raise ValueError('at least one transient intact/intermediate state required')
    sub=g[transient][:,transient]
    flux=sparse.csr_matrix(np.vstack([np.asarray(g[b][:,transient].sum(axis=0)).ravel() for b in sets]))
    augmented=sparse.bmat([[sub,None],[flux,sparse.csr_matrix((len(names),len(names)))]],format='csc')
    check_generator(augmented)
    return dict(generator=augmented,transient_generator=sub,flux_matrix=flux,
        transient_indices=transient,cause_names=names,
        all_states_can_reach_B=_all_can_reach(g,occupied),physical_clock=None,
        status='specified-basin first passage only; no material calibration')


def evolve_first_passage(model,initial_transient_mass,times):
    p=np.asarray(initial_transient_mass,float);times=np.asarray(times,float)
    n=len(model['transient_indices']);m=len(model['cause_names'])
    if (p.shape!=(n,) or np.any(p<0) or not np.all(np.isfinite(p))
            or not np.isclose(p.sum(),1,rtol=0,atol=1e-13)
            or times.ndim!=1 or not len(times) or not np.all(np.isfinite(times))
            or times[0]<0 or np.any(np.diff(times)<0)):
        raise ValueError('normalized nonnegative initial intact mass and nondecreasing times required')
    state=np.r_[p,np.zeros(m)];last=0.;states=[];flux=[]
    for t in times:
        if t>last:state=expm_multiply((t-last)*model['generator'],state)
        states.append(state.copy());flux.append(model['flux_matrix']@state[:n]);last=t
    states=np.asarray(states);cumulative=states[:,n:]
    return dict(times=times,transient_mass=states[:,:n],survival=states[:,:n].sum(axis=1),
        cause_cumulative=cumulative,cause_flux=np.asarray(flux),mass_residual=states.sum(axis=1)-1,
        minimum_probability=float(states.min()),physical_clock=None)


def basin_committor(generator,*,intact_core,initiated_basins):
    """Probability to hit B before returning to A; not the first-passage CDF."""
    g=check_generator(generator);n=g.shape[0]
    a,b=_indices(intact_core,n),_indices(initiated_basins,n)
    if np.intersect1d(a,b).size:raise ValueError('A and B must be disjoint')
    if not _all_can_reach(g,np.r_[a,b]):
        raise ValueError('states unable to reach A or B: committor Dirichlet problem undefined')
    t=np.setdiff1d(np.arange(n),np.r_[a,b]);h=np.zeros(n);h[b]=1
    if len(t):
        backward=g.T.toarray();block=backward[np.ix_(t,t)]
        h[t]=_reliable_solve(block,-backward[np.ix_(t,b)].sum(axis=1))
    if np.min(h)<-1e-10 or np.max(h)>1+1e-10:
        raise ValueError('committor solve violates probability bounds')
    residual=np.asarray(g.T@h)[t]
    return dict(committor=h,interior_residual=residual,
        meaning='P(hit specified B before specified A); basins supplied externally')


def first_passage_moments(model):
    """Mean first-entry time and eventual cause probabilities in model units.

    Closed nonabsorbing classes give infinite MFPT; reject instead of adding
    regularization that would manufacture a finite rate.
    """
    if not model['all_states_can_reach_B']:
        raise ValueError('nonabsorbing closed class; finite MFPT is unavailable')
    q=model['transient_generator'].toarray();b=model['flux_matrix'].toarray()
    mean=_reliable_solve(q.T,-np.ones(len(q)))
    causes=_reliable_solve(q.T,-b.T).T
    if (np.any(mean<0) or np.min(causes)<-1e-10 or np.max(causes)>1+1e-10
            or np.max(abs(causes.sum(axis=0)-1))>1e-9):
        raise ValueError('first-passage solve violates time or probability bounds')
    return dict(mean_first_passage_time=mean,eventual_cause_probability=causes,
        mean_residual=q.T@mean+1,cause_residual=causes@q+b,physical_clock=None)
