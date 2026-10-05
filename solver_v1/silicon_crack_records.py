"""Joint first-formation and qualified-persistence records, research only.

Physical state labels are supplied externally. Neither a geometric bond cutoff
nor a persistence duration is inferred here. A duration-dependent qualification
requires age/history in the supplied state space, not future-looking labels.
The original dynamics, including crack closure/healing, remain unchanged.
"""
from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from solver_v1.silicon_initiation_probability import check_generator


def _states(values, size):
    raw = np.asarray(values)
    if (raw.ndim != 1 or not raw.size or not np.issubdtype(raw.dtype, np.integer)
            or np.any(raw < 0) or np.any(raw >= size)
            or len(np.unique(raw)) != len(raw)):
        raise ValueError('distinct nonempty in-range integer states required')
    return raw.astype(int)


def joint_crack_records(generator, *, formation_states, qualified_states):
    """Lift a column generator into three monotone history labels.

    Flags: 0=no formation yet, 1=formed but not yet qualified, 2=qualified.
    Qualified states must be a subset of crack states, including closed cracks
    where appropriate. Entering a state updates flags; physical exits are kept.
    Qualified means confirmed under an EXTERNALLY specified finite protocol,
    never a proof that future healing has probability zero.
    """
    g = check_generator(generator)
    size = g.shape[0]
    formation = _states(formation_states, size)
    qualified = _states(qualified_states, size)
    if np.setdiff1d(qualified, formation).size:
        raise ValueError('qualified states must also be formation states')
    is_formation = np.isin(np.arange(size), formation)
    is_qualified = np.isin(np.arange(size), qualified)
    off = g.tocoo()
    rows, cols, data = [], [], []
    first_flux = np.zeros((2, 3*size))
    for flag in range(3):
        for destination, source, rate in zip(off.row, off.col, off.data):
            if source == destination:
                rows.append(flag*size+source)
                cols.append(flag*size+source)
                data.append(rate)
                continue
            updated = max(flag, 2 if is_qualified[destination]
                          else 1 if is_formation[destination] else 0)
            rows.append(updated*size+destination)
            cols.append(flag*size+source)
            data.append(rate)
            if flag == 0 and updated >= 1:
                first_flux[0, flag*size+source] += rate
            if flag < 2 and updated == 2:
                first_flux[1, flag*size+source] += rate
    lifted = sparse.csc_matrix((data, (rows, cols)), shape=(3*size, 3*size))
    check_generator(lifted)
    return dict(generator=lifted, original_generator=g, state_count=size,
                formation_mask=is_formation, qualified_mask=is_qualified,
                first_event_flux=sparse.csr_matrix(first_flux),
                physical_clock=None, material_labels_verified=False,
                status='joint history records with externally supplied qualification')


def evolve_joint_records(model, initial_mass, times):
    """CDFs of first formation and first qualification; do NOT add them.

    The CDF difference is formed-but-not-yet-qualified history probability,
    not the probability of a permanently recoverable crack. Present crack
    occupancy is reported independently of either cumulative history.
    """
    size = model['state_count']
    p, time = np.asarray(initial_mass, float), np.asarray(times, float)
    if (p.shape != (size,) or not np.isfinite(p).all() or np.any(p < 0)
            or not np.isclose(p.sum(), 1., rtol=0, atol=1e-13)
            or time.ndim != 1 or not len(time) or not np.isfinite(time).all()
            or time[0] < 0 or np.any(np.diff(time) < 0)):
        raise ValueError('normalized finite mass and nondecreasing times required')
    state = np.zeros((3, size))
    for i, value in enumerate(p):
        flag = (2 if model['qualified_mask'][i]
                else 1 if model['formation_mask'][i] else 0)
        state[flag, i] = value
    state = state.ravel()
    saved, fluxes, previous = [], [], 0.
    for t in time:
        if t > previous:
            state = expm_multiply((t-previous)*model['generator'], state)
        saved.append(state.reshape(3, size).copy())
        fluxes.append(model['first_event_flux']@state)
        previous = t
    mass = np.asarray(saved)
    flag_mass = mass.sum(axis=2)
    physical = mass.sum(axis=1)
    return dict(times=time, joint_mass=mass, physical_mass=physical,
                formation_cumulative=flag_mass[:, 1:].sum(axis=1),
                qualified_cumulative=flag_mass[:, 2],
                formed_not_yet_qualified=flag_mass[:, 1],
                present_crack_probability=physical[:, model['formation_mask']].sum(axis=1),
                first_event_flux=np.asarray(fluxes),
                mass_residual=mass.sum(axis=(1, 2))-1.,
                minimum_mass=float(mass.min()), physical_clock=None)
