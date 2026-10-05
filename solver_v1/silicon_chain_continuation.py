"""Exact auxiliary HMC continuation with one-file committed checkpoints.

An interrupted proposal changes neither the retained state nor its RNG state.
These coordinates, momenta and indices have no physical time interpretation.
"""
from __future__ import annotations
from copy import deepcopy
import json
from pathlib import Path
import numpy as np
from .silicon_thermal_research import harmonic_split_proposal


def advance(evaluate, state, *, step, steps):
    """Return a new state; do not mutate the caller's arrays or RNG history."""
    if not np.isfinite(step) or step <= 0 or len(steps) != 2 or not 1 <= steps[0] <= steps[1]:
        raise ValueError('positive finite step and ordered integer trajectory range required')
    rng = np.random.default_rng()
    rng.bit_generator.state = deepcopy(state['rng_state'])
    u = state['u'].copy()
    p = rng.standard_normal(len(u))
    length = int(rng.integers(steps[0], steps[1]+1))
    current = (state['phi'], state['gradient'].copy(), state['observation'].copy())
    h0 = float(.5*(u@u+p@p)+current[0])
    proposal = harmonic_split_proposal(evaluate, u, p, step, length, current)
    uniform = difference = None
    accept = False
    if proposal is not None:
        difference = float(.5*(proposal[0]@proposal[0]+proposal[1]@proposal[1])+proposal[2]-h0)
        if not np.isfinite(difference):
            raise FloatingPointError('nonfinite Hamiltonian difference')
        uniform = float(rng.random())
        accept = bool(np.log(uniform) < min(0., -difference))
        if accept:
            u = proposal[0].copy()
            current = proposal[2:]
    result = dict(u=u, phi=float(current[0]), gradient=current[1].copy(),
                  observation=current[2].copy(), rng_state=deepcopy(rng.bit_generator.state),
                  completed=state['completed']+1, retained_evaluation=state['retained_evaluation'])
    transaction = dict(steps=length, delta_h=difference, uniform=uniform,
                       accepted=accept, domain_rejected=proposal is None)
    return result, transaction


def save_checkpoint(path: Path, state, transaction=None):
    """Commit all retained state and RNG bytes together by atomic replacement."""
    temporary = path.with_suffix(path.suffix+'.tmp')
    with temporary.open('wb') as stream:
        np.savez_compressed(stream, u=state['u'], phi=state['phi'], gradient=state['gradient'],
                            observation=state['observation'], completed=state['completed'],
                            retained_evaluation=state['retained_evaluation'],
                            rng_json=json.dumps(state['rng_state'], allow_nan=False),
                            transaction_json=json.dumps(transaction, allow_nan=False))
    temporary.replace(path)


def load_checkpoint(path: Path):
    with np.load(path, allow_pickle=False) as z:
        state = dict(u=z['u'].copy(), phi=float(z['phi']), gradient=z['gradient'].copy(),
                     observation=z['observation'].copy(), completed=int(z['completed']),
                     retained_evaluation=int(z['retained_evaluation']),
                     rng_state=json.loads(str(z['rng_json'])))
        transaction = json.loads(str(z['transaction_json']))
    if (state['u'].ndim != 1 or state['gradient'].shape != state['u'].shape
            or not np.isfinite(state['u']).all() or not np.isfinite(state['gradient']).all()
            or not np.isfinite(state['phi']) or not np.isfinite(state['observation']).all()
            or state['completed'] < 0 or state['retained_evaluation'] < 0):
        raise ValueError('invalid committed state')
    return state, transaction


def parent_state(folder: Path):
    """Restart from the retained state, including a rejected final proposal."""
    meta = json.loads((folder/'summary.json').read_text(encoding='utf-8'))
    records = json.loads((folder/'transactions.json').read_text(encoding='utf-8'))
    if not meta['complete'] or len(records) != meta['completed_proposals']:
        raise ValueError('complete verified parent required')
    with np.load(folder/'chain.npz', allow_pickle=False) as chain:
        index = int(chain['retained_evaluation'][-1])
        if index != records[-1]['retained_evaluation']:
            raise ValueError('parent retained history mismatch')
        u, observation, phi = chain['samples'][-1].copy(), chain['observations'][-1].copy(), float(chain['potentials'][-1])
    with np.load(folder/f'evaluation_{index:04d}.npz', allow_pickle=False) as record:
        if not bool(record['inside']):
            raise ValueError('parent retained state outside domain')
        for left, right in ((u, record['u']), (observation, record['observation']), (phi, record['correction'])):
            if not np.array_equal(left, right):
                raise ValueError('parent retained bytes mismatch')
        gradient = record['gradient'].copy()
    return dict(u=u, phi=phi, gradient=gradient, observation=observation,
                rng_state=deepcopy(meta['final_rng_state']), completed=0, retained_evaluation=index)
