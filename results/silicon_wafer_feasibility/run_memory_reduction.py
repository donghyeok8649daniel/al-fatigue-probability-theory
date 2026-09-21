"""Exact-Hessian spectral compression of memory, with withheld time windows."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from solver_v1.silicon_conditional_research import harmonic_release
from solver_v1.silicon_memory_reduction import (
    positive_memory_quadrature,quadrature_memory,quadrature_release,
)
from .run_local_crack_audit import save_json


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--references',type=Path,default=Path('.cache/si-thermal-v4/front4'))
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    if any(args.output.iterdir()):
        parser.error('use an empty output')
    start=time.perf_counter()
    v3=Path('results/silicon_conditional_v3')
    spectrum=np.load(v3/'memory_spectra.npz',allow_pickle=False)
    stationary=json.loads((v3/'stationary_gaussian.json').read_text(encoding='utf-8'))
    times=np.arange(20001)*.001
    records=[]
    for name in ['initial','saddle','opened_minimum']:
        values=spectrum[name+'_lambda_eV_A2'];weights=spectrum[name+'_weights_eV_A2']
        exact=quadrature_memory(values,weights,times)
        histories=dict(time_ps=times,exact_memory_eV_A2=exact)
        for rank in [1,2,4,8,16,32,64,128,256]:
            reduced=positive_memory_quadrature(values,weights,rank)
            memory=quadrature_memory(reduced['nodes_eV_A2'],reduced['weights_eV_A2'],times)
            histories[f'memory_rank{rank}']=memory
            np.savez_compressed(args.output/f'{name}_rank{rank}.npz',**reduced)
            errors={str(window):dict(max_relative_to_Gamma0=float(np.max(abs(memory[times<=window]-exact[times<=window]))/weights.sum()),
                rms_relative_to_Gamma0=float(np.sqrt(np.mean((memory[times<=window]-exact[times<=window])**2))/weights.sum()))
                for window in [.1,.5,1.,5.,20.]}
            records.append(dict(state=name,rank=rank,kernel_errors=errors,orthogonality_error=reduced['orthogonality_error']))
        np.savez_compressed(args.output/(name+'_histories.npz'),**histories)
        print(name,'kernel ranks complete',flush=True)
    # Independent Cartesian response, not the spectral quadrature itself.
    reference=np.load(args.references/'initial.npz',allow_pickle=False)
    free=np.flatnonzero(~reference['fixed']);bond=reference['bond']
    from solver_v1.silicon_crack_research import RelaxedCoordinates
    coordinates=RelaxedCoordinates(reference['reference_positions'],reference['fixed'],bond=bond)
    dofs=(free[:,None]*3+np.arange(3)).ravel()
    b=coordinates.tangent_matrix()[dofs]
    observation=np.zeros(len(dofs))
    for atom,sign in zip(bond,[-1.,1.]): observation[3*np.searchsorted(free,atom)+1]=sign
    d=.5*observation
    displacement=d-b@reference['response']
    expected=harmonic_release(reference['full_hessian'],displacement,observation,times)
    responses=dict(time_ps=times,cartesian_response=expected)
    for row in [r for r in records if r['state']=='initial']:
        rank=row['rank'];reduced=np.load(args.output/f'initial_rank{rank}.npz',allow_pickle=False)
        response=quadrature_release(reduced['nodes_eV_A2'],reduced['weights_eV_A2'],
            stationary['initial']['relaxed_curvature_a_eV_A2'],stationary['initial']['q_mass_eV_ps2_A2'],times)
        responses[f'gap_rank{rank}']=response
        row['response_max_errors']={str(window):float(np.max(abs(response[times<=window]-expected[times<=window])))
                                     for window in [.1,.5,1.,5.,20.]}
    np.savez_compressed(args.output/'release.npz',**responses)
    save_json(args.output/'summary.json',dict(records=records,
        elapsed_seconds=time.perf_counter()-start,coefficient_fitting=False,
        material_parameters_added=0,constant_friction_added=False,
        exact_scope='positive spectral quadrature of the unchanged finite harmonic bath',
        production_PDE_enabled=False,production_t0_seconds=None,
        code_sha256={name:hashlib.sha256(Path(name).read_bytes()).hexdigest() for name in (
            'solver_v1/silicon_memory_reduction.py','results/silicon_wafer_feasibility/run_memory_reduction.py')},
        source_spectrum_sha256=hashlib.sha256((v3/'memory_spectra.npz').read_bytes()).hexdigest()))
    for row in [r for r in records if r['state']=='initial']:
        print(row['rank'],row['response_max_errors'],flush=True)


if __name__=='__main__':
    main()
