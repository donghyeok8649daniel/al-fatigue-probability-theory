"""Track the lower conditional branch found by the thermal search, without fitting."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.linalg import cholesky
from scipy.sparse.linalg import eigsh

from solver_v1.silicon_crack_research import RelaxedCoordinates,SpatialSW,relax_atoms
from .run_local_crack_audit import save_json,polish,atomic_hessian
from .run_static_probe import source_parameters


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed',type=Path,required=True)
    parser.add_argument('--references',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    if any(args.output.iterdir()):
        parser.error('use an empty output')
    start=time.perf_counter()
    meta=json.loads((args.references/'summary.json').read_text(encoding='utf-8'))
    seed=np.load(args.seed,allow_pickle=False)
    reference=np.load(args.references/'initial.npz',allow_pickle=False)
    coords=RelaxedCoordinates(reference['reference_positions'],reference['fixed'],bond=reference['bond'])
    p,_=source_parameters();model=SpatialSW(p,front_period=meta['geometry']['front_period_A'])
    rows=[];positions_by_name={};cache={float(reference['gap_A']):seed['quenched_positions']}
    # The warm start is always the nearest ALREADY FOUND reconstructed branch,
    # so the original crystalline source cannot silently replace this branch.
    points=sorted(meta['rows'],key=lambda row:abs(row['gap_A']-float(reference['gap_A'])))
    for point in points:
        gap=point['gap_A'];nearest=min(cache,key=lambda value:abs(value-gap))
        positions,info=relax_atoms(model,coords,initial=cache[nearest],values=[gap],tolerance=2e-6)
        positions,correction=polish(model,coords,positions,reference['fixed'],values=[gap])
        result=model.evaluate(positions);gradient,reaction=coords.pullback(result.gradient)
        h,_,_=atomic_hessian(model,coords,positions,reference['fixed'])
        factor=cholesky(h.toarray(),lower=True)
        eigenvalue=eigsh(h,k=1,which='SA',return_eigenvectors=False,tol=1e-9,
            v0=np.sin(np.arange(h.shape[0])+.7))[0]
        row=dict(name=point['name'],gap_A=gap,
            energy_from_original_initial_eV=float(np.sum(result.site_energy-reference['sites'])),
            energy_below_original_branch_eV=float(np.sum(result.site_energy-reference['sites'])-point['energy_difference_eV']),
            logdet_bath=float(2*np.log(np.diag(factor)).sum()),
            conditional_lowest_eigenvalue_eV_A2=float(eigenvalue),
            reaction_eV_A=float(reaction[0]),force_residual_eV_A=float(np.max(abs(gradient))),
            max_bath_coordinate_A=float(abs(coords.encode(positions)).max()),
            relaxation=info,correction=correction)
        rows.append(row);positions_by_name[point['name']]=positions;cache[gap]=positions
        save_json(args.output/'running.json',rows)
        print(row['name'],gap,row['energy_below_original_branch_eV'],row['conditional_lowest_eigenvalue_eV_A2'],flush=True)
    # Exact same-q line is only a path upper bound, not a minimum-energy path.
    original=reference['positions'];reconstructed=positions_by_name['initial']
    line=[]
    for fraction in np.linspace(0,1,65):
        r=(1-fraction)*original+fraction*reconstructed
        line.append([float(fraction),float(np.sum(model.evaluate(r).site_energy-reference['sites']))])
    np.savez_compressed(args.output/'positions.npz',**positions_by_name,
        reference=reference['reference_positions'],fixed=reference['fixed'],bond=reference['bond'])
    save_json(args.output/'summary.json',dict(rows=rows,fixed_gap_linear_path=line,
        path_interpretation='sampled straight path, an upper-bound diagnostic; no saddle/MEP certification',
        elapsed_seconds=time.perf_counter()-start,source_seed_sha256=hashlib.sha256(args.seed.read_bytes()).hexdigest(),
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        real_Si_surface_reconstruction_validated=False,finite_T_global_PMF_certified=False))


if __name__=='__main__':
    main()
