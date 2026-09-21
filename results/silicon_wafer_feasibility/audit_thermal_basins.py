"""Identify conditional stationary branches reached by actual thermal samplers."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from solver_v1.silicon_crack_research import RelaxedCoordinates,SpatialSW,relax_atoms
from .run_local_crack_audit import save_json,polish
from .run_static_probe import source_parameters


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ensembles',type=Path,required=True)
    parser.add_argument('--references',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--labels',nargs='+',required=True)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    if any(args.output.iterdir()):
        parser.error('choose an empty output')
    start=time.perf_counter();rows=[]
    meta=json.loads((args.references/'summary.json').read_text(encoding='utf-8'))
    p,_=source_parameters()
    model=SpatialSW(p,front_period=meta['geometry']['front_period_A'])
    for label in args.labels:
        record=json.loads((args.ensembles/(label+'.json')).read_text(encoding='utf-8'))
        data=np.load(args.ensembles/(label+'.npz'),allow_pickle=False)
        ref=np.load(args.references/(record['reference']+'.npz'),allow_pickle=False)
        coords=RelaxedCoordinates(ref['reference_positions'],ref['fixed'],bond=ref['bond'])
        source=data['final_positions_A']
        positions,info=relax_atoms(model,coords,initial=source,values=[record['gap_A']],tolerance=2e-6,maxiter=3000)
        correction=None
        if not info['force_converged']:
            try:
                positions,correction=polish(model,coords,positions,ref['fixed'],values=[record['gap_A']])
            except RuntimeError as error:
                correction={'failed':str(error)}
        result=model.evaluate(positions)
        gradient,reaction=coords.pullback(result.gradient)
        change=positions-ref['positions']
        atom_norm=np.linalg.norm(change,axis=1)
        largest=np.argsort(atom_norm)[-12:][::-1]
        bath=coords.encode(positions)
        row=dict(label=label,temperature_K=record['temperature_K'],seed=record['seed'],
            reference=record['reference'],conditional_gap_A=record['gap_A'],relaxation=info,polish=correction,
            energy_difference_from_original_branch_eV=float(np.sum(result.site_energy-ref['sites'])),
            remaining_bath_force_eV_A=float(np.max(abs(gradient))),reaction_eV_A=float(reaction[0]),
            rms_atom_displacement_A=float(np.sqrt(np.mean(atom_norm[~ref['fixed']]**2))),
            maximum_atom_displacement_A=float(atom_norm.max()),
            max_bath_coordinate_A=float(abs(bath).max()),
            inside_sampling_box=bool(np.all(abs(bath)<record['bath_box_halfwidth_A'])),
            largest_displaced_atom_ids=largest.tolist(),largest_displacements_A=atom_norm[largest].tolist(),
            interpretation='unbounded fixed-gap quench identifies minima; it does not change or certify the finite-box thermal measure',
            input_sha256=hashlib.sha256((args.ensembles/(label+'.npz')).read_bytes()).hexdigest())
        np.savez_compressed(args.output/(label+'.npz'),thermal_positions=source,quenched_positions=positions,
            original_branch=ref['positions'],fixed=ref['fixed'],bond=ref['bond'])
        rows.append(row);save_json(args.output/'running.json',rows)
        print(label,row['energy_difference_from_original_branch_eV'],row['maximum_atom_displacement_A'],
              row['remaining_bath_force_eV_A'],flush=True)
    save_json(args.output/'summary.json',dict(rows=rows,elapsed_seconds=time.perf_counter()-start,
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),global_PMF_certified=False))


if __name__=='__main__':
    main()
