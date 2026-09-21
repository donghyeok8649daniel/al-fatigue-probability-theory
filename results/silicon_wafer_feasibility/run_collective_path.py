"""Relax other atoms at fixed crack gap and retained reconstruction displacement."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.linalg import cho_factor,cho_solve

from solver_v1.silicon_crack_research import RelaxedCoordinates,SpatialSW,relax_atoms
from solver_v1.silicon_collective_research import RetainedBathDirection
from .run_local_crack_audit import save_json,atomic_hessian
from .run_static_probe import source_parameters


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reconstructed',type=Path,required=True)
    parser.add_argument('--references',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--intervals',type=int,default=32)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    if any(args.output.iterdir()):
        parser.error('choose an empty output')
    start=time.perf_counter()
    reference=np.load(args.references/'initial.npz',allow_pickle=False)
    reconstructed=np.load(args.reconstructed,allow_pickle=False)['initial']
    meta=json.loads((args.references/'summary.json').read_text(encoding='utf-8'))
    base=RelaxedCoordinates(reference['reference_positions'],reference['fixed'],bond=reference['bond'])
    direction=base.encode(reconstructed)
    collective=RetainedBathDirection(base,direction)
    p,_=source_parameters();model=SpatialSW(p,front_period=meta['geometry']['front_period_A'])
    gap=float(reference['gap_A']);rows=[];positions_by_label={}
    for direction_name,indices,initial in [
        ('forward',range(args.intervals+1),reference['positions']),
        ('reverse',range(args.intervals,-1,-1),reconstructed)]:
        previous=initial
        for index in indices:
            fraction=index/args.intervals;retained=fraction*collective.path_length
            positions,info=relax_atoms(model,collective,initial=previous,values=[gap,retained],tolerance=2e-6)
            result=model.evaluate(positions);gradient,reactions=collective.pullback(result.gradient)
            correction=[]
            # Use the same explicit Hessian restriction if L-BFGS stops on
            # roundoff, retaining all nonzero/negative curvature in the record.
            for _ in range(3):
                if np.max(abs(gradient))<2e-6: break
                h,_,_=atomic_hessian(model,base,positions,reference['fixed'])
                restricted=collective.restricted_hessian(h.toarray())
                try:
                    step=cho_solve(cho_factor(restricted,lower=True),-gradient)
                except np.linalg.LinAlgError:
                    correction.append({'nonpositive_conditional_hessian':True});break
                if np.max(abs(step))>.05:
                    correction.append({'large_newton_step_A':float(np.max(abs(step)))});break
                positions=collective.positions(collective.encode(positions)+step,[gap,retained])
                result=model.evaluate(positions);gradient,reactions=collective.pullback(result.gradient)
                correction.append({'remaining_force_eV_A':float(np.max(abs(gradient)))})
            row=dict(direction=direction_name,index=index,fraction=fraction,collective_length_A=retained,
                energy_from_original_eV=float(np.sum(result.site_energy-reference['sites'])),
                q_reaction_eV_A=float(reactions[0]),collective_reaction_eV_A=float(reactions[1]),
                free_force_residual_eV_A=float(np.max(abs(gradient))),
                force_converged=bool(np.max(abs(gradient))<2e-6),relaxation=info,correction=correction,
                maximum_bath_coordinate_A=float(abs(base.encode(positions)).max()))
            label=f'{direction_name}_{index:03d}';positions_by_label[label]=positions
            rows.append(row);previous=positions
            save_json(args.output/'running.json',rows)
            print(direction_name,index,row['energy_from_original_eV'],row['free_force_residual_eV_A'],flush=True)
    pairs={name:{row['index']:row for row in rows if row['direction']==name} for name in ['forward','reverse']}
    difference=max(abs(pairs['forward'][i]['energy_from_original_eV']-pairs['reverse'][i]['energy_from_original_eV'])
                   for i in range(args.intervals+1))
    np.savez_compressed(args.output/'positions.npz',**positions_by_label,bath_direction=collective.direction,
        fixed=reference['fixed'],bond=reference['bond'],reference=reference['reference_positions'])
    save_json(args.output/'summary.json',dict(rows=rows,path_length_A=collective.path_length,
        maximum_forward_reverse_energy_difference_eV=difference,
        maximum_energy_forward_eV=max(row['energy_from_original_eV'] for row in rows if row['direction']=='forward'),
        all_relaxations_force_converged=all(row['force_converged'] for row in rows),
        stationary_saddle_certified=False,global_minimum_path_certified=False,
        interpretation='conditional path at a fixed first gap; collective coordinate in Angstrom and all other atoms relaxed',
        elapsed_seconds=time.perf_counter()-start,
        input_sha256=hashlib.sha256(args.reconstructed.read_bytes()).hexdigest(),
        code_sha256={name:hashlib.sha256(Path(name).read_bytes()).hexdigest() for name in (
            'solver_v1/silicon_collective_research.py','results/silicon_wafer_feasibility/run_collective_path.py')}))


if __name__=='__main__':
    main()
