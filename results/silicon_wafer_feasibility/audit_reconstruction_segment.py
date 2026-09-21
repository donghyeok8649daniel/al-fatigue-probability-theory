"""Two-direction constrained continuation between specified conditional minima.

Sampled maxima are never barriers by themselves. Candidate stationary points
are polished in the full fixed-q bath and audited for index and connectivity.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from solver_v1.silicon_crack_research import RelaxedCoordinates,SpatialSW,relax_atoms
from solver_v1.silicon_collective_research import RetainedBathDirection
from .run_local_crack_audit import save_json,polish,atomic_hessian,spectrum
from .run_static_probe import source_parameters


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--references',type=Path,required=True)
    parser.add_argument('--minima',type=Path,required=True)
    parser.add_argument('--left',default='state04')
    parser.add_argument('--right',default='state01')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    if any(args.output.iterdir()):parser.error('use an empty output')
    ref=np.load(args.references/'initial.npz',allow_pickle=False)
    states=np.load(args.minima,allow_pickle=False)
    metadata=json.loads((args.references/'summary.json').read_text(encoding='utf-8'))
    p,_=source_parameters();model=SpatialSW(p,front_period=metadata['geometry']['front_period_A'])
    base=RelaxedCoordinates(ref['reference_positions'],ref['fixed'],bond=ref['bond'])
    left,right=states[args.left],states[args.right]
    coordinates=RetainedBathDirection(base,base.encode(right)-base.encode(left))
    first,last=[float(coordinates.direction@base.encode(r)) for r in [left,right]]
    gap=float(ref['gap_A']);points={};rows=[]
    for direction,indices,start in [('forward',range(49),left),('reverse',range(48,-1,-1),right)]:
        previous=start
        for i in indices:
            value=first+(last-first)*i/48
            r,info=relax_atoms(model,coordinates,initial=previous,values=[gap,value],tolerance=1e-6)
            if not info['force_converged']:raise RuntimeError('constrained continuation did not converge')
            gradient,reaction=coordinates.pullback(model.evaluate(r).gradient)
            label=f'{direction}_{i:03d}';points[label]=r;previous=r
            rows.append(dict(label=label,direction=direction,index=i,fraction=i/48,
                energy_from_original_eV=float(np.sum(model.evaluate(r).site_energy-ref['sites'])),
                reaction_eV_A=float(reaction[1]),force_residual_eV_A=float(np.max(abs(gradient))),
                relaxation=info))
        save_json(args.output/'running.json',rows)
        print(direction,'complete',flush=True)
    audits=[]
    for direction in ['forward','reverse']:
        curve=sorted([r for r in rows if r['direction']==direction],key=lambda r:r['index'])
        for i in range(1,len(curve)-1):
            row=curve[i]
            if row['energy_from_original_eV']<max(curve[i-1]['energy_from_original_eV'],curve[i+1]['energy_from_original_eV']):
                continue
            label=row['label'];audit=dict(source=label)
            try:
                r,correction=polish(model,base,points[label],ref['fixed'],values=[gap])
                h,_,_=atomic_hessian(model,base,r,ref['fixed']);spec,mode=spectrum(h,8)
                audit.update(stationary_correction=correction,conditional_spectrum=spec,
                    energy_from_original_eV=float(np.sum(model.evaluate(r).site_energy-ref['sites'])))
                saved=dict(stationary=r,unstable_mode=mode)
                if spec['negative_index']==1:
                    ends=[]
                    for sign in [-1.,1.]:
                        seed=base.positions(base.encode(r)+sign*.03*mode,[gap])
                        minimum,info=relax_atoms(model,base,initial=seed,values=[gap],tolerance=1e-6)
                        minimum,fix=polish(model,base,minimum,ref['fixed'],values=[gap])
                        hm,_,_=atomic_hessian(model,base,minimum,ref['fixed']);ms,_=spectrum(hm,5)
                        ends.append(dict(sign=sign,energy_from_original_eV=float(np.sum(model.evaluate(minimum).site_energy-ref['sites'])),
                            distances_A={name:float(np.linalg.norm(minimum-states[name])) for name in states.files if name.startswith('state')},
                            conditional_spectrum=ms,relaxation=info,correction=fix))
                        saved[f'descent_{int(sign):+d}']=minimum
                    audit['descents']=ends
                np.savez_compressed(args.output/(label+'_stationary.npz'),**saved)
            except (RuntimeError,np.linalg.LinAlgError) as error:
                audit['failure']=str(error)
            audits.append(audit);print(audit,flush=True)
    np.savez_compressed(args.output/'positions.npz',**points,direction=coordinates.direction)
    save_json(args.output/'summary.json',dict(rows=rows,stationary_audits=audits,left=args.left,right=args.right,
        collective_length_A=coordinates.path_length,
        maximum_directional_disagreement_eV=max(abs(next(r['energy_from_original_eV'] for r in rows if r['label']==f'forward_{i:03d}')-
                                                   next(r['energy_from_original_eV'] for r in rows if r['label']==f'reverse_{i:03d}')) for i in range(49)),
        source_minima_sha256=hashlib.sha256(args.minima.read_bytes()).hexdigest(),
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        global_path_certified=False,physical_transition_rate=None))


if __name__=='__main__':main()
