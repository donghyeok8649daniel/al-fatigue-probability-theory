"""Resolve a continuous fixed-gap reconstruction path and audit its saddle."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from solver_v1.silicon_crack_research import RelaxedCoordinates,SpatialSW,relax_atoms
from solver_v1.silicon_path_research import reparameterize_path,relax_neb
from .run_local_crack_audit import save_json,polish,atomic_hessian,spectrum
from .run_static_probe import source_parameters


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--references',type=Path,required=True)
    parser.add_argument('--path',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--images',type=int,default=25)
    parser.add_argument('--steps',type=int,default=2200)
    parser.add_argument('--spring',type=float,default=1.)
    parser.add_argument('--fire-step',type=float,default=.02)
    parser.add_argument('--fire-max-step',type=float,default=.15)
    parser.add_argument('--climb-after',type=int,default=200)
    parser.add_argument('--max-image-move',type=float,default=.08)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    if any(args.output.iterdir()):
        parser.error('use an empty output')
    start=time.perf_counter()
    reference=np.load(args.references/'initial.npz',allow_pickle=False)
    meta=json.loads((args.references/'summary.json').read_text(encoding='utf-8'))
    path=np.load(args.path,allow_pickle=False)
    coordinates=RelaxedCoordinates(reference['reference_positions'],reference['fixed'],bond=reference['bond'])
    labels=sorted(name for name in path.files if name.startswith('forward_'))
    if not labels:
        raise ValueError('expected saved forward path positions')
    initial=reparameterize_path(np.array([coordinates.encode(path[name]) for name in labels]),args.images)
    p,_=source_parameters()
    models=[SpatialSW(p,front_period=meta['geometry']['front_period_A']) for _ in range(args.images)]
    gap=float(reference['gap_A']);sites=reference['sites']
    def evaluate(index,x):
        result=models[index].evaluate(coordinates.positions(x,[gap]))
        return float(np.sum(result.site_energy-sites)),coordinates.pullback(result.gradient)[0]
    def progress(record,x,energies):
        save_json(args.output/'progress.json',record)
        if record['iteration']%100==0:
            np.savez_compressed(args.output/'checkpoint.npz',bath_coordinates=x,energies=energies)
            print(record,flush=True)
    result=relax_neb(evaluate,initial,spring=args.spring,max_steps=args.steps,progress=progress,
        step=args.fire_step,max_step=args.fire_max_step,climb_after=args.climb_after,
        max_image_move=args.max_image_move)
    images=np.array([coordinates.positions(x,[gap]) for x in result['images']])
    np.savez_compressed(args.output/'band.npz',positions=images,bath_coordinates=result['images'],
        energies=result['energies'],gradients=result['gradients'],initial_bath_coordinates=initial)
    audit={}
    saddle=images[result['highest_image']]
    try:
        saddle,correction=polish(models[0],coordinates,saddle,reference['fixed'],values=[gap])
        h,_,b=atomic_hessian(models[0],coordinates,saddle,reference['fixed'])
        spec,unstable=spectrum(h,8)
        energy=float(np.sum(models[0].evaluate(saddle).site_energy-sites))
        audit=dict(stationary_correction=correction,conditional_spectrum=spec,energy_from_initial_eV=energy,
                   conditional_saddle_index_one=spec['negative_index']==1)
        endpoints=[];endpoint_positions={}
        for sign in [-1.,1.]:
            seed=coordinates.positions(coordinates.encode(saddle)+sign*.03*unstable,[gap])
            relaxed,info=relax_atoms(models[0],coordinates,initial=seed,values=[gap],tolerance=2e-6)
            relaxed,minimum_correction=polish(models[0],coordinates,relaxed,reference['fixed'],values=[gap])
            end_h,_,_=atomic_hessian(models[0],coordinates,relaxed,reference['fixed'])
            end_spec,_=spectrum(end_h,5)
            endpoint_positions[f'descent_{int(sign):+d}']=relaxed
            endpoints.append(dict(sign=sign,energy_from_initial_eV=float(np.sum(models[0].evaluate(relaxed).site_energy-sites)),
                distance_from_original_A=float(np.linalg.norm(relaxed-reference['positions'])),
                distance_from_reconstructed_A=float(np.linalg.norm(relaxed-images[-1])),
                spectrum=end_spec,relaxation=info,correction=minimum_correction))
        audit['descents']=endpoints
        np.savez_compressed(args.output/'stationary.npz',saddle=saddle,unstable_mode=unstable,**endpoint_positions)
    except (RuntimeError,np.linalg.LinAlgError) as error:
        audit['stationary_audit_failure']=str(error)
    save_json(args.output/'summary.json',dict(
        images=args.images,spring=args.spring,iteration_budget=args.steps,history=result['history'],
        fire_step=args.fire_step,fire_max_step=args.fire_max_step,climb_after=args.climb_after,
        max_image_move=args.max_image_move,
        band_converged=result['converged'],iterations=result['iterations'],neb_force_max=result['neb_force_max'],
        maximum_sampled_band_energy_eV=float(result['energies'].max()),highest_image=result['highest_image'],
        saddle_audit=audit,elapsed_seconds=time.perf_counter()-start,
        code_sha256={name:hashlib.sha256(Path(name).read_bytes()).hexdigest() for name in (
            'solver_v1/silicon_path_research.py','results/silicon_wafer_feasibility/run_reconstruction_neb.py')},
        source_path_sha256=hashlib.sha256(args.path.read_bytes()).hexdigest(),
        global_minimum_path_certified=False,physical_transition_rate=None,real_wafer_calibrated=False))


if __name__=='__main__':
    main()
