"""Quench saved collective-path states into distinct minima and short NEB seeds."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from solver_v1.silicon_crack_research import RelaxedCoordinates,SpatialSW,relax_atoms
from .run_local_crack_audit import save_json,polish,atomic_hessian,spectrum
from .run_static_probe import source_parameters


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--references',type=Path,required=True)
    parser.add_argument('--path',type=Path,required=True)
    parser.add_argument('--resolved',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    if any(args.output.iterdir()):parser.error('use an empty output')
    start=time.perf_counter()
    ref=np.load(args.references/'initial.npz',allow_pickle=False)
    meta=json.loads((args.references/'summary.json').read_text(encoding='utf-8'))
    path=np.load(args.path,allow_pickle=False)
    base=RelaxedCoordinates(ref['reference_positions'],ref['fixed'],bond=ref['bond'])
    p,_=source_parameters();model=SpatialSW(p,front_period=meta['geometry']['front_period_A'])
    gap=float(ref['gap_A']);sites=ref['sites'];states={};rows=[];mapping={}
    candidates=[(name,path[name]) for name in sorted(path.files) if name.startswith(('forward_','reverse_'))
                and int(name.split('_')[1])%4==0]
    for sample in sorted(args.resolved.glob('sign*.npz')):
        data=np.load(sample,allow_pickle=False)
        if 'stationary' in data.files:candidates.append((sample.stem,data['stationary']))
    for label,seed in candidates:
        r,info=relax_atoms(model,base,initial=seed,values=[gap],tolerance=2e-6)
        r,correction=polish(model,base,r,ref['fixed'],values=[gap])
        match=next((name for name,old in states.items() if np.linalg.norm(r-old)<1e-5),None)
        if match is None:
            name=f'state{len(states):02d}';h,_,_=atomic_hessian(model,base,r,ref['fixed'])
            spec,_=spectrum(h,5)
            if spec['negative_index']!=0:raise RuntimeError('quench did not give a stable conditional minimum')
            states[name]=r
            rows.append(dict(name=name,energy_from_original_eV=float(np.sum(model.evaluate(r).site_energy-sites)),
                source_label=label,conditional_spectrum=spec,relaxation=info,correction=correction))
            print('new',rows[-1]['name'],rows[-1]['energy_from_original_eV'],label,flush=True)
            match=name
        mapping[label]=match
    edges=set()
    for direction in ['forward','reverse']:
        labels=sorted((key for key in mapping if key.startswith(direction)),key=lambda name:int(name.split('_')[1]))
        for left,right in zip(labels,labels[1:]):
            a,b=mapping[left],mapping[right]
            if a!=b:edges.add(tuple(sorted((a,b))))
    # Include each nearby resolved minimum against the original state. The
    # resulting band remains a candidate until its index/connection checks pass.
    origin=mapping['forward_000']
    for label,state in mapping.items():
        if label.startswith('sign') and state!=origin:edges.add(tuple(sorted((origin,state))))
    edge_rows=[]
    for a,b in sorted(edges):
        label=a+'_'+b
        np.savez_compressed(args.output/(label+'.npz'),forward_000=states[a],forward_001=states[b])
        line=[]
        for fraction in np.linspace(0,1,33):
            r=(1-fraction)*states[a]+fraction*states[b]
            line.append([float(fraction),float(np.sum(model.evaluate(r).site_energy-sites))])
        edge_rows.append(dict(label=label,initial=a,final=b,linear_path=line,
            endpoint_distance_A=float(np.linalg.norm(states[a]-states[b]))))
    np.savez_compressed(args.output/'minima.npz',**states,fixed=ref['fixed'],bond=ref['bond'],reference=ref['reference_positions'])
    save_json(args.output/'summary.json',dict(states=rows,source_to_state=mapping,edges=edge_rows,
        elapsed_seconds=time.perf_counter()-start,global_state_space_exhausted=False,
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))


if __name__=='__main__':main()
