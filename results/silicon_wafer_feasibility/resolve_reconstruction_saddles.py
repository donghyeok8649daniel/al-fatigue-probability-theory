"""Break a transverse instability of an index-two conditional stationary point."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.sparse.linalg import eigsh

from solver_v1.silicon_crack_research import RelaxedCoordinates,SpatialSW,relax_atoms
from solver_v1.silicon_collective_research import RetainedBathDirection
from .run_local_crack_audit import save_json,polish,atomic_hessian,spectrum
from .run_static_probe import source_parameters


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--references',type=Path,required=True)
    parser.add_argument('--stationary',type=Path,required=True)
    parser.add_argument('--collective-path',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    if any(args.output.iterdir()):parser.error('use an empty output')
    start=time.perf_counter()
    reference=np.load(args.references/'initial.npz',allow_pickle=False)
    meta=json.loads((args.references/'summary.json').read_text(encoding='utf-8'))
    source=np.load(args.stationary,allow_pickle=False)
    path=np.load(args.collective_path,allow_pickle=False)
    base=RelaxedCoordinates(reference['reference_positions'],reference['fixed'],bond=reference['bond'])
    p,_=source_parameters();model=SpatialSW(p,front_period=meta['geometry']['front_period_A'])
    fixed=reference['fixed'];gap=float(reference['gap_A']);sites=reference['sites']
    center=source['stationary'];x0=base.encode(center)
    h,_,_=atomic_hessian(model,base,center,fixed)
    values,vectors=eigsh(h,k=5,which='SA',tol=1e-10,v0=np.sin(np.arange(h.shape[0])+.71))
    negative=vectors[:,values<0]
    if negative.shape[1]!=2:raise ValueError('this audit requires exactly two unstable directions')
    coefficients=negative.T@path['bath_direction']
    coefficients/=np.linalg.norm(coefficients)
    retained=negative@coefficients
    transverse=negative@np.array([coefficients[1],-coefficients[0]])
    coordinates=RetainedBathDirection(base,retained)
    progress0=float(retained@x0)
    # Check each negative direction using independent finite differences of forces.
    derivative_errors=[]
    for mode in negative.T:
        step=1e-5
        plus=base.pullback(model.evaluate(base.positions(x0+step*mode,[gap])).gradient)[0]
        minus=base.pullback(model.evaluate(base.positions(x0-step*mode,[gap])).gradient)[0]
        derivative_errors.append(float(np.max(abs((plus-minus)/(2*step)-h@mode))))
    rows=[]
    for sign in [-1.,1.]:
        for amplitude in [.05,.2]:
            for offset in [0.,-.1,.1]:
                label=f'sign{int(sign):+d}_amp{amplitude:.2f}_offset{offset:+.2f}'
                seed=base.positions(x0+sign*amplitude*transverse+offset*retained,[gap])
                r,info=relax_atoms(model,coordinates,initial=seed,values=[gap,progress0+offset],tolerance=2e-6)
                row=dict(label=label,sign=sign,amplitude_A=amplitude,progress_offset_A=offset,
                    constrained_relaxation=info,constrained_energy_eV=float(np.sum(model.evaluate(r).site_energy-sites)))
                saved={'constrained_positions':r}
                try:
                    stationary,correction=polish(model,base,r,fixed,values=[gap])
                    hstat,_,_=atomic_hessian(model,base,stationary,fixed)
                    spec,unstable=spectrum(hstat,8)
                    energy=float(np.sum(model.evaluate(stationary).site_energy-sites))
                    row.update(energy_from_original_eV=energy,conditional_spectrum=spec,
                        stationary_correction=correction,index=spec['negative_index'])
                    saved.update(stationary=stationary,unstable_mode=unstable)
                    if spec['negative_index']==1:
                        descents=[]
                        for downhill in [-1.,1.]:
                            local_seed=base.positions(base.encode(stationary)+downhill*.03*unstable,[gap])
                            minimum,minimum_info=relax_atoms(model,base,initial=local_seed,values=[gap],tolerance=2e-6)
                            minimum,minimum_correction=polish(model,base,minimum,fixed,values=[gap])
                            hmin,_,_=atomic_hessian(model,base,minimum,fixed)
                            minspec,_=spectrum(hmin,5)
                            saved[f'descent_{int(downhill):+d}']=minimum
                            descents.append(dict(sign=downhill,
                                energy_from_original_eV=float(np.sum(model.evaluate(minimum).site_energy-sites)),
                                distance_from_original_A=float(np.linalg.norm(minimum-reference['positions'])),
                                distance_from_fully_reconstructed_A=float(np.linalg.norm(minimum-path['forward_032'])),
                                conditional_spectrum=minspec,relaxation=minimum_info,correction=minimum_correction))
                        row['descents']=descents
                except (RuntimeError,np.linalg.LinAlgError) as error:
                    row['stationary_failure']=str(error)
                rows.append(row);save_json(args.output/'running.json',rows)
                np.savez_compressed(args.output/(label+'.npz'),**saved)
                print(label,row.get('energy_from_original_eV'),row.get('index'),row.get('stationary_failure'),flush=True)
    save_json(args.output/'summary.json',dict(rows=rows,source_negative_eigenvalues=values[:2].tolist(),
        independent_hessian_vector_errors=derivative_errors,elapsed_seconds=time.perf_counter()-start,
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        source_stationary_sha256=hashlib.sha256(args.stationary.read_bytes()).hexdigest(),
        scope='stationary branches at fixed selected gap; index/connection certification does not fit a physical transition rate',
        physical_transition_rate=None,global_path_certified=False))


if __name__=='__main__':main()
