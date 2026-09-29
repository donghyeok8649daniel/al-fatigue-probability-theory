"""Actual MP-0b3 evaluations of full-dimensional Gaussian proposals, not MD.

Antithetic pairs and repeated temperatures are dependent diagnostic probes.
Weights diagnose overlap; they do not certify a partition function or a basin.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.linalg import eigh, solve
from scipy.spatial.distance import pdist
from scipy.special import logsumexp
from solver_v1.silicon_thermal_validity import KB, validate_fixed_cartesian
from results.silicon_wafer_feasibility.mace_force_only_v9 import force_only
from results.silicon_wafer_feasibility.run_mace_boron_v9 import MODEL_SHA256


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,d): p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def main(a):
    if a.output.exists(): raise ValueError('fresh output required')
    if a.max_seconds<=0: raise ValueError('positive budget required')
    if sha(a.model)!=MODEL_SHA256: raise ValueError('model hash mismatch')
    geometry=a.root/'results/silicon_initiation_v11/intact_prism_360_linesearch_v2/geometry.npz'
    free=np.load(geometry)['free'];sources=[]
    for name in ('loading8','return8'):
        path=a.root/('results/silicon_initiation_v13/dense_'+name+'/raw_hessian.npz')
        with np.load(path) as data: s={key:data[key].copy() for key in data.files}
        s['hessian'],_=validate_fixed_cartesian(s['hessian'],s['basis'],free)
        lam,v=eigh(s['hessian'])
        if lam[0]<=0: raise ValueError('positive Gaussian modes required')
        s.update(name=name,eigenvalues=lam,eigenvectors=v,sha256=sha(path),
                 center=-solve(s['hessian'],s['gradient'],assume_a='pos'))
        sources.append(s)
    a.output.mkdir(parents=True)
    protocol=dict(model_sha256=MODEL_SHA256,source_sha256={s['name']:s['sha256'] for s in sources},
                  geometry_sha256=sha(geometry),runner_sha256=sha(Path(__file__)),
                  temperatures_K=[50.,300.],seeds=[271828,314159],normal_vectors_per_seed=4,
                  signs=[-1,1],planned_points=64,planned_source_replays=2,
                  proposal='d=-H^-1 g + V sqrt(kBT/lambda) z, z~N(0,I648)',
                  proposed_points_are_not_canonical_MACE_samples=True,
                  antithetic_pairs_and_temperature_pairs_are_not_independent=True,
                  no_displacement_clipping=True,collision_guard_minimum_pair_A=1.0,
                  guard_is_computational_not_a_crack_criterion=True,
                  allocation_guard_available_bytes=6000000000,
                  allocation_scope='360-atom pure Si; not the failed 321-1200-atom oxygen batch',
                  max_seconds=a.max_seconds,new_MD=0,new_DFT=0)
    dump(a.output/'protocol.json',protocol)
    started=time.perf_counter();rows=[];baselines=[];calls=0;status='partial';error=None
    try:
        import psutil
        import torch
        from ase import Atoms
        from mace.calculators import MACECalculator
        if psutil.virtual_memory().available<6000000000: raise MemoryError('360-atom allocation guard')
        torch.set_num_threads(2);torch.set_num_interop_threads(1)
        calc=MACECalculator(model_paths=str(a.model),device='cpu',default_dtype='float64')
        for s in sources:
            atom=Atoms(numbers=s['numbers'],positions=s['positions'],cell=s['cell'],pbc=False)
            e0,f0=force_only(calc,atom);calls+=1
            if abs(e0-float(s['energy']))>1e-8 or np.max(abs(f0-s['forces']))>1e-8:
                raise ValueError('source replay changed')
            baselines.append(dict(state=s['name'],energy_error_eV=e0-float(s['energy']),force_error_eV_A=float(np.max(abs(f0-s['forces'])))))
            for seed in protocol['seeds']:
                normals=np.random.default_rng(seed).standard_normal((4,len(s['eigenvalues'])))
                for sample,z in enumerate(normals):
                    for temperature in protocol['temperatures_K']:
                        fluct=s['eigenvectors']@(np.sqrt(KB*temperature/s['eigenvalues'])*z)
                        for sign in protocol['signs']:
                            if time.perf_counter()-started>a.max_seconds: raise TimeoutError('bounded full-dimensional probe budget')
                            d=s['center']+sign*fluct
                            r=s['positions']+(s['basis']@d).reshape(-1,3)
                            minpair=float(pdist(r).min())
                            if minpair<1.0: raise ValueError('proposal collision guard; no clipping/substitution')
                            if not np.array_equal(r[~free],s['positions'][~free]): raise ValueError('fixed grip moved')
                            atom.positions[:]=r;e,f=force_only(calc,atom);calls+=1
                            if not np.isfinite(e) or not np.isfinite(f).all(): raise ValueError('nonfinite actual model evaluation')
                            predicted=float(s['gradient']@d+.5*d@s['hessian']@d)
                            err=float(e-e0-predicted)
                            g=-s['basis'].T@f.ravel();linear=s['gradient']+s['hessian']@d
                            row=dict(state=s['name'],seed=seed,sample=sample,temperature_K=temperature,sign=sign,
                                     actual_delta_energy_eV=float(e-e0),quadratic_delta_energy_eV=predicted,
                                     energy_residual_eV=err,energy_residual_over_kBT=err/(KB*temperature),
                                     log_importance_weight=-err/(KB*temperature),minimum_pair_A=minpair,
                                     free_atom_RMS_displacement_A=float(np.sqrt(np.mean(np.sum((r[free]-s['positions'][free])**2,axis=1)))),
                                     max_atom_displacement_A=float(np.linalg.norm(r-s['positions'],axis=1).max()),
                                     full_gradient_residual_norm_eV_A=float(np.linalg.norm(g-linear)),
                                     linear_gradient_norm_eV_A=float(np.linalg.norm(linear)))
                            index=len(rows);rows.append(row)
                            np.savez_compressed(a.output/f'point_{index:03d}.npz',positions=r,forces=f,energy=e,
                                                reduced_displacement=d,normal_vector=z)
                            dump(a.output/'points.json',rows)
                            if len(rows)%8==0: print('POINTS',len(rows),'/64',round(time.perf_counter()-started,1),'s',flush=True)
        status='complete'
    except Exception as exc:
        error=type(exc).__name__+': '+str(exc)
    groups=[]
    for name in ('loading8','return8'):
        for temperature in protocol['temperatures_K']:
            r=[v for v in rows if v['state']==name and v['temperature_K']==temperature]
            if not r: continue
            logs=np.array([v['log_importance_weight'] for v in r]);w=np.exp(logs-logsumexp(logs))
            groups.append(dict(state=name,temperature_K=temperature,points=len(r),
                               deltaU_over_kBT_min=float(-logs.max()),deltaU_over_kBT_max=float(-logs.min()),
                               weight_ESS_diagnostic=float(1/np.dot(w,w)),maximum_normalized_weight=float(w.max()),
                               max_gradient_residual_over_linear=max(v['full_gradient_residual_norm_eV_A']/v['linear_gradient_norm_eV_A'] for v in r),
                               independent_draws=len(r)//2,statistical_ESS_certified=False))
    result=dict(complete=status=='complete',error=error,points=len(rows),new_model_calls=calls,
                source_replays=baselines,groups=groups,elapsed_seconds=time.perf_counter()-started,
                new_MD=0,new_DFT=0,physical_clock=None,material_approved=False,
                canonical_sampling=False,basin_membership_certified=False,free_energy_estimated=False,
                interpretation='small full-dimensional Gaussian-proposal overlap diagnostic with actual model evaluations')
    dump(a.output/'summary.json',result)
    dump(a.output/'artifact_manifest.json',{'sha256':{p.name:sha(p) for p in a.output.iterdir() if p.is_file()}})
    print(json.dumps(result,indent=2),flush=True)
    if status!='complete': raise SystemExit(2)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path('.'))
    p.add_argument('--model',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--max-seconds',type=float,default=300.)
    main(p.parse_args())
