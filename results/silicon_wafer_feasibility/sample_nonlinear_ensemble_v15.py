"""Bounded exact-energy pCN pilot at 300 K; MC indices are NOT physical time.

Two different proposal references share ONE Cartesian target domain. A short
run cannot establish equilibration, global mixing, a PMF, or crack initiation.
The pure-Si 6 GB allocation guard is preserved before importing torch/MACE.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from solver_v1.silicon_device_ensemble import FixedGripTarget
from solver_v1.silicon_thermal_research import pcn_sample

MODEL_SHA='2f2be696351ac9e94fbe01cdfb6f017679acdbd2db7645209ef55fec9826b012'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,d): p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n',encoding='utf-8')

def main(a):
    if a.output.exists() or a.max_seconds<=0: raise ValueError('fresh output and positive budget required')
    a.root=a.root.resolve(); a.output.mkdir(parents=True)
    if sha(a.model)!=MODEL_SHA: raise ValueError('model hash mismatch')
    sources=[]; paths=[]
    for name in ('loading8','return8'):
        path=a.root/('results/silicon_initiation_v13/dense_'+name+'/raw_hessian.npz')
        with np.load(path) as z: sources.append({k:z[k].copy() for k in z.files})
        paths.append(path)
    if not np.array_equal(sources[0]['basis'],sources[1]['basis']): raise ValueError('same coordinate basis required')
    common=(sources[0]['positions']+sources[1]['positions'])/2
    targets=[FixedGripTarget(s['positions'],s['basis'],s['hessian'],s['gradient'],300.,common) for s in sources]
    if any(not t.inside(s['positions']) for t,s in zip(targets,sources)): raise ValueError('source outside shared domain')
    protocol=dict(started_utc=datetime.now(timezone.utc).isoformat(),
        source_sha256={str(p.relative_to(a.root)).replace('\\','/'):sha(p) for p in paths},
        model_sha256=MODEL_SHA, runner_sha256=sha(Path(__file__)),
        target_module_sha256=sha(a.root/'solver_v1/silicon_device_ensemble.py'),
        sampler_module_sha256=sha(a.root/'solver_v1/silicon_thermal_research.py'),
        seeds=[151,953], temperature_K=300., proposal_scale=.35, warmup=16, draws=48,
        planned_proposals=128, planned_initial_evaluations=2,
        target='exp(-U_MP0b3/kBT) dr_free on the common finite Cartesian domain',
        domain='each free Cartesian component within 3 A of midpoint of two source positions; all pair distances >= 1 A',
        domain_is_not_an_intact_basin=True, domain_sensitivity_not_verified=True,
        fixed_grip_count=int((~targets[0].free).sum()), free_dimension=targets[0].basis.shape[1],
        correction='Phi=(U-U_source)/kBT - dot(u,u)/2; constant coordinate Jacobian',
        accept='min(1,exp(Phi_current-Phi_proposed)); retain repeated states after rejection',
        source_states_share_same_fixed_grip_positions=True,
        gaussian_references_may_differ_but_target_domain_and_energy_are_same=True,
        no_displacement_clipping=True, allocation_guard_available_bytes=6000000000,
        max_seconds=a.max_seconds, new_MD=0,new_DFT=0, physical_clock=None)
    dump(a.output/'protocol.json',protocol)
    np.savez_compressed(a.output/'common_domain.npz',center=common,free=targets[0].free)
    start=time.perf_counter(); calls=0; rows=[]; completed=[]; error=None; phase='preallocation'
    try:
        import psutil
        available=psutil.virtual_memory().available
        dump(a.output/'allocation_check.json',dict(available_bytes=available,required_bytes=6000000000,
             passed=available>=6000000000,checked_utc=datetime.now(timezone.utc).isoformat()))
        if available<6000000000: raise MemoryError('360-atom pure-Si 6 GB allocation guard; no model loaded')
        import torch
        from ase import Atoms
        from mace.calculators import MACECalculator
        from results.silicon_wafer_feasibility.mace_force_only_v9 import force_only
        torch.set_num_threads(2);torch.set_num_interop_threads(1)
        calc=MACECalculator(model_paths=str(a.model),device='cpu',default_dtype='float64')
        phase='sampling'
        for chain,(s,target,seed) in enumerate(zip(sources,targets,protocol['seeds'])):
            folder=a.output/f'chain_{chain}';folder.mkdir()
            atom=Atoms(numbers=s['numbers'],positions=s['positions'],cell=s['cell'],pbc=False)
            count=0
            def evaluate(u):
                nonlocal count,calls
                if time.perf_counter()-start>a.max_seconds: raise TimeoutError('pilot wall-time budget')
                r=target.coordinates(u); inside=target.inside(r)
                if inside:
                    atom.positions[:]=r; e,f=force_only(calc,atom); calls+=1
                    if not np.isfinite(e) or not np.isfinite(f).all(): raise FloatingPointError('nonfinite model output')
                    phi=target.correction(u,e,float(s['energy']))
                    if count==0 and (abs(e-float(s['energy']))>1e-8 or np.max(abs(f-s['forces']))>1e-8):
                        raise ValueError('initial source force/energy replay changed')
                    obs=np.array([e,float(np.sqrt(np.mean(np.sum((r[target.free]-s['positions'][target.free])**2,axis=1)))),
                        float(np.max(abs(r[target.free]-common[target.free])))])
                else: e,f,phi,obs=np.nan,np.full_like(r,np.nan),np.inf,np.full(3,np.nan)
                np.savez_compressed(folder/f'evaluation_{count:03d}.npz',u=u,positions=r,energy=e,forces=f,
                                    correction=phi,inside=inside,observation=obs)
                rows.append(dict(chain=chain,evaluation=count,inside=bool(inside),model_calls_total=calls))
                count+=1;dump(a.output/'progress.json',dict(evaluations=rows,elapsed_seconds=time.perf_counter()-start))
                if count%8==0: print('CHAIN',chain,'EVALUATIONS',count,'CALLS',calls,flush=True)
                return phi,obs
            result=pcn_sample(evaluate,target.reference.whiten(np.zeros(target.basis.shape[1])),
                proposal_scale=.35,draws=48,warmup=16,rng=np.random.default_rng(seed))
            arrays={k:result.pop(k) for k in ('samples','observations','potentials','accepted_at_saved_step')}
            np.savez_compressed(folder/'chain.npz',**arrays)
            dump(folder/'summary.json',result);completed.append(chain)
        phase='complete'
    except Exception as exc: error=type(exc).__name__+': '+str(exc)
    dump(a.output/'summary.json',dict(complete=phase=='complete',phase=phase,error=error,
        completed_chains=completed,evaluations_recorded=len(rows),new_model_calls=calls,
        elapsed_seconds=time.perf_counter()-start,new_MD=0,new_DFT=0,
        equilibrium_certified=False,free_energy_estimated=False,crack_probability_estimated=False,
        physical_clock=None,production_changed=False))
    print((a.output/'summary.json').read_text(),flush=True)
    if phase!='complete': raise SystemExit(2)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path('.'));p.add_argument('--model',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--max-seconds',type=float,default=1200.)
    main(p.parse_args())
