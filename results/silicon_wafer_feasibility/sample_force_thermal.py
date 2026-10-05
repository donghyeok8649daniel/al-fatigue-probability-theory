"""Bounded force-aware sampling pilot; auxiliary HMC coordinates are not time.

Preserve every force evaluation and rejected-state repeat. The shared finite
Cartesian domain is a numerical integration domain, not a crack definition.
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
from solver_v1.silicon_sampling_score import nonlinear_correction_gradient
from solver_v1.silicon_thermal_research import harmonic_split_proposal

MODEL_SHA = '2f2be696351ac9e94fbe01cdfb6f017679acdbd2db7645209ef55fec9826b012'

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def dump(path, value):
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    temporary.replace(path)

def main(a):
    a.root = a.root.resolve()
    if a.output.exists() or a.max_seconds <= 0:
        raise ValueError('fresh output and positive time budget required')
    deadline = datetime.fromisoformat(a.deadline_utc).astimezone(timezone.utc)
    if datetime.now(timezone.utc) >= deadline:
        raise TimeoutError('research deadline already passed')
    if sha(a.model) != MODEL_SHA:
        raise ValueError('model hash mismatch')
    a.output.mkdir(parents=True)
    paths = [a.root/f'results/silicon_initiation_v13/dense_{name}/raw_hessian.npz'
             for name in ('loading8', 'return8')]
    sources = []
    for path in paths:
        with np.load(path) as data:
            sources.append({key:data[key].copy() for key in data.files})
    if not np.array_equal(sources[0]['basis'], sources[1]['basis']):
        raise ValueError('identical free-coordinate basis required')
    common = (sources[0]['positions']+sources[1]['positions'])/2
    direction = (sources[1]['positions']-sources[0]['positions']).ravel()
    direction /= np.linalg.norm(direction)
    targets = [FixedGripTarget(s['positions'],s['basis'],s['hessian'],s['gradient'],300.,common)
               for s in sources]
    source_files = paths+[a.root/'solver_v1'/name for name in (
        'silicon_device_ensemble.py','silicon_sampling_score.py','silicon_thermal_research.py',
        'silicon_thermal_identity.py')]+[a.root/'results/silicon_wafer_feasibility/mace_force_only_v9.py']
    protocol = dict(started_utc=datetime.now(timezone.utc).isoformat(), deadline_utc=deadline.isoformat(),
        runner_sha256=sha(Path(__file__)), source_sha256={str(p.relative_to(a.root)).replace('\\','/'):sha(p) for p in source_files},
        model_sha256=MODEL_SHA, temperature_K=300., seeds=[1151,1953], warmup=24, draws=32,
        step=.3, steps=[2,4], planned_proposals=112, max_seconds=a.max_seconds,
        allocation_guard_available_bytes=6000000000, free_dimension=targets[0].basis.shape[1],
        fixed_grip_count=int((~targets[0].free).sum()), halfwidth_A=3., minimum_pair_A=1.,
        target='exp(-U/kBT) dr_free on ONE common finite Cartesian domain',
        domain_is_not_intact_or_crack_basin=True, domain_sensitivity_verified=False,
        score='-L^-1 B.T F_total / sqrt(kBT) - u',
        hamiltonian='(u.u+p.p)/2+Phi(u); full auxiliary momentum refresh; Metropolis correction',
        step_and_momentum_are_not_physical=True, no_adaptation=True, no_clipping=True,
        new_MD=0,new_DFT=0,physical_clock=None,material_approved=False)
    dump(a.output/'protocol.json',protocol)
    np.savez_compressed(a.output/'common_domain.npz',center=common,free=targets[0].free,direction=direction)
    start=time.perf_counter();calls=0;total_evaluations=0;completed=[];error=None;phase='preallocation'
    try:
        import psutil
        available=psutil.virtual_memory().available
        dump(a.output/'allocation_check.json',dict(available_bytes=available,required_bytes=6000000000,
             passed=available>=6000000000,checked_utc=datetime.now(timezone.utc).isoformat()))
        if available<6000000000:
            raise MemoryError('pure-Si 6 GB allocation guard; model not loaded')
        import torch
        from ase import Atoms
        from mace.calculators import MACECalculator
        from results.silicon_wafer_feasibility.mace_force_only_v9 import force_only
        torch.set_num_threads(2);torch.set_num_interop_threads(1)
        calc=MACECalculator(model_paths=str(a.model),device='cpu',default_dtype='float64')
        phase='preflight_and_sampling'
        for chain,(source,target,seed) in enumerate(zip(sources,targets,protocol['seeds'])):
            folder=a.output/f'chain_{chain}';folder.mkdir()
            atom=Atoms(numbers=source['numbers'],positions=source['positions'],cell=source['cell'],pbc=False)
            count=0;last_index=None;context='preflight';iteration=-1
            def evaluate(u):
                nonlocal calls,total_evaluations,count,last_index
                if time.perf_counter()-start>a.max_seconds or datetime.now(timezone.utc)>=deadline:
                    raise TimeoutError('absolute or elapsed research budget reached')
                r=target.coordinates(u);inside=target.inside(r)
                if inside:
                    atom.positions[:]=r;e,f=force_only(calc,atom);calls+=1
                    phi=target.correction(u,e,float(source['energy']))
                    gradient=nonlinear_correction_gradient(target,u,f)
                    obs=np.array([e,float((r-common).ravel()@direction),
                        float(np.sqrt(np.mean(np.sum((r[target.free]-common[target.free])**2,axis=1))))])
                    if not np.isfinite([e,phi]).all() or not np.isfinite(f).all():
                        raise FloatingPointError('nonfinite actual energy or force')
                else:
                    e=np.nan;f=np.full_like(r,np.nan);phi=np.inf
                    gradient=np.full_like(u,np.nan);obs=np.full(3,np.nan)
                last_index=count
                np.savez_compressed(folder/f'evaluation_{count:04d}.npz',u=u,positions=r,energy=e,forces=f,
                    correction=phi,gradient=gradient,inside=inside,observation=obs,
                    context=context,iteration=iteration)
                count+=1;total_evaluations+=1
                dump(a.output/'progress.json',dict(chain=chain,phase=context,iteration=iteration,
                    chain_evaluations=count,evaluations_recorded=total_evaluations,new_model_calls=calls,
                    elapsed_seconds=time.perf_counter()-start,updated_utc=datetime.now(timezone.utc).isoformat()))
                if count%12==0:
                    print('CHAIN',chain,'PHASE',context,'ITERATION',iteration,'CALLS',calls,flush=True)
                return phi,gradient,obs
            initial=target.reference.whiten(np.zeros(target.basis.shape[1]))
            current=evaluate(initial)
            with np.load(folder/'evaluation_0000.npz') as record:
                energy_error=abs(float(record['energy'])-float(source['energy']))
                force_error=float(np.max(abs(record['forces']-source['forces'])))
            if energy_error>1e-8 or force_error>1e-8:
                raise ValueError('source energy/force replay mismatch')
            check_rng=np.random.default_rng(711+chain)
            point=initial+.08*check_rng.standard_normal(len(initial))
            vector=check_rng.standard_normal(len(initial));vector/=np.linalg.norm(vector)
            center=evaluate(point);epsilon=1e-4
            plus=evaluate(point+epsilon*vector);minus=evaluate(point-epsilon*vector)
            finite=(plus[0]-minus[0])/(2*epsilon);analytic=float(center[1]@vector)
            preflight=dict(source_energy_error_eV=energy_error,source_force_error_eV_A=force_error,
                directional_difference_epsilon=epsilon,finite_difference=finite,analytic=analytic,
                absolute_score_error=abs(finite-analytic),score_tolerance=2e-4)
            if not np.isfinite(finite) or abs(finite-analytic)>2e-4:
                dump(folder/'preflight.json',preflight)
                raise ValueError('actual force/energy directional derivative mismatch')
            if chain==0:
                momentum=check_rng.standard_normal(len(initial))
                proposal=harmonic_split_proposal(evaluate,point,momentum,.3,2,center)
                if proposal is None:raise ValueError('preflight path left domain')
                reverse=harmonic_split_proposal(evaluate,proposal[0],-proposal[1],.3,2,proposal[2:])
                if reverse is None:raise ValueError('reverse preflight left domain')
                preflight['actual_reverse_position_error']=float(np.max(abs(reverse[0]-point)))
                preflight['actual_reverse_momentum_error']=float(np.max(abs(reverse[1]+momentum)))
                if max(preflight['actual_reverse_position_error'],preflight['actual_reverse_momentum_error'])>1e-8:
                    raise ValueError('actual-force reverse path mismatch')
            dump(folder/'preflight.json',preflight)
            # Reset to the declared source; preflight calls are not chain draws.
            u=initial.copy();current_index=0;rng=np.random.default_rng(seed)
            samples=[];observations=[];potentials=[];accepted=[];retained=[];errors=[];records=[]
            context='sampling';accepted_warmup=0;outside=0
            for iteration in range(protocol['warmup']+protocol['draws']):
                momentum=rng.standard_normal(len(u))
                length=int(rng.integers(2,5));start_index=current_index;first=count
                h0=float(.5*(u@u+momentum@momentum)+current[0])
                proposal=harmonic_split_proposal(evaluate,u,momentum,.3,length,current)
                if proposal is None:
                    accept=False;difference=None;uniform=None;outside+=1
                else:
                    difference=float(.5*(proposal[0]@proposal[0]+proposal[1]@proposal[1])+proposal[2]-h0)
                    uniform=float(rng.random());accept=bool(np.log(uniform)<min(0.,-difference))
                    if accept:
                        u=proposal[0].copy();current=proposal[2:];current_index=last_index
                records.append(dict(iteration=iteration,steps=length,start_evaluation=start_index,
                    first_evaluation=first,last_evaluation=last_index,retained_evaluation=current_index,
                    delta_h=difference,uniform=uniform,accepted=accept,domain_rejected=proposal is None))
                if iteration<protocol['warmup']:
                    accepted_warmup+=int(accept)
                else:
                    samples.append(u.copy());observations.append(current[2].copy());potentials.append(current[0])
                    accepted.append(accept);retained.append(current_index);errors.append(np.inf if difference is None else difference)
                np.savez_compressed(folder/'chain.npz',samples=np.asarray(samples),observations=np.asarray(observations),
                    potentials=np.asarray(potentials),accepted_at_saved_step=np.asarray(accepted),
                    retained_evaluation=np.asarray(retained),hamiltonian_errors=np.asarray(errors))
                dump(folder/'transactions.json',records)
                dump(folder/'summary.json',dict(complete=iteration==55,completed_proposals=iteration+1,
                    acceptance_fraction=float(np.mean(accepted)) if accepted else None,
                    warmup_acceptance_fraction=accepted_warmup/min(iteration+1,24),
                    domain_rejected_trajectories=outside,final_rng_state=rng.bit_generator.state,physical_time_ps=None))
            completed.append(chain)
        phase='complete'
    except Exception as exc:
        error=type(exc).__name__+': '+str(exc)
    dump(a.output/'summary.json',dict(complete=phase=='complete',phase=phase,error=error,completed_chains=completed,
        evaluations_recorded=total_evaluations,new_model_calls=calls,elapsed_seconds=time.perf_counter()-start,
        finished_utc=datetime.now(timezone.utc).isoformat(),new_MD=0,new_DFT=0,equilibrium_certified=False,
        free_energy_estimated=False,crack_probability_estimated=False,physical_clock=None,production_changed=False))
    print((a.output/'summary.json').read_text(),flush=True)
    if phase!='complete':raise SystemExit(2)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True);parser.add_argument('--model',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--max-seconds',type=float,default=2350.)
    parser.add_argument('--deadline-utc',required=True)
    main(parser.parse_args())
