"""New bounded continuation of the two verified fixed-grip Si HMC chains."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib, json, time
from pathlib import Path
import numpy as np
from solver_v1.silicon_device_ensemble import FixedGripTarget
from solver_v1.silicon_sampling_score import nonlinear_correction_gradient
from solver_v1.silicon_chain_continuation import advance, parent_state, save_checkpoint

PARENT_MANIFEST_SHA='d3209e42a7022ae5649beaa675b3250e5e0ce441477abbeba725e96db96de4a1'
MODEL_SHA='2f2be696351ac9e94fbe01cdfb6f017679acdbd2db7645209ef55fec9826b012'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def arrays(path):
    with np.load(path,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
def dump(path,value):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    temporary.replace(path)
def verify_parent(root,parent):
    manifest=parent/'package_manifest.json'
    if sha(manifest)!=PARENT_MANIFEST_SHA:raise ValueError('parent manifest fingerprint changed')
    content=json.loads(manifest.read_text(encoding='utf-8'))
    for relative,digest in content['files'].items():
        if sha(root/relative)!=digest:raise ValueError('parent bytes changed: '+relative)
    pilot=parent/'atomistic_pilot'
    p=json.loads((pilot/'protocol.json').read_text(encoding='utf-8'))
    if not json.loads((pilot/'summary.json').read_text(encoding='utf-8'))['complete']:
        raise ValueError('parent run incomplete')
    for relative,digest in p['source_sha256'].items():
        if sha(root/relative)!=digest:raise ValueError('parent source changed: '+relative)
    return p

def main(a):
    a.root=a.root.resolve()
    if a.output.exists() or a.draws<=0 or a.max_seconds<=0:
        raise ValueError('fresh output and positive bounded plan required')
    parent=a.root/'results/silicon_force_thermal_comparison'
    original=verify_parent(a.root,parent)
    if sha(a.model)!=MODEL_SHA:raise ValueError('model changed')
    paths=[a.root/f'results/silicon_initiation_v13/dense_{n}/raw_hessian.npz' for n in ('loading8','return8')]
    sources=[arrays(p) for p in paths]
    domain=arrays(parent/'atomistic_pilot/common_domain.npz')
    common,direction=domain['center'],domain['direction']
    targets=[FixedGripTarget(s['positions'],s['basis'],s['hessian'],s['gradient'],300.,common) for s in sources]
    states=[parent_state(parent/f'atomistic_pilot/chain_{i}') for i in range(2)]
    a.output.mkdir(parents=True)
    protocol=dict(started_utc=datetime.now(timezone.utc).isoformat(),
        parent_package='results/silicon_force_thermal_comparison',parent_manifest_sha256=PARENT_MANIFEST_SHA,
        parent_final_retained_evaluations=[s['retained_evaluation'] for s in states],
        parent_completed_proposals=56,parent_saved_draws=32,new_warmup=0,additional_draws_per_chain=a.draws,
        step=original['step'],steps=original['steps'],temperature_K=300.,halfwidth_A=3.,minimum_pair_A=1.,
        max_seconds=a.max_seconds,allocation_guard_available_bytes=6000000000,
        runner_sha256=sha(Path(__file__)),continuation_core_sha256=sha(a.root/'solver_v1/silicon_chain_continuation.py'),
        inherited_source_sha256=original['source_sha256'],model_sha256=MODEL_SHA,
        exact_parent_rng_restored=True,no_reseed=True,no_adaptation=True,
        checkpoint='atomic single NPZ after each complete proposal, RNG included; incomplete paths retained as orphans',
        sequential_interleaving='one proposal per chain, prevents a one-chain-only budget allocation',
        target=original['target'],domain_is_not_crack_basin=True,domain_sensitivity_verified=False,
        physical_time=None,new_MD=0,new_DFT=0,kBT_fitted=False,production_changed=False)
    dump(a.output/'protocol.json',protocol)
    start=time.perf_counter();calls=0;counts=[0,0];error=None;phase='allocation';folders=[];atoms=[]
    active_chain=0;context='restart_verification';iteration=-1
    try:
        import psutil
        available=psutil.virtual_memory().available
        dump(a.output/'allocation_check.json',dict(available_bytes=available,required_bytes=6000000000,
            passed=available>=6000000000,checked_utc=datetime.now(timezone.utc).isoformat()))
        if available<6000000000:raise MemoryError('6 GB guard failed before loading torch/model')
        import torch
        from ase import Atoms
        from mace.calculators import MACECalculator
        from results.silicon_wafer_feasibility.mace_force_only_v9 import force_only
        torch.set_num_threads(2);torch.set_num_interop_threads(1)
        calc=MACECalculator(model_paths=str(a.model),device='cpu',default_dtype='float64')
        phase='verification'
        for i,s in enumerate(sources):
            folder=a.output/f'chain_{i}';folder.mkdir();folders.append(folder)
            atoms.append(Atoms(numbers=s['numbers'],positions=s['positions'],cell=s['cell'],pbc=False))
        def evaluate(u):
            nonlocal calls
            if time.perf_counter()-start>a.max_seconds:raise TimeoutError('bounded elapsed calculation budget reached')
            i=active_chain;target=targets[i];s=sources[i];r=target.coordinates(u);inside=target.inside(r)
            if inside:
                atoms[i].positions[:]=r;e,f=force_only(calc,atoms[i]);calls+=1
                phi=target.correction(u,e,float(s['energy']))
                gradient=nonlinear_correction_gradient(target,u,f)
                obs=np.array([e,float((r-common).ravel()@direction),
                    float(np.sqrt(np.mean(np.sum((r[target.free]-common[target.free])**2,axis=1))))])
                if not np.isfinite(f).all() or not np.isfinite([e,phi]).all() or not np.isfinite(gradient).all():
                    raise FloatingPointError('nonfinite actual evaluation')
            else:
                e=np.nan;f=np.full_like(r,np.nan);phi=np.inf;gradient=np.full_like(u,np.nan);obs=np.full(3,np.nan)
            index=counts[i]
            np.savez_compressed(folders[i]/f'evaluation_{index:04d}.npz',u=u,positions=r,energy=e,
                forces=f,correction=phi,gradient=gradient,inside=inside,observation=obs,context=context,iteration=iteration)
            counts[i]+=1
            dump(a.output/'progress.json',dict(chain=i,context=context,iteration=iteration,new_model_calls=calls,
                evaluations_recorded=sum(counts),committed_proposals=[s['completed'] for s in states],
                elapsed_seconds=time.perf_counter()-start,updated_utc=datetime.now(timezone.utc).isoformat()))
            if calls%16==0:print('CALLS',calls,'COMMITTED',[s['completed'] for s in states],flush=True)
            return phi,gradient,obs
        for active_chain,state in enumerate(states):
            current=evaluate(state['u'])
            parent_record=arrays(parent/f"atomistic_pilot/chain_{active_chain}/evaluation_{state['retained_evaluation']:04d}.npz")
            fresh=arrays(folders[active_chain]/'evaluation_0000.npz')
            errors=dict(energy_error_eV=abs(float(fresh['energy'])-float(parent_record['energy'])),
                force_error_eV_A=float(np.max(abs(fresh['forces']-parent_record['forces']))),
                correction_error=abs(current[0]-state['phi']),score_error=float(np.max(abs(current[1]-state['gradient']))))
            dump(folders[active_chain]/'restart_check.json',errors)
            if max(errors.values())>1e-8:raise ValueError('actual retained state replay mismatch')
            # Preserve the exact parent score/phi for the first kick; do not quietly refresh them.
            state['retained_evaluation']=0
            save_checkpoint(folders[active_chain]/'checkpoint_0000.npz',state)
        phase='sampling';context='sampling'
        for iteration in range(a.draws):
            for active_chain in range(2):
                old=states[active_chain];first=counts[active_chain]
                new,tx=advance(evaluate,old,step=protocol['step'],steps=protocol['steps'])
                new['retained_evaluation']=counts[active_chain]-1 if tx['accepted'] else old['retained_evaluation']
                tx.update(iteration=iteration,global_proposal=56+iteration,
                    start_evaluation=old['retained_evaluation'],first_evaluation=first,
                    last_evaluation=counts[active_chain]-1,retained_evaluation=new['retained_evaluation'])
                save_checkpoint(folders[active_chain]/f"checkpoint_{new['completed']:04d}.npz",new,tx)
                states[active_chain]=new
        phase='complete'
    except Exception as exc:
        error=type(exc).__name__+': '+str(exc)
    dump(a.output/'summary.json',dict(complete=phase=='complete',phase=phase,error=error,
        committed_proposals=[s['completed'] for s in states],evaluations_recorded=sum(counts),new_model_calls=calls,
        elapsed_seconds=time.perf_counter()-start,finished_utc=datetime.now(timezone.utc).isoformat(),
        incomplete_proposal_evaluations_preserved=True,new_MD=0,new_DFT=0,
        equilibrium_certified=False,free_energy_estimated=False,crack_probability_estimated=False,
        physical_clock=None,production_changed=False))
    print((a.output/'summary.json').read_text(encoding='utf-8'),flush=True)
    if phase!='complete':raise SystemExit(2)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True);p.add_argument('--model',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--draws',type=int,default=64)
    p.add_argument('--max-seconds',type=float,default=2400.);main(p.parse_args())
