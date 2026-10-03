"""Recompute only four failed shape features; independently check cached winner."""
from pathlib import Path
import sys,json,hashlib,time
import numpy as np
cache=Path(__file__).resolve().parent
root=cache.parents[2]/'aft-silicon-wafer'
sys.path.insert(0,str(root))
from results.silicon_wafer_feasibility.run_dft_material_audit import read_frames
from results.silicon_wafer_feasibility.run_static_probe import source_parameters
from solver_v1.run_silicon_sw_guarded_fit import design
from solver_v1.silicon_sw_material_fit import force_features,fit_with_bulk_energy_guard
out=cache/'dual_repair_results'
assert not out.exists();out.mkdir()
start=time.perf_counter();frames,_=read_frames(root/'.cache/si-atomistic-v5/Si_PRX_GAP.zip')
fit=json.loads((cache/'results/fit.json').read_text());training=fit['training_frames'];anchor=fit['anchor_frame']
old=np.load(cache/'results/features.npz');e0=old['energy_basis'];f0=old['force_basis']
offsets=old['offsets'];ref_e=old['reference_energy'];ref_f=old['reference_force']
_,_,b0,t0=design(e0,f0,training,frames,offsets,ref_e,ref_f,anchor)
bound=float(np.sum((b0@np.ones(3)-t0)**2));p,_=source_parameters();rc=p['sigma']*p['a']
profiles=json.loads((cache/'guarded_results/profiles.json').read_text())
chosen=json.loads((cache/'guarded_results/selected_fit.json').read_text())
targets=[r for r in profiles if r['status']!='converged_feasible']+[chosen]
results=[];evaluations=0
for row in targets:
    shape=dict(p,sigma=row['sigma_A'],a=rc/row['sigma_A'],gamma=row['gamma'])
    reused=row['profile']==chosen['profile']
    if reused:
        saved=np.load(cache/'guarded_results/selected_features.npz')
        e,f=saved['energy_basis'],saved['force_basis']
    else:
        e,f=np.zeros_like(e0),np.zeros_like(f0)
        for i in training:
            a=frames[i];lo,hi=offsets[i:i+2]
            e[i],f[lo:hi]=force_features(a.positions,a.cell.array,a.pbc,shape)
            evaluations+=1
    x,y,b,t=design(e,f,training,frames,offsets,ref_e,ref_f,anchor)
    coefficients,record=fit_with_bulk_energy_guard(x,y,b,t,bound)
    np.savez_compressed(out/f"profile_{row['profile']}_design.npz",design=x,target=y,
                        bulk_design=b,bulk_target=t,maximum_bulk_sq=bound)
    result=dict(profile=row['profile'],sigma_A=row['sigma_A'],gamma=row['gamma'],
                previous_status=row['status'],previous_reason=row.get('reason'),
                cached_features=reused,fit=record)
    results.append(result)
    print(json.dumps(result),flush=True)
    (out/'profiles.json').write_text(json.dumps(results,indent=2)+'\n')
summary=dict(complete=True,success=len(results)==5,new_feature_evaluations=evaluations,
             recovered_failed_profiles=4,cached_selected_profile=chosen['profile'],
             old_selected_objective=chosen['fit']['objective_squared'],
             all_recovered_objectives_exceed_previous_best=all(r['fit']['objective_squared']>
                    chosen['fit']['objective_squared'] for r in results if not r['cached_features']),
             elapsed_seconds=time.perf_counter()-start,material_approved=False)
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
sources=[dict(path='solver_v1/silicon_sw_material_fit.py',sha256=hashlib.sha256(
            (root/'solver_v1/silicon_sw_material_fit.py').read_bytes()).hexdigest()),
         dict(path='solver_v1/run_silicon_sw_guarded_fit.py',sha256=hashlib.sha256(
            (root/'solver_v1/run_silicon_sw_guarded_fit.py').read_bytes()).hexdigest()),
         dict(path='repair_failed_profiles.py',sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())]
(out/'source_manifest.json').write_text(json.dumps(sources,indent=2)+'\n')
print(json.dumps(summary))
