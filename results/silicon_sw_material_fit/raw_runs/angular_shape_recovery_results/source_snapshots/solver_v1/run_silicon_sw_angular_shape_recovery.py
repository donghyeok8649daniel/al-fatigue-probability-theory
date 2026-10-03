"""Independent recovery after final NumPy-bool JSON summary serialization failed.

Do not rerun completed profile/validation/interface calculations. Verify source
snapshots, all polynomial designs, locked selection, full raw validation arrays,
and direct/Bessel recorded jets. Reconstruct counters from completed control
flow explicitly; the lost original perf_counter elapsed value is unavailable.
"""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from .run_silicon_sw_material_fit import dump
from .run_silicon_sw_guarded_fit import design
from .run_silicon_sw_stress_fit import source_stress
from results.silicon_wafer_feasibility.run_dft_material_audit import read_frames


def main(args):
    if args.output.exists():raise ValueError('fresh recovery output required')
    args.output.mkdir(parents=True);start=time.perf_counter();rows=json.loads((args.run/'profiles.json').read_text())
    selected=json.loads((args.run/'selected_fit.json').read_text());lock=json.loads((args.run/'selection_lock.json').read_text())
    if hashlib.sha256((args.run/'selected_fit.json').read_bytes()).hexdigest()!=lock['sha256']:raise RuntimeError('locked selection changed')
    sources=json.loads((args.run/'source_manifest.json').read_text())
    for r in sources:
        if hashlib.sha256((args.run/'source_snapshots'/r['path']).read_bytes()).hexdigest()!=r['sha256']:
            raise RuntimeError('calculation source snapshot differs')
    if len(rows)!=36 or any(r['status'] not in ['completed_feasible','no_feasible_beta'] for r in rows):
        raise RuntimeError('cannot reconstruct counters across unresolved profile failures')
    verified=[]
    for row in rows:
        if row['fourth_beta_maximum_error']>3e-9:raise RuntimeError('fourth-beta identity not passed')
        raw=np.load(args.run/f"profile_{row['profile']}_polynomial_design.npz");alpha=raw['amplitudes']
        if row['status']=='completed_feasible':
            beta=row['selected']['beta'];x=raw['design'];b=raw['bulk_design'];e=raw['EOS_error_basis']
            residual=(x[0]+beta*x[1]+beta**2*x[2])@alpha-raw['target']
            er=e[0]+beta*e[1]+beta**2*e[2];objective=float(residual@residual+np.mean((er/.02)**2))
            br=(b[0]+beta*b[1]+beta**2*b[2])@alpha-raw['bulk_target'];guard=float(br@br)
            if abs(objective-row['selected']['objective_squared'])>2e-8 or abs(guard-row['selected']['bulk_energy_squared'])>2e-9:
                raise RuntimeError('independent design objective/guard differs')
            verified.append(dict(profile=row['profile'],objective_squared=objective,bulk_energy_squared=guard))
    winner=min([r for r in rows if r['status']=='completed_feasible'],key=lambda r:r['selected']['objective_squared'])
    if selected['profile']!=winner['profile']:raise RuntimeError('saved winner differs from training-only candidates')
    features=np.load(args.run/'selected_features.npz');frames,_=read_frames(args.archive)
    fit=json.loads((args.baseline/'fit.json').read_text());training=fit['training_frames'];anchor=fit['anchor_frame']
    e,f,offsets,re,rf=[features[k] for k in ['energy_basis','force_basis','offsets','reference_energy','reference_force']]
    if len(e)!=len(frames) or int(offsets[-1])!=sum(len(a) for a in frames):raise RuntimeError('raw validation geometry count differs')
    alpha=np.array(selected['amplitudes']);ea=e[anchor]/len(frames[anchor]);ra=re[anchor]/len(frames[anchor])
    metrics=json.loads((args.run/'frame_metrics.json').read_text());errors=[]
    for record in metrics:
        i=record['frame'];lo,hi=offsets[i:i+2];diff=np.einsum('nck,k->nc',f[lo:hi],alpha)-rf[lo:hi]
        error=abs(float(np.sum(diff**2))-record['force_squared_error']);errors.append(error)
        if error>2e-8:raise RuntimeError('raw force metrics differ')
        if record['relative_energy_error_eV_atom'] is not None:
            energy=float((e[i]/len(frames[i])-ea)@alpha-(re[i]/len(frames[i])-ra))
            if abs(energy-record['relative_energy_error_eV_atom'])>2e-10:raise RuntimeError('raw relative energy metrics differ')
    if len(metrics)!=len(frames) or len(set(r['frame'] for r in metrics))!=len(frames):raise RuntimeError('validation metrics incomplete')
    audit=json.loads((args.run/'bessel_audit.json').read_text());consistent=True
    if len(audit)!=8:raise RuntimeError('interface audit incomplete')
    for r in audit:
        b,d=r['bessel'],r['direct'];ev=abs(b['energy']-d['energy']);g=float(np.max(abs(np.array(b['gradient'])-d['gradient'])));h=float(np.max(abs(np.array(b['hessian'])-d['hessian'])))
        if max(abs(ev-r['energy_error_eV']),abs(g-r['gradient_error_eV_A']),abs(h-r['hessian_error_eV_A2']))>1e-14:raise RuntimeError('recorded Bessel jet errors differ')
        consistent&=ev<2e-8 and g<2e-7 and h<3e-6
    prior=json.loads((args.prior/'profiles.json').read_text());cached=0
    for row in rows:
        for beta in [-.375,0.,.375]:
            cached+=int(any(r['shape']['sigma']==row['sigma_A'] and r['shape']['gamma']==row['gamma'] and r['angle_beta']==beta for r in prior))
    new_controls=len(rows)*3-cached
    counters=dict(hessian=len(rows)*5,feature=new_controls*len(training)+len(frames),eos=new_controls*12,
                  identity_controls=len(rows)*12,cached_designs=cached)
    dump(args.output/'verified_profiles.json',verified)
    dump(args.output/'summary.json',dict(complete=True,original_exit_code=1,original_failure='only final summary JSON rejected NumPy bool; all raw profiles and full validation saved',
        recovered_by='independent raw array/source/selection/jet audit; no repeated geometry calculation',profiles=len(rows),
        selected_profile=selected['profile'],selected_beta=selected['selected']['beta'],
        reconstructed_counts=counters,counter_provenance='completed fixed loops and cached branch identity; not recovered live counter memory',
        original_elapsed_seconds=None,recovery_elapsed_seconds=time.perf_counter()-start,
        source_snapshots_verified=True,locked_selection_verified=True,validation_frames_verified=len(metrics),
        maximum_force_squared_metric_difference=max(errors),bessel_consistency_passed=bool(consistent),
        material_approved=False,first_initiation_validated=False,physical_clock_validated=False,new_DFT=0,new_MD=0))
    raw=Path(__file__).read_bytes();out=args.output/'source_snapshots/solver_v1/run_silicon_sw_angular_shape_recovery.py';out.parent.mkdir(parents=True);out.write_bytes(raw)
    dump(args.output/'source_manifest.json',[dict(path='solver_v1/run_silicon_sw_angular_shape_recovery.py',sha256=hashlib.sha256(raw).hexdigest())])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['run','prior','baseline','archive','output']:p.add_argument('--'+name,required=True,type=Path)
    main(p.parse_args())
