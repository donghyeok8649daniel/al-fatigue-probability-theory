"""A fourth static anchor through SW preferred cosine, before force validation.

Motivation is the known four-neighbor harmonic identity and the independently
declared C44 calibration target, not the post-selection external force set.
This changes an existing SW parameter, not the Bessel representation. Full
nonzero angular bulk backgrounds and their derivatives must be kept.
"""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from results.silicon_wafer_feasibility.run_static_probe import source_parameters
from .run_silicon_sw_extended_bulk import angular_anchor
from .run_silicon_sw_material_fit import dump


def main(args):
    if args.output.exists():raise ValueError('fresh output required')
    args.output.mkdir(parents=True);start=time.perf_counter();root=Path(__file__).resolve().parents[1]
    bindings=[]
    for name in ['run_silicon_sw_angle_preference.py','run_silicon_sw_extended_bulk.py',
        'silicon_sw_anharmonic.py','silicon_sw_material_fit.py','silicon_environment_research.py']:
        relative='solver_v1/'+name;raw=(root/relative).read_bytes();target=args.output/'source_snapshots'/relative
        target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(raw)
        bindings.append(dict(path=relative,sha256=hashlib.sha256(raw).hexdigest()))
    dump(args.output/'source_manifest.json',bindings)
    ref=json.loads(args.reference.read_text());lattice=ref['diamond_a0'];c11,c12,c44=[ref['diamond_'+k] for k in ['c11','c12','c44']]
    p,_=source_parameters();rc=p['sigma']*p['a'];counts=dict(hessian=0)
    sigmas=args.sigmas;gammas=args.gammas;cosines=np.linspace(-.6,-.05,23)
    dump(args.output/'protocol.json',dict(reason='known C44 calibration mismatch of tetrahedral four-neighbor harmonic identity',
        cutoff_A=rc,sigma_values_A=sigmas,gamma_values=gammas,preferred_cosine_interval=[-.6,-.05],
        interval_is_research_bound_not_physical_prior=True,initial_cosine_samples=23,
        angle_betas=args.betas if not args.positive_boundary else None,
        beta_policy='beta=.5/(1-c0), monotone positive endpoint' if args.positive_boundary else 'fixed declared beta grid',
        anchors='a0 stationarity,C11,C12 at each cosine; first sign-change roots of relaxed C44-target',
        target_role='known published finite-strain calibration, not exact DFT zero-strain Hessian certification',
        keep_nonzero_angular_bulk_background=True,global_root_completeness_certified=False,
        external_force_set_used=False,force_geometries_evaluated=False,material_approved=False,new_DFT=0,new_MD=0))
    scans=[];accepted=[];profile=0
    for sigma in sigmas:
        for gamma,beta in [(gamma,beta) for gamma in gammas for beta in ([None] if args.positive_boundary else args.betas)]:
            shape=dict(p,sigma=sigma,a=rc/sigma,gamma=gamma);cache={};samples=[]
            def evaluate(c):
                key=float(c)
                if key not in cache:
                    actual_beta=.5/(1-key) if args.positive_boundary else beta
                    parameters,alpha,anchor=angular_anchor(dict(shape,costheta0=key),lattice,c11,c12,actual_beta,counts)
                    if np.any(alpha<=0):raise ValueError('nonpositive amplitude')
                    cache[key]=(parameters,alpha,anchor)
                return cache[key]
            for cosine in cosines:
                row=dict(sigma_A=sigma,gamma=gamma,angle_beta=.5/(1-float(cosine)) if args.positive_boundary else beta,
                         preferred_cosine=float(cosine))
                try:
                    _,alpha,anchor=evaluate(cosine)
                    row.update(status='admissible_local_bulk',C44_GPa=anchor['C44_GPa'],
                               residual_GPa=anchor['C44_GPa']-c44,amplitudes=alpha.tolist())
                except (ValueError,RuntimeError) as error:row.update(status='failed_or_inadmissible',reason=str(error))
                samples.append(row);scans.append(row)
            brackets=[(a['preferred_cosine'],b['preferred_cosine']) for a,b in zip(samples,samples[1:])
                if a['status']==b['status']=='admissible_local_bulk' and a['residual_GPa']*b['residual_GPa']<0]
            for lo,hi in brackets:
                cosine=brentq(lambda c:evaluate(c)[2]['C44_GPa']-c44,lo,hi,xtol=2e-13)
                parameters,alpha,anchor=evaluate(cosine)
                if abs(anchor['C44_GPa']-c44)>2e-6:raise RuntimeError('fourth static anchor unresolved')
                accepted.append(dict(profile=profile,status='stable_bulk_anchor',shape=dict(shape,costheta0=cosine),
                    parameters=parameters,amplitudes=alpha.tolist(),angle_beta=.5/(1-cosine) if args.positive_boundary else beta,
                    preferred_cosine=float(cosine),
                    sigma_A=sigma,gamma=gamma,cutoff_A=rc,anchor=anchor,elastic_gate_passed=True,root_bracket=[lo,hi]))
                profile+=1
            dump(args.output/'cosine_scan.json',scans);dump(args.output/'admissible_profiles.json',accepted)
            print(json.dumps(dict(sigma=sigma,gamma=gamma,beta=beta,sign_change_brackets=len(brackets),admissible_roots=len(accepted))),flush=True)
    dump(args.output/'summary.json',dict(complete=True,sampled_cosines=len(scans),admissible_root_profiles=len(accepted),
        new_Hessian_evaluations=counts['hessian'],elapsed_seconds=time.perf_counter()-start,
        force_dataset_evaluated=False,material_approved=False,first_initiation_validated=False,
        physical_clock_validated=False,new_DFT=0,new_MD=0))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--reference',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--betas',type=float,nargs='+',default=[0.])
    parser.add_argument('--positive-boundary',action='store_true')
    parser.add_argument('--sigmas',type=float,nargs='+',default=[1.6,2.1,2.8])
    parser.add_argument('--gammas',type=float,nargs='+',default=[.9,1.2,1.5])
    main(parser.parse_args())
