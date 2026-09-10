"""Budgeted deterministic LOCAL joint interface/source-core shape probe.

No new energy term; no material/default adoption. Normalize each data block by
its unchanged candidate's squared residual, so weight=1 minimizes their sum of
fractional errors. This is an explicit Pareto diagnostic, NOT uncertainty or a
unique Al fit. Preserve seven exact anchors and all supplied bulk inequalities.
"""
import argparse
import hashlib
import json
from dataclasses import replace
from functools import lru_cache
import time
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from .coordination_screening import CoordinationScreenedInterface
from .frozen_core_force_target import FrozenCoreForceTarget
from .interface_tangent_calibration import TangentCalibrationProblem,impose_tangents,NONNEGATIVE
from .report_current_material_core import restore_case
from .run_current_material_core import ROOT,load_current_material
from .run_source_core_reference import load_source_material,build_source_core
from .tail_constrained_material import spectral_profile
from .vector_material_calibration import MaterialObservation
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json


class ProbeBudgetReached(RuntimeError):
    pass


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--max-seconds',type=float,default=1800.)
    parser.add_argument('--max-nfev',type=int,default=8)
    parser.add_argument('--resume-probe',type=Path)
    parser.add_argument('--log-search-radius',type=float,default=.2)
    parser.add_argument('--expand-search-domain',action='store_true')
    args=parser.parse_args()
    if args.out.exists(): raise FileExistsError('fresh shape-probe output required')
    if args.max_seconds<=0 or not np.isfinite(args.max_seconds) or args.max_nfev<1:
        raise ValueError('positive declared calculation budget required')
    if not np.isfinite(args.log_search_radius) or args.log_search_radius<=0:
        raise ValueError('positive declared log-shape search radius required')
    began=time.perf_counter()
    model,_,binding=load_current_material()
    source,tensor,source_binding=load_source_material()
    source_path=ROOT/'results/current_material_core_v22/source_reference/escaped_plus_R8r5'
    core,field,*_=restore_case(source_path,source,tensor,source_binding,core_builder=build_source_core)
    frozen=FrozenCoreForceTarget(core,field,radius=2.)
    problem=TangentCalibrationProblem(ROOT/'results/fcc111_active_interface/coordination_screening_v19/joint_refinement')
    original=impose_tangents(problem.observations,problem.source_tangents)
    selected=[i for i,o in enumerate(original) if o.role=='fit']
    shape0=np.asarray(model.screened_shape); x0=np.log(shape0[:5])
    base_matrix=problem.matrix(tuple(shape0)); core0=frozen.coefficient_matrix(model)
    static_residual=(base_matrix@model.coefficients-np.array([o.target for o in original]))/np.array([o.scale for o in original])
    scale_static=float(np.linalg.norm(static_residual[selected]))
    scale_core=float(np.linalg.norm(core0['design']@model.coefficients-core0['target']))
    if min(scale_static,scale_core)<=1e-12: raise ValueError('nonzero baseline blocks required for fractional-error scaling')
    observations=[replace(o,scale=o.scale*scale_static) if o.role=='fit' else o for o in original]
    observations += [MaterialObservation(f'frozen_core_gradient_{i}',float(value),scale_core,'eV/L0','fit')
                     for i,value in enumerate(core0['target'])]
    lower=np.maximum(x0-args.log_search_radius,np.log([.7,2.,2.,1e-5,2.]))
    upper=np.minimum(x0+args.log_search_radius,np.log([8.,12.,12.,1e6,12.]))
    start=x0;resume=None
    if args.resume_probe:
        prior_definition=json.loads((args.resume_probe/'definition.json').read_bytes())
        prior_raw=(args.resume_probe/'completion.json').read_bytes()
        prior=json.loads(prior_raw)
        old_lower,old_upper=np.asarray(prior_definition['log_shape_bounds'])
        same_bounds=np.allclose([old_lower,old_upper],[lower,upper],atol=1e-13,rtol=0)
        expanded=(args.expand_search_domain and np.all(lower<=old_lower+1e-13)
                  and np.all(upper>=old_upper-1e-13))
        if (prior_definition['candidate_sha256']!=binding['parameter_sha256']
                or prior_definition['source_sha256']!=source_binding['parameter_sha256']
                or not (same_bounds or expanded)
                or not np.allclose([prior_definition['static_baseline_squared_loss'],prior_definition['core_baseline_squared_loss']],
                                   [scale_static**2,scale_core**2],rtol=1e-12,atol=0)):
            raise ValueError('continuation must preserve data, objective scales, material and original bounds')
        start=np.asarray(prior['optimizer']['last_accepted_coordinates'],float)
        if start.shape!=(5,) or np.any(~np.isfinite(start)) or np.any(start<lower) or np.any(start>upper):
            raise ValueError('valid previous accepted coordinates required, not an incidental FD trial')
        resume=dict(source=args.resume_probe.resolve().relative_to(ROOT).as_posix(),
            completion_sha256=hashlib.sha256(prior_raw).hexdigest(),
            accepted_coordinates=start,new_bounds_centered_on_previous_point=False,
            explicit_superset_search_box=bool(expanded and not same_bounds),
            previous_log_bounds=[old_lower,old_upper])
    save_json(args.out/'definition.json',dict(candidate_sha256=binding['parameter_sha256'],
        source_sha256=source_binding['parameter_sha256'],source_state=source_path.relative_to(ROOT).as_posix(),
        initial_shape=shape0,baseline_coefficients=model.coefficients,
        target_radius_over_L0=2.,source_free_radius_over_L0=8.,candidate_ring=7,
        exact_observations=[o.name for o in original if o.role=='exact'],
        static_baseline_squared_loss=scale_static**2,core_baseline_squared_loss=scale_core**2,
        objective='static_loss/static_baseline_loss + core_loss/core_baseline_loss',
        block_scales_are_not_uncertainties=True,core_block_weight=1.,
        log_shape_bounds=[lower,upper],screening_exponent_fixed=0.,
        log_search_radius=args.log_search_radius,
        local_probe_not_global_fit=True,max_nfev=args.max_nfev,max_seconds=args.max_seconds,
        actual_start_log_coordinates=start,resume=resume,
        no_new_energy_term=True,physical_yield_used=False,material_accepted=False,
        core_columns_are_energy_gradients=True,production_changed=False))
    profiles=[]
    @lru_cache(maxsize=96)
    def evaluate(coordinates):
        if time.perf_counter()-began>args.max_seconds: raise ProbeBudgetReached('declared wall budget reached')
        tick=time.perf_counter(); shape=np.r_[np.exp(coordinates),0.]
        probe_model=CoordinationScreenedInterface(shape,np.ones(10),law='power',tolerance=2e-12)
        static=problem.matrix(tuple(shape)); target=frozen.coefficient_matrix(probe_model)
        M=np.vstack([static,target['design']]); operators,tails=problem.operators(tuple(shape))
        fit=spectral_profile(M,observations,operators,tails,nonnegative=NONNEGATIVE)
        c=fit['coefficients']
        raw=(static@c-np.array([o.target for o in original]))/np.array([o.scale for o in original])
        error=target['design']@c-target['target']
        fit.update(shape=shape,static_squared_loss=float(raw[selected]@raw[selected]),
            core_force_rms_eV_L0=float(np.sqrt(np.mean(error**2))),
            evaluation_seconds=time.perf_counter()-tick,material_accepted=False,
            relaxed_core_recomputed=False)
        profiles.append(fit)
        save_json(args.out/'checkpoint.json',dict(completed=False,profiles=profiles))
        print(f'profile {len(profiles)}: joint={fit["squared_loss"]:.7g}, '
              f'static={fit["static_squared_loss"]:.7g}, coreRMS={fit["core_force_rms_eV_L0"]:.7g}, '
              f'LJpositive={fit["strictly_positive_LJ"]}, {fit["evaluation_seconds"]:.2f}s',flush=True)
        return fit
    def fun(x):
        fit=evaluate(tuple(x))
        return np.asarray(fit['residuals'])[fit['selected_rows']]
    try:
        result=least_squares(fun,start,bounds=(lower,upper),jac='2-point',diff_step=1e-3,
            max_nfev=args.max_nfev,ftol=1e-7,xtol=1e-7,gtol=2e-6)
        status=dict(success=bool(result.success),message=str(result.message),nfev=result.nfev,
            njev=result.njev,last_accepted_coordinates=result.x,last_accepted_loss=float(result.fun@result.fun),
            optimality=result.optimality,active_mask=result.active_mask)
    except (ArithmeticError,ValueError,ProbeBudgetReached) as error:
        status=dict(success=False,error_type=type(error).__name__,message=str(error))
    if not profiles:
        save_json(args.out/'completion.json',dict(completed=False,optimizer=status,material_accepted=False))
        return
    eligible=[p for p in profiles if p['strictly_positive_LJ']]
    best=min(eligible,key=lambda p:p['squared_loss']) if eligible else None
    endpoint=False
    if best is not None and 'last_accepted_coordinates' in status:
        endpoint=bool(np.allclose(best['shape'][:5],np.exp(status['last_accepted_coordinates']),rtol=1e-12,atol=0))
    write_csv(args.out/'profile_summary.csv',[dict(index=i,joint_loss=p['squared_loss'],
        static_loss=p['static_squared_loss'],core_RMS=p['core_force_rms_eV_L0'],positive_LJ=p['strictly_positive_LJ'],
        exact_residual=p['exact_residual'],kkt_residual=p['kkt_residual'],
        minimum_robust_margin=p['minimum_robust_margin'],seconds=p['evaluation_seconds']) for i,p in enumerate(profiles)])
    save_json(args.out/'completion.json',dict(completed=True,optimizer=status,best=best,
        selected_profile_is_optimizer_endpoint=endpoint,profiles_evaluated=len(profiles),
        material_accepted=False,full_material_calibration_complete=False,relaxed_core_recomputed=False,
        local_probe_only=True,production_changed=False,physical_yield_validated=False,
        physical_Hz=False,elapsed_seconds=time.perf_counter()-began))


if __name__=='__main__': main()
