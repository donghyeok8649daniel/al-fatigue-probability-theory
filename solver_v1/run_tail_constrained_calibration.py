"""Actual same-family deterministic calibration with analytic spectral tails."""
import argparse
from dataclasses import asdict, replace
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import least_squares
from scipy.special import expit, logit

from .range_resolved_material import RangeObservationCache
from .run_vector_material_calibration import source_and_targets
from .run_vector_registry_audit import save_json
from .run_low_stress_cyclic_diagnostic import write_csv
from .tail_constrained_material import TailBulkCoefficientBasis, augmented_wavepoints, spectral_profile
from .yield_elastic_metric import cubic_metric_problem, cubic_to_mode_matrix


ROOT=Path(__file__).resolve().parents[1]
DEFAULT_OUT=ROOT/'results/fcc111_active_interface/tail_calibration_v14/baseline'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=DEFAULT_OUT)
    parser.add_argument('--max-nfev',type=int,default=28)
    parser.add_argument('--radius',type=float,default=12.)
    parser.add_argument('--wavepoint-step',type=float,default=.25,
        help='stability-constraint grid spacing in the same corrected cubic reciprocal coordinates')
    parser.add_argument('--additional-wavepoints',type=Path,
        help='JSON {points,provenance}: explicitly recorded independent stability counterexamples')
    parser.add_argument('--stability-stretches',type=float,nargs='+',default=[1.],
        help='fixed-material isotropic geometry dilations for stability constraints; not temperature/target fitting')
    parser.add_argument('--upper-decays',type=float,nargs=3,default=[8.,12.,12.],
                        help='explicit search bounds for a recorded range-continuation test')
    parser.add_argument('--symmetry-channel',action='store_true',help='separate one-amplitude nested Eg test')
    parser.add_argument('--exact-bulk',action='store_true',help='staged bulk-constrained interface calibration, not zero target uncertainty')
    parser.add_argument('--cubic-embedding',action='store_true',help='separate one-coefficient density-shape extension after measured failure')
    parser.add_argument('--quartic-angular',action='store_true',help='one per-site squared angular-norm extension, not density cubic')
    parser.add_argument('--saturating-angular',action='store_true',help='nested nonnegative rational angular shape parameter')
    parser.add_argument('--starts',type=Path,help='explicit recorded deterministic start JSON list')
    parser.add_argument('--difference-scheme',choices=('2-point','3-point'),default='2-point',
        help='declared profiled-residual derivative scheme, not a material change')
    parser.add_argument('--difference-step',type=float,default=2e-4,
        help='relative step in log shape coordinates; independently refine to audit stationarity')
    parser.add_argument('--interface-development',action='store_true',
        help='v15: reclassify inspected shape failures as development; new off-grid tests stay held out')
    parser.add_argument('--mixed-density',action='store_true',
        help='v15 nested positive two-exponential site density; fourth shape is a mixing weight')
    parser.add_argument('--static-bloch-targets',action='store_true',
        help='source spatial Hessian magnitudes, NOT phonon Hz; separate declared fit/heldout q')
    parser.add_argument('--density-angular-cross',action='store_true',
        help='one density/rank3 mixed invariant with an explicit convex nonlinear-sector constraint')
    args=parser.parse_args(); out=args.out
    if not np.isfinite(args.difference_step) or args.difference_step<=0:
        raise ValueError('finite positive declared numerical derivative step required')
    stretches=sorted(set([1.]+args.stability_stretches))
    if np.any(~np.isfinite(stretches)) or min(stretches)<=0:
        raise ValueError('finite positive explicit stability stretches required')
    if len(stretches)>1 and (not args.quartic_angular or args.mixed_density or args.cubic_embedding):
        raise ValueError('fixed-material stretched stability currently supports quartic/rational/cross families only')
    if np.any(~np.isfinite(args.upper_decays)) or np.any(np.asarray(args.upper_decays)<=np.array([.7,2.,2.])):
        raise ValueError('finite upper decays greater than lower bounds required')
    if out.exists(): raise FileExistsError('a fresh result directory is required')
    started=time.perf_counter(); source,observations,states=source_and_targets()
    observation_provenance=None
    if args.interface_development:
        from .interface_development_targets import development_observations
        observations,observation_provenance=development_observations(source,observations,states)
    bloch_targets=()
    if args.static_bloch_targets:
        from .static_bloch_targets import source_bloch_targets,append_bloch_problem
        bloch_targets=source_bloch_targets(source.reference)
    starts=[[2.579101374767287,5.996534839345317,4.962329189737055],
            [5.3707449112586065,5.6320174463985,3.2556291351974864],
            [3.,8.,7.]]
    if args.starts: starts=json.loads(args.starts.read_bytes())
    if args.mixed_density:
        if not (args.interface_development and args.quartic_angular and args.symmetry_channel) or args.saturating_angular or args.cubic_embedding:
            raise ValueError('density-mixture ablation requires the v15 quartic development family alone')
        if not args.starts:
            starts=[[2.5482578561,4.6782110749,11.0487420063,w] for w in (.2,.7)]
        if any(len(s)!=4 or not 0<s[3]<1 for s in starts):
            raise ValueError('mixed starts require three positive decays and an interior weight')
    if args.saturating_angular:
        if not args.quartic_angular: raise ValueError('rational ablation requires the quartic angular family')
        if not args.starts:
            # Three fixed starts around the already completed quartic basin;
            # alpha is a dimensionless moment-shape parameter, never an area.
            starts=[[2.782954,5.646841,11.999,alpha] for alpha in (1.e3,1.e5,1.e7)]
    wavepoint_provenance=None;additional=()
    if args.additional_wavepoints:
        raw=args.additional_wavepoints.read_bytes();record=json.loads(raw)
        if not record.get('provenance'):raise ValueError('additional wavepoint provenance required')
        additional=record['points']
        wavepoint_provenance=dict(record=record,sha256=hashlib.sha256(raw).hexdigest())
    points=augmented_wavepoints(args.wavepoint_step,additional);cache=RangeObservationCache(observations)
    if args.symmetry_channel:
        from .symmetry_resolved_material import SymmetryObservationCache, SymmetryTailBulkBasis
        cache=SymmetryObservationCache(observations)
        bulk_type=SymmetryTailBulkBasis
        nonnegative=(0,1,2,4,5,8)
    else:
        bulk_type=TailBulkCoefficientBasis
        nonnegative=(0,1,2,4,5)
    if args.cubic_embedding:
        if not args.symmetry_channel: raise ValueError('this runner requires the separately tested Eg family first')
        from .symmetry_resolved_material import CubicSymmetryObservationCache, CubicSymmetryTailBulkBasis
        cache=CubicSymmetryObservationCache(observations);bulk_type=CubicSymmetryTailBulkBasis
    if args.quartic_angular:
        if not args.symmetry_channel or args.cubic_embedding:
            raise ValueError('separate nested angular ablation; no simultaneous cubic-density extension')
        from .quartic_angular_material import QuarticSymmetryObservationCache, QuarticSymmetryTailBulkBasis
        cache=QuarticSymmetryObservationCache(observations);bulk_type=QuarticSymmetryTailBulkBasis
        nonnegative=(0,1,2,4,5,8,9)
    if args.mixed_density:
        from .mixed_density_material import MixedDensityObservationCache, MixedDensityTailBulkBasis
        cache=MixedDensityObservationCache(observations);bulk_type=MixedDensityTailBulkBasis
    profile_function=spectral_profile
    if args.density_angular_cross:
        if not (args.interface_development and args.quartic_angular) or args.mixed_density or args.saturating_angular or args.cubic_embedding:
            raise ValueError('cross-invariant ablation requires the unchanged quartic v15 family alone')
        from .density_angular_cross_material import CrossDensityAngularCache,CrossDensityAngularBulk,cross_sector_profile
        cache=CrossDensityAngularCache(observations);bulk_type=CrossDensityAngularBulk
        profile_function=cross_sector_profile
    def problem(decays,basis=None):
        raw=cache.matrix(decays)
        matrix,obs=cubic_metric_problem(raw[:,:8],observations)
        if args.symmetry_channel:
            extra=raw[:,8:].copy();extra[2:5]=np.linalg.solve(cubic_to_mode_matrix(),extra[2:5])
            matrix=np.column_stack([matrix,extra])
        if args.exact_bulk:
            obs=[replace(o,role='exact') if i<5 else o for i,o in enumerate(obs)]
        if bloch_targets:
            if basis is None:basis=bulk_type(decays if args.mixed_density else decays[:3],radius=args.radius)
            matrix,obs,_=append_bloch_problem(matrix,obs,basis,bloch_targets)
        return matrix,obs
    profiles=[]; stops=[]; failed_profiles=[]
    save_json(out/'definition.json',dict(source_sha256=source.reference.sha256,
        observations=[asdict(o) for o in observations],source_states=states,starts=starts,
        interface_development=bool(args.interface_development),
        observation_provenance=observation_provenance,
        static_bloch_targets=[asdict(t) for t in bloch_targets],
        static_bloch_target_scale='10% per directional restoring stiffness; model discrepancy, NOT measured uncertainty',
        wavepoints_cubic=points,radius_over_L0=args.radius,
        wavepoint_grid_step=args.wavepoint_step,additional_wavepoint_provenance=wavepoint_provenance,
        stability_stretches=stretches,stability_row_order='stretch-major, then declared wavepoints',
        fixed_material_gauge_during_strain=True,
        criterion='min eig(H_R) >= sum |c_j| analytic FCC Voronoi-tail envelope_j',
        whole_zone_proof=False,new_energy_terms=bool(args.symmetry_channel),
        extra_amplitude='D_E>=0' if args.symmetry_channel else None,
        exact_bulk_stage=bool(args.exact_bulk),
        cubic_density_extension=bool(args.cubic_embedding),
        quartic_angular_extension=bool(args.quartic_angular),
        rational_angular_extension=bool(args.saturating_angular),
        density_mixture_extension=bool(args.mixed_density),
        density_angular_cross_extension=bool(args.density_angular_cross),
        environment_sector_constraint='[[C,E_xI/2],[E_xI/2,K3]] positive semidefinite' if args.density_angular_cross else None,
        fourth_shape_coordinate='logit density weight' if args.mixed_density else 'log shape if present',
        density_weight_search_bounds=expit(np.array([-9.,9.])) if args.mixed_density else None,
        lower_decays=[.7,2.,2.],upper_decays=args.upper_decays,
        numerical_derivative_scheme=args.difference_scheme,numerical_relative_step=args.difference_step,
        max_nfev_per_start=args.max_nfev,
        physical_kinetics_calibrated=False))
    @lru_cache(maxsize=192)
    def evaluate(logs):
        before=time.perf_counter(); decays=np.exp(logs)
        if args.mixed_density:decays[3]=expit(logs[3])
        basis=bulk_type(decays if args.mixed_density else decays[:3],radius=args.radius)
        matrix,obs=problem(decays,basis)
        computed=[]
        for stretch in stretches:
            state_basis=basis
            if stretch!=1.:
                from .isotropic_bulk_validation import IsotropicBulkBasis
                state_basis=IsotropicBulkBasis(decays[:3],stretch=stretch,radius=args.radius,
                    include_cross=args.density_angular_cross)
            computed.extend(state_basis.evaluate(q) for q in points)
        columns,tails=map(np.asarray,zip(*computed))
        try:
            fit=profile_function(matrix,obs,columns,tails,nonnegative=nonnegative)
        except (ArithmeticError,ValueError) as error:
            failed_profiles.append(dict(log_shape=logs,decays=decays,message=str(error),
                elapsed_seconds=time.perf_counter()-before))
            save_json(out/'failed_profiles.json',dict(completed=False,failures=failed_profiles))
            raise
        hold=[i for i,o in enumerate(obs) if o.role=='heldout']
        row=dict(decays=decays,**fit,cubic_GPa=(matrix@fit['coefficients'])[2:5],
            heldout_normalized_rms=float(np.sqrt(np.mean(fit['residuals'][hold]**2))),
            seconds=time.perf_counter()-before)
        profiles.append(row)
        print(f"profile {len(profiles)} loss={row['squared_loss']:.8g} C={row['cubic_GPa']} "
              f"cuts={row['cut_iterations']} margin={row['minimum_robust_margin']:.3g} "
              f"{row['seconds']:.2f}s",flush=True)
        save_json(out/'progress.json',dict(completed=False,profiles=len(profiles),
            best_loss=min(p['squared_loss'] for p in profiles),elapsed_seconds=time.perf_counter()-started))
        save_json(out/'checkpoint.json',dict(completed=False,profiles=profiles,optimizers=stops))
        return row
    for start in starts:
        def fun(logs):
            row=evaluate(tuple(logs)); return row['residuals'][row['selected_rows']]
        try:
            bounds=[[.7,2.,2.],args.upper_decays]
            if args.saturating_angular: bounds=[bounds[0]+[1e-3],bounds[1]+[1e9]]
            initial=np.log(start);log_bounds=np.log(bounds)
            if args.mixed_density:
                initial[3]=logit(start[3]);log_bounds=np.column_stack([log_bounds,[-9.,9.]])
            run=least_squares(fun,initial,jac=args.difference_scheme,diff_step=args.difference_step,
                bounds=log_bounds,max_nfev=args.max_nfev,
                ftol=1e-7,xtol=1e-7,gtol=2e-6)
            final=np.exp(run.x)
            if args.mixed_density:final[3]=expit(run.x[3])
            stops.append(dict(start=start,success=bool(run.success),message=str(run.message),
                nfev=run.nfev,njev=run.njev,optimality=float(run.optimality),
                final_decays=final,final_loss=float(run.fun@run.fun)))
        except (ArithmeticError,ValueError) as exc:
            stops.append(dict(start=start,success=False,message=str(exc),numerical_stop=True))
            print(f'REJECTED numerical profile in start {start}: {exc}',flush=True)
        save_json(out/'checkpoint.json',dict(completed=False,profiles=profiles,optimizers=stops))
    best=min((r for r in profiles if r['strictly_positive_LJ']),key=lambda r:r['squared_loss'])
    matrix,obs=problem(best['decays'])
    write_csv(out/'residuals.csv',[dict(observable=o.name,role=o.role,units=o.units,
        target=o.target,prediction=p,scale=o.scale,normalized_residual=r)
        for o,p,r in zip(obs,best['predictions'],best['residuals'])])
    save_json(out/'calibration.json',dict(completed=True,best=best,profiles=profiles,
        optimizers=stops,failed_profiles=failed_profiles,elapsed_seconds=time.perf_counter()-started,
        source_sha256=source.reference.sha256,material_accepted=False))
    save_json(out/'failed_profiles.json',dict(completed=True,failures=failed_profiles))
    save_json(out/'progress.json',dict(completed=True,profiles=len(profiles),best_loss=best['squared_loss']))
    print('completed actual calibration',best['squared_loss'],flush=True)


if __name__=='__main__': main()
