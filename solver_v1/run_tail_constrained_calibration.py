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
from .material_calibration_controls import append_fixed_pair, select_profile_result, feasible_simplex_search


ROOT=Path(__file__).resolve().parents[1]
DEFAULT_OUT=ROOT/'results/fcc111_active_interface/tail_calibration_v14/baseline'


def finalize_profile_run(out, profiles, stops, failed_profiles, observations, *, elapsed_seconds, source_sha256):
    """Save actual finished attempts, including the no-admissible-LJ outcome."""
    best, admissible = select_profile_result(profiles)
    write_csv(out/'residuals.csv',[dict(observable=o.name,role=o.role,units=o.units,
        target=o.target,prediction=p,scale=o.scale,normalized_residual=r)
        for o,p,r in zip(observations,best['predictions'],best['residuals'])])
    save_json(out/'calibration.json',dict(completed=True,best=best,profiles=profiles,
        optimizers=stops,failed_profiles=failed_profiles,elapsed_seconds=elapsed_seconds,
        source_sha256=source_sha256,admissible_solution_found=admissible,
        best_is_diagnostic_closure_only=not admissible,material_accepted=False))
    save_json(out/'failed_profiles.json',dict(completed=True,failures=failed_profiles))
    save_json(out/'progress.json',dict(completed=True,profiles=len(profiles),best_loss=best['squared_loss'],
        admissible_solution_found=admissible))
    print('completed actual calibration',best['squared_loss'],'admissible LJ:',admissible,flush=True)


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
    parser.add_argument('--even-development',action='store_true',
        help='v16: inspected v15 holdouts become development, with new heldout states')
    parser.add_argument('--quadrupole-saturation',action='store_true',
        help='v16: one shear-gauge-normalized shape for the two rank2 site invariants')
    parser.add_argument('--fixed-pair-from',type=Path,
        help='explicit embedding-only ablation: inherit two positive LJ coefficients from a completed fit')
    parser.add_argument('--optimizer',choices=('least-squares','feasible-simplex'),default='least-squares',
        help='simplex rejects undefined profiles with +inf; no unstable finite penalty fit')
    parser.add_argument('--profile-only',action='store_true',
        help='re-solve coefficients at each exact declared shape; no shape optimization')
    parser.add_argument('--normal-development',action='store_true',
        help='v17: inspected v16 states become development, new normal/vector validation')
    parser.add_argument('--rank-one-range',action='store_true',
        help='v17: one independent rank1 exponential range, not a new force law')
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
    fixed_pair=None
    if args.fixed_pair_from:
        raw=(args.fixed_pair_from/'calibration.json').read_bytes();inherited=json.loads(raw)
        if not inherited['completed'] or not inherited['best']['strictly_positive_LJ']:
            raise ValueError('completed positive-LJ reference required for pair control')
        if inherited['source_sha256']!=source.reference.sha256:
            raise ValueError('pair reference target-source mismatch')
        fixed_pair=dict(coefficients=inherited['best']['coefficients'][:2],
            calibration_sha256=hashlib.sha256(raw).hexdigest(),
            reference_run=args.fixed_pair_from.name,
            interpretation='inherited LJ control, NOT two new observed Al data')
    observation_provenance=None
    if args.interface_development:
        from .interface_development_targets import development_observations
        observations,observation_provenance=development_observations(source,observations,states)
    if args.even_development:
        if args.interface_development:
            raise ValueError('select one target generation, not both v15 and v16')
        from .interface_even_development_targets import even_development_observations
        observations,observation_provenance=even_development_observations(source,observations,states)
    if args.normal_development:
        if args.even_development or args.interface_development:
            raise ValueError('choose exactly one generation of development observations')
        from .interface_normal_development_targets import normal_development_observations
        observations,observation_provenance=normal_development_observations(source,observations,states)
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
    if args.quadrupole_saturation:
        if (not (args.even_development or args.normal_development) or not args.quartic_angular or not args.symmetry_channel
                or args.saturating_angular or args.mixed_density or args.density_angular_cross
                or args.cubic_embedding):
            raise ValueError('even saturation is a separate v16 quartic-family ablation')
        if not args.starts:
            starts=[[2.5482578561,4.6782110749,11.0487420063,a] for a in (1.,10.)]
        if any(len(s)!=(5 if args.rank_one_range else 4) or s[3]<=0 for s in starts):
            raise ValueError('three decays and positive saturation starts required')
    if args.rank_one_range:
        if not (args.normal_development and args.quadrupole_saturation and args.starts):
            raise ValueError('independent rank1 range requires v17 saturation and explicit five-component starts')
        if any(s[4]<=0 for s in starts):
            raise ValueError('positive rank1 microscopic decay required')
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
    if args.quadrupole_saturation:
        from .quadrupole_saturation import SaturatedQuadrupoleCache,SaturatedQuadrupoleBulk
        cache=SaturatedQuadrupoleCache(observations);bulk_type=SaturatedQuadrupoleBulk
    if args.rank_one_range:
        from .rank_one_range_material import RankOneRangeCache, RankOneRangeBulk
        cache=RankOneRangeCache(observations);bulk_type=RankOneRangeBulk
    def problem(decays,basis=None):
        raw=cache.matrix(decays)
        matrix,obs=cubic_metric_problem(raw[:,:8],observations)
        if args.symmetry_channel:
            extra=raw[:,8:].copy();extra[2:5]=np.linalg.solve(cubic_to_mode_matrix(),extra[2:5])
            matrix=np.column_stack([matrix,extra])
        if args.exact_bulk:
            obs=[replace(o,role='exact') if i<5 else o for i,o in enumerate(obs)]
        if bloch_targets:
            if basis is None:basis=bulk_type(decays if (args.mixed_density or args.rank_one_range) else decays[:3],radius=args.radius)
            matrix,obs,_=append_bloch_problem(matrix,obs,basis,bloch_targets)
        matrix,obs=append_fixed_pair(matrix,obs,fixed_pair)
        return matrix,obs
    profiles=[]; stops=[]; failed_profiles=[]
    save_json(out/'definition.json',dict(source_sha256=source.reference.sha256,
        observations=[asdict(o) for o in observations],source_states=states,starts=starts,
        interface_development=bool(args.interface_development),
        even_development=bool(args.even_development),
        normal_development=bool(args.normal_development),
        rank_one_range_extension=bool(args.rank_one_range),
        quadrupole_saturation_extension=bool(args.quadrupole_saturation),
        fixed_pair_control=fixed_pair,
        optimizer_backend=args.optimizer,
        profile_only=bool(args.profile_only),
        quadrupole_saturation_gauge='fixed reference 2||d_gamma Q2||^2; Eg gauge=1' if args.quadrupole_saturation else None,
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
        rank1_decay_bounds=[2.,12.] if args.rank_one_range else None,
        physical_kinetics_calibrated=False))
    @lru_cache(maxsize=192)
    def evaluate(logs):
        before=time.perf_counter(); decays=np.exp(logs)
        if args.mixed_density:decays[3]=expit(logs[3])
        basis=bulk_type(decays if (args.mixed_density or args.rank_one_range) else decays[:3],radius=args.radius)
        matrix,obs=problem(decays,basis)
        computed=[]
        for stretch in stretches:
            state_basis=basis
            if stretch!=1.:
                from .isotropic_bulk_validation import IsotropicBulkBasis
                state_basis=IsotropicBulkBasis(decays[:3],stretch=stretch,radius=args.radius,
                    include_cross=args.density_angular_cross,
                    rank1_decay=decays[4] if args.rank_one_range else None)
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
            if args.quadrupole_saturation: bounds=[bounds[0]+[1e-5],bounds[1]+[1e6]]
            if args.rank_one_range: bounds=[bounds[0]+[2.],bounds[1]+[12.]]
            initial=np.log(start);log_bounds=np.log(bounds)
            if args.mixed_density:
                initial[3]=logit(start[3]);log_bounds=np.column_stack([log_bounds,[-9.,9.]])
            if args.profile_only:
                if np.any(initial<log_bounds[0]) or np.any(initial>log_bounds[1]):
                    raise ValueError('declared profile-only shape is outside the recorded bounds')
                row=evaluate(tuple(initial))
                stops.append(dict(start=start,success=True,message='coefficient profile verified; shape not optimized',
                    backend='coefficient-profile-only',nfev=1,njev=0,optimality=None,
                    shape_optimization_performed=False,final_decays=row['decays'],final_loss=row['squared_loss']))
                save_json(out/'checkpoint.json',dict(completed=False,profiles=profiles,optimizers=stops))
                continue
            if args.optimizer=='feasible-simplex':
                run=feasible_simplex_search(fun,initial,log_bounds,max_evaluations=args.max_nfev)
            else:
                run=least_squares(fun,initial,jac=args.difference_scheme,diff_step=args.difference_step,
                    bounds=log_bounds,max_nfev=args.max_nfev,
                    ftol=1e-7,xtol=1e-7,gtol=2e-6)
            final=np.exp(run.x)
            if args.mixed_density:final[3]=expit(run.x[3])
            stops.append(dict(start=start,success=bool(run.success),message=str(run.message),
                backend=args.optimizer,nfev=run.nfev,njev=getattr(run,'njev',None),
                optimality=float(run.optimality) if hasattr(run,'optimality') else None,
                rejected_profiles=getattr(run,'rejected_profiles',[]),
                final_decays=final,final_loss=float(run.fun if args.optimizer=='feasible-simplex' else run.fun@run.fun)))
        except (ArithmeticError,ValueError) as exc:
            stops.append(dict(start=start,success=False,message=str(exc),numerical_stop=True))
            print(f'REJECTED numerical profile in start {start}: {exc}',flush=True)
        save_json(out/'checkpoint.json',dict(completed=False,profiles=profiles,optimizers=stops))
    best,_=select_profile_result(profiles)
    matrix,obs=problem(best['decays'])
    finalize_profile_run(out,profiles,stops,failed_profiles,obs,
        elapsed_seconds=time.perf_counter()-started,source_sha256=source.reference.sha256)


if __name__=='__main__': main()
