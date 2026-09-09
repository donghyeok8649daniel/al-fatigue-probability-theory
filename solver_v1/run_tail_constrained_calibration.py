"""Actual same-family deterministic calibration with analytic spectral tails."""
import argparse
from dataclasses import asdict, replace
from functools import lru_cache
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import least_squares

from .range_resolved_material import RangeObservationCache
from .run_vector_material_calibration import source_and_targets
from .run_vector_registry_audit import save_json
from .run_low_stress_cyclic_diagnostic import write_csv
from .tail_constrained_material import TailBulkCoefficientBasis, declared_wavepoints, spectral_profile
from .yield_elastic_metric import cubic_metric_problem, cubic_to_mode_matrix


ROOT=Path(__file__).resolve().parents[1]
DEFAULT_OUT=ROOT/'results/fcc111_active_interface/tail_calibration_v14/baseline'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=DEFAULT_OUT)
    parser.add_argument('--max-nfev',type=int,default=28)
    parser.add_argument('--radius',type=float,default=12.)
    parser.add_argument('--upper-decays',type=float,nargs=3,default=[8.,12.,12.],
                        help='explicit search bounds for a recorded range-continuation test')
    parser.add_argument('--symmetry-channel',action='store_true',help='separate one-amplitude nested Eg test')
    parser.add_argument('--exact-bulk',action='store_true',help='staged bulk-constrained interface calibration, not zero target uncertainty')
    parser.add_argument('--cubic-embedding',action='store_true',help='separate one-coefficient density-shape extension after measured failure')
    parser.add_argument('--quartic-angular',action='store_true',help='one per-site squared angular-norm extension, not density cubic')
    parser.add_argument('--saturating-angular',action='store_true',help='nested nonnegative rational angular shape parameter')
    parser.add_argument('--starts',type=Path,help='explicit recorded deterministic start JSON list')
    args=parser.parse_args(); out=args.out
    if np.any(~np.isfinite(args.upper_decays)) or np.any(np.asarray(args.upper_decays)<=np.array([.7,2.,2.])):
        raise ValueError('finite upper decays greater than lower bounds required')
    if out.exists(): raise FileExistsError('a fresh result directory is required')
    started=time.perf_counter(); source,observations,states=source_and_targets()
    starts=[[2.579101374767287,5.996534839345317,4.962329189737055],
            [5.3707449112586065,5.6320174463985,3.2556291351974864],
            [3.,8.,7.]]
    if args.starts: starts=json.loads(args.starts.read_bytes())
    if args.saturating_angular:
        if not args.quartic_angular: raise ValueError('rational ablation requires the quartic angular family')
        if not args.starts:
            # Three fixed starts around the already completed quartic basin;
            # alpha is a dimensionless moment-shape parameter, never an area.
            starts=[[2.782954,5.646841,11.999,alpha] for alpha in (1.e3,1.e5,1.e7)]
    points=declared_wavepoints(); cache=RangeObservationCache(observations)
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
    def problem(decays):
        raw=cache.matrix(decays)
        matrix,obs=cubic_metric_problem(raw[:,:8],observations)
        if args.symmetry_channel:
            extra=raw[:,8:].copy();extra[2:5]=np.linalg.solve(cubic_to_mode_matrix(),extra[2:5])
            matrix=np.column_stack([matrix,extra])
        if args.exact_bulk:
            obs=[replace(o,role='exact') if i<5 else o for i,o in enumerate(obs)]
        return matrix,obs
    profiles=[]; stops=[]
    save_json(out/'definition.json',dict(source_sha256=source.reference.sha256,
        observations=[asdict(o) for o in observations],source_states=states,starts=starts,
        wavepoints_cubic=points,radius_over_L0=args.radius,
        criterion='min eig(H_R) >= sum |c_j| analytic FCC Voronoi-tail envelope_j',
        whole_zone_proof=False,new_energy_terms=bool(args.symmetry_channel),
        extra_amplitude='D_E>=0' if args.symmetry_channel else None,
        exact_bulk_stage=bool(args.exact_bulk),
        cubic_density_extension=bool(args.cubic_embedding),
        quartic_angular_extension=bool(args.quartic_angular),
        rational_angular_extension=bool(args.saturating_angular),
        lower_decays=[.7,2.,2.],upper_decays=args.upper_decays,
        physical_kinetics_calibrated=False))
    @lru_cache(maxsize=192)
    def evaluate(logs):
        before=time.perf_counter(); decays=np.exp(logs)
        matrix,obs=problem(decays)
        basis=bulk_type(decays[:3],radius=args.radius)
        computed=[basis.evaluate(q) for q in points]
        columns,tails=map(np.asarray,zip(*computed))
        fit=spectral_profile(matrix,obs,columns,tails,nonnegative=nonnegative)
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
            run=least_squares(fun,np.log(start),jac='2-point',diff_step=2e-4,
                bounds=np.log(bounds),max_nfev=args.max_nfev,
                ftol=1e-7,xtol=1e-7,gtol=2e-6)
            stops.append(dict(start=start,success=bool(run.success),message=str(run.message),
                nfev=run.nfev,njev=run.njev,optimality=float(run.optimality),
                final_decays=np.exp(run.x),final_loss=float(run.fun@run.fun)))
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
        optimizers=stops,elapsed_seconds=time.perf_counter()-started,
        source_sha256=source.reference.sha256,material_accepted=False))
    save_json(out/'progress.json',dict(completed=True,profiles=len(profiles),best_loss=best['squared_loss']))
    print('completed actual calibration',best['squared_loss'],flush=True)


if __name__=='__main__': main()
