"""Separate coefficient, stability and radial-search errors on the SAME loss.

Unconstrained coefficient projections are lower bounds at fixed radial shapes,
not admissible potentials. Held-out data never enter any fitted projection.
"""
import argparse
from dataclasses import replace
import json
from pathlib import Path

import numpy as np
from scipy.linalg import null_space

from .interface_development_targets import development_observations
from .quartic_angular_material import QuarticSymmetryObservationCache
from .run_vector_material_calibration import source_and_targets
from .run_vector_registry_audit import save_json
from .run_low_stress_cyclic_diagnostic import write_csv
from .tail_constrained_material import convex_profile
from .yield_elastic_metric import cubic_metric_problem, cubic_to_mode_matrix


def coefficient_diagnostic(matrix, observations):
    y=np.array([o.target for o in observations]);s=np.array([o.scale for o in observations])
    exact=[i for i,o in enumerate(observations) if o.role=='exact']
    fit=[i for i,o in enumerate(observations) if o.role=='fit']
    E=matrix[exact]/s[exact,None];rhs=y[exact]/s[exact]
    offset=np.linalg.lstsq(E,rhs,rcond=None)[0];Z=null_space(E)
    design=matrix[fit]@Z/s[fit,None]
    free=np.linalg.lstsq(design,(y[fit]-matrix[fit]@offset)/s[fit],rcond=None)[0]
    c=offset+Z@free
    residual=(matrix@c-y)/s
    constrained=convex_profile(matrix,observations,nonnegative=(0,1,2,4,5,8,9))
    assert np.max(abs(E@c-rhs))<1e-8
    assert residual[fit]@residual[fit]<=constrained['squared_loss']+1e-7
    return dict(unrestricted_fixed_shape_loss=float(residual[fit]@residual[fit]),
        unrestricted_coefficients=c,unrestricted_residuals=residual,
        unrestricted_is_accepted=False,
        sign_constrained_loss=constrained['squared_loss'],
        sign_constrained_coefficients=constrained['coefficients'],
        sign_constrained_residuals=constrained['residuals'],
        sign_constrained_kkt=constrained['kkt_residual'],
        exact_rank=int(np.linalg.matrix_rank(E)),
        free_coefficient_singular_values=np.linalg.svd(design,compute_uv=False),
        interpretation='fixed-shape least-squares floor only; neither radial global bound nor material validation')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--starts',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--continue-on-numerical-failure',action='store_true')
    args=p.parse_args()
    if args.out.exists():raise FileExistsError('preserve previous diagnostics')
    source,old,states=source_and_targets()
    observations,_=development_observations(source,old,states)
    cache=QuarticSymmetryObservationCache(observations)
    cases=[];residual_rows=[]
    for index,decays in enumerate(json.loads(args.starts.read_bytes())):
        raw=cache.matrix(decays)
        matrix,obs=cubic_metric_problem(raw[:,:8],observations)
        extra=raw[:,8:].copy();extra[2:5]=np.linalg.solve(cubic_to_mode_matrix(),extra[2:5])
        matrix=np.column_stack([matrix,extra])
        obs=[replace(o,role='exact') if i<5 else o for i,o in enumerate(obs)]
        try:
            result=coefficient_diagnostic(matrix,obs)
        except (ArithmeticError,ValueError) as exc:
            if not args.continue_on_numerical_failure:raise
            cases.append(dict(decays=decays,numerical_failure=str(exc),accepted=False))
            save_json(args.out/'checkpoint.json',dict(completed=False,cases=cases))
            print(index,'NUMERICAL FAILURE',str(exc),flush=True)
            continue
        cases.append(dict(decays=decays,**result))
        for j,o in enumerate(obs):
            residual_rows.append(dict(case=index,observable=o.name,role=o.role,scale=o.scale,
                unrestricted=result['unrestricted_residuals'][j],
                sign_constrained=result['sign_constrained_residuals'][j]))
        print(index,result['unrestricted_fixed_shape_loss'],result['sign_constrained_loss'],flush=True)
        save_json(args.out/'checkpoint.json',dict(completed=False,cases=cases))
    save_json(args.out/'coefficient_diagnosis.json',dict(completed=True,source_sha256=source.reference.sha256,cases=cases))
    write_csv(args.out/'residual_comparison.csv',residual_rows)


if __name__=='__main__':main()
