"""Reproduce a declared finite-difference probe that failed spectral profiling."""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import time

import numpy as np

from .interface_development_targets import development_observations
from .quartic_angular_material import QuarticSymmetryObservationCache,QuarticSymmetryTailBulkBasis
from .run_vector_material_calibration import source_and_targets
from .run_vector_registry_audit import save_json
from .tail_constrained_material import spectral_profile,declared_wavepoints
from .yield_elastic_metric import cubic_metric_problem,cubic_to_mode_matrix


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('calibration',type=Path)
    p.add_argument('--profile-index',type=int,required=True)
    p.add_argument('--shape-index',type=int,required=True)
    p.add_argument('--log-relative-step',type=float,default=2e-4)
    p.add_argument('--max-cuts',type=int,default=32)
    p.add_argument('--out',type=Path,required=True)
    args=p.parse_args()
    if args.out.exists():raise FileExistsError('preserve previous diagnostic')
    began=time.perf_counter()
    data=json.loads(args.calibration.read_bytes())
    decays=np.array(data['profiles'][args.profile_index]['decays'])
    base=decays.copy();decays[args.shape_index]=np.exp(np.log(decays[args.shape_index])*(1+args.log_relative_step))
    source,obs,states=source_and_targets();obs,_=development_observations(source,obs,states)
    cache=QuarticSymmetryObservationCache(obs);raw=cache.matrix(decays)
    M,obs=cubic_metric_problem(raw[:,:8],obs)
    extra=raw[:,8:].copy();extra[2:5]=np.linalg.solve(cubic_to_mode_matrix(),extra[2:5])
    M=np.column_stack([M,extra]);obs=[replace(o,role='exact') if i<5 else o for i,o in enumerate(obs)]
    basis=QuarticSymmetryTailBulkBasis(decays,radius=12.)
    values=[basis.evaluate(q) for q in declared_wavepoints()]
    columns,tails=map(np.asarray,zip(*values))
    result=dict(source_sha256=source.reference.sha256,profile_index=args.profile_index,
        base_decays=base,shape_index=args.shape_index,probe_decays=decays,
        log_relative_step=args.log_relative_step,max_cuts=args.max_cuts,
        no_material_or_tolerance_change=True,accepted_as_material=False)
    try:
        fit=spectral_profile(M,obs,columns,tails,max_cuts=args.max_cuts,nonnegative=(0,1,2,4,5,8,9))
        result.update(numerically_verified=True,fit=fit)
    except (ArithmeticError,ValueError) as error:
        result.update(numerically_verified=False,numerical_failure=str(error))
    save_json(args.out/'precision_audit.json',dict(result,completed=True,elapsed_seconds=time.perf_counter()-began))
    print(result['numerically_verified'],result.get('numerical_failure'),flush=True)


if __name__=='__main__':main()
