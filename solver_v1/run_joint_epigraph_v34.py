"""Deterministic existing-family constrained calibration; research only."""
import argparse
import json
from pathlib import Path
from functools import lru_cache
import time
import numpy as np
from scipy.optimize import minimize
from .run_profile_shape_v34 import ShapeProfile
from .profile_shape_sensitivity import difference_stencil
from .joint_shape_epigraph import constraints_and_jacobian
from .interface_tangent_calibration import NONNEGATIVE
from .run_vector_registry_audit import save_json


def run(parent, out, maxiter, resume=None, angular_upper=12.):
    if maxiter < 1: raise ValueError('positive iteration budget required')
    study=ShapeProfile(parent,out)
    if not np.isfinite(angular_upper) or angular_upper < 12:
        raise ValueError('angular upper must be at least original 12')
    study.upper[[1,2,4]]=np.log(angular_upper)
    if resume is not None:
        previous=json.loads(Path(resume).read_bytes())
        if not previous['completed']:raise ValueError('completed study required')
        study.start=np.array(previous['final_shape_coordinates'])
    save_json(study.out/'joint_definition.json',dict(resume=str(resume) if resume else None,
        actual_start=study.start,lower=study.lower,upper=study.upper,
        angular_upper=angular_upper,interpretation='deterministic search-domain sensitivity, not material calibration',
        target_scales_changed=False,new_energy_terms=False))
    p=study.profile(tuple(study.start))
    c=np.asarray(p['coefficients'])
    D=np.maximum(abs(c),np.asarray(p['coefficient_scales']))
    targets=np.array([o.target for o in study.obs])
    scales=np.array([o.scale for o in study.obs])
    exact=p['exact_rows']; rows=p['tested_rows']; nx=len(study.start)
    initial=np.r_[study.start,c/D,p['minimax_normalized_error']]
    records=[]

    @lru_cache(maxsize=4)
    def derivative(x):
        return np.array([sum(w*study.matrix(tuple(v)) for v,w in
            difference_stencil(np.array(x),i,study.lower,study.upper,1e-4)) for i in range(nx)])

    def evaluate(z, jac=False):
        x=tuple(z[:nx]); M=study.matrix(x)
        dm=derivative(x) if jac else np.zeros((nx,*M.shape))
        return constraints_and_jacobian(M,dm,targets,scales,exact,rows,z[nx:-1]*D,D,z[-1])

    def checkpoint(z):
        eq,ineq,*_=evaluate(z)
        record=dict(z=z,eta=float(z[-1]),exact_residual=float(np.max(abs(eq))),
                    inequality_violation=float(max(0,-min(ineq))),
                    elapsed_seconds=time.perf_counter()-study.began)
        records.append(record)
        save_json(study.out/'joint_checkpoint.json',dict(completed=False,iterations=records))
        print('joint',len(records),'eta',z[-1],'eq',record['exact_residual'],flush=True)

    checkpoint(initial)
    objective_jac=np.r_[np.zeros(len(initial)-1),1.]
    bounds=list(zip(study.lower,study.upper))+[(0,None) if i in NONNEGATIVE else (None,None)
        for i in range(len(c))]+[(0,None)]
    result=minimize(lambda z:z[-1],initial,jac=lambda z:objective_jac,
        method='SLSQP',bounds=bounds,
        constraints=[dict(type='eq',fun=lambda z:evaluate(z)[0],jac=lambda z:evaluate(z,True)[2]),
                     dict(type='ineq',fun=lambda z:evaluate(z)[1],jac=lambda z:evaluate(z,True)[3])],
        callback=checkpoint,options=dict(maxiter=maxiter,ftol=1e-8))
    # LP independently restores exact anchors and certifies coefficient optimality
    # at the proposed shape. This is not a global shape convergence certificate.
    final=study.profile(tuple(result.x[:nx]))
    save_json(study.out/'joint_summary.json',dict(completed=True,optimizer_success=bool(result.success),
        message=str(result.message),iterations=int(result.nit),initial_eta=float(initial[-1]),
        final_profile=final,final_shape_coordinates=result.x[:nx],joint_iterations=records,
        coefficient_scales=D,elapsed_seconds=time.perf_counter()-study.began,
        material_accepted=False,physical_time_calibrated=False,production_changed=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(__doc__)
    p.add_argument('--parent',required=True);p.add_argument('--out',required=True)
    p.add_argument('--maxiter',type=int,default=6)
    p.add_argument('--resume');p.add_argument('--angular-upper',type=float,default=12.)
    a=p.parse_args();run(a.parent,a.out,a.maxiter,a.resume,a.angular_upper)
