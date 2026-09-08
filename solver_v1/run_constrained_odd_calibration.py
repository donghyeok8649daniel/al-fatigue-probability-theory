"""Hierarchical calibration with EXACT target equilibrium and cohesion.

An approximate-force loss permits residual pre-stress at the target lattice.
This separate study first satisfies the two bulk reference constraints exactly,
then optimizes elasticity/interface observations. No data weight is retuned to
obtain fatigue or slip; the eliminated pair coefficients must remain positive.
"""
import time
import numpy as np
from scipy.optimize import Bounds, LinearConstraint, minimize

from .odd_moment_calibration import odd_columns, odd_targets
from .run_odd_moment_quadratic_probe import scalar_basis
from .run_matched_interface_study import ROOT, save_json


def constrained_fit(decay,angular_decay,target,scales):
    matrix=np.column_stack((scalar_basis(float(decay)),odd_columns(float(angular_decay))))
    values=constrained_matrix_fit(matrix,target,scales,
        coefficient_order=["u","v","A","B","C","D3","D1"],
        lower_rest=[0.,0.,0.,0.,-np.inf])
    return dict(scalar_decay=float(decay),angular_decay=float(angular_decay),**values)


def constrained_matrix_fit(matrix,target,scales,*,coefficient_order,lower_rest,nonnegative_observation_rows=None):
    """Exact equilibrium/cohesion elimination; the remaining fit is convex QP.

    Matrix columns, their physical meaning, and lower bounds are explicit.
    Mathematical feasibility never certifies a physical material candidate.
    """
    # c=[u,v,A,B,C,D3,D1]. Eliminate u,v by exact linear constraints:
    # normal strain force=0 and atomization energy=3.36 eV/atom.
    block=matrix[:2,:2]
    offset=np.linalg.solve(block,target[:2])
    mapping=-np.linalg.solve(block,matrix[:2,2:])
    constant=matrix[2:,:2]@offset
    reduced=matrix[2:,2:]+matrix[2:,:2]@mapping
    design=reduced/scales[2:,None]; rhs=(target[2:]-constant)/scales[2:]
    norm=np.linalg.norm(design,axis=0)
    scaled=design/norm
    objective=lambda z:float(np.sum((scaled@z-rhs)**2))
    jac=lambda z:2*scaled.T@(scaled@z-rhs)
    # Pure LJ constrained solution is a deterministic feasible starting point.
    if np.any(offset<=0):
        raise ValueError("pure-LJ equilibrium/cohesion constrained start not positive")
    constraints=[LinearConstraint(mapping/norm,-offset,np.full(2,np.inf))]
    extra=None
    if nonnegative_observation_rows is not None:
        extra=np.asarray(nonnegative_observation_rows,dtype=float)
        row_norm=np.linalg.norm(extra,axis=1)
        if np.any(row_norm==0):
            raise ValueError("inequality rows must be nonzero")
        reduced_extra=(extra[:,2:]+extra[:,:2]@mapping)/norm/row_norm[:,None]
        lower_extra=-(extra[:,:2]@offset)/row_norm
        constraints.append(LinearConstraint(reduced_extra,lower_extra,np.full(len(extra),np.inf)))
    solved=minimize(objective,np.zeros(len(lower_rest)),jac=jac,method="SLSQP",
        bounds=Bounds(lower_rest,np.full(len(lower_rest),np.inf)),
        constraints=constraints,
        options={"maxiter":180,"ftol":1e-10})
    rest=solved.x/norm; pair=offset+mapping@rest
    coefficients=np.r_[pair,rest]
    prediction=matrix@coefficients; residual=(prediction-target)/scales
    # Positivity of a 1e-16 residual after coefficient elimination is not a
    # resolved nonzero attractive LJ term. Propagate cancellation roundoff;
    # retain raw coefficients and distinguish feasibility from physical use.
    eps=np.finfo(float).eps
    roundoff=(32*eps/(1-32*eps))*np.linalg.cond(block)*(abs(offset)+abs(mapping)@abs(rest))
    resolved_positive=bool(np.all(pair>roundoff))
    valid=bool(solved.success and np.all(pair>=-roundoff) and np.max(abs(residual[:2]))<1e-8)
    extra_values=extra@coefficients if extra is not None else np.array([])
    if extra is not None:
        valid=valid and bool(np.min(extra_values/row_norm)>=-1e-9)
    return dict(coefficients=coefficients,
        coefficient_order=list(coefficient_order),predictions=prediction,
        residuals=residual,squared_loss=float(residual@residual),success=bool(solved.success),
        admissible=valid,inner_message=str(solved.message),inner_nit=int(solved.nit),
        exact_constraints=["normal_force_at_target=0","cohesion=3.36 eV/atom"],
        pair_elimination_roundoff=roundoff,
        strictly_positive_LJ_resolved=resolved_positive,
        extra_inequality_observations=extra_values,
        physical_candidate_accepted=False,
        total_stability_checked=False)


def main():
    started=time.perf_counter(); target,scales=odd_targets()
    rows=[constrained_fit(d,k,target,scales) for d in np.geomspace(.7,8.,13) for k in np.geomspace(2.,12.,9)]
    eligible=[row for row in rows if row["admissible"]]
    if not eligible:
        raise ValueError("no admissible coefficient fit; no artificial penalty/clipping")
    logs=[]
    for start in sorted(eligible,key=lambda row:row["squared_loss"])[:2]:
        def objective(x):
            row=constrained_fit(x[0],x[1],target,scales); rows.append(row)
            # A failed inner solve is recorded and is never eligible as a fit.
            # No surrogate physical energy/force is constructed from it.
            return row["squared_loss"]
        result=minimize(objective,[start["scalar_decay"],start["angular_decay"]],method="Powell",
            bounds=((.7,8.),(2.,12.)),options={"maxfev":120,"maxiter":25,"ftol":1e-7,"xtol":1e-4})
        logs.append(dict(start=[start["scalar_decay"],start["angular_decay"]],success=bool(result.success),
            nfev=int(result.nfev),message=str(result.message),squared_loss=float(result.fun)))
    eligible=[row for row in rows if row["admissible"]]
    best=min(eligible,key=lambda row:row["squared_loss"])
    save_json(ROOT/"constrained_odd_fit.json",dict(best=best,profile=rows,local_optimization=logs,
        target=target,scales=scales,elapsed_seconds=time.perf_counter()-started,
        status="exact reference geometry/cohesion; still requires stability and held-out validation"))
    print(best,flush=True)


if __name__=="__main__":
    main()
