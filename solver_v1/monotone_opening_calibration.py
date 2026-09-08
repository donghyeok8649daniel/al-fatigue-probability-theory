"""Research inequality test after observed unphysical-to-reference overshoot.

The fit targets and their weights are UNCHANGED. Require nonnegative opening
traction on a declared discrete path, then independently refine that path.
This is a reference-compatibility constraint, not a universal theorem that no
material surface can ever have an oscillating disjoining pressure.
"""
from functools import lru_cache
import json
import time
import numpy as np
from scipy.optimize import minimize

from .angular_environment_reference import AngularInterfaceInvariant
from .fcc111_active_interface import FCC111ActiveInterface
from .joint_fcc_interface_calibration import IDEAL_H, JointFit, build_joint_model
from .even_moment_calibration import even_basis
from .odd_moment_calibration import odd_targets
from .run_constrained_odd_calibration import constrained_matrix_fit


OPENING_RATIOS=(1.1,1.2,1.5,2.,2.5,3.,4.,5.,8.)


def _probe(decay):
    bulk=build_joint_model(JointFit("basis_only",float(decay),np.array([.1,.15,1.,0.,0.]),0.,None,None))
    bulk.a0=IDEAL_H
    return FCC111ActiveInterface(bulk,tolerance=2e-11)


@lru_cache(maxsize=256)
def scalar_force_rows(decay,ratios=OPENING_RATIOS):
    face=_probe(decay); rows=[]
    for ratio in ratios:
        a=ratio*IDEAL_H
        repulsive=face._power_cross_difference(a,0.,6.)[0][1]
        attractive=-face._power_cross_difference(a,0.,3.)[0][1]
        depths=np.array(face._density_changes(a,0.)[0])/face.rho_bulk
        dx,xa=depths[:,0],depths[:,1]
        rows.append([repulsive,attractive,-np.sum(xa/np.sqrt(1+dx)),2*np.sum(xa),4*np.sum(dx*xa)])
    return np.asarray(rows)


@lru_cache(maxsize=256)
def angular_force_rows(decay,ratios=OPENING_RATIOS):
    face=_probe(decay)
    invariants=[AngularInterfaceInvariant(face,rank=rank,tolerance=2e-11) for rank in (3,1,2)]
    return np.array([[invariant.evaluate(ratio*IDEAL_H,0.)[0][1] for invariant in invariants]
                     for ratio in ratios])


def opening_force_basis(decay,angular_decay,ratios=OPENING_RATIOS):
    return np.column_stack((scalar_force_rows(float(decay),tuple(ratios)),
                            angular_force_rows(float(angular_decay),tuple(ratios))))


def constrained_opening_fit(decay,angular_decay,target,scales,*,ratios=OPENING_RATIOS):
    result=constrained_matrix_fit(even_basis(decay,angular_decay,convex=True),target,scales,
        coefficient_order=["u","v","A","B","C","D3","D1","D2"],
        lower_rest=[0.,-np.inf,0.,0.,-np.inf,-np.inf],
        nonnegative_observation_rows=opening_force_basis(decay,angular_decay,tuple(ratios)))
    return dict(scalar_decay=float(decay),angular_decay=float(angular_decay),
                opening_inequality_a_over_h=tuple(ratios),**result)


def main():
    from .run_matched_interface_study import ROOT,save_json
    started=time.perf_counter(); target,scales=odd_targets()
    # Systematic predefined grid, not random starts or fatigue optimization.
    rows=[constrained_opening_fit(d,k,target,scales)
          for d in np.geomspace(.7,8.,9) for k in np.geomspace(2.,12.,7)]
    logs=[]
    eligible=[r for r in rows if r["admissible"]]
    if not eligible:
        raise ValueError("no feasible opening-constrained coefficient study")
    start=min(eligible,key=lambda r:r["squared_loss"])
    def objective(x):
        trial=constrained_opening_fit(*x,target,scales); rows.append(trial)
        return trial["squared_loss"]
    solved=minimize(objective,[start["scalar_decay"],start["angular_decay"]],method="Powell",
        bounds=((.7,8.),(2.,12.)),options={"maxfev":100,"xtol":1e-4,"ftol":1e-7})
    logs.append(dict(start=[start["scalar_decay"],start["angular_decay"]],success=bool(solved.success),
                     nfev=int(solved.nfev),message=str(solved.message)))
    eligible=[r for r in rows if r["admissible"]]
    best=min(eligible,key=lambda r:r["squared_loss"])
    save_json(ROOT/"monotone_opening_fit.json",dict(best=best,profile=rows,local_optimization=logs,
        target=target,scales=scales,elapsed_seconds=time.perf_counter()-started,
        status="sampled opening inequalities only; dense-path/material stability not certified"))
    print(best,flush=True)


if __name__=="__main__":
    main()
