"""Bounded nested-family check: the existing convex C term with I1+I3.

No new functional form or production model is added. All ingredients already
have analytic derivatives. This script records failure/success of the nested
coefficient family rather than attributing a failure to a omitted C term.
"""
from functools import lru_cache
import time
import numpy as np
from scipy.optimize import lsq_linear, minimize

from .odd_moment_calibration import odd_columns, odd_targets, scalar_curvature_row
from .joint_fcc_interface_calibration import (
    JointFit, build_joint_model, joint_basis, IDEAL_H,
)
from .fcc111_active_interface import FCC111ActiveInterface
from .run_matched_interface_study import ROOT, save_json


@lru_cache(maxsize=256)
def scalar_basis(decay):
    bulk=build_joint_model(JointFit("probe",float(decay),np.array([.1,.15,1.,0.,0.]),0.,None,None))
    bulk.a0=IDEAL_H; face=FCC111ActiveInterface(bulk,tolerance=2e-11)
    changes=face._density_changes(IDEAL_H,0)[0]
    normal=4*sum((r[1]/face.rho_bulk)**2 for r in changes)
    return np.vstack((joint_basis(float(decay)),np.r_[scalar_curvature_row(float(decay)),normal]))


def fit(decay,angular_decay,target,scales):
    matrix=np.column_stack((scalar_basis(float(decay)),odd_columns(float(angular_decay))))
    normalized=matrix/scales[:,None]; norm=np.linalg.norm(normalized,axis=0)
    result=lsq_linear(normalized/norm,target/scales,
        bounds=([0.]*6+[-np.inf],np.full(7,np.inf)),tol=2e-10,max_iter=250)
    coefficients=result.x/norm; prediction=matrix@coefficients
    residual=(prediction-target)/scales
    return dict(scalar_decay=float(decay),angular_decay=float(angular_decay),coefficients=coefficients,
                coefficient_order=["u","v","A","B","C","D3","D1"],predictions=prediction,
                residuals=residual,squared_loss=float(residual@residual),success=bool(result.success))


def main():
    started=time.perf_counter(); target,scales=odd_targets()
    rows=[fit(d,k,target,scales) for d in np.geomspace(.7,8.,13) for k in np.geomspace(2.,12.,9)]
    logs=[]
    for start in sorted(rows,key=lambda row:row["squared_loss"])[:2]:
        def objective(x):
            row=fit(x[0],x[1],target,scales); rows.append(row)
            return row["squared_loss"]
        optimum=minimize(objective,[start["scalar_decay"],start["angular_decay"]],
            method="Powell",bounds=((.7,8.),(2.,12.)),
            options={"maxfev":120,"maxiter":25,"xtol":1e-4,"ftol":1e-7})
        logs.append(dict(start=[start["scalar_decay"],start["angular_decay"]],
                         success=bool(optimum.success),nfev=int(optimum.nfev),
                         message=str(optimum.message),squared_loss=float(optimum.fun)))
    best=min(rows,key=lambda row:row["squared_loss"])
    save_json(ROOT/"odd_moment_quadratic_probe.json",dict(best=best,profile=rows,local_optimization=logs,
        target=target,scales=scales,elapsed_seconds=time.perf_counter()-started,
        status="nested-family feasibility only; needs actual bulk/interface/stress validation"))
    print(best,flush=True)


if __name__=="__main__":
    main()
