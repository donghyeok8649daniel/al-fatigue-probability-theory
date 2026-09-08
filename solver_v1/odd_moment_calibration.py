"""Scenario-driven minimal I1+I3 test after I3 over-stiffens the interface.

The new normal-interface curvature is an atomistic target, not a desired
fatigue curve. The scalar LJ/EAM and the original PDE remain unchanged.
Signed I1 energy requires TOTAL stability checks; it is not assumed physical.
"""
from functools import lru_cache
import math
import numpy as np
from scipy.optimize import lsq_linear, minimize

from .angular_environment_reference import AngularInterfaceInvariant
from .angular_interface_calibration import angular_observation_column
from .fcc111_active_interface import FCC111ActiveInterface
from .fcc111_geometry import DIRECT_110
from .joint_fcc_interface_calibration import (
    IDEAL_H, INTERFACE_STATES, JointFit, build_joint_model, joint_basis, matched_targets,
)
from .reference_eam_targets import MishinRigidFCCReference


@lru_cache(maxsize=256)
def scalar_curvature_row(decay):
    bulk=build_joint_model(JointFit("probe",float(decay),np.array([.1,.15,1.,0.,0.]),0.,None,None))
    bulk.a0=IDEAL_H
    face=FCC111ActiveInterface(bulk,tolerance=2e-11)
    s6=face._power_cross_difference(IDEAL_H,0.,6.)[0][3]
    ms3=-face._power_cross_difference(IDEAL_H,0.,3.)[0][3]
    depths=face._density_changes(IDEAL_H,0.)[0]
    first=np.array([r[1] for r in depths])/face.rho_bulk
    second=np.array([r[3] for r in depths])/face.rho_bulk
    return np.array([s6,ms3,2*np.sum(.25*first**2-.5*second),2*np.sum(second)])


@lru_cache(maxsize=256)
def odd_columns(decay):
    bulk=build_joint_model(JointFit("probe",float(decay),np.array([.1,.15,1.,0.,0.]),0.,None,None))
    bulk.a0=IDEAL_H; vector_values=[]
    for _,opening,slip,path in INTERFACE_STATES:
        face=FCC111ActiveInterface(bulk,path_id=path,tolerance=2e-11)
        vector=AngularInterfaceInvariant(face,tolerance=2e-11,rank=1)
        a=400*IDEAL_H if math.isinf(opening) else opening*IDEAL_H
        vector_values.append(vector.evaluate(a,slip)[0][0])
    face=FCC111ActiveInterface(bulk,path_id=DIRECT_110,tolerance=2e-11)
    curvature1=AngularInterfaceInvariant(face,tolerance=2e-11,rank=1).evaluate(IDEAL_H,0.)[0][3]
    curvature3=AngularInterfaceInvariant(face,tolerance=2e-11).evaluate(IDEAL_H,0.)[0][3]
    return np.column_stack((np.r_[angular_observation_column(float(decay)),curvature3],
                            np.r_[np.zeros(5),vector_values,curvature1]))


def odd_targets(reference=None):
    reference=reference or MishinRigidFCCReference()
    target,scale=matched_targets(reference)
    length_angstrom=4.05/math.sqrt(2)
    curvature=reference.interface_derivatives(reference.h,0.)[3]*length_angstrom**2
    # 10% local-interface curvature allowance, separate from bulk elasticity.
    return np.r_[target,curvature],np.r_[scale,.1*curvature]


def odd_basis(scalar_decay,angular_decay):
    scalar=np.vstack((joint_basis(float(scalar_decay))[:,:4],scalar_curvature_row(float(scalar_decay))))
    return np.column_stack((scalar,odd_columns(float(angular_decay))))


def odd_fit(scalar_decay,angular_decay,target,scales):
    matrix=odd_basis(scalar_decay,angular_decay)
    scaled=matrix/scales[:,None]
    norm=np.linalg.norm(scaled,axis=0)
    result=lsq_linear(scaled/norm,target/scales,
        bounds=([0.,0.,0.,0.,0.,-np.inf],np.full(6,np.inf)),tol=2e-10,max_iter=200)
    coeff=result.x/norm
    prediction=matrix@coeff; residual=(prediction-target)/scales
    return dict(scalar_decay=float(scalar_decay),angular_decay=float(angular_decay),
        coefficients=coeff,coefficient_order=["u","v","A","B","D3","D1"],
        predictions=prediction,residuals=residual,squared_loss=float(residual@residual),
        linear_solve_success=bool(result.success))


def odd_profile(target,scales,*,max_nfev=150):
    scalar=np.geomspace(.7,8.,13); angular=np.geomspace(2.,12.,9)
    rows=[odd_fit(d,k,target,scales) for d in scalar for k in angular]
    logs=[]
    for start in sorted(rows,key=lambda r:r["squared_loss"])[:2]:
        def objective(x):
            trial=odd_fit(x[0],x[1],target,scales); rows.append(trial)
            return trial["squared_loss"]
        result=minimize(objective,[start["scalar_decay"],start["angular_decay"]],method="Powell",
            bounds=((scalar[0],scalar[-1]),(angular[0],angular[-1])),
            options={"maxfev":max_nfev,"maxiter":30,"xtol":1e-4,"ftol":1e-7})
        logs.append(dict(start=[start["scalar_decay"],start["angular_decay"]],
            success=bool(result.success),message=str(result.message),nfev=int(result.nfev),
            squared_loss=float(result.fun)))
    return min(rows,key=lambda r:r["squared_loss"]),rows,logs


def main():
    import time
    from .run_matched_interface_study import ROOT,save_json
    started=time.perf_counter(); target,scales=odd_targets()
    best,profile,logs=odd_profile(target,scales)
    save_json(ROOT/"odd_moment_fit.json",dict(best=best,profile=profile,local_optimization=logs,
        target=target,scales=scales,elapsed_seconds=time.perf_counter()-started,
        status="total stability and held-out stress/curves required; no production adoption"))
    print(best,flush=True)


if __name__=="__main__":
    main()
