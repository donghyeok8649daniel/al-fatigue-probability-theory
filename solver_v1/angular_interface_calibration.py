"""Observable-first test of one angular invariant; no production promotion."""
from functools import lru_cache
import math
import numpy as np
from scipy.optimize import minimize, nnls

from .angular_environment_reference import AngularInterfaceInvariant
from .fcc111_active_interface import FCC111ActiveInterface
from .joint_fcc_interface_calibration import (
    IDEAL_H, INTERFACE_STATES, JointFit, build_joint_model, joint_basis,
)


@lru_cache(maxsize=256)
def angular_observation_column(decay):
    bulk=build_joint_model(JointFit("probe",float(decay),np.array([.1,.15,1.,0.,0.]),0.,None,None))
    bulk.a0=IDEAL_H
    values=[]
    for _,opening,slip,path in INTERFACE_STATES:
        face=FCC111ActiveInterface(bulk,path_id=path,tolerance=2e-11)
        moment=AngularInterfaceInvariant(face,tolerance=2e-11)
        a=400*IDEAL_H if math.isinf(opening) else opening*IDEAL_H
        values.append(moment.evaluate(a,slip)[0][0])
    # Odd tensor is exactly zero in affine centrosymmetric bulk and isolated
    # atoms. It contributes neither cohesion nor homogeneous elastic constants.
    return np.r_[np.zeros(5),values]


def angular_fit(scalar_decay,angular_decay,target,scales,*,quadratic=False):
    matrix=np.column_stack((joint_basis(float(scalar_decay))[:,:5 if quadratic else 4],
                            angular_observation_column(float(angular_decay))))
    coeff,_=nnls(matrix/scales[:,None],target/scales,maxiter=500)
    prediction=matrix@coeff; residual=(prediction-target)/scales
    return dict(scalar_decay=float(scalar_decay),angular_decay=float(angular_decay),
                coefficients=coeff,predictions=prediction,residuals=residual,
                squared_loss=float(residual@residual),quadratic=bool(quadratic))


def angular_profile(target,scales,*,tied=True,scalar_grid=None,angular_grid=None,quadratic=False,max_nfev=None):
    scalar_grid=np.geomspace(.5,8.,17) if scalar_grid is None else scalar_grid
    angular_grid=np.geomspace(2.,16.,9) if angular_grid is None else angular_grid
    if tied:
        rows=[angular_fit(d,d,target,scales,quadratic=quadratic) for d in scalar_grid]
    else:
        rows=[angular_fit(d,k,target,scales,quadratic=quadratic) for d in scalar_grid for k in angular_grid]
    logs=[]
    for row in sorted(rows,key=lambda row:row["squared_loss"])[:2]:
        start=[row["scalar_decay"]] if tied else [row["scalar_decay"],row["angular_decay"]]
        bounds=[(min(scalar_grid),max(scalar_grid))]
        if not tied:
            bounds.append((min(angular_grid),max(angular_grid)))
        def evaluate(x):
            trial=angular_fit(x[0],x[0] if tied else x[1],target,scales,quadratic=quadratic)
            rows.append(trial)
            return trial
        options={"xtol":1e-4,"ftol":1e-7,"maxiter":40}
        if max_nfev is not None:
            options["maxfev"]=int(max_nfev)
        result=minimize(lambda x:evaluate(x)["squared_loss"],start,method="Powell",
                        bounds=bounds,options=options)
        rows.append(evaluate(result.x))
        logs.append(dict(start=start,success=bool(result.success),nfev=int(result.nfev),
                         message=str(result.message),squared_loss=float(result.fun)))
    return min(rows,key=lambda row:row["squared_loss"]),rows,logs
