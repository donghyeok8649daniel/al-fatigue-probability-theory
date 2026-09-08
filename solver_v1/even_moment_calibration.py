"""Nested I2 test: independent bulk anisotropy without changing LJ or PDE.

The derived quadrupole has one shared radial range with I1/I3. Unlike the odd
moments it responds to affine anisotropic strain, but not isotropic expansion.
All coefficients remain energy coefficients of per-atom environmental sums.
"""
from functools import lru_cache
import math
import time
import numpy as np
from scipy.optimize import minimize

from .angular_environment_reference import AngularInterfaceInvariant
from .fcc111_active_interface import FCC111ActiveInterface
from .joint_fcc_interface_calibration import IDEAL_H, INTERFACE_STATES, JointFit, build_joint_model
from .odd_moment_calibration import odd_basis, odd_targets
from .run_odd_moment_quadratic_probe import scalar_basis
from .run_constrained_odd_calibration import constrained_matrix_fit


def quadrupole_bulk_curvatures(invariant):
    """Hydro/normal/shear Hessian per unit D2 at isotropically scaled FCC.

    Same-plane atoms do not move under these normal/simple-shear modes.
    Q_bulk=0 by cubic symmetry; its first variations are infinite plane sums.
    No finite differences or finite neighbor cutoff in this implementation.
    """
    h=invariant.h
    derivatives=np.zeros((2,3,3)); small=0
    for layer in range(1,invariant.max_layers+1):
        v,_,_=invariant._plane(layer*h,invariant.interface._baseline_delta(layer))
        term=2*layer*h*v[1:3]
        derivatives+=term
        if np.max(np.abs(term))<invariant.tolerance:
            small+=1
            if small>=4:
                return np.r_[0.,2*np.sum(derivatives**2,axis=(1,2))]
        else:
            small=0
    raise RuntimeError("bulk quadrupole derivative layer sum failed convergence")


@lru_cache(maxsize=256)
def even_column(decay):
    bulk=build_joint_model(JointFit("probe_only",float(decay),np.array([.1,.15,1.,0.,0.]),0.,None,None))
    bulk.a0=IDEAL_H
    face=FCC111ActiveInterface(bulk,tolerance=2e-11)
    invariant=AngularInterfaceInvariant(face,tolerance=2e-11,rank=2)
    column=[0.,0.,*quadrupole_bulk_curvatures(invariant)]
    for _,opening,slip,path in INTERFACE_STATES:
        face=FCC111ActiveInterface(bulk,path_id=path,tolerance=2e-11)
        invariant=AngularInterfaceInvariant(face,tolerance=2e-11,rank=2)
        a=400*IDEAL_H if math.isinf(opening) else opening*IDEAL_H
        column.append(invariant.evaluate(a,slip)[0][0])
    face=FCC111ActiveInterface(bulk,tolerance=2e-11)
    column.append(AngularInterfaceInvariant(face,tolerance=2e-11,rank=2).evaluate(IDEAL_H,0.)[0][3])
    return np.asarray(column)


def even_basis(decay,angular_decay,*,convex=False):
    odd=odd_basis(float(decay),float(angular_decay))
    if convex:
        odd=np.column_stack((scalar_basis(float(decay)),odd[:,4:]))
    return np.column_stack((odd,even_column(float(angular_decay))))


def even_fit(decay,angular_decay,target,scales,*,convex=False,signed=False,free_linear=False):
    names=["u","v","A","B"]+(["C"] if convex else [])+["D3","D1","D2"]
    lower=[0.,-np.inf if free_linear else 0.]+([0.] if convex else [])+[0.,-np.inf,-np.inf if signed else 0.]
    result=constrained_matrix_fit(even_basis(decay,angular_decay,convex=convex),target,scales,
                                 coefficient_order=names,lower_rest=lower)
    return dict(scalar_decay=float(decay),angular_decay=float(angular_decay),
                convex_embedding=convex,signed_quadrupole=signed,signed_linear_embedding=free_linear,**result)


def profile(target,scales,*,convex=False,signed=False,free_linear=False,max_nfev=100):
    scalar=np.geomspace(.7,8.,11); angular=np.geomspace(2.,12.,8)
    rows=[even_fit(d,k,target,scales,convex=convex,signed=signed,free_linear=free_linear) for d in scalar for k in angular]
    eligible=[r for r in rows if r["admissible"]]
    logs=[]
    for start in sorted(eligible,key=lambda r:r["squared_loss"])[:2]:
        def objective(x):
            trial=even_fit(*x,target,scales,convex=convex,signed=signed,free_linear=free_linear); rows.append(trial)
            return trial["squared_loss"]
        result=minimize(objective,[start["scalar_decay"],start["angular_decay"]],method="Powell",
            bounds=((.7,8.),(2.,12.)),options={"maxfev":max_nfev,"xtol":1e-4,"ftol":1e-7})
        logs.append(dict(start=[start["scalar_decay"],start["angular_decay"]],success=bool(result.success),
                         nfev=int(result.nfev),message=str(result.message),squared_loss=float(result.fun)))
    eligible=[r for r in rows if r["admissible"]]
    return min(eligible,key=lambda r:r["squared_loss"]),rows,logs


def main():
    import argparse
    from .run_matched_interface_study import ROOT,save_json
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--convex",action="store_true")
    parser.add_argument("--signed",action="store_true")
    parser.add_argument("--free-linear",action="store_true",help="permit signed B; separate coefficient-sector audit")
    args=parser.parse_args()
    started=time.perf_counter(); target,scales=odd_targets()
    best,rows,logs=profile(target,scales,convex=args.convex,signed=args.signed,free_linear=args.free_linear)
    name="even_moment"+("_convex" if args.convex else "")+("_signed" if args.signed else "")+("_free_B" if args.free_linear else "")
    save_json(ROOT/(name+"_fit.json"),dict(best=best,profile=rows,local_optimization=logs,
        target=target,scales=scales,elapsed_seconds=time.perf_counter()-started,
        status="exact reference geometry/cohesion; candidate not accepted without physical validation"))
    print(best,flush=True)


if __name__=="__main__":
    main()
