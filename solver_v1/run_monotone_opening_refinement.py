"""Refine opening inequalities at ACTUAL force extrema, without clipping.

Fixed radial parameters; deterministic coefficient QP constraint exchange.
The coarse inequality grid can miss a narrow negative-traction interval. Add
the actual minimum found by analytic W_aa roots, then repeat and independently
refine root brackets. This still is not a global interval-arithmetic proof.
"""
import json
import time
import numpy as np
from scipy.optimize import brentq

from .angular_environment_reference import AngularInterfaceResearchSurface
from .fcc111_active_interface import FCC111ActiveInterface
from .joint_fcc_interface_calibration import JointFit,IDEAL_H,build_joint_model
from .monotone_opening_calibration import constrained_opening_fit,opening_force_basis
from .odd_moment_calibration import odd_targets
from .run_matched_interface_study import ROOT,save_json


def surface_from_fit(fit,tolerance=2e-11):
    c=np.asarray(fit["coefficients"])
    if not fit["strictly_positive_LJ_resolved"]:
        raise ValueError("unresolved or absent LJ attraction; no material surface")
    bulk=build_joint_model(JointFit("opening_constraint_research",fit["scalar_decay"],c[:5],
                                  fit["squared_loss"],None,None)); bulk.a0=IDEAL_H
    return AngularInterfaceResearchSurface(FCC111ActiveInterface(bulk,tolerance=tolerance),c[5],
        tolerance=tolerance,angular_decay=fit["angular_decay"],vector_amplitude_ev=c[6],
        quadrupole_amplitude_ev=c[7])


def opening_force_extrema(surface,count=80):
    h=surface.h; grid=np.linspace(1.001,12.,count)
    values=[surface.packed(float(r*h),0.) for r in grid]
    ratios=[]
    for left,right,vl,vr in zip(grid[:-1],grid[1:],values[:-1],values[1:]):
        if vl[3]*vr[3]<0:
            ratios.append(brentq(lambda r:surface.packed(float(r*h),0.)[3],left,right,xtol=2e-11))
    ratios=sorted([1.001,*ratios,12.])
    return [dict(a_over_h=r,force_ev_coordinate=float(surface.packed(r*h,0.)[1]),
                 normal_curvature=float(surface.packed(r*h,0.)[3])) for r in ratios]


def main():
    started=time.perf_counter()
    fit=json.loads((ROOT/"monotone_opening_fit.json").read_text(encoding="utf-8"))["best"]
    target,scales=odd_targets(); ratios=list(fit["opening_inequality_a_over_h"])
    history=[]; stopped=False
    for iteration in range(12):
        surface=surface_from_fit(fit)
        extrema=opening_force_extrema(surface)
        minimum=min(extrema,key=lambda row:row["force_ev_coordinate"])
        point=minimum["a_over_h"]
        fine=surface_from_fit(fit,tolerance=2e-12).packed(point*IDEAL_H,0.)[1]
        coarse=surface_from_fit(fit,tolerance=2e-9).packed(point*IDEAL_H,0.)[1]
        row=opening_force_basis(fit["scalar_decay"],fit["angular_decay"],(point,))[0]
        c=np.asarray(fit["coefficients"])
        arithmetic=64*np.finfo(float).eps*float(abs(row)@abs(c))
        floor=float(max(abs(fine-coarse),abs(fine-row@c),arithmetic))
        history.append(dict(iteration=iteration,fit=fit,extrema=extrema,
            minimum_force_ev_coordinate=minimum["force_ev_coordinate"],
            force_refinement_floor=floor,minimum_ratio=point))
        if minimum["force_ev_coordinate"]>=-4*floor:
            stopped=True; break
        ratios.append(point)
        fit=constrained_opening_fit(fit["scalar_decay"],fit["angular_decay"],target,scales,ratios=ratios)
        if not fit["admissible"] or not fit["strictly_positive_LJ_resolved"]:
            break
    refined=opening_force_extrema(surface_from_fit(fit),160) if fit["strictly_positive_LJ_resolved"] else []
    save_json(ROOT/"monotone_opening_refined.json",dict(best=fit,exchange_iterations=history,
        refined_extrema=refined,stopped_at_force_resolution=stopped,
        elapsed_seconds=time.perf_counter()-started,
        material_accepted=False,global_monotonicity_proven=False,
        scope="1.001<=a/h<=12; analytic curvature roots on 80/160 brackets, numerical forces retained"))
    print(dict(squared_loss=fit["squared_loss"],coefficients=fit["coefficients"],
        iterations=len(history),stopped_at_resolution=stopped,refined_extrema=refined),flush=True)


if __name__=="__main__":
    main()
