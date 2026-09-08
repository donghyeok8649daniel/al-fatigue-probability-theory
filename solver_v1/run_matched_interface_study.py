"""Reproducible calibration THEN physical/static stress-scenario execution.

Outputs are isolated from old calibration files and production defaults.
External EAM interpolation supplies comparison targets only; every candidate
uses the analytic LJ infinite sums. No dynamics or physical time is inferred.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
import json
import math
from pathlib import Path
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import lsq_linear, minimize_scalar

from .aluminum_calibration import EV_J
from .analytic_density_feasibility import (
    deterministic_mixture_search, intrinsic_fault_upper_bound, squared_envelope_profile,
)
from .angular_interface_calibration import angular_profile
from .angular_environment_reference import AngularInterfaceResearchSurface
from .fcc111_active_interface import FCC111ActiveInterface
from .fcc111_geometry import DIRECT_110, SHOCKLEY_112, registry_path
from .full_fcc_calibration_audit import cubic_constants_gpa, relaxed_bulk_observables
from .interface_static_scenarios import (
    InterfaceUnits, continue_static_branch, hessian_from_packed, packed_evaluation,
)
from .joint_fcc_interface_calibration import (
    INTERFACE_STATES, JOINT_NAMES, JointFit, joint_basis, joint_equilibrium,
    joint_profile, matched_targets,
)
from .reference_eam_targets import MishinRigidFCCReference, SOURCE_SHA256, SOURCE_URL
from .run_full_fcc_calibration_audit import interface_refinement, write_csv, json_safe
from .run_full_fcc_calibration_audit import fixed_registry_opening_barriers


ROOT=Path(__file__).resolve().parents[1]/"results/fcc111_active_interface/matched_v3"


def save_json(path,value):
    path.write_text(json.dumps(value,indent=2,default=json_safe,allow_nan=False)+"\n",encoding="utf-8")


def record_fit(name,fit,target,scales):
    return [dict(family=name,target=key,reference=float(y),prediction=float(p),
                 absolute_error=float(p-y),relative_error=float((p-y)/y) if y else None,
                 normalization_scale=float(scale),normalized_residual=float(r),used_for_fit=True)
            for key,y,p,scale,r in zip(JOINT_NAMES,target,fit["predictions"],scales,fit["residuals"])]


def signed_quadratic_study(target,scales):
    """Do not confuse a convex-sector failure with unrestricted C failure."""
    def fit(decay):
        matrix=joint_basis(float(decay))
        result=lsq_linear(matrix/scales[:,None],target/scales,
            bounds=([0.,0.,0.,-np.inf,-np.inf],np.full(5,np.inf)),tol=1e-11,max_iter=300)
        prediction=matrix@result.x; residual=(prediction-target)/scales
        return dict(decay=float(decay),coefficients=result.x,predictions=prediction,
                    residuals=residual,squared_loss=float(residual@residual),
                    success=bool(result.success),message=result.message)
    grid=np.geomspace(.35,8.,33); rows=[fit(d) for d in grid]
    for i in range(1,len(grid)-1):
        if rows[i]["squared_loss"]<=min(rows[i-1]["squared_loss"],rows[i+1]["squared_loss"]):
            optimum=minimize_scalar(lambda d:fit(d)["squared_loss"],bounds=grid[i-1:i+2:2],
                                    method="bounded",options={"xatol":2e-6})
            rows.append(fit(optimum.x))
    return min(rows,key=lambda r:r["squared_loss"]),rows


def calibrate(out):
    started=time.perf_counter(); source=MishinRigidFCCReference()
    target,scales=matched_targets(source)
    targets=[]
    for i,(name,y,scale) in enumerate(zip(JOINT_NAMES,target,scales)):
        units=("eV/atom per unit normal strain" if i==0 else "eV/atom" if i==1 else
               "eV/atom per squared unit strain" if i<5 else "eV/interface_cell")
        targets.append(dict(quantity=name,value=float(y),units=units,
            scale=float(scale),scale_meaning="declared model discrepancy, not measurement uncertainty",
            method="Mishin1999 EAM rounded bulk values" if i<5 else "recomputed NIST setfl on rigid halves",
            temperature_K=0,lattice_angstrom=4.05,orientation="FCC; interface (111)",
            relaxation="affine bulk; rigid half-crystals for interface",source_doi="10.1103/PhysRevB.59.3393",
            source_url=SOURCE_URL,source_sha256=SOURCE_SHA256))
    write_csv(out/"matched_targets.csv",targets)
    save_json(out/"reference_provenance.json",dict(source_url=SOURCE_URL,sha256=source.sha256,
        source_role="independent target generator ONLY; not a canonical potential",
        source_cutoff_angstrom=source.cutoff,interpolation="cubic; linear cross-check",
        lattice_angstrom=4.05,atomic_cell_area_angstrom2=source.geometry.atomic_cell_area,
        atomic_area_is_not_statistical_correlation_area=True,
        cohesion_ev_atom=source.cohesion(),work_separation_J_m2=source.work_separation()*16.02176634/source.geometry.atomic_cell_area))
    all_fits={}; residuals=[]
    for family in ("sqrt","linear","quadratic"):
        result,profile=joint_profile(target,scales,family=family)
        data=asdict(result); all_fits[family]=data
        save_json(out/(family+"_profile.json"),[asdict(row) for row in profile])
        residuals.extend(record_fit(family,data,target,scales))
        save_json(out/"fitted_candidates.json",all_fits)
        print(f"{family}: actual fit loss={result.squared_loss:.9g}",flush=True)
    signed,profiles=signed_quadratic_study(target,scales)
    all_fits["signed_quadratic_exploratory"]=signed
    save_json(out/"signed_quadratic_profile.json",profiles)
    residuals.extend(record_fit("signed_quadratic_exploratory",signed,target,scales))
    for name,runner in (("positive_two_exponential",deterministic_mixture_search),
                        ("squared_envelope",lambda y,z:squared_envelope_profile(y,z,grid=np.geomspace(.5,12.,13))),
                        ("angular_tied",lambda y,z:angular_profile(y,z,tied=True)),
                        ("angular_independent_range",lambda y,z:angular_profile(y,z,tied=False))):
        best,profile,logs=runner(target,scales)
        all_fits[name]=best
        save_json(out/(name+"_profile.json"),dict(rows=profile,local_optimization=logs))
        save_json(out/"fitted_candidates.json",all_fits)
        residuals.extend(record_fit(name,best,target,scales))
        print(f"{name}: actual fit loss={best['squared_loss']:.9g}",flush=True)
    write_csv(out/"fit_residuals.csv",residuals)
    feasibility=[intrinsic_fault_upper_bound(d,target,scales) for d in np.geomspace(.35,8.,49)]
    save_json(out/"fixed_decay_feasibility.json",feasibility)
    # Local log sensitivity in the independent coefficient/gauge coordinates.
    fit=all_fits["linear"]; coeff=np.asarray(fit["coefficients"]); decay=fit["decay"]
    columns=[joint_basis(decay)[:,i]*coeff[i]/scales for i in range(4)]
    step=2e-4
    columns.append((joint_basis(decay*math.exp(step))@coeff-joint_basis(decay*math.exp(-step))@coeff)/(2*step*scales))
    jac=np.column_stack(columns); singular=np.linalg.svd(jac,compute_uv=False)
    write_csv(out/"sensitivity_matrix.csv",[dict(target=name,**{p:float(v) for p,v in zip(("log_u","log_v","log_A","log_B","log_decay"),row)}) for name,row in zip(JOINT_NAMES,jac)])
    write_csv(out/"singular_values.csv",[dict(index=i,singular_value=float(value)) for i,value in enumerate(singular)])
    save_json(out/"calibration_execution.json",dict(actual_optimization_executed=True,
        elapsed_seconds=time.perf_counter()-started,linear_sensitivity_condition=float(singular[0]/singular[-1]),
        fit_targets_count=len(target),heldout="remaining registry/opening curves and stress responses",
        source_is_not_independent_experiment=True,production_changed=False))


def surfaces_from_results(out):
    data=json.loads((out/"fitted_candidates.json").read_text(encoding="utf-8"))
    length=4.05/math.sqrt(2)*1e-10
    results=[]; bulk_rows=[]
    for name in ("linear","angular_tied","angular_independent_range","angular_vector_and_third",
                 "angular_even_odd","angular_monotone_opening"):
        if name not in data:
            continue
        fit=data[name]
        D1=D2=0.
        if name=="linear":
            coefficients=np.asarray(fit["coefficients"]); decay=fit["decay"]; D=None
        elif name in ("angular_even_odd","angular_monotone_opening"):
            params=dict(zip(fit["coefficient_order"],fit["coefficients"]))
            coefficients=np.array([params.get(key,0.) for key in ("u","v","A","B","C")])
            decay=fit["scalar_decay"]; D=params["D3"]; D1=params["D1"]; D2=params["D2"]
            if not fit.get("strictly_positive_LJ_resolved",False):
                bulk_rows.append(dict(candidate=name,status="unresolved positive LJ coefficient; no surface constructed"))
                continue
        else:
            coefficients=np.r_[np.asarray(fit["coefficients"])[:4],0.]
            decay=fit["scalar_decay"]; D=fit["coefficients"][4]
            if name=="angular_vector_and_third":
                D1=fit["coefficients"][5]
        if np.any(coefficients[:2]<=0):
            bulk_rows.append(dict(candidate=name,status="invalid nonpositive LJ coefficient")); continue
        qfit=JointFit(name,decay,coefficients,fit["squared_loss"],None,None)
        try:
            bulk,stretch=joint_equilibrium(qfit)
        except ValueError as exc:
            bulk_rows.append(dict(candidate=name,status=str(exc))); continue
        obs=relaxed_bulk_observables(bulk)
        if D2:
            from .angular_environment_reference import AngularInterfaceInvariant
            from .even_moment_calibration import quadrupole_bulk_curvatures
            even=AngularInterfaceInvariant(FCC111ActiveInterface(bulk,tolerance=2e-11),
                rank=2,angular_decay=fit["angular_decay"],tolerance=2e-11)
            obs[2:]+=D2*quadrupole_bulk_curvatures(even)
        constants=cubic_constants_gpa(obs,volume_scale=stretch**3)
        area=bulk.geometry.atomic_cell_area*length**2
        bulk_rows.append(dict(candidate=name,status="stable FCC root; interface not yet accepted",
            lattice_angstrom=4.05*stretch,cohesion_ev_atom=obs[1],normal_force_ev=obs[0],
            epsilon_LJ_ev=bulk.p.epsilon,sigma_LJ_angstrom=bulk.p.sigma_lj*length/1e-10,
            **constants))
        for path in (DIRECT_110,SHOCKLEY_112):
            face=FCC111ActiveInterface(bulk,path_id=path,tolerance=2e-11)
            if D is None:
                evaluate=lambda a,s,face=face:packed_evaluation(face,a,s)
            else:
                surface=AngularInterfaceResearchSurface(face,D,angular_decay=fit["angular_decay"],
                    vector_amplitude_ev=D1,quadrupole_amplitude_ev=D2)
                evaluate=surface.packed
            results.append(dict(name=name,path=path,evaluate=evaluate,h=face.h,
                                period=face.period,units=InterfaceUnits(length,area),bulk=bulk))
    source=MishinRigidFCCReference()
    for path in (DIRECT_110,SHOCKLEY_112):
        evaluate=lambda a,s,path=path:source.interface_derivatives(a,s,path_id=path)
        results.append(dict(name="Mishin_source_rigid",path=path,evaluate=evaluate,h=source.h,
            period=registry_path(path).period_over_b*source.geometry.b,
            units=InterfaceUnits(1e-10,source.geometry.atomic_cell_area*1e-20),bulk=None))
    write_csv(out/"bulk_scenario_summary.csv",bulk_rows)
    return results


def scenarios(out):
    started=time.perf_counter(); surfaces=surfaces_from_results(out)
    registry=[]; opening=[]; stress=[]; snapshots=[]; heldout=[]; refinements=[]
    for surface in surfaces:
        name,path=surface["name"],surface["path"]
        ev=surface["evaluate"]; h=surface["h"]; period=surface["period"]; units=surface["units"]
        # These full curves were not all included in the fit; fit points remain
        # marked. Reference and candidates have their own zero-pressure lattice.
        stop=1. if path==DIRECT_110 else 1/3
        for fraction in np.linspace(0,stop,49):
            v=ev(h,float(fraction*period))
            row=dict(candidate=name,path=path,s_over_period=fraction,
                s_angstrom=fraction*period*units.length_scale_m/1e-10,
                energy_J_m2=float(units.energy_to_surface(v[0])),energy_ev_cell=v[0],
                normal_traction_GPa=float(units.force_to_traction(v[1])),
                shear_traction_GPa=float(units.force_to_traction(v[2])),
                min_hessian_eigenvalue=float(np.linalg.eigvalsh(hessian_from_packed(v))[0]))
            registry.append(row)
        if path==DIRECT_110:
            for ratio in np.r_[np.linspace(.95,2.,36),3.,5.,20.,100.]:
                v=ev(float(ratio*h),0.)
                opening.append(dict(candidate=name,a_over_h=ratio,
                    opening_angstrom=(ratio-1)*h*units.length_scale_m/1e-10,
                    energy_J_m2=float(units.energy_to_surface(v[0])),
                    normal_traction_GPa=float(units.force_to_traction(v[1])),
                    registry_curvature=v[5]))
        # Explicit loads: no undocumented mapping from nominal axial stress.
        tensile=np.array([[t,0.] for t in (0.,.05,.15,.4,1.,2.,4.,6.,8.,10.)])
        compress=np.array([[t,0.] for t in (0.,-.1,-.5,-.9,-2.)])
        shear=np.array([[0.,t] for t in (0.,.05,.15,.4,.8,1.2,2.,3.,4.)])
        combined=np.array([[t/2,t/2] for t in (0.,.1,.3,.8,2.,4.,8.,12.)])
        unload=np.array([[0.,0.],[.5,.2],[1.,.4],[.5,.2],[0.,0.]])
        for protocol,loads in (("normal_tension",tensile),("compression_control",compress),
                               ("pure_resolved_shear",shear),("explicit_45degree_axis",combined),
                               ("quasistatic_unload",unload)):
            rows=continue_static_branch(ev,h=h,period=period,units=units,loads_gpa=loads,
                                        max_traction_step=.2)
            stress.extend(dict(candidate=name,path=path,protocol=protocol,**row) for row in rows)
        for label,a,s in (("pristine",h,0.),("compressed",.98*h,0.),("opened",1.2*h,0.),
                          ("partial_registry",h,period/3 if path==SHOCKLEY_112 else .5*period),
                          ("coupled",1.05*h,.12*period)):
            v=ev(a,s); step=2e-5*min(h,period)
            va=(ev(a+step,s)-ev(a-step,s))/(2*step)
            vs=(ev(a,s+step)-ev(a,s-step))/(2*step)
            snapshots.append(dict(candidate=name,path=path,scenario=label,
                energy_ev_cell=v[0],energy_J_m2=float(units.energy_to_surface(v[0])),
                normal_traction_GPa=float(units.force_to_traction(v[1])),
                shear_traction_GPa=float(units.force_to_traction(v[2])),
                min_hessian_eigenvalue=float(np.linalg.eigvalsh(hessian_from_packed(v))[0]),
                gradient_FD_error=float(max(abs(va[0]-v[1]),abs(vs[0]-v[2]))),
                hessian_FD_error=float(max(abs(va[1]-v[3]),abs(va[2]-v[4]),abs(vs[1]-v[4]),abs(vs[2]-v[5])))))
        write_csv(out/"stress_scenarios.csv",stress)
        write_csv(out/"scenario_validation.csv",snapshots)
        print(f"{name}/{path}: actual static stress and unload scenarios complete",flush=True)
    for path in (DIRECT_110,SHOCKLEY_112):
        reference=[r for r in registry if r["candidate"]=="Mishin_source_rigid" and r["path"]==path]
        for name in dict.fromkeys(r["name"] for r in surfaces if r["name"]!="Mishin_source_rigid"):
            candidate=[r for r in registry if r["candidate"]==name and r["path"]==path]
            if len(candidate)!=len(reference):
                continue
            fit_fractions=(0.,.5,1.) if path==DIRECT_110 else (0.,1/6,1/3)
            heldout_indices=[i for i,row in enumerate(reference)
                             if not any(abs(row["s_over_period"]-f)<1e-10 for f in fit_fractions)]
            error=np.array([candidate[i]["energy_J_m2"]-reference[i]["energy_J_m2"] for i in heldout_indices])
            heldout.append(dict(candidate=name,path=path,curve_rmse_J_m2=float(np.sqrt(np.mean(error**2))),
                maximum_absolute_error_J_m2=float(np.max(abs(error))),
                endpoint_prediction_J_m2=candidate[-1]["energy_J_m2"],
                endpoint_reference_J_m2=reference[-1]["energy_J_m2"],
                heldout_count=len(heldout_indices),fitted_fractions_excluded=str(fit_fractions),
                endpoint_is_fit_or_periodic_reference=True,
                convention="rigid halves; candidate's own relaxed bulk lattice; fitted fractions excluded from RMSE"))
    # Independent a/s reciprocal/direct refinement of unchanged scalar core.
    linear=next((s for s in surfaces if s["name"]=="linear" and s["path"]==DIRECT_110),None)
    refinements=interface_refinement("linear_joint_fit",linear["bulk"]) if linear is not None else []
    for filename,rows in (("registry_curves.csv",registry),("opening_curves.csv",opening),
                          ("heldout_curve_validation.csv",heldout),("convergence_summary.csv",refinements)):
        write_csv(out/filename,rows)
    fig,axes=plt.subplots(1,3,figsize=(13,4))
    for name in dict.fromkeys(r["name"] for r in surfaces):
        for axis,path in zip(axes[:2],(DIRECT_110,SHOCKLEY_112)):
            rows=[r for r in registry if r["candidate"]==name and r["path"]==path]
            axis.plot([r["s_over_period"] for r in rows],[r["energy_J_m2"] for r in rows],label=name)
        rows=[r for r in opening if r["candidate"]==name and r["a_over_h"]<=3]
        axes[2].plot([r["opening_angstrom"] for r in rows],[r["energy_J_m2"] for r in rows],label=name)
    for axis,title in zip(axes,("Direct <110> registry","Shockley partial registry","Normal opening")):
        axis.set_title(title); axis.set_ylabel("Energy [J/m2]"); axis.grid(alpha=.25)
    axes[0].set_xlabel("s / full repeat"); axes[1].set_xlabel("s / full repeat")
    axes[2].set_xlabel("Opening displacement [angstrom]"); axes[0].legend(fontsize=6)
    fig.suptitle("Static research scenarios; not accepted Al / no physical kinetic time")
    fig.tight_layout(); fig.savefig(out/"scenario_comparison.png",dpi=150); plt.close(fig)
    save_json(out/"scenario_execution.json",dict(elapsed_seconds=time.perf_counter()-started,
        source_geometry="rigid halves",stress_units="GPa; explicit normal and resolved shear tractions",
        unload_protocol="quasistatic equilibrium only, NOT a dynamic hold or residual-plasticity test",
        kinetics="unavailable",physical_Hz_enabled=False,production_PDE_changed=False,
        no_branch_failure_is_automatically_a_spinodal=True))


def summarize_saved_curves(out):
    """Postprocess existing numerical curves; exclude training points exactly."""
    with (out/"registry_curves.csv").open(encoding="utf-8") as stream:
        rows=list(csv.DictReader(stream))
    summaries=[]
    for path in (DIRECT_110,SHOCKLEY_112):
        source=[r for r in rows if r["candidate"]=="Mishin_source_rigid" and r["path"]==path]
        fractions=(0.,.5,1.) if path==DIRECT_110 else (0.,1/6,1/3)
        selected=[i for i,r in enumerate(source)
                  if not any(abs(float(r["s_over_period"])-f)<1e-10 for f in fractions)]
        for name in dict.fromkeys(r["candidate"] for r in rows if r["candidate"]!="Mishin_source_rigid"):
            candidate=[r for r in rows if r["candidate"]==name and r["path"]==path]
            if len(source)!=len(candidate):
                raise ValueError("candidate/source curve sampling mismatch")
            error=np.array([float(candidate[i]["energy_J_m2"])-float(source[i]["energy_J_m2"]) for i in selected])
            summaries.append(dict(candidate=name,path=path,heldout_count=len(selected),
                curve_rmse_J_m2=float(np.sqrt(np.mean(error**2))),
                maximum_absolute_error_J_m2=float(np.max(abs(error))),
                endpoint_prediction_J_m2=float(candidate[-1]["energy_J_m2"]),
                endpoint_reference_J_m2=float(source[-1]["energy_J_m2"]),
                fitted_fractions_excluded=str(fractions),endpoint_is_fit_or_periodic_reference=True,
                convention="rigid halves; candidate's own relaxed bulk lattice; training points excluded from RMSE"))
    write_csv(out/"heldout_curve_validation.csv",summaries)


def refine_scenarios(out):
    started=time.perf_counter(); surfaces=surfaces_from_results(out)
    continuation=[]; barriers=[]; local_stiffness=[]
    for surface in surfaces:
        if surface["path"]!=DIRECT_110:
            continue
        name=surface["name"]; evaluate=surface["evaluate"]
        h=surface["h"]; period=surface["period"]; units=surface["units"]
        for protocol,loads in (("normal",[[0,0],[.05,0],[.15,0],[.4,0],[1.,0],[2.,0],[4.,0]]),
                               ("mixed",[[0,0],[.5,.2],[1.,.4],[0.,0.]])):
            histories=[]
            for step in (.2,.1,.05):
                history=continue_static_branch(evaluate,h=h,period=period,units=units,
                    loads_gpa=loads,max_traction_step=step)
                histories.append(history)
            for step,history in zip((.2,.1,.05),histories):
                for row,reference in zip(history,histories[-1]):
                    valid=row["a"] is not None and reference["a"] is not None
                    continuation.append(dict(candidate=name,protocol=protocol,traction_step_GPa=step,
                        **row,a_error_angstrom=(row["a"]-reference["a"])*units.length_scale_m/1e-10 if valid else None,
                        s_error_angstrom=(row["s"]-reference["s"])*units.length_scale_m/1e-10 if valid else None))
        class ConstrainedReference:
            # STATIC diagnostic adapter, not registered as a production surface.
            def __init__(self):
                self.h=h
            def energy(self,a,s):
                return evaluate(a,s)[0]
            def grad(self,a,s):
                return evaluate(a,s)[1:3]
            def hessian(self,a,s):
                return hessian_from_packed(evaluate(a,s))
        local_barriers=fixed_registry_opening_barriers(name,ConstrainedReference(),
            units.length_scale_m,units.atomic_cell_area_m2)
        barriers.extend(dict(**row,coordinate_length_scale_m=units.length_scale_m,
                             atomic_cell_area_m2=units.atomic_cell_area_m2) for row in local_barriers)
        v=evaluate(h,0)
        local_stiffness.append(dict(candidate=name,coordinate_length_scale_m=units.length_scale_m,
            atomic_cell_area_m2=units.atomic_cell_area_m2,normal_tangent_GPa=float(h*units.force_to_traction(v[3])),
            registry_tangent_GPa=float(h*units.force_to_traction(v[5])),
            coupling_GPa=float(h*units.force_to_traction(v[4])),
            definition="traction per opening strain delta_a/h or shear s/h, fixed area; not bulk Young modulus"))
        print(f"{name}: continuation refinement and constrained opening branch computed",flush=True)
    write_csv(out/"stress_continuation_refinement.csv",continuation)
    write_csv(out/"constrained_opening_barriers.csv",barriers)
    write_csv(out/"local_interface_stiffness.csv",local_stiffness)
    summarize_saved_curves(out)
    save_json(out/"refinement_execution.json",dict(elapsed_seconds=time.perf_counter()-started,
        refinement="traction continuation step .2/.1/.05 GPa; not a dynamical dt",
        barrier_status="fixed s=0 constrained opening, NOT the full coupled spinodal or dynamic event ordering"))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase",choices=("calibration","scenarios","refinement","summaries","angular-extension",
        "odd-scenarios","even-scenarios","monotone-scenarios","all"),default="all")
    args=parser.parse_args(); ROOT.mkdir(parents=True,exist_ok=True)
    if args.phase in ("calibration","all"):
        calibrate(ROOT)
    if args.phase in ("scenarios","all"):
        scenarios(ROOT)
    if args.phase in ("refinement","all"):
        refine_scenarios(ROOT)
    if args.phase=="summaries":
        summarize_saved_curves(ROOT)
    if args.phase=="angular-extension":
        target,scales=matched_targets()
        started=time.perf_counter()
        best,rows,logs=angular_profile(target,scales,tied=False,quadratic=True,max_nfev=220)
        save_json(ROOT/"angular_quadratic_extension.json",dict(best=best,profile=rows,local_optimization=logs,
            elapsed_seconds=time.perf_counter()-started,
            classification="extra complexity not accepted; compare every residual and search-bound sensitivity"))
        print(best,flush=True)
    if args.phase=="odd-scenarios":
        fit=json.loads((ROOT/"odd_moment_fit.json").read_text(encoding="utf-8"))["best"]
        folder=ROOT/"odd_moment"; folder.mkdir(parents=True,exist_ok=True)
        save_json(folder/"fitted_candidates.json",dict(angular_vector_and_third=fit))
        scenarios(folder)
        refine_scenarios(folder)
    if args.phase=="even-scenarios":
        fit=json.loads((ROOT/"even_moment_convex_signed_free_B_fit.json").read_text(encoding="utf-8"))["best"]
        folder=ROOT/"even_moment"; folder.mkdir(parents=True,exist_ok=True)
        save_json(folder/"fitted_candidates.json",dict(angular_even_odd=fit))
        scenarios(folder)
        refine_scenarios(folder)
    if args.phase=="monotone-scenarios":
        fit=json.loads((ROOT/"monotone_opening_refined.json").read_text(encoding="utf-8"))["best"]
        folder=ROOT/"monotone_opening"; folder.mkdir(parents=True,exist_ok=True)
        save_json(folder/"fitted_candidates.json",dict(angular_monotone_opening=fit))
        scenarios(folder)
        refine_scenarios(folder)


if __name__=="__main__":
    main()
