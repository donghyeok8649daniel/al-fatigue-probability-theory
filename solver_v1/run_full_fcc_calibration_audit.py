"""Actually rerun deterministic calibration and static physical-readiness gates."""
from __future__ import annotations

import csv
from dataclasses import asdict, replace
import json
from pathlib import Path
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import root_scalar

from .aluminum_calibration import EV_J, cell_ev_to_energy_density, mishin_lu_wei_targets
from .aluminum_full_fcc_calibration import FullFCCParameters
from .fcc111_active_interface import FCC111ActiveInterface
from .fcc111_full_energy import FullFCC111StackEnergy
from .fcc111_geometry import DIRECT_110, SHOCKLEY_112
from .full_fcc_calibration_audit import (
    IDEAL_H, NAMES, bounded_independent_fit, cubic_constants_gpa,
    deterministic_profile, independent_bulk_targets, independent_observables,
    independent_sensitivity, profile_point,
    isotropic_bulk_equilibrium, relaxed_bulk_observables,
)
from .run_aluminum_full_fcc_calibration import LINEAR_CANDIDATES

RESULT_ROOT = Path(__file__).resolve().parents[1]/"results"
FIELDS = ("energy","d_da","d_ds","d2_daa","d2_das","d2_dss")


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w",newline="",encoding="utf-8") as stream:
        keys=list(dict.fromkeys(key for row in rows for key in row))
        writer=csv.DictWriter(stream,fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def json_safe(value):
    if isinstance(value,np.ndarray):
        return value.tolist()
    if isinstance(value,(np.bool_,np.integer,np.floating)):
        return value.item()
    raise TypeError(type(value).__name__)


def interface_refinement(name, bulk):
    a,s=1.05*bulk.a0,.17
    cfg=bulk.stack_config
    accurate=FullFCC111StackEnergy(
        bulk.p,density_params=bulk.density_params,embedding=bulk.embedding,
        stack_config=replace(cfg,tol=2e-12,
            reciprocal=replace(cfg.reciprocal,tol=2e-14)),
        find_equilibrium=False)
    accurate.a0=bulk.a0
    reference=FCC111ActiveInterface(accurate,tolerance=2e-12)
    rv=reference.evaluate(a,s)
    rows=[]
    for label,model,tol,consecutive in (
        ("reference",bulk,2e-8,4),("density_layer_refined",bulk,2e-10,6),
        ("reciprocal_and_layer_refined",accurate,2e-12,6)):
        interface=FCC111ActiveInterface(model,tolerance=tol,
                                       consecutive_small_layers=consecutive)
        value=interface.evaluate(a,s)
        for field in FIELDS:
            rows.append({"set":name,"study":label,"field":field,
                "value":getattr(value,field),"reference":getattr(rv,field),
                "absolute_error":abs(getattr(value,field)-getattr(rv,field)),
                "tolerance":tol,"layers_used":value.layers_used,
                "shells_used":value.maximum_plane_shells,
                "tail_status":"empirical proxy; refinement differences authoritative"})
    for radius,layers in ((24,24),(48,24),(24,48),(48,48),(80,80)):
        value=reference.direct_reference(a,s,radial_index=radius,layers=layers)
        for field in FIELDS:
            rows.append({"set":name,"study":"independent_real_space","field":field,
                "value":getattr(value,field),"reference":getattr(rv,field),
                "absolute_error":abs(getattr(value,field)-getattr(rv,field)),
                "direct_radius":radius,"direct_layers":layers})
    return rows


def fixed_registry_opening_barriers(name, interface, length_m, area_m2):
    """Follow the first stable s=0 branch, not a unimodal global-peak guess.

    A narrow positive traction peak can precede a negative-traction region.
    A single bounded minimization over the whole interval can miss it.  Track
    the first positive-to-negative W_aa crossing from the intact reference.
    This remains a constrained normal spinodal, never a coupled stability claim.
    """
    curvature=lambda a:float(interface.hessian(float(a),0)[0,0])
    if curvature(interface.h)<=0:
        return [{"set":name,"status":"no locally stable normal reference"}]
    scan=interface.h*(1+np.r_[0.,np.geomspace(1e-5,3.,96)])
    left=scan[0]; c_left=curvature(left); peak_a=None
    for right in scan[1:]:
        c_right=curvature(right)
        if c_left>0>=c_right:
            peak_a=float(root_scalar(curvature,bracket=(left,right),
                                     method="brentq",xtol=1e-11).root)
            break
        left,c_left=right,c_right
    if peak_a is None:
        return [{"set":name,"status":"first normal spinodal not bracketed in a/h<=4"}]
    maximum=float(interface.grad(peak_a,0)[0])
    if maximum<=0:
        return [{"set":name,"status":"first normal spinodal has nonpositive traction"}]
    rows=[]
    for fraction in (.25,.5,.75,.95):
        force=fraction*maximum
        derivative=lambda a:float(interface.grad(float(a),0)[0]-force)
        minimum=root_scalar(derivative,
            bracket=(.97*interface.h,peak_a),method="brentq").root
        # First unstable crossing adjacent to that branch.  Do not jump over
        # intervening extrema or claim a global event order.
        saddle_bracket=None
        left=peak_a
        for right in peak_a+interface.h*np.geomspace(1e-6,100.,128):
            if derivative(right)<=0:
                saddle_bracket=(left,float(right))
                break
            left=float(right)
        if saddle_bracket is None:
            rows.append({"set":name,"fraction_of_fixed_s_peak":fraction,
                         "status":"saddle not bracketed"})
            continue
        saddle=root_scalar(derivative,bracket=saddle_bracket,
                           method="brentq").root
        barrier=interface.energy(saddle,0)-interface.energy(minimum,0)-force*(saddle-minimum)
        rows.append({"set":name,"fraction_of_fixed_s_peak":fraction,
            "normal_force_eV_per_reduced_a":force,
            "normal_traction_GPa":force*EV_J/(length_m*area_m2)/1e9,
            "bound_a_reduced_L0":minimum,"saddle_a_reduced_L0":saddle,
            "opening_barrier_eV_cell":barrier,
            "fixed_s_spinodal_a_reduced_L0":peak_a,
            "fixed_s_spinodal_force_eV_per_reduced_a":maximum,
            "registry_curvature_at_minimum":interface.hessian(minimum,0)[1,1],
            "status":"fixed s=0; zero resolved shear; not coupled spinodal"})
    return rows


def main():
    started=time.perf_counter()
    out=RESULT_ROOT/"aluminum_full_fcc_calibration"/"audited_v2"
    interface_out=RESULT_ROOT/"fcc111_active_interface"/"audited_v2"
    out.mkdir(parents=True,exist_ok=True)
    interface_out.mkdir(parents=True,exist_ok=True)
    target,scales=independent_bulk_targets()
    geometry,_=mishin_lu_wei_targets()
    area=geometry.atomic_cell_area_angstrom2
    length_m=geometry.burgers_angstrom*1e-10
    write_csv(out/"independent_targets.csv",[
        {"observable":name,"target":float(y),"scale":float(scale),
         "units":"eV/atom; derivatives with respect to dimensionless strain",
         "temperature_K":0,"method":"Mishin 1999 EAM reference",
         "source":"doi:10.1103/PhysRevB.59.3393","used_for_fit":True}
        for name,y,scale in zip(NAMES,target,scales)])
    grid=np.geomspace(.35,8.,49)
    profiles=[]; selected={}
    for extended in (False,True):
        family="linear" if extended else "square"
        for indices,stage in (((0,1),"equilibrium_cohesion"),
                              ((0,1,2,3),"normal_hydro"),
                              ((0,1,2,3,4),"complete_elastic")):
            if stage=="equilibrium_cohesion":
                rows=[profile_point(float(d),extended=extended,
                                    observation_indices=indices) for d in grid]
            else:
                best,rows=deterministic_profile(extended=extended,decay_grid=grid,
                                                observation_indices=indices)
                if stage=="complete_elastic":
                    selected[f"{family}_profile"]=best.parameters
                    print(f"{family}: loss={best.squared_loss:.8g}, decay={best.decay:.8g}",flush=True)
            for row in rows:
                profiles.append({"family":family,"stage":stage,"decay":row.decay,
                    "squared_loss":row.squared_loss,"admissible":row.admissible,
                    "reason":row.reason,"u":row.coefficients[0],"v":row.coefficients[1],
                    "A":row.coefficients[2],"B":row.coefficients[3],
                    "independent_observation_count":len(indices),
                    "parameter_count":5 if extended else 4})
    write_csv(out/"staged_profile.csv",profiles)
    optimization=[]; bounded=[]
    starts=(*LINEAR_CANDIDATES,FullFCCParameters(1.,.75,5.,4.,10.))
    for i,start in enumerate(starts):
        fitted,record=bounded_independent_fit(start)
        optimization.append({"index":i,**record})
        bounded.append((fitted,record["squared_loss"]))
        print(f"bounded start {i}: loss={record['squared_loss']:.8g}",flush=True)
    selected["linear_historical_box"]=min(bounded,key=lambda row:row[1])[0]
    (out/"optimization_log.json").write_text(
        json.dumps(optimization,indent=2,default=json_safe)+"\n",encoding="utf-8")
    parameters=[]; residuals=[]; sensitivities=[]; svds=[]; convergence=[]; models={}
    for name,p in {**selected,"legacy_unidentified":LINEAR_CANDIDATES[0]}.items():
        y=independent_observables(p)
        model,stretch=isotropic_bulk_equilibrium(p)
        models[name]=model
        relaxed=relaxed_bulk_observables(model)
        # Target-geometry decomposition is evaluated on the undeformed basis,
        # not at a non-FCC point of the newly relaxed b.
        from .aluminum_full_fcc_calibration import build_full_fcc_calibration_model
        value=build_full_fcc_calibration_model(p).evaluate(IDEAL_H,0)
        relative_embedding=-p.embedding_sqrt_ev+p.embedding_linear_ev
        cancellation=(abs(value.pair.value)+abs(relative_embedding))/abs(y[1])
        parameters.append({"set":name,**asdict(p),"a0_reduced_L0":model.a0,
            "equilibrium_lattice_angstrom":stretch*geometry.lattice_angstrom,
            "isotropic_lattice_scale":stretch,
            "atomic_cell_area_angstrom2":area*stretch**2,
            "a0_target":IDEAL_H,"normal_force_eV":y[0],
            "pair_energy_ev_atom":value.pair.value,
            "embedding_relative_atomized_ev":relative_embedding,
            "cancellation_ratio":cancellation,
            "zero_pressure_interpretation":abs(y[0])<1e-5,
            **{f"target_geometry_{key}":v for key,v in cubic_constants_gpa(y).items()},
            **cubic_constants_gpa(relaxed,volume_scale=stretch**3),
            "relaxed_cohesion_ev_atom":relaxed[1],
            "squared_loss":float(np.sum(((y-target)/scales)**2)),
            "production_valid":False})
        for obs,pred,ref,scale in zip(NAMES,y,target,scales):
            residuals.append({"set":name,"observable":obs,"prediction":pred,
                "reference":ref,"absolute_error":abs(pred-ref),
                "relative_error":None if ref==0 else (pred-ref)/ref,
                "normalized_residual":(pred-ref)/scale,"used_for_fit":name!="legacy_unidentified"})
        for step in (4e-4,2e-4,1e-4):
            jac=independent_sensitivity(p,extended=p.embedding_linear_ev>0,log_step=step)
            _,singular,vt=np.linalg.svd(jac,full_matrices=False)
            for mode,sv in enumerate(singular):
                svds.append({"set":name,"log_step":step,"mode":mode+1,
                    "singular_value":sv,"condition_number":singular[0]/singular[-1],
                    **{f"right_vector_{i}":x for i,x in enumerate(vt[mode])}})
            if step==2e-4:
                for i,obs in enumerate(NAMES):
                    for j,v in enumerate(jac[i]):
                        sensitivities.append({"set":name,"observable":obs,
                            "log_parameter_index":j,"sensitivity":v})
        for step in (1.6e-3,8e-4,4e-4):
            for obs,v in zip(NAMES,independent_observables(p,strain_step=step)):
                convergence.append({"set":name,"strain_step":step,"observable":obs,"value":v})
    for filename,rows in (("parameter_sets.csv",parameters),("fit_residuals.csv",residuals),
        ("sensitivity_matrix.csv",sensitivities),("identifiability_svd.csv",svds),
        ("derivative_refinement.csv",convergence)):
        write_csv(out/filename,rows)
    print("bulk fits, units, independent-mode audits saved",flush=True)

    heldout=[]; gsf=[]; opening=[]; hessians=[]; refinements=[]; barriers=[]; grids=[]
    for name in ("linear_profile","linear_historical_box","legacy_unidentified"):
        bulk=models[name]
        area=bulk.geometry.atomic_cell_area*(length_m/1e-10)**2
        direct=FCC111ActiveInterface(bulk,tolerance=2e-10)
        work_sep=direct.separated_limit()
        for path in (DIRECT_110,SHOCKLEY_112):
            interface=FCC111ActiveInterface(bulk,path_id=path,tolerance=2e-10)
            endpoint=interface.period if path==DIRECT_110 else bulk.p.b/np.sqrt(3)
            energies=[]
            for s in np.linspace(0,endpoint,61):
                value=interface.energy(interface.h,float(s)); energies.append(value)
                gsf.append({"set":name,"path":path,"s_over_b":s/bulk.p.b,
                    "s_angstrom":s*length_m/1e-10,"energy_ev_cell":value,
                    "energy_J_m2":cell_ev_to_energy_density(value,area),
                    "relaxation":"rigid half crystals; fixed interface spacing"})
            heldout.append({"set":name,"observable":f"{path}_usf",
                "prediction_J_m2":cell_ev_to_energy_density(max(energies),area),
                "reference_J_m2":.250 if path==DIRECT_110 else .224,
                "used_for_fit":False,"source":"doi:10.1103/PhysRevB.62.3099",
                "comparison_status":"relaxed DFT/different lattice vs rigid model; context only"})
            if path==SHOCKLEY_112:
                heldout.append({"set":name,"observable":"intrinsic_fault",
                    "prediction_J_m2":cell_ev_to_energy_density(energies[-1],area),
                    "reference_J_m2":.164,"used_for_fit":False,
                    "source":"doi:10.1103/PhysRevB.62.3099",
                    "comparison_status":"relaxed DFT vs rigid model; context only"})
        for ref,source in ((1.74,"doi:10.1103/PhysRevB.59.3393"),
                           (2.12,"doi:10.1039/C6RA08958E")):
            heldout.append({"set":name,"observable":"work_separation",
                "prediction_J_m2":cell_ev_to_energy_density(work_sep,area),
                "reference_J_m2":ref,"used_for_fit":False,"source":source,
                "comparison_status":"unrelaxed model vs source surface convention"})
        for a in np.r_[np.linspace(direct.h,3*direct.h,41),
                       [10*direct.h,30*direct.h,100*direct.h]]:
            value=direct.evaluate(float(a),0)
            opening.append({"set":name,"a_reduced_L0":a,
                "opening_angstrom":(a-direct.h)*length_m/1e-10,
                "energy_ev_cell":value.energy,"energy_J_m2":cell_ev_to_energy_density(value.energy,area),
                "traction_GPa":value.d_da*EV_J/(length_m*area*1e-20)/1e9,
                "work_separation_limit_J_m2":cell_ev_to_energy_density(work_sep,area)})
        h=direct.hessian(direct.h,0)
        hessians.append({"set":name,"a0_reduced_L0":direct.h,"Waa":h[0,0],"Was":h[0,1],
            "Wss":h[1,1],"eig_min":float(np.min(np.linalg.eigvalsh(h))),
            "normal_force_ev":direct.grad(direct.h,0)[0]})
        refinements.extend(interface_refinement(name,bulk))
        barriers.extend(fixed_registry_opening_barriers(name,direct,length_m,area*1e-20))
        if name=="linear_profile":
            for a in np.linspace(.95*direct.h,2*direct.h,19):
                for s in np.linspace(0,direct.period,25):
                    grids.append({"a_reduced_L0":a,"s_reduced_L0":s,
                                  "a_over_b":a/bulk.p.b,"s_over_b":s/bulk.p.b,
                                  "energy_ev_cell":direct.energy(float(a),float(s))})
        print(f"{name}: static held-outs and independent refinement done",flush=True)
    for row in heldout:
        row["relative_error"]=row["prediction_J_m2"]/row["reference_J_m2"]-1
    for filename,rows in (("heldout_validation.csv",heldout),("gsf_curve.csv",gsf),
        ("opening_curve.csv",opening),("hessian.csv",hessians),
        ("convergence_summary.csv",refinements),("opening_barriers.csv",barriers),
        ("active_interface_energy_grid.csv",grids)):
        write_csv(interface_out/filename,rows)
    fig,axes=plt.subplots(1,3,figsize=(13,4))
    for name in ("linear_profile","linear_historical_box","legacy_unidentified"):
        for ax,path in zip(axes[:2],(DIRECT_110,SHOCKLEY_112)):
            rows=[r for r in gsf if r["set"]==name and r["path"]==path]
            ax.plot([r["s_over_b"] for r in rows],[r["energy_J_m2"] for r in rows],label=name)
        rows=[r for r in opening if r["set"]==name and r["a_reduced_L0"]<3]
        axes[2].plot([r["opening_angstrom"] for r in rows],[r["energy_J_m2"] for r in rows],label=name)
    for ax,title in zip(axes,("Direct <110>","Shockley <112>","Local opening")):
        ax.set_title(title); ax.set_ylabel("J/m2"); ax.grid(alpha=.25)
    axes[0].set_xlabel("s/b"); axes[1].set_xlabel("s/b")
    axes[2].set_xlabel("Opening displacement [angstrom]")
    axes[0].legend(fontsize=6)
    fig.suptitle("STATIC research candidates - not validated Al/PDE inputs")
    fig.tight_layout(); fig.savefig(interface_out/"physical_gate_comparison.png",dpi=160)
    plt.close(fig)
    summary={
        "legacy_rank_upper_bound":4,"legacy_rank5_claim":"rejected",
        "missing_elastic_mode_added":"homogeneous simple shear gamma=s/h111",
        "staged_optimization_actually_executed":True,
        "bulk_fit":"independent targets fitted; high cancellation or search-bound dependence",
        "interface_status":"held-outs not validated; relaxation mapping incomplete",
        "fully_validated_Al_solver":False,"ready_for_production_PDE":False,
        "ready_for_specimen_UI":False,"next_UI_phase_requires_user_confirmation":True,
        "kinetic_calibration":"unavailable; physical seconds/Hz disabled",
        "specimen_spatial_probability_field":"not supplied by a single local P(a,s,t)",
        "A_c_in_atomistic_calibration":False,"elapsed_seconds":time.perf_counter()-started}
    (out/"readiness.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(summary,indent=2),flush=True)


if __name__=="__main__":
    main()
