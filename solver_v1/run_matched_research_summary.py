"""Auditable summaries of SAVED fits plus fresh exact separation evaluation.

This does not re-optimize and does not certify an aluminum potential. Sources,
normalizations, numerical boundary solutions, and material failure are kept
separate. No default calibration or probability result is overwritten.
"""
import csv
import hashlib
import json
import numpy as np

from .fcc111_active_interface import FCC111ActiveInterface
from .joint_fcc_interface_calibration import JOINT_NAMES
from .run_matched_interface_study import ROOT,save_json,surfaces_from_results
from .run_full_fcc_calibration_audit import write_csv
from .reference_eam_targets import SOURCE_URL,SOURCE_SHA256


def read_csv(path):
    with path.open(encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def target_units(index):
    if index==0:
        return "eV/atom per unit normal strain"
    if index==1:
        return "eV/atom"
    if index<5:
        return "eV/atom per squared unit strain"
    if index<10:
        return "eV/interface_cell"
    return "eV/interface_cell per squared reduced length"


def main():
    targets=read_csv(ROOT/"matched_targets.csv")
    for i,row in enumerate(targets):
        row["units"]=target_units(i)
    write_csv(ROOT/"matched_targets.csv",targets)
    target=np.array([float(r["value"]) for r in targets])
    scales=np.array([float(r["scale"]) for r in targets])
    data=json.loads((ROOT/"fitted_candidates.json").read_text(encoding="utf-8"))
    files=("angular_quadratic_extension","odd_moment_fit","odd_moment_quadratic_probe",
           "constrained_odd_fit","even_moment_fit","even_moment_convex_signed_free_B_fit",
           "monotone_opening_fit","monotone_opening_refined")
    for name in files:
        path=ROOT/(name+".json")
        if path.exists():
            record=json.loads(path.read_text(encoding="utf-8"))
            data[name]=record.get("best_positivity_audit",record["best"])
    normal=json.loads((ROOT/"odd_moment_fit.json").read_text(encoding="utf-8"))
    extended_target=np.r_[target,normal["target"][-1]]
    extended_scales=np.r_[scales,normal["scales"][-1]]
    names=(*JOINT_NAMES,"local_normal_interface_curvature")
    residuals=[]; overview=[]
    for name,fit in data.items():
        prediction=np.asarray(fit["predictions"])
        y=extended_target[:len(prediction)]; scale=extended_scales[:len(prediction)]
        residual=(prediction-y)/scale
        np.testing.assert_allclose(residual,fit["residuals"],atol=2e-8,rtol=2e-9)
        coefficients=np.asarray(fit["coefficients"])
        for i,(p,reference,r) in enumerate(zip(prediction,y,residual)):
            residuals.append(dict(candidate=name,target=names[i],reference=reference,prediction=p,
                units=target_units(i),absolute_error=p-reference,
                relative_error=(p-reference)/reference if reference else None,
                normalization_scale=scale[i],normalized_residual=r,used_for_fit=True,
                geometry="fixed target FCC; not an error silently relaxed away",
                source_doi="10.1103/PhysRevB.59.3393",source_url=SOURCE_URL,
                reference_method="0 K EAM rigid interface / affine bulk; not experiment"))
        overview.append(dict(candidate=name,observation_count=len(y),
            squared_loss=float(residual@residual),max_absolute_normalized_residual=float(max(abs(residual))),
            each_fit_residual_within_one_declared_scale=bool(max(abs(residual))<=1),
            strictly_positive_LJ_resolved=fit.get("strictly_positive_LJ_resolved",None),
            u=coefficients[0],v=coefficients[1],physical_candidate_accepted=False,
            status="research only; actual branch/held-out/mobility validation separate"))
    write_csv(ROOT/"all_fit_residuals.csv",residuals)
    write_csv(ROOT/"nested_fit_overview.csv",overview)
    comparisons=[]; extrema=[]; unload=[]; refinements=[]
    for relative in (".","odd_moment","even_moment","monotone_opening"):
        folder=ROOT/relative
        if not (folder/"scenario_execution.json").exists():
            continue
        surfaces=surfaces_from_results(folder)
        bulk=read_csv(folder/"bulk_scenario_summary.csv")
        registry=read_csv(folder/"registry_curves.csv")
        opening=read_csv(folder/"opening_curves.csv")
        stress=read_csv(folder/"stress_scenarios.csv")
        local=read_csv(folder/"local_interface_stiffness.csv")
        for row in bulk:
            name=row["candidate"]
            surface=next(s for s in surfaces if s["name"]==name and s["path"]=="direct_110")
            face=FCC111ActiveInterface(surface["bulk"],tolerance=2e-11)
            separated=face.separated_limit()
            bound=getattr(surface["evaluate"],"__self__",None)
            if bound is not None:
                separated+=bound.amplitude_ev*bound.angular.evaluate(400*face.h,0.)[0][0]
                if bound.vector is not None:
                    separated+=bound.vector_amplitude_ev*bound.vector.evaluate(400*face.h,0.)[0][0]
                if bound.quadrupole is not None:
                    separated+=bound.quadrupole_amplitude_ev*bound.quadrupole.evaluate(400*face.h,0.)[0][0]
            gsf=[r for r in registry if r["candidate"]==name]
            normal=[r for r in opening if r["candidate"]==name]
            values=[r for r in local if r["candidate"]==name][0]
            comparisons.append(dict(result_folder=relative,**row,
                normal_interface_tangent_GPa=values["normal_tangent_GPa"],
                registry_interface_tangent_GPa=values["registry_tangent_GPa"],
                direct_max_J_m2=max(float(r["energy_J_m2"]) for r in gsf if r["path"]=="direct_110"),
                shockley_max_J_m2=max(float(r["energy_J_m2"]) for r in gsf if r["path"]=="shockley_112"),
                intrinsic_J_m2=[float(r["energy_J_m2"]) for r in gsf if r["path"]=="shockley_112"][-1],
                separation_J_m2=float(surface["units"].energy_to_surface(separated)),
                separation_method="analytic LJ limit + exponentially saturated per-site environment",
                minimum_sampled_opening_traction_GPa=min(float(r["normal_traction_GPa"]) for r in normal if float(r["a_over_h"])>1.),
                material_adopted=False))
            unload.extend(dict(result_folder=relative,**r) for r in stress
                          if r["candidate"]==name and r["protocol"]=="quasistatic_unload")
        # Make historical source-coordinate L0 explicit in every barrier row.
        path=folder/"constrained_opening_barriers.csv"
        bars=read_csv(path)
        for row in bars:
            scale=next(r for r in local if r["candidate"]==row["set"])
            row["coordinate_length_scale_m"]=scale["coordinate_length_scale_m"]
            row["atomic_cell_area_m2"]=scale["atomic_cell_area_m2"]
        write_csv(path,bars)
    write_csv(ROOT/"physical_model_comparison.csv",comparisons)
    write_csv(ROOT/"all_quasistatic_unload.csv",unload)
    provenance={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()
                for path in sorted(ROOT.rglob("*_fit.json"))}
    save_json(ROOT/"physical_readiness.json",dict(
        status="not accepted for quantitative Al or production interface PDE",
        mathematical_research_extensions_implemented=True,
        remaining_gates=["large independent elastic mismatch","parameter sensitivity",
                        "held-out curve mismatch","broader material/kinetic validation absent"],
        historical_failures_not_hidden=["large pair/embedding cancellation in prior candidates",
                                       "unconstrained intermediate opening overshoot"],
        monotone_constraints_are_not_global_proof=True,
        all_family_impossibility_proven=False,production_defaults_changed=False,
        kinetic_mobility_calibrated=False,physical_seconds_enabled=False,physical_Hz_enabled=False,
        specimen_correlation_area_calibrated=False,source_sha256=SOURCE_SHA256,
        fit_input_hashes=provenance))
    print("Saved residuals, physical comparisons, provenance and non-adoption status.",flush=True)


if __name__=="__main__":
    main()
