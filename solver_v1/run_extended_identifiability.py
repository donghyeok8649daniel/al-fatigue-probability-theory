"""Signed-coefficient local sensitivities of the new odd/even research fits.

Analytic coefficient columns plus deterministic radial-range differences.
Log-absolute coordinates retain each coefficient's sign; no log(negative).
This is inverse-curvature information, not confidence in a calibrated material.
"""
import json
import math
import numpy as np

from .even_moment_calibration import even_basis
from .odd_moment_calibration import odd_basis,odd_targets
from .joint_fcc_interface_calibration import JOINT_NAMES
from .run_matched_interface_study import ROOT,save_json
from .run_full_fcc_calibration_audit import write_csv


def main():
    target,scales=odd_targets(); summaries=[]; rows=[]; correlations=[]
    for name,filename,basis in (
        ("odd_I1_I3","odd_moment_fit.json",odd_basis),
        ("even_I1_I2_I3_C","even_moment_convex_signed_free_B_fit.json",
         lambda d,k:even_basis(d,k,convex=True)),
        ("monotone_I1_I2_I3_C","monotone_opening_refined.json",
         lambda d,k:even_basis(d,k,convex=True))):
        if not (ROOT/filename).exists():
            continue
        fit=json.loads((ROOT/filename).read_text(encoding="utf-8"))["best"]
        c=np.asarray(fit["coefficients"]); d=fit["scalar_decay"]; k=fit["angular_decay"]
        if np.any(c==0):
            raise ValueError("zero boundary coefficient requires one-sided, not logarithmic, sensitivity")
        names=["log_abs_"+x for x in fit["coefficient_order"]]+["log_scalar_decay","log_angular_decay"]
        matrices=[]
        for step in (4e-4,2e-4,1e-4):
            coefficient=basis(d,k)*c[None,:]/scales[:,None]
            radial=np.column_stack((
                (basis(d*math.exp(step),k)@c-basis(d*math.exp(-step),k)@c)/(2*step*scales),
                (basis(d,k*math.exp(step))@c-basis(d,k*math.exp(-step))@c)/(2*step*scales)))
            matrix=np.column_stack((coefficient,radial)); matrices.append(matrix)
            u,s,vh=np.linalg.svd(matrix,full_matrices=False)
            rows.extend(dict(candidate=name,log_step=step,index=i,singular_value=float(value),
                condition_number=float(s[0]/s[-1]),rank_relative_1e8=int(np.sum(s>s[0]*1e-8)),
                parameter_count=len(names),scale_status="declared model discrepancy, not measured uncertainty")
                for i,value in enumerate(s))
        covariance=(vh.T/s**2)@vh
        std=np.sqrt(np.diag(covariance)); corr=covariance/std[:,None]/std[None,:]
        for i in range(len(names)):
            for j in range(i+1,len(names)):
                correlations.append(dict(candidate=name,first=names[i],second=names[j],
                    local_inverse_curvature_correlation=float(corr[i,j])))
        write_csv(ROOT/(name+"_sensitivity.csv"),[dict(target=key,**dict(zip(names,row)))
            for key,row in zip((*JOINT_NAMES,"local_normal_curvature_ev_L0sq"),matrices[-1])])
        # Hard equilibrium/cohesion constraints restrict perturbations to the
        # nullspace of the first two rows. Report this tangent, not 10 free
        # parameters as though the equalities were absent.
        _,constraints_s,cvh=np.linalg.svd(matrices[-1][:2],full_matrices=True)
        tangent=cvh[2:].T
        reduced=np.linalg.svd(matrices[-1][2:]@tangent,compute_uv=False)
        summaries.append(dict(candidate=name,coefficient_signs=np.sign(c),
            singular_values=s,condition_number=float(s[0]/s[-1]),
            max_stencil_change=float(np.max(abs(matrices[-1]-matrices[-2]))),
            hard_constraints_used_in_fit=(name!="odd_I1_I3"),
            tangent_covers_two_bulk_equalities_only=True,
            additional_active_opening_inequalities=(name=="monotone_I1_I2_I3_C"),
            exact_constraint_tangent_singular_values=reduced,
            exact_constraint_tangent_condition=float(reduced[0]/reduced[-1]),
            confidence_intervals_available=False,material_accepted=False))
    write_csv(ROOT/"extended_identifiability_svd.csv",rows)
    write_csv(ROOT/"extended_parameter_correlations.csv",correlations)
    save_json(ROOT/"extended_identifiability_summary.json",summaries)
    print(summaries,flush=True)


if __name__=="__main__":
    main()
