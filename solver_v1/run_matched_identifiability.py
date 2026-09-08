"""Local identifiability and stencil refinement of every evaluated extension."""
import json
import math
import numpy as np

from .angular_interface_calibration import angular_observation_column
from .joint_fcc_interface_calibration import joint_basis, matched_targets, JOINT_NAMES
from .run_full_fcc_calibration_audit import write_csv
from .run_matched_interface_study import ROOT, save_json


def main():
    data=json.loads((ROOT/"fitted_candidates.json").read_text(encoding="utf-8"))
    extra=ROOT/"angular_quadratic_extension.json"
    if extra.exists():
        data["angular_quadratic"]=json.loads(extra.read_text(encoding="utf-8"))["best"]
    _,scales=matched_targets(); rows=[]; correlations=[]; parameters=[]
    for name in ("linear","angular_tied","angular_independent_range","angular_quadratic"):
        if name not in data:
            continue
        fit=data[name]; angular=name!="linear"; quadratic=fit.get("quadratic",False)
        count=5 if quadratic else 4
        coeff=np.asarray(fit["coefficients"])
        scalar=coeff[:count]
        if np.any(scalar<=0):
            raise ValueError("boundary coefficients need one-sided analysis, not logarithms")
        d=fit.get("scalar_decay",fit.get("decay"))
        amplitude=float(coeff[-1]) if angular else 0.
        k=fit.get("angular_decay",d)
        tied=name=="angular_tied"
        keys=["log_u","log_v","log_A","log_B"]+(["log_C"] if quadratic else [])
        logvalues=list(np.log(scalar))
        if angular:
            keys.append("log_D"); logvalues.append(math.log(amplitude))
        keys.append("log_scalar_decay"); logvalues.append(math.log(d))
        if angular and not tied:
            keys.append("log_angular_decay"); logvalues.append(math.log(k))
        x=np.asarray(logvalues)
        def predict(logparams):
            values=np.exp(logparams)
            c=values[:count]; index=count
            D=values[index] if angular else 0.
            index+=int(angular)
            scalar_decay=values[index]; index+=1
            angular_decay=scalar_decay if tied else values[index] if angular else scalar_decay
            y=joint_basis(float(scalar_decay))[:,:count]@c
            if angular:
                y=y+D*angular_observation_column(float(angular_decay))
            return y
        matrices=[]
        for step in (4e-4,2e-4,1e-4):
            jac=[]
            for i in range(len(x)):
                dx=np.zeros(len(x)); dx[i]=step
                jac.append((predict(x+dx)-predict(x-dx))/(2*step*scales))
            matrix=np.column_stack(jac); matrices.append(matrix)
            u,s,vh=np.linalg.svd(matrix,full_matrices=False)
            for index,value in enumerate(s):
                rows.append(dict(candidate=name,log_step=step,index=index,singular_value=float(value),
                    condition_number=float(s[0]/s[-1]),rank_at_relative_1e8=int(np.sum(s>s[0]*1e-8)),
                    parameter_count=len(keys),confidence_status="local inverse curvature, not statistical confidence"))
        matrix=matrices[-1]; u,s,vh=np.linalg.svd(matrix,full_matrices=False)
        covariance=(vh.T/s**2)@vh
        diagonal=np.sqrt(np.diag(covariance))
        correlation=covariance/diagonal[:,None]/diagonal[None,:]
        for i in range(len(keys)):
            for j in range(i+1,len(keys)):
                correlations.append(dict(candidate=name,first=keys[i],second=keys[j],
                    local_inverse_curvature_correlation=float(correlation[i,j])))
        write_csv(ROOT/(name+"_sensitivity.csv"),[
            dict(target=target,**{key:float(value) for key,value in zip(keys,row)})
            for target,row in zip(JOINT_NAMES,matrix)])
        eps=scalar[1]**2/(4*scalar[0]); sigma=(scalar[0]/scalar[1])**(1/6)
        parameters.append(dict(candidate=name,epsilon_LJ_eV=float(eps),sigma_LJ_over_L0=float(sigma),
            A_eV=float(scalar[2]),B_eV=float(scalar[3]),C_eV=float(scalar[4]) if quadratic else 0.,
            scalar_decay_L0=float(d),angular_D_eV=amplitude,angular_decay_L0=float(k) if angular else None,
            L0_angstrom=4.05/math.sqrt(2),density_gauge="fixed target-FCC density",
            Jacobian_stencil_absolute_change=float(np.max(abs(matrices[-1]-matrices[-2]))),
            condition_number=float(s[0]/s[-1]),status="research; not accepted Al",
            boundary_note="scalar decay at declared lower search bound" if quadratic and d<.501 else "none detected in this record"))
    write_csv(ROOT/"all_identifiability_svd.csv",rows)
    write_csv(ROOT/"parameter_correlations.csv",correlations)
    write_csv(ROOT/"physical_parameter_sets.csv",parameters)
    save_json(ROOT/"identifiability_status.json",dict(
        observable_rank_is_not_material_validation=True,
        intervals_are_not_statistical_confidence=True,
        model_discrepancy_scales_not_experimental_uncertainties=True,
        numerical_stencil_refinement_performed=True,
        quadratic_boundary_sensitive=True))
    print(parameters)


if __name__=="__main__":
    main()
