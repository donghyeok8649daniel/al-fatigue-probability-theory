"""Full replay and excluded-q validation of a completed v31 radial search."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.linalg import null_space
from .run_finite_q_compatibility_v31 import points
from .coordination_screening import CoordinationScreenedBulk
from .interface_tangent_calibration import TangentCalibrationProblem,impose_tangents
from .periodic_plane_covariance import SourceEAMBlochHessian
from .run_current_material_core import ROOT
from .run_source_core_reference import load_source_material
from .vector_material_calibration import LENGTH_M
from .run_vector_registry_audit import save_json
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_low_frequency_forcing_v29 import sha


def run(study,out):
    study,out=Path(study),Path(out)
    if out.exists():raise FileExistsError('fresh validation directory required')
    report=json.loads((study/'summary.json').read_bytes())
    if not report['completed'] or not report['best']:raise ValueError('completed search required')
    best=report['best'];shape=np.asarray(best['shape']);coeff=np.asarray(best['coefficients'])
    problem=TangentCalibrationProblem(ROOT/'results/fcc111_active_interface/coordination_screening_v19/joint_refinement')
    obs=impose_tangents(problem.observations,problem.source_tangents)
    base=problem.matrix(tuple(shape));pred=base@coeff
    selected=[i for i,o in enumerate(obs) if o.role!='exact' and o.units=='eV/L0^2']
    exact=[i for i,o in enumerate(obs) if o.role=='exact']
    residuals=[dict(name=o.name,target=o.target,prediction=float(p),scale=o.scale,
        normalized_residual=float((p-o.target)/o.scale),units=o.units,
        used_in_objective=i in selected,exact_constraint=i in exact,
        prior_inspected=True) for i,(o,p) in enumerate(zip(obs,pred))]
    source,_,binding=load_source_material();operator=SourceEAMBlochHessian(source.source)
    bulk=CoordinationScreenedBulk(shape,radius=16.,law='power')
    low=CoordinationScreenedBulk(shape,radius=12.,law='power')
    components=[(0,0),(1,1),(2,2),(0,1),(0,2),(1,2)]
    rows=[];extra=[];scales=[];qnorm=[]
    for index,(role,q) in enumerate(points()):
        reference=operator.evaluate(np.asarray(q)*2*np.pi/source.source.geometry.lattice_constant)*(LENGTH_M/1e-10)**2
        columns,tails=bulk.evaluate(q);columns12,_=low.evaluate(q)
        H=np.einsum('c,cij->ij',coeff,columns);scale=.05*np.linalg.norm(reference)
        error=H-reference
        features=np.array([error[i,j]*(1 if i==j else np.sqrt(2)) for i,j in components])/scale
        if role=='fit':
            extra.extend(columns[:,i,j]*(1 if i==j else np.sqrt(2)) for i,j in components)
            scales.extend([scale]*6);qnorm.extend(features)
        rows.append(dict(point_index=index,role=role,q=q,
            source_matrix_eV_L02=reference,predicted_matrix_eV_L02=H,
            relative_Frobenius_error=float(np.linalg.norm(error)/np.linalg.norm(reference)),
            max_normalized_component_error=float(np.max(abs(features))),
            minimum_eigenvalue=float(np.linalg.eigvalsh(H).min()),
            tail_bound=float(tails@abs(coeff)),
            radius_change_bound=float(np.linalg.norm(columns-columns12,axis=(1,2))@abs(coeff)),
            production_admissible=False))
    eta=max(max(abs(residuals[i]['normalized_residual']) for i in selected),max(abs(np.array(qnorm))))
    if abs(eta-best['value'])>2e-5:raise ArithmeticError('full matrix replay differs from search objective')
    target_scales=np.array([o.scale for o in obs])
    E=base[exact]/target_scales[exact,None]
    J=np.vstack([base[selected]/target_scales[selected,None],np.array(extra)/np.array(scales)[:,None]])
    D=1/np.linalg.norm(np.vstack([E,J]),axis=0)
    N=null_space(E*D);singular=np.linalg.svd((J*D)@N,compute_uv=False)
    out.mkdir(parents=True);write_csv(out/'target_residuals.csv',residuals)
    # Matrices stay in JSON; scalar diagnostics remain directly usable CSV.
    write_csv(out/'excluded_q_validation.csv',[{k:v for k,v in r.items() if not k.endswith('_eV_L02')} for r in rows])
    save_json(out/'summary.json',dict(completed=True,search_sha256=sha(study/'summary.json'),
        source=binding,shape=shape,coefficients=coeff,replayed_eta=eta,
        maximum_exact_normalized_error=max(abs(residuals[i]['normalized_residual']) for i in exact),
        conditional_joint_singular_values=singular,remaining_coefficient_directions=N.shape[1],
        q_validation=rows,optimizer_success=report['optimizer_success'],
        material_accepted=False,full_Brillouin_stability_certified=False,
        actual_yield_validated=False,production_changed=False))
    print(json.dumps(dict(eta=eta,excluded_max_relative=max(r['relative_Frobenius_error'] for r in rows if r['role']=='excluded'),
        singular_values=singular.tolist(),minimum_tested_eigenvalue=min(r['minimum_eigenvalue'] for r in rows))))


if __name__=='__main__':
    p=argparse.ArgumentParser(__doc__);p.add_argument('--study',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();run(a.study,a.out)
