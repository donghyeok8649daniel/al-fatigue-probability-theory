"""Joint finite-q/interface identifiability of the unchanged LJ/Bessel family."""
import argparse
from dataclasses import asdict
from pathlib import Path
import time
import numpy as np
from scipy.linalg import null_space
from .coordination_screening import CoordinationScreenedBulk
from .core_interface_compatibility import minimax_compatibility
from .interface_tangent_calibration import TangentCalibrationProblem,impose_tangents,NONNEGATIVE
from .periodic_plane_covariance import SourceEAMBlochHessian
from .run_current_material_core import ROOT,load_current_material
from .run_source_core_reference import load_source_material
from .vector_material_calibration import LENGTH_M,MaterialObservation
from .run_vector_registry_audit import save_json
from .run_low_stress_cyclic_diagnostic import write_csv


def points():
    return [('fit',(t,t,0.)) for t in (.125,.375)]+[('fit',(t,t,t)) for t in (.125,.375)]+[
        ('excluded',(t,t,0.)) for t in (.25,.625)]+[('excluded',(t,t,t)) for t in (.25,.45)]


def run(out):
    out=Path(out)
    if out.exists():raise FileExistsError('fresh compatibility audit required')
    started=time.perf_counter();source,_,binding=load_source_material();operator=SourceEAMBlochHessian(source.source)
    problem=TangentCalibrationProblem(ROOT/'results/fcc111_active_interface/coordination_screening_v19/joint_refinement')
    baseobs=impose_tangents(problem.observations,problem.source_tangents)
    models=[('original',None),('wider',ROOT/'results/current_material_core_v22/wider_probe_validation/research_candidate_snapshot.json'),
            ('radial',ROOT/'results/core_interface_compatibility_v23/radial_full_validation/joint_LS/research_candidate_snapshot.json')]
    out.mkdir(parents=True)
    save_json(out/'definition.json',dict(points=points(),source=binding,coefficient_signs=NONNEGATIVE,
        full_symmetric_matrix_features='H00,H11,H22,sqrt(2)H01,sqrt(2)H02,sqrt(2)H12',
        component_scale='.05*Frobenius norm of source matrix; discrepancy choice, not experimental uncertainty',
        units='eV/L0^2',exact_anchors=7,old_inspected_jets='development, not blind heldout',
        production_changed=False,physical_yield_used=False,atomic_mass_used=False))
    rows=[];reports=[];components=[(0,0),(1,1),(2,2),(0,1),(0,2),(1,2)]
    for name,path in models:
        model,_,modelbinding=load_current_material() if path is None else load_current_material(path)
        base=problem.matrix(tuple(model.screened_shape));obs=list(baseobs);extra=[];tail_rows=[];targets=[]
        bulk=CoordinationScreenedBulk(model.screened_shape,radius=16.,law='power')
        bulk12=CoordinationScreenedBulk(model.screened_shape,radius=12.,law='power')
        for k,(role,point) in enumerate(points()):
            target=operator.evaluate(np.asarray(point)*2*np.pi/source.source.geometry.lattice_constant)*(LENGTH_M/1e-10)**2
            columns,tails=bulk.evaluate(point);low,_=bulk12.evaluate(point);targets.append(target)
            for i,j in components:
                w=1. if i==j else np.sqrt(2.)
                extra.append(w*columns[:,i,j]);obs.append(MaterialObservation(
                    f'v31_q{k}_H{i}{j}',float(w*target[i,j]),float(.05*np.linalg.norm(target)),
                    'eV/L0^2','heldout' if role=='excluded' else role))
            tail_rows.append(dict(point=point,coefficient_tail=tails,columns_radius_change=np.linalg.norm(columns-low,axis=(1,2))))
        M=np.vstack([base,extra]);exact=[i for i,o in enumerate(obs) if o.role=='exact']
        qfit=[i for i in range(len(base),len(obs)) if obs[i].role=='fit']
        old=[i for i,o in enumerate(baseobs) if o.role!='exact' and o.units=='eV/L0^2']
        scales=np.array([o.scale for o in obs]);D=1/np.linalg.norm(M/scales[:,None],axis=0)
        N=null_space(M[exact]*D/scales[exact,None]);S=(M[qfit]*D/scales[qfit,None])@N
        sv=np.linalg.svd(S,compute_uv=False)
        result=dict(model=name,binding=modelbinding,exact_rank=len(D)-N.shape[1],null_dimension=N.shape[1],
                    finite_q_null_singular_values=sv,coefficient_scale=D,tail_checks=tail_rows,profiles=[])
        for label,selected in [('finite_q',qfit),('interface',old),('joint',old+qfit)]:
            fit=minimax_compatibility(M,obs,selected,nonnegative=NONNEGATIVE)
            result['profiles'].append(dict(label=label,**fit))
            if fit['completed']:
                coeff=fit['coefficients']
                for k,(role,point) in enumerate(points()):
                    predicted=np.array(fit['predictions'][len(base)+6*k:len(base)+6*(k+1)])
                    target=np.array([obs[len(base)+6*k+j].target for j in range(6)])
                    rows.append(dict(model=name,profile=label,point_index=k,role=role,q=point,
                        matrix_relative_error=float(np.linalg.norm(predicted-target)/np.linalg.norm(target)),
                        minimax=fit['minimax_normalized_error'],strictly_positive_LJ=fit['strictly_positive_LJ'],
                        tail_bound=float(tail_rows[k]['coefficient_tail']@abs(coeff)),
                        radius_change_bound=float(tail_rows[k]['columns_radius_change']@abs(coeff)),
                        material_accepted=False))
            print(name,label,fit.get('minimax_normalized_error'),flush=True)
        result['observations']=[asdict(o) for o in obs]
        save_json(out/(name+'.json'),result);reports.append(result)
    write_csv(out/'comparison.csv',rows)
    save_json(out/'summary.json',dict(completed=True,elapsed_seconds=time.perf_counter()-started,
        models=[dict(model=r['model'],exact_rank=r['exact_rank'],null_dimension=r['null_dimension'],
                     singular_values=r['finite_q_null_singular_values'],
                     profiles=[dict(label=f['label'],completed=f['completed'],eta=f.get('minimax_normalized_error'),
                                    positive_LJ=f.get('strictly_positive_LJ')) for f in r['profiles']]) for r in reports],
        family_impossibility_proved=False,material_accepted=False,actual_yield_validated=False,production_changed=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(__doc__);p.add_argument('--out',type=Path,required=True);run(p.parse_args().out)
