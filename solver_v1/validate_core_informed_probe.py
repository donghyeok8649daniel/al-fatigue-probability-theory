"""Excluded-state and local-rank audit of a completed bounded shape probe.

Snapshot is research-only, explicitly not an adopted Al or kinetic model.
Previously inspected source data are labeled retrospective, not blind tests.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from .coordination_screening import CoordinationScreenedInterface
from .frozen_core_force_target import FrozenCoreForceTarget
from .interface_tangent_calibration import TangentCalibrationProblem,impose_tangents,NONNEGATIVE
from .report_current_material_core import restore_case
from .run_current_material_core import ROOT,load_current_material
from .run_source_core_reference import load_source_material,build_source_core
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--probe',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    if args.out.exists(): raise FileExistsError('fresh probe validation output required')
    raw=(args.probe/'completion.json').read_bytes(); data=json.loads(raw)
    best=data.get('best')
    if not data['completed'] or best is None or not best['strictly_positive_LJ']:
        raise ValueError('completed positive-LJ diagnostic candidate required')
    old,_,binding=load_current_material()
    shape=np.asarray(best['shape']);c=np.asarray(best['coefficients'])
    model=CoordinationScreenedInterface(shape,c,law='power',tolerance=2e-12)
    problem=TangentCalibrationProblem(ROOT/'results/fcc111_active_interface/coordination_screening_v19/joint_refinement')
    observations=impose_tangents(problem.observations,problem.source_tangents)
    M=problem.matrix(tuple(shape)); original_M=problem.matrix(tuple(old.screened_shape))
    if np.max(abs(M@c-np.asarray(best['predictions'])[:len(M)]))>2e-8:
        raise ArithmeticError('selected static profile did not replay')
    rows=[]
    for o,pred,before in zip(observations,M@c,original_M@old.coefficients):
        rows.append(dict(observable=o.name,role=o.role,target=o.target,scale=o.scale,units=o.units,
            baseline=before,prediction=pred,normalized_residual=(pred-o.target)/o.scale,
            excluded_from_current_loss=o.role=='heldout',blind_validation=False))
    write_csv(args.out/'static_observable_audit.csv',rows)
    source,tensor,source_binding=load_source_material();checks=[]
    cases=[('R8_fit','escaped_plus_R8r5',0.),('R10_excluded','stable_plus_R10r5',0.),
           ('R16_excluded','stable_plus_R16r5',0.),
           ('R16_perturbed_plus','stable_plus_R16r5',.002),
           ('R16_perturbed_minus','stable_plus_R16r5',-.002)]
    for name,path,amplitude in cases:
        source_path=ROOT/'results/current_material_core_v22/source_reference'/path
        core,field,*_=restore_case(source_path,source,tensor,source_binding,core_builder=build_source_core)
        field=field+amplitude*np.cos(np.arange(field.size)*.37).reshape(field.shape)
        target=FrozenCoreForceTarget(core,field,radius=2.)
        current=target.coefficient_matrix(model);before=target.coefficient_matrix(old)
        error=current['design']@c-current['target']
        previous=before['design']@old.coefficients-before['target']
        checks.append(dict(case=name,rows=len(current['logical_rows']),
            original_RMS=float(np.sqrt(np.mean(previous**2))),new_RMS=float(np.sqrt(np.mean(error**2))),
            maximum_error=float(np.max(abs(error))),source_gradient_RMS=float(np.sqrt(np.mean(current['target']**2))),
            units='eV/L0',included_in_loss=name=='R8_fit',blind_validation=False,
            state_sha256=hashlib.sha256((source_path/'state.csv').read_bytes()).hexdigest()))
        print(f'{name}: source-core RMS {checks[-1]["original_RMS"]:.6g} -> {checks[-1]["new_RMS"]:.6g}',flush=True)
    write_csv(args.out/'excluded_core_force_audit.csv',checks)
    # Derive the Jacobian from ACTUAL saved difference profiles at the last
    # accepted optimizer state, not the lowest incidental difference trial.
    profiles=json.loads((args.probe/'checkpoint.json').read_bytes())['profiles']
    opt=data['optimizer'];rank=dict(available=False)
    if 'last_accepted_coordinates' in opt:
        x=np.asarray(opt['last_accepted_coordinates'])
        base=next((p for p in profiles if np.allclose(np.log(p['shape'][:5]),x,rtol=0,atol=2e-13)),None)
        if base is not None:
            selected=base['selected_rows'];r=np.asarray(base['residuals'])[selected];columns=[];steps=[]
            for j in range(5):
                candidates=[]
                for p in profiles:
                    dx=np.log(p['shape'][:5])-x
                    if abs(dx[j])>1e-7 and np.max(abs(np.delete(dx,j)))<2e-13:
                        candidates.append((abs(dx[j]),dx[j],p))
                if not candidates: break
                _,step,p=min(candidates,key=lambda item:item[0])
                columns.append((np.asarray(p['residuals'])[selected]-r)/step);steps.append(step)
            if len(columns)==5:
                J=np.array(columns).T;sv=np.linalg.svd(J,compute_uv=False)
                normalized=J/np.linalg.norm(J,axis=0)
                rank=dict(available=True,method='saved actual one-sided log-shape profile differences',
                    log_steps=steps,singular_values=sv,condition_number=float(sv[0]/sv[-1]),
                    numerical_rank=int(np.sum(sv>max(J.shape)*np.finfo(float).eps*sv[0])),
                    column_correlations=normalized.T@normalized,
                    structural_identifiability_proved=False,positive_LJ_endpoint=base['strictly_positive_LJ'])
                write_csv(args.out/'shape_sensitivity_matrix.csv',[dict(row=i,**{f'log_shape_{j}':v for j,v in enumerate(row)}) for i,row in enumerate(J)])
    save_json(args.out/'local_identifiability.json',rank)
    # Explicit snapshot format accepted by the EXISTING research runner; the
    # user/default material file is never modified by this validator.
    save_json(args.out/'research_candidate_snapshot.json',dict(completed=True,best=best,
        screening_law='power',source_probe_sha256=hashlib.sha256(raw).hexdigest(),
        selected_profile_is_optimizer_endpoint=data['selected_profile_is_optimizer_endpoint'],
        optimizer=opt,material_accepted=False,for_static_research_only=True,
        production_changed=False,kinetic_calibration=False))
    save_json(args.out/'completion.json',dict(completed=True,current_candidate_sha256=binding['parameter_sha256'],
        source_sha256=source_binding['parameter_sha256'],static_profile_replayed=True,
        actual_excluded_force_evaluations=True,finite_q_scope='unchanged sampled points and analytic tail inequalities only',
        signs_satisfied=bool(np.all(c[list(NONNEGATIVE)]>=0)),
        material_accepted=False,physical_yield_validated=False,physical_Hz=False,
        files_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in args.out.glob('*.csv')}))


if __name__=='__main__': main()
