"""Replay augmented fixed/force-ensemble matrices and physical unit conversions.

Uses saved raw matrices/forces, scipy eigensystems and a positive-block solve.
Stored finite-difference check records are inspected, not recomputed MACE calls.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh, cho_factor, cho_solve


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def main(args):
    if args.output.exists(): raise ValueError('Fresh audit output required')
    full=json.loads((args.augmented/'summary.json').read_text(encoding='utf-8'))
    protocol=json.loads((args.augmented/'protocol.json').read_text(encoding='utf-8'))
    if full.get('complete') is not True: raise ValueError('Completed augmentation required')
    if protocol['geometry_sha256'] != sha(args.geometry): raise ValueError('Geometry hash differs')
    with np.load(args.geometry) as file:
        ref=file['reference'];free=file['free'];lower=file['lower'];upper=file['upper']
        area=float(file['area_A2']);gauge=float(file['gauge_length_A'])
    if area<=0 or gauge<=0 or np.any(lower&upper) or not np.array_equal(free,~(lower|upper)):
        raise ValueError('Invalid declared specimen geometry')
    norm=np.sqrt(float(np.count_nonzero(lower|upper))/4.)
    free_components=np.flatnonzero(np.repeat(free,3));dim=len(free_components)+1
    expected_basis=np.zeros((3*len(free),dim));expected_basis[free_components,np.arange(dim-1)]=1.
    grip_vector=np.zeros(ref.shape);grip_vector[lower,2]=-.5/norm;grip_vector[upper,2]=.5/norm
    expected_basis[:,-1]=grip_vector.ravel()
    metrics={}
    def check(name,actual,expected,tol):
        actual,expected=np.asarray(actual,float),np.asarray(expected,float)
        if actual.shape!=expected.shape or not np.all(np.isfinite(actual)) or not np.all(np.isfinite(expected)):
            raise ValueError('Shape/nonfinite mismatch: '+name)
        error=float(np.max(abs(actual-expected)))
        if error>tol: raise ValueError(name+' mismatch: '+str(error))
        metrics[name]=error
    if [x['state'] for x in full['states']]!=['loading8','return8','loading10']:
        raise ValueError('Exact three-state order required')
    check('protocol_grip_norm',protocol['normalized_grip_coordinate_norm'],norm,1e-12)
    states=[];pins=[]
    for state in full['states']:
        name=state['state'];folder=args.augmented/name;source=args.results/('dense_'+name)
        saved=json.loads((folder/'summary.json').read_text(encoding='utf-8'))
        if saved!=state or not saved.get('complete'): raise ValueError('Aggregate and per-state summaries differ')
        if sha(source/'raw_hessian.npz')!=state['source_raw_hessian_sha256']: raise ValueError('Fixed matrix input changed')
        with np.load(source/'raw_hessian.npz') as file: original={k:file[k].copy() for k in file.files}
        with np.load(folder/'raw_augmented_hessian.npz') as file: augmented={k:file[k].copy() for k in file.files}
        h=augmented['hessian'];basis=augmented['basis'];r=original['positions'];f=original['forces']
        check(name+'_basis',basis,expected_basis,1e-12)
        check(name+'_positions',augmented['positions'],r,1e-12)
        check(name+'_fixed_basis',basis[:,:-1],original['basis'],1e-12)
        check(name+'_fixed_block',h[:-1,:-1],(original['hessian']+original['hessian'].T)/2,1e-12)
        check(name+'_grip_row',h[-1],augmented['grip_autograd_row'],1e-12)
        check(name+'_symmetry',h,h.T,1e-12)
        # Independent physical axial reaction from lower/upper force sums.
        reaction=float((f[lower,2].sum()-f[upper,2].sum())/2)
        gradient=-expected_basis.T@f.ravel();gradient[-1]-=reaction/norm
        check(name+'_gradient',augmented['gradient'],gradient,1e-8)
        check(name+'_matched_force',reaction,state['matched_dead_force_eV_A'],1e-8)
        stress=reaction/area*160.2176634
        check(name+'_matched_stress',stress,state['matched_nominal_stress_GPa'],1e-8)
        values=eigh(h,eigvals_only=True,driver='evr',check_finite=True)
        check(name+'_spectrum',values,augmented['eigenvalues'],1e-9)
        check(name+'_minimum',values[0],state['force_ensemble_minimum_curvature_eV_A2'],1e-9)
        if int(np.sum(values<0))!=state['force_ensemble_negative_modes']: raise ValueError('Negative mode count differs')
        block=h[:-1,:-1];mixed=h[:-1,-1]
        fixed_values=eigh(block,eigvals_only=True,driver='evr',check_finite=True)
        if fixed_values[0]<=0: raise ValueError('Positive free block required for this Schur audit')
        response=cho_solve(cho_factor(block,check_finite=True),mixed,check_finite=True)
        schur=float(h[-1,-1]-mixed@response);physical=schur*norm**2
        extension=float(np.mean(r[upper,2]-ref[upper,2])-np.mean(r[lower,2]-ref[lower,2]))
        tangent=physical/area*160.2176634*(gauge+extension)
        check(name+'_schur',schur,state['schur_curvature_eV_A2'],1e-9)
        check(name+'_physical_stiffness',physical,state['relaxed_axial_stiffness_eV_A2'],1e-8)
        check(name+'_tangent',tangent,state['finite_prism_tangent_GPa'],1e-6)
        if int(np.sum(values<0))!=int(schur<0): raise ValueError('Positive-block Schur inertia does not agree')
        force_checks=state['independent_mixed_column_checks'];energy_checks=state['grip_energy_checks']
        if len(force_checks)!=2 or len(energy_checks)!=2: raise ValueError('Missing recorded derivative controls')
        check(name+'_force_steps',[x['step_A'] for x in force_checks],[2e-4,1e-4],0)
        check(name+'_energy_steps',[x['step_A'] for x in energy_checks],[2e-3,1e-3],0)
        derivative_errors=np.array([x['mixed_column_error_norm_eV_A2'] for x in force_checks])
        if not np.all(np.isfinite(derivative_errors)) or np.min(derivative_errors)<0 or np.max(derivative_errors)>1e-5:
            raise ValueError('Recorded force derivative checks fail')
        check(name+'_energy_error_definition',[x['curvature_eV_A2']-h[-1,-1] for x in energy_checks],
              [x['error_eV_A2'] for x in energy_checks],1e-12)
        states.append(dict(state=name,matched_nominal_stress_GPa=stress,
            fixed_minimum_curvature_eV_A2=float(fixed_values[0]),force_minimum_curvature_eV_A2=float(values[0]),
            computed_negative_modes=int(np.sum(values<0)),schur_curvature_eV_A2=schur,
            finite_prism_tangent_GPa=tangent,free_gradient_max_eV_A=float(np.max(abs(gradient[:-1]))),
            independent_spectrum_max_error_eV_A2=metrics[name+'_spectrum'],
            recorded_mixed_column_max_error_eV_A2=float(derivative_errors.max()),
            recorded_energy_curvature_max_error_eV_A2=max(abs(x['error_eV_A2']) for x in energy_checks)))
        pins.append(dict(state=name,augmented_raw_sha256=sha(folder/'raw_augmented_hessian.npz'),
                         fixed_raw_sha256=sha(source/'raw_hessian.npz'),state_summary_sha256=sha(folder/'summary.json')))
    if full['new_grip_autograd_rows']!=3 or full['reused_fixed_grip_rows']!=3*(dim-1) or full['new_graph_evaluations']!=27:
        raise ValueError('Completed three-state call accounting differs from declared runner schedule')
    report=dict(complete=True,geometry_sha256=sha(args.geometry),aggregate_sha256=sha(args.augmented/'summary.json'),
        runner_sha256=sha(Path(__file__)),sources=pins,states=states,per_comparison_errors=metrics,
        checked_scalar_or_array_comparisons=len(metrics),new_model_calls=0,new_DFT=0,new_MD=0,
        actual_force_differences_recomputed=False,
        method='saved raw forces plus independent Cartesian basis, SciPy evr spectra, Cholesky Schur complement and physical unit conversion',
        scope='each state own matched dead load; approximate stationarity; finite rigid-grip prism; no bulk/material/first-crack/clock certification')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('per_comparison_errors','sources')},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    for name in ('results','augmented','geometry','output'):parser.add_argument('--'+name,type=Path,required=True)
    main(parser.parse_args())
