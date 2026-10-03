"""Independent raw prediction/Gram/source/hash verification; no geometry fitting."""
import argparse,hashlib,json,re
from pathlib import Path
from collections import defaultdict
import numpy as np


def read(path):return json.loads(path.read_text(encoding='utf-8'))


def main(args):
    root=args.package;catalog=read(root/'raw_cache_catalog.json');verified=0
    for record in catalog:
        if record['exported_path'] is None:continue
        path=root/record['exported_path'];raw=path.read_bytes()
        if len(raw)!=record['bytes'] or hashlib.sha256(raw).hexdigest()!=record['sha256']:
            raise RuntimeError('exported raw bytes differ: '+record['cache_relative_path'])
        verified+=1
    source_count=0
    for manifest in (root/'raw_runs').glob('*/source_manifest.json'):
        folder=manifest.parent
        for binding in read(manifest):
            if not binding['path'].endswith('.py') and not binding['path'].endswith('.sw'):continue
            path=folder/'source_snapshots'/binding['path']
            if not path.exists():
                # The dual repair kept its one changed helper under this name.
                alternative=folder/'calculation_fit_source.py'
                if alternative.exists() and hashlib.sha256(alternative.read_bytes()).hexdigest()==binding['sha256']:
                    source_count+=1;continue
                if args.repo is not None:
                    current=args.repo/binding['path']
                    if current.exists() and hashlib.sha256(current.read_bytes()).hexdigest()==binding['sha256']:
                        source_count+=1;continue
                raise RuntimeError('calculation source is unavailable: '+str(path))
            if hashlib.sha256(path.read_bytes()).hexdigest()!=binding['sha256']:
                raise RuntimeError('calculation source hash differs: '+str(path))
            source_count+=1
    predictions=np.load(root/'selected_predictions.npz');meta=read(root/'frame_metadata.json');cases=read(root/'prediction_cases.json')
    offsets=predictions['offsets'];reference=predictions['reference_force'];reference_e=predictions['reference_energy'];forces=predictions['model_force'];energy=predictions['model_energy']
    if len(meta)!=2475 or len(reference)!=171815 or len(cases['cases'])!=forces.shape[0]:raise RuntimeError('raw prediction sizes differ')
    if [r['frame'] for r in meta]!=list(range(len(meta))) or any(offsets[i+1]-offsets[i]!=r['atoms'] for i,r in enumerate(meta)):
        raise RuntimeError('frame metadata and atom offsets differ')
    grouped=read(root/'reaggregated_group_metrics.json');maximum_force_error=0.;maximum_energy_error=0.
    for row in grouped:
        case=next(i for i,r in enumerate(cases['cases']) if r['label']==row['case'])
        members=[r for r in meta if (r['config_type'],r['declared_xc'],r['training'])==(row['config_type'],row['declared_xc'],row['training'])]
        selected_indices=np.concatenate([np.arange(offsets[r['frame']],offsets[r['frame']+1]) for r in members])
        force_rmse=float(np.linalg.norm((forces[case]-reference)[selected_indices].ravel())/np.sqrt(3*len(selected_indices)))
        maximum_force_error=max(maximum_force_error,abs(force_rmse-row['force_RMSE_eV_A']))
        if row['relative_energy_RMSE_eV_atom'] is not None:
            indices=np.array([r['frame'] for r in members]);atoms=np.array([r['atoms'] for r in members]);anchor=cases['anchor_frame'];an=meta[anchor]['atoms']
            residual=(energy[case,indices]-reference_e[indices])/atoms-(energy[case,anchor]-reference_e[anchor])/an
            rmse=float(np.linalg.norm(residual)/np.sqrt(len(members)))
            maximum_energy_error=max(maximum_energy_error,abs(rmse-row['relative_energy_RMSE_eV_atom']))
    if maximum_force_error>2e-12 or maximum_energy_error>2e-11:raise RuntimeError('independent group aggregation differs')
    if args.lammps_control is not None:
        control=np.load(args.lammps_control)
        if not np.array_equal(control['offsets'],offsets) or not np.array_equal(control['reference_force'],reference):raise RuntimeError('reference archive control differs')
        lmp_force=float(np.max(abs(control['original_sw_1985_force']-forces[0])))
        lmp_energy=float(np.max(abs(control['original_sw_1985_energy']-energy[0])))
        if lmp_force>1e-9 or lmp_energy>1e-9:raise RuntimeError('archived independent LAMMPS SW control differs')
    else:lmp_force=lmp_energy=None
    replayed=0;max_normalized_loss_error=0.
    for record in read(root/'design_gram_replay.json'):
        rows=read(root/'raw_runs'/record['run']/'profiles.json');index=int(re.search(r'profile_(\d+)',record['raw_file']).group(1));row=next(r for r in rows if r['profile']==index)
        fit=row.get('fit',row)
        if record.get('beta_quadratic'):
            if 'selected' not in row:continue
            beta=row['selected']['beta'];v=np.array([1.,beta,beta**2]);alpha=np.array(row['amplitudes'])
            gram=np.einsum('i,j,ijab->ab',v,v,np.array(record['xTx']));linear=v@np.array(record['xTy'])
            bg=np.einsum('i,j,ijab->ab',v,v,np.array(record['bTb']));bl=v@np.array(record['bTt'])
            expected=row['selected'];eos=v@np.array(record['EOS_error_basis']);extra=float(np.mean((eos/.02)**2))
        else:
            if 'amplitudes' not in fit:continue
            alpha=np.array(fit['amplitudes']);gram=np.array(record['xTx']);linear=np.array(record['xTy']);bg=np.array(record['bTb']);bl=np.array(record['bTt']);expected=fit
            extra=(row['EOS_RMSE_eV_atom']/.02)**2 if 'training_objective_squared' in row else 0.
        loss=float(alpha@gram@alpha-2*alpha@linear+record['yTy']+extra)
        bulk=float(alpha@bg@alpha-2*alpha@bl+record['tTt'])
        target=expected['objective_squared'];terms=abs(float(alpha@gram@alpha))+2*abs(float(alpha@linear))+abs(record['yTy'])+abs(extra)
        tolerance=512*np.finfo(float).eps*max(1.,terms)
        error=abs(loss-target);max_normalized_loss_error=max(max_normalized_loss_error,error/max(1.,terms))
        if error>max(2e-8,tolerance):raise RuntimeError('compact loss replay differs: '+record['run']+'/'+record['raw_file'])
        if abs(bulk-expected['bulk_energy_squared'])>max(2e-9,512*np.finfo(float).eps*(abs(float(alpha@bg@alpha))+2*abs(float(alpha@bl))+abs(record['tTt']))):
            raise RuntimeError('compact bulk guard replay differs')
        replayed+=1
    for prefix in ['angular_results','angular_shape_results','extended_fit_results','angle_preference_fit_results',
                   'joint_angle_fit_results','angular_boundary_fit_results','refined_boundary_fit_results',
                   'local_boundary_fit_results','guard_boundary_fit_results']:
        path=root/'raw_runs'/prefix/'bessel_audit.json'
        records=read(path)
        if len(records)!=8:raise RuntimeError('Bessel audit has missing states')
        for r in records:
            b,d=r['bessel'],r['direct'];e=abs(b['energy']-d['energy']);g=float(np.max(abs(np.array(b['gradient'])-d['gradient'])));h=float(np.max(abs(np.array(b['hessian'])-d['hessian'])))
            if not (e<2e-8 and g<2e-7 and h<3e-6):raise RuntimeError('recorded Bessel audit failed')
    final_branch=read(root/'raw_runs/final_selected_interface_results/branch_roots.json')
    branch=read(root/'raw_runs/interface_recovery_results/branch_roots.json')+read(root/'raw_runs/new_candidate_interface_results/branch_roots.json')+final_branch
    for row in branch:
        si_traction=row['jet']['gradient_eV_A'][0]*1.602176634e-19/1e-10/(row['atomic_cell_area_A2']*1e-20)/1e9
        if abs(si_traction-row['normal_traction_GPa'])>1e-10:raise RuntimeError('independent SI traction conversion differs')
    checks=read(root/'raw_runs/final_selected_interface_results/finite_difference_checks.json')
    maximum_fd=max(c['maximum_hessian_error_eV_A2'] for r in checks for c in r['checks'])
    if maximum_fd>1e-8:raise RuntimeError('recorded independent force-Hessian check failed')
    new_branch=final_branch
    root_mesh_error=max(abs(next(r['opening_A'] for r in new_branch if r['cut_kind']==cut and r['mesh_points']==81)-
                            next(r['opening_A'] for r in new_branch if r['cut_kind']==cut and r['mesh_points']==161)) for cut in ['shuffle','glide'])
    if root_mesh_error>1e-8:raise RuntimeError('sampled branch root brackets disagree')
    resolution=read(root/'raw_runs/final_selected_convergence_results/resolution_checks.json')
    if len(resolution)!=20:raise RuntimeError('resolution audit is incomplete')
    for r in resolution:
        b,d=r['bessel'],r['direct']
        errors=[abs(b['energy']-d['energy']),float(np.max(abs(np.array(b['gradient'])-d['gradient']))),
                float(np.max(abs(np.array(b['hessian'])-d['hessian'])))]
        if r['shell_index']>=96 and r['nodes']>=512 and any(e>=t for e,t in zip(errors,[2e-8,2e-7,3e-6])):
            raise RuntimeError('nominal or refined Bessel resolution fails recorded tolerance')
    fits=read(root/'raw_runs/final_selected_finite_strain_results/fits.json')
    coupled=sorted((r for r in fits if r['pattern']=='cubic_exx_plus_yz'),key=lambda r:r['increment'],reverse=True)
    errors=[np.abs(np.array([r['C11_GPa'],r['C12_GPa'],r['C44_GPa']])-r['exact_local_constants_GPa']) for r in coupled]
    if len(errors)!=3 or not (np.all(errors[0]>errors[1]) and np.all(errors[1]>errors[2])):
        raise RuntimeError('independent conservative stress fits do not approach exact anchors')
    selection=read(root/'final_selection.json');score,name,profile=min((r['training_plus_EOS_objective_squared'],r['run'],r['profile'])
        for r in selection['profiles'] if r['eligible'])
    selected_path=root/'raw_runs'/name/'selected_fit.json'
    if (name,profile,score)!=(selection['selected_run'],selection['selected_profile'],selection['selected_training_plus_EOS_objective_squared']):
        raise RuntimeError('locked final training selection differs')
    if hashlib.sha256(selected_path.read_bytes()).hexdigest()!=selection['selected_fit_sha256']:
        raise RuntimeError('final selected parameter bytes differ')
    if len(selection['profiles'])!=30 or sum(r['eligible'] for r in selection['profiles'])!=12:
        raise RuntimeError('positive-boundary selection census differs')
    span=read(root/'raw_runs/force_span_diagnostic_recovery_results/group_force_bounds.json')
    for r in span:
        gram=np.array(r['Gram']);linear=np.array(r['linear']);n=r['components']
        for key,field in [('locked_amplitudes','locked_force_RMSE_eV_A'),
                          ('unconstrained_force_only_amplitudes','unconstrained_force_only_RMSE_eV_A'),
                          ('nonnegative_force_only_amplitudes','nonnegative_force_only_RMSE_eV_A')]:
            a=np.array(r[key]);squared=float(a@gram@a-2*a@linear+r['target_squared'])
            if squared<0 or abs(np.sqrt(squared/n)-r[field])>2e-12:raise RuntimeError('fixed-shape force bound replay differs')
        a=np.array(r['unconstrained_force_only_amplitudes'])
        if np.linalg.norm(gram@a-linear,np.inf)>1e-8:raise RuntimeError('force-only least-squares stationarity differs')
        expected=next(row['force_RMSE_eV_A'] for row in grouped if row['case']=='joint_static_candidate' and
            row['config_type']==r['config_type'] and row['declared_xc']==r['declared_xc'])
        if abs(expected-r['locked_force_RMSE_eV_A'])>2e-12:raise RuntimeError('fixed-shape diagnostic uses different final candidate')
    report=dict(complete=True,exported_raw_files_verified=verified,exact_source_bindings_verified=source_count,
        raw_prediction_frames_verified=len(meta),raw_prediction_cases_verified=len(cases['cases']),group_metrics_verified=len(grouped),
        maximum_force_group_difference_eV_A=maximum_force_error,maximum_energy_group_difference_eV_atom=maximum_energy_error,
        archived_LAMMPS_force_maximum_difference_eV_A=lmp_force,archived_LAMMPS_energy_maximum_difference_eV=lmp_energy,
        compact_Gram_profiles_replayed=replayed,maximum_loss_error_over_quadratic_term_scale=max_normalized_loss_error,
        actual_source_jets_audited=True,rigid_interface_SI_units_verified=True,
        new_candidate_force_Hessian_maximum_difference_eV_A2=maximum_fd,
        new_candidate_branch_mesh_maximum_difference_A=root_mesh_error,new_geometry_evaluations=0,
        recorded_Bessel_resolution_states_verified=len(resolution),
        new_candidate_finite_strain_smallest_step_error_GPa=errors[-1].tolist(),
        locked_final_training_selection_verified=True,
        fixed_shape_force_diagnostic_groups_verified=len(span),
        material_approved=False,first_initiation_validated=False,physical_clock_validated=False)
    if args.output is not None:args.output.write_bytes((json.dumps(report,indent=2)+'\n').encode())
    print(json.dumps(report))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--package',type=Path,default=Path(__file__).resolve().parent)
    p.add_argument('--repo',type=Path);p.add_argument('--lammps-control',type=Path);p.add_argument('--output',type=Path)
    main(p.parse_args())
