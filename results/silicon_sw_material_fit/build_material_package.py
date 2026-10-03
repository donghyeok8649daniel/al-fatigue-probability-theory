"""Curate exact sources, compact replay statistics, and raw model predictions.

Large repeated design arrays stay in the original local cache, with exact SHA,
size and source bindings. The public package preserves every small raw report,
historical source snapshots, selected predictions, and quadratic Gram statistics.
No research evaluations are performed by this packaging operation.
"""
import argparse,hashlib,json,shutil
from pathlib import Path
from collections import defaultdict
import numpy as np


def save(path,data):
    path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes((json.dumps(data,indent=2,allow_nan=False)+'\n').encode())


def main(args):
    import sys
    sys.path.insert(0,str(args.repo.resolve()))
    from results.silicon_wafer_feasibility.run_dft_material_audit import read_frames
    if args.output.exists():raise ValueError('fresh package directory required')
    args.output.mkdir(parents=True);catalog=[];gram=[]
    names=['results','guarded_results','dual_repair_results','stress_results','elastic_anchor_results',
        'density_results','coordination_results','eos_reaudit_results','angular_results',
        'angular_continuum_results','angular_continuum_corrected_results','angular_shape_results',
        'angular_shape_recovery_results','phonon_audit_results','external_audit_results',
        'extended_bulk_results','extended_fit_results','angle_preference_results','angle_preference_fit_results',
        'interface_audit_results','interface_recovery_results','finite_strain_results','final_phonon_results','independent_elastic_force_root_results',
        'independent_elastic_results','elastic_guard_results','anharmonic_first_test_failure',
        'joint_angle_preflight_results','joint_angle_fit_results','angular_boundary_preflight_results',
        'angular_boundary_fit_results','refined_boundary_preflight_results','refined_boundary_fit_results',
        'new_candidate_external_results','new_candidate_phonon_results','new_candidate_interface_results',
        'new_candidate_finite_strain_results','final_convergence_results',
        'local_boundary_preflight_results','local_boundary_fit_results',
        'guard_boundary_preflight_results','guard_boundary_fit_results','final_selection_sources',
        'final_selected_external_results','final_selected_phonon_results','final_selected_interface_results',
        'final_selected_finite_strain_results','final_selected_convergence_results','final_selected_phonon_recovery_results',
        'force_span_diagnostic_results','force_span_diagnostic_recovery_results']
    for name in names:
        folder=args.cache/name
        for source in sorted(folder.rglob('*')):
            if not source.is_file():continue
            raw=source.read_bytes();relative=source.relative_to(args.cache).as_posix();suffix=source.suffix.lower()
            copied=suffix in ['.json','.csv','.py','.sw']
            if name in ['phonon_audit_results','final_phonon_results','external_audit_results',
                        'new_candidate_external_results','new_candidate_phonon_results',
                        'final_selected_external_results','final_selected_phonon_results',
                        'final_selected_phonon_recovery_results'] and suffix=='.npz':copied=True
            if copied:
                target=args.output/'raw_runs'/relative;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(raw)
            catalog.append(dict(cache_relative_path=relative,sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw),
                exported_path=('raw_runs/'+relative) if copied else None,
                large_raw_storage='original working cache retained; regenerate with exact source snapshot and archived DFT input' if not copied else None))
            if suffix=='.npz' and ('design' in source.stem):
                d=np.load(source)
                if 'design' not in d:continue
                x,y=d['design'],d['target'];b,t=d['bulk_design'],d['bulk_target']
                if x.ndim==2:
                    rec=dict(run=name,raw_file=source.name,source_sha256=hashlib.sha256(raw).hexdigest(),
                        design_shape=list(x.shape),xTx=(x.T@x).tolist(),xTy=(x.T@y).tolist(),yTy=float(y@y),
                        bTb=(b.T@b).tolist(),bTt=(b.T@t).tolist(),tTt=float(t@t),
                        maximum_bulk_sq=float(d['maximum_bulk_sq']))
                else:
                    rec=dict(run=name,raw_file=source.name,source_sha256=hashlib.sha256(raw).hexdigest(),design_shape=list(x.shape),
                        beta_quadratic=True,xTx=np.einsum('ira,jrb->ijab',x,x).tolist(),xTy=np.einsum('ira,r->ia',x,y).tolist(),
                        yTy=float(y@y),bTb=np.einsum('ira,jrb->ijab',b,b).tolist(),bTt=np.einsum('ira,r->ia',b,t).tolist(),
                        tTt=float(t@t),EOS_error_basis=d['EOS_error_basis'].tolist(),maximum_bulk_sq=float(d['maximum_bulk_sq']))
                gram.append(rec)
    # Earlier private audit runners recorded their helper path without nesting
    # a source snapshot. Resolve only an exactly matching recorded SHA.
    supplemental=[]
    for manifest in (args.output/'raw_runs').glob('*/source_manifest.json'):
        folder=manifest.parent
        for binding in json.loads(manifest.read_text()):
            if not binding['path'].endswith(('.py','.sw')):continue
            target=folder/'source_snapshots'/binding['path']
            if target.exists():
                if hashlib.sha256(target.read_bytes()).hexdigest()!=binding['sha256']:raise RuntimeError('copied snapshot mismatch')
                continue
            providers=[args.cache/folder.name/'calculation_fit_source.py',args.cache/binding['path'],args.repo/binding['path']]
            source=next((f for f in providers if f.exists() and hashlib.sha256(f.read_bytes()).hexdigest()==binding['sha256']),None)
            if source is None:raise RuntimeError('historical source bytes unavailable: '+folder.name+'/'+binding['path'])
            target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(source.read_bytes())
            supplemental.append(dict(path=target.relative_to(args.output).as_posix(),sha256=binding['sha256'],
                exact_recorded_source=True,source_role='recorded helper recovered from original exact bytes'))
    save(args.output/'supplemental_source_bindings.json',supplemental)
    save(args.output/'raw_cache_catalog.json',catalog);save(args.output/'design_gram_replay.json',gram)
    baseline=np.load(args.cache/'results/features.npz');offsets=baseline['offsets'];ref_e=baseline['reference_energy'];ref_f=baseline['reference_force']
    first=json.loads((args.cache/'results/fit.json').read_text())
    cases=[('original_SW','results/features.npz',np.ones(3)),('amplitude_fit','results/features.npz',np.array(first['primary']['amplitudes']))]
    final_selection=json.loads((args.cache/'final_selection.json').read_text())
    for name,label in [('guarded_results','bulk_energy_guard'),('angular_results','angular_grid'),('angular_shape_results','angular_shape'),
                       ('extended_fit_results','extended_support_diagnostic'),('angle_preference_fit_results','angle_preference_diagnostic'),
                       ('joint_angle_fit_results','joint_angle_diagnostic'),('angular_boundary_fit_results','angular_boundary_diagnostic'),
                       ('refined_boundary_fit_results','early_joint_candidate'),('local_boundary_fit_results','local_joint_candidate'),
                       ('guard_boundary_fit_results','guard_joint_candidate')]:
        if name==final_selection['selected_run']:label='joint_static_candidate'
        selected=json.loads((args.cache/name/'selected_fit.json').read_text());alpha=selected.get('amplitudes',selected.get('fit',{}).get('amplitudes'))
        cases.append((label,name+'/selected_features.npz',np.array(alpha)))
    energies=[];forces=[];cases_manifest=[]
    for label,path,alpha in cases:
        d=np.load(args.cache/path)
        if not np.array_equal(d['offsets'],offsets) or not np.array_equal(d['reference_force'],ref_f) or not np.array_equal(d['reference_energy'],ref_e):
            raise RuntimeError('prediction sources disagree on original references')
        energy=d['energy_basis']@alpha;force=np.einsum('nck,k->nc',d['force_basis'],alpha)
        if label=='original_SW':energy=d['energy_basis'].sum(axis=1);force=d['force_basis'].sum(axis=2)
        energies.append(energy);forces.append(force)
        cases_manifest.append(dict(label=label,raw_cache_relative_path=path,source_sha256=hashlib.sha256((args.cache/path).read_bytes()).hexdigest(),
            amplitudes=alpha.tolist(),production_enabled=False,material_approved=False))
    np.savez_compressed(args.output/'selected_predictions.npz',offsets=offsets,reference_energy=ref_e,reference_force=ref_f,
                         model_energy=np.stack(energies),model_force=np.stack(forces))
    frames,archive=read_frames(args.archive);training=set(first['training_frames']);af=first['anchor_frame'];metadata=[]
    if len(frames)!=len(ref_e) or sum(len(a) for a in frames)!=len(ref_f):raise RuntimeError('prediction frame count differs from archived input')
    for i,a in enumerate(frames):metadata.append(dict(frame=i,atoms=len(a),config_type=str(a.info['config_type']),
        declared_xc=str(a.info.get('xc_functional','UNSPECIFIED')),training=i in training))
    save(args.output/'frame_metadata.json',metadata);save(args.output/'prediction_cases.json',dict(cases=cases_manifest,anchor_frame=af,
        energy_anchor='same PW91 thermal frame per atom, fixed and not relaxed DFT ground state',
        archive_sha256=hashlib.sha256(args.archive.read_bytes()).hexdigest(),frames=len(frames),atoms=len(ref_f),new_DFT=0,new_MD=0))
    metrics=[];groups=defaultdict(list)
    for record in metadata:groups[(record['config_type'],record['declared_xc'],record['training'])].append(record)
    for case_index,(label,_,_) in enumerate(cases):
        for (kind,xc,trained),g in sorted(groups.items()):
            squared=0.;ee=[]
            for r in g:
                i=r['frame'];lo,hi=offsets[i:i+2];diff=forces[case_index][lo:hi]-ref_f[lo:hi];squared+=float(np.sum(diff**2))
                if xc=='PW91':ee.append(float(energies[case_index][i]/r['atoms']-energies[case_index][af]/len(frames[af])-(ref_e[i]/r['atoms']-ref_e[af]/len(frames[af]))))
            metrics.append(dict(case=label,config_type=kind,declared_xc=xc,training=trained,frames=len(g),atoms=sum(r['atoms'] for r in g),
                force_RMSE_eV_A=float(np.sqrt(squared/(3*sum(r['atoms'] for r in g)))),relative_energy_RMSE_eV_atom=float(np.sqrt(np.mean(np.square(ee)))) if ee else None))
    save(args.output/'reaggregated_group_metrics.json',metrics)
    # Copy only factual primary tables and the small external raw structure set.
    primary=['test-results/model-CASTEP_ASE-test-bulk_diamond-properties.json','test-results/model-SW-test-bulk_diamond-properties.json',
             'test-results/model-CASTEP_ASE-test-phonon_diamond-properties.json','test-results/model-SW-test-phonon_diamond-properties.json',
             'tests/testing_database_force_error/testing_database.xyz','manifest.json','method_manifest.json','additional_manifest.json',
             'historical_elasticity/manifest.json']
    for relative in primary:
        source=args.cache/'framework_sources'/relative;target=args.output/'primary_sources'/relative
        target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(source.read_bytes())
    save(args.output/'packaging_summary.json',dict(complete=True,large_raw_files_retained_locally=True,cache_files=len(catalog),
        raw_files_exported=sum(r['exported_path'] is not None for r in catalog),designs_with_Gram_statistics=len(gram),
        core_prediction_cases=len(cases),raw_prediction_frames=len(frames),raw_prediction_atoms=len(ref_f),
        no_new_research_geometry_evaluations=True,material_approved=False,first_initiation_validated=False,physical_clock_validated=False,
        known_invalid_results=['angular_continuum_results EOS sum normalization; corrected run is authoritative',
                               'interface_audit_results GPa display ten times too small; interface_recovery_results is authoritative'],
        corrections=['density EOS reuse audited separately; angular shape final JSON error recovered by raw audit; elastic anchor call count audited separately']))
    for name in ['elastic_structure_audit.json','elastic_structure_audit.py','elastic_bulk_feasibility.py']:
        source=args.cache/name;target=args.output/'structural_audits'/name;target.parent.mkdir(parents=True,exist_ok=True)
        target.write_bytes(source.read_bytes())
    raw=Path(__file__).read_bytes();(args.output/'build_material_package.py').write_bytes(raw)
    (args.output/'final_selection.json').write_bytes((args.cache/'final_selection.json').read_bytes())
    (args.output/'select_final_shape.py').write_bytes((args.cache/'select_final_shape.py').read_bytes())
    (args.output/'diagnose_selected_force_span.py').write_bytes((args.cache/'diagnose_selected_force_span.py').read_bytes())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['cache','repo','archive','output']:p.add_argument('--'+name,required=True,type=Path)
    main(p.parse_args())
