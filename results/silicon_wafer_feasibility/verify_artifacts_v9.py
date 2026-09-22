"""Read-only audit of saved v9 scientific artifacts; optional new manifest.

This does not rerun MACE, DFT, minimizations, or pytest. Hashes certify file
identity, while the explicitly linked numerical controls establish their scope.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
from scipy.linalg import null_space


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8'))


def main(args):
    root=args.root.resolve();repo=root.parents[1]
    assert '<!-- FINAL_' not in (root/'COMPLETED_SUMMARY.md').read_text(encoding='utf-8')
    comparisons=list(csv.DictReader((root/'mace_si_dft_validation/same_geometry_comparison.csv').open(encoding='utf-8')))
    source=read_json(root/'boron_source_audit/source_audit.json')
    assert digest(root/'sources/durham_boron_2025.zip')==source['zip_sha256']
    dft=read_json(root/'mace_si_dft/summary.json')
    validation=read_json(root/'mace_si_dft_validation/raw_validation.json')
    assert dft['frames']==validation['raw_frames_checked']==2475
    assert dft['atoms']==171815 and validation['source_arrays_exactly_equal_to_v5']
    assert sum(int(r['frames']) for r in comparisons)==2475 and len(comparisons)==35
    assert not dft['material_calibrated'] and not dft['kinetic_calibrated']
    hessians=[]
    for directory in ['mace_boron_hessian','mace_boron_follow_hessian','mace_boron_neutral_substitutional_hessian']:
        for path in sorted((root/directory).glob('*_hessian.json')):
            r=read_json(path);raw=np.load(path.with_suffix('.npz'))
            values=raw['Cartesian_projected_eigenvalues']
            assert float(values[0])==r['Cartesian_projected_min_eigenvalue_eV_A2']
            assert np.all(np.isfinite(raw['hessian']))
            assert raw['hessian'].shape==(3*r['atoms'],3*r['atoms'])
            translations=np.tile(np.eye(3),(r['atoms'],1))/np.sqrt(r['atoms'])
            basis=null_space(translations.T)
            h=.5*(raw['hessian']+raw['hessian'].T)
            recomputed=np.linalg.eigvalsh(basis.T@h@basis)
            spectrum_error=float(np.max(abs(recomputed-values)))
            assert spectrum_error<1e-10
            hessians.append(dict(file=path.relative_to(root).as_posix(),structure=r['structure'],
                min_curvature_eV_A2=float(values[0]),force_max_eV_A=r['max_force_eV_A'],
                raw_matrix_spectrum_recompute_error_eV_A2=spectrum_error,
                force_evaluations=r['force_evaluations']))
    assert len(hessians)==9
    unstable=[r['structure'] for r in hessians if r['min_curvature_eV_A2']<0]
    assert unstable==['3boron_int_cluster_confi_2']
    endpoints=read_json(root/'boron_endpoint_validation.json')
    assert len(endpoints['corrected_endpoints'])==2
    assert all(r['fixed_cell_minimum_verified'] for r in endpoints['corrected_endpoints'])
    for name in ['mace_boron_hessian_refinement','mace_boron_follow_hessian_refinement',
                 'mace_boron_neutral_substitutional_hessian_refinement']:
        r=read_json(root/name/'mode_refinement_summary.json')
        assert r['finest_max_curvature_error']<1e-4
    mode_summary=read_json(root/'mace_boron_largest_mode/summary.json')
    assert mode_summary['complete'] and len(mode_summary['cases'])==3
    for record,directory in zip(mode_summary['cases'],['mace_boron','mace_boron_216host','mace_boron_512host']):
        assert digest(root/directory/record['state_file'])==record['state_sha256']
        assert max(c['eigen_residual_L2_eV_A2_u'] for c in record['independent_step_checks'])<1e-6
        assert not record['full_mechanical_stability_checked']
        raw=np.load(root/'mace_boron_largest_mode'/(record['label']+'_largest_mode.npz'))
        assert abs(np.linalg.norm(raw['mass_weighted_eigenvector'])-1)<1e-10
    assert mode_summary['cases'][0]['full_Hessian_control']['absolute_eigenvalue_difference']<1e-6
    bulk_first=read_json(root/'mace_bulk_elastic/summary.json')
    bulk_counted=read_json(root/'mace_bulk_elastic_counted/summary.json')
    bulk_independent=read_json(root/'mace_bulk_elastic_independent/validation.json')
    assert bulk_first['finest_constants']==bulk_counted['finest_constants']
    assert bulk_counted['states']==36 and bulk_counted['tangents']==18
    assert bulk_counted['standard_ASE_calculator_calls']==25
    assert bulk_counted['standard_ASE_state_assessments']==37
    assert bulk_independent['reference_sha256']==digest(root/'mace_bulk_elastic/summary.json')
    assert max(abs(v) for v in bulk_independent['differences_from_hydro_tetragonal_xy_GPa'].values())<.001
    for name in ['states.csv','tangents.csv','elastic_constants.csv']:
        assert (root/'mace_bulk_elastic'/name).read_bytes()==(root/'mace_bulk_elastic_counted'/name).read_bytes()
    suite=read_json(root/'verification/full_solver.json')
    assert suite['exit_code']==0
    suite_log=(root/'verification/full_solver.log').read_text(encoding='utf-8')
    assert '946 passed' in suite_log and '6 subtests passed' in suite_log
    replay=[]
    if args.replay:
        for directory in ['dynamics','memory','global_charge','boron_source_audit']:
            for path in sorted((root/directory).iterdir()):
                if path.suffix not in {'.csv','.json'}:continue
                other=args.replay/path.relative_to(root)
                equal=path.read_bytes()==other.read_bytes()
                assert equal,'replay mismatch: '+str(path)
                replay.append(dict(file=path.relative_to(root).as_posix(),byte_equal=equal,sha256=digest(path)))
        assert len(replay)==19
    if args.write:
        if replay:
            (root/'verification/deterministic_replay.json').write_text(json.dumps(replay,indent=2)+'\n',encoding='utf-8')
        report=dict(status='completed saved-artifact consistency audit; not a rerun of the heavy calculations',
            DFT_reference_geometries=2475,DFT_reference_atoms=171815,source_groups=35,
            MACE_lower_RMSE_than_both_baselines_groups=sum(float(r['MACE_component_rmse_eV_A'])<min(float(r['SW_component_rmse_eV_A']),float(r['Tersoff_component_rmse_eV_A'])) for r in comparisons),
            source_groups_tied_at_zero=[r['config_type'] for r in comparisons if all(float(r[n+'_component_rmse_eV_A'])==0 for n in ['MACE','SW','Tersoff'])],
            checked_full_Hessians=hessians,unstable_original_states_preserved=unstable,
            matrix_free_highest_mode_cases=len(mode_summary['cases']),
            bulk_elastic_replay_CSV_files=3,bulk_elastic_independent_modes_checked=True,
            bulk_counter_correction=dict(initial_standard_ASE_evaluations_label_meant_state_assessments=37,
                counted_replay_actual_standard_calculator_calls=25,force_only_calls=34,
                numerical_CSV_changed=False),
            deterministic_replay_files=len(replay),full_solver_exit_code=suite['exit_code'],
            heavy_calculations_byte_replayed=False,new_DFT_runs=0,new_MD_runs=0,
            material_calibrated=False,physical_clock_calibrated=False)
        (root/'verification/saved_artifact_audit.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        paths=[p for p in root.rglob('*') if p.is_file() and p.name!='artifact_manifest.json']
        paths.extend(repo/'solver_v1'/name for name in ['silicon_charge_dynamics.py','test_silicon_charge_dynamics.py','SILICON_CHARGE_DYNAMICS_V9.md'])
        paths.extend((repo/'results/silicon_wafer_feasibility').glob('*_v9.py'))
        records=[dict(file=p.relative_to(repo).as_posix(),size_bytes=p.stat().st_size,sha256=digest(p)) for p in sorted(set(paths))]
        manifest=dict(scope='v9 results and source code; excludes this manifest and mutable project handoff',files=records)
        (root/'artifact_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    manifest=read_json(root/'artifact_manifest.json')
    for record in manifest['files']:
        path=repo/record['file']
        assert path.stat().st_size==record['size_bytes'] and digest(path)==record['sha256'],record['file']
    if args.git_index:
        # Git can otherwise normalize CRLF and silently invalidate raw hashes.
        # Check the actual staged blob bytes in one bounded batch process.
        records=manifest['files']+[dict(file=(root/'artifact_manifest.json').relative_to(repo).as_posix(),
            size_bytes=(root/'artifact_manifest.json').stat().st_size,sha256=digest(root/'artifact_manifest.json'))]
        request=''.join(':'+r['file']+'\n' for r in records).encode('utf-8')
        batch=subprocess.run(['git','cat-file','--batch'],cwd=repo,input=request,capture_output=True,check=True)
        cursor=0
        for record in records:
            end=batch.stdout.index(b'\n',cursor);header=batch.stdout[cursor:end].split();cursor=end+1
            assert len(header)==3 and header[1]==b'blob','missing staged blob: '+record['file']
            size=int(header[2]);data=batch.stdout[cursor:cursor+size];cursor+=size+1
            assert size==record['size_bytes'] and hashlib.sha256(data).hexdigest()==record['sha256'],record['file']
        assert cursor==len(batch.stdout)
    print(json.dumps(dict(manifest_files_verified=len(manifest['files']),full_Hessians_checked=len(hessians),
        deterministic_replay_files=len(replay),source_groups=len(comparisons),git_index_verified=args.git_index,status='passed')))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    p.add_argument('--replay',type=Path);p.add_argument('--write',action='store_true')
    p.add_argument('--git-index',action='store_true')
    main(p.parse_args())
