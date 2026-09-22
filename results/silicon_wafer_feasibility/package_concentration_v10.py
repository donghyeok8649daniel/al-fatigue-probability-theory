"""Hash completed research evidence and verify the actual staged Git bytes."""
import argparse
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import platform
import subprocess


REPO=Path(__file__).resolve().parents[2]
DATA=REPO/'results/silicon_concentration_v10'
CODES=[
    'solver_v1/silicon_concentration_research.py',
    'solver_v1/test_silicon_concentration_research.py',
    'solver_v1/SILICON_CONCENTRATION_V10.md',
    *['results/silicon_wafer_feasibility/'+name+'.py' for name in [
        'run_concentration_cleavage_v10','run_loaded_concentration_v10',
        'validate_concentration_cleavage_v10','check_concentration_stationarity_v10','validate_concentration_modes_v10',
        'summarize_concentration_v10','audit_concentration_coverage_v10',
        'package_concentration_v10','mace_force_only_v9','run_mace_boron_v9']]]


def sha(data):return hashlib.sha256(data).hexdigest()


def main(args):
    manifest=DATA/'artifact_manifest.json'
    if args.write:
        for name,number in [('fixed_cell_v2',25),('loaded_pure_100MPa',1),('loaded_dopants_100MPa',24)]:
            r=json.loads((DATA/name/'summary.json').read_text())
            if not r['complete'] or len(r['cases'])!=number:raise AssertionError('incomplete execution: '+name)
        for name in ['independent_pure','independent_B8','independent_P8']:
            if not json.loads((DATA/name/'summary.json').read_text())['passed']:
                raise AssertionError('independent control failed: '+name)
        for name in ['lowest_B8_cluster','lowest_P8_cluster_v2']:
            r=json.loads((DATA/name/'summary.json').read_text())
            if not r['solver_converged'] or len(r['eigenvalues_eV_A2'])!=2:
                raise AssertionError('selected curvature search incomplete: '+name)
        for name in ['lowest_B8_energy_control_v2','lowest_P8_energy_control']:
            if not json.loads((DATA/name/'summary.json').read_text())['passed']:
                raise AssertionError('independent energy curvature failed: '+name)
        for species,mode_folder,energy_folder in [('B','lowest_B8_cluster','lowest_B8_energy_control_v2'),
                                                  ('P','lowest_P8_cluster_v2','lowest_P8_energy_control')]:
            mode=json.loads((DATA/mode_folder/'summary.json').read_text())
            energy=json.loads((DATA/energy_folder/'summary.json').read_text())
            state=DATA/'fixed_cell_v2'/f'{species}_08_interface_cluster'/'opening_0.000.npz'
            if mode['state_sha256']!=sha(state.read_bytes()):raise AssertionError('curvature source identity')
            if energy['modes_sha256']!=sha((DATA/mode_folder/'modes.npz').read_bytes()):
                raise AssertionError('energy curvature mode identity')
        replay=json.loads((DATA/'verification/postprocessing_replay.json').read_text())
        if not replay['passed'] or not all(r['byte_equal'] for r in replay['files']):
            raise AssertionError('postprocessing replay')
        analysis=json.loads((DATA/'analysis/summary.json').read_text())
        if not analysis['complete'] or analysis['raw_states']!=425 or analysis['loaded_raw_states']!=444:
            raise AssertionError('raw state reconstruction')
        tests=json.loads((DATA/'verification/final_related_tests.json').read_text())
        if tests['exit_code']!=0:raise AssertionError('related tests')
        runtime=DATA/'runtime.json'
        if not runtime.exists():
            runtime.write_text(json.dumps(dict(python=platform.python_version(),
                packages={key:version(key) for key in ['numpy','scipy','ase','mace-torch','torch','e3nn','matplotlib']},
                device='CPU',model_dtype='float64',threads_per_MACE_process=2,
                model_weights_committed=False,physical_clock_calibrated=False),indent=2)+'\n',encoding='utf-8')
        files=sorted({p for p in DATA.rglob('*') if p.is_file() and p!=manifest}
                     | {REPO/name for name in CODES})
        entries=[]
        for p in files:
            data=p.read_bytes();entries.append(dict(path=p.relative_to(REPO).as_posix(),bytes=len(data),sha256=sha(data)))
        manifest.write_text(json.dumps(dict(schema=1,files=entries,file_count=len(entries),
            scope='immutable v10 evidence and exact research sources; mutable project handoff excluded'),indent=2)+'\n',encoding='utf-8')
    expected=json.loads(manifest.read_text())
    for r in expected['files']:
        p=REPO/r['path'];p.resolve().relative_to(REPO.resolve())
        data=p.read_bytes()
        if len(data)!=r['bytes'] or sha(data)!=r['sha256']:raise AssertionError('working file mismatch: '+r['path'])
    if args.git_index:
        request=''.join(':'+r['path']+'\n' for r in expected['files']).encode()
        output=subprocess.run(['git','cat-file','--batch'],cwd=REPO,input=request,stdout=subprocess.PIPE,check=True).stdout
        start=0
        for r in expected['files']:
            end=output.index(b'\n',start);header=output[start:end].decode();fields=header.split()
            if len(fields)!=3 or fields[1]!='blob':raise AssertionError('missing staged blob: '+r['path'])
            size=int(fields[2]);data=output[end+1:end+1+size]
            if size!=r['bytes'] or sha(data)!=r['sha256']:raise AssertionError('Git byte mismatch: '+r['path'])
            start=end+2+size
        if start!=len(output):raise AssertionError('unexpected cat-file trailing output')
    print(json.dumps(dict(files_verified=len(expected['files']),git_index_verified=args.git_index,status='passed')),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');p.add_argument('--git-index',action='store_true')
    main(p.parse_args())
