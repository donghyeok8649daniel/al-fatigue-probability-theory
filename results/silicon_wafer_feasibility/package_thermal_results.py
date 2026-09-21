"""Validate artifact integrity and record honest completion/verification status."""
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import subprocess
import zipfile

from .run_local_crack_audit import save_json


def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))


def main():
    root=Path('results/silicon_thermal_v4')
    sampling=[];dynamics=[]
    for name in ['ensemble_pilot','hmc_pilot','mean_force_profile','front8_ensembles','multibasin_box15']:
        path=root/name/'summary.json'
        if not path.exists():path=root/name/'interrupted_summary.json'
        meta=read(path)
        sampling.append(dict(name=name,chains=len(meta['records']),
            saved_draws=sum(r['draws'] for r in meta['records']),complete=meta.get('computation_completed',False)))
    for name in ['thermal_dynamics','front8_dynamics','timestep_half','reconstruction_dynamics']:
        meta=read(root/name/'summary.json')
        dynamics.append(dict(name=name,trajectories=len(meta['records']),
            total_ps=meta['total_new_trajectory_ps'],complete=meta['computation_completed']))
    raw=read(root/'raw_validation.json')
    npz_count=0
    for path in root.rglob('*.npz'):
        with zipfile.ZipFile(path) as archive:
            failure=archive.testzip()
            if failure is not None:raise RuntimeError(f'corrupt archive member: {path.name}/{failure}')
        npz_count+=1
    for path in (root/'source_snapshots').glob('*.py'):
        if hashlib.sha256(path.read_bytes()).hexdigest()!=path.stem:raise RuntimeError('historical source checksum differs')
    report=dict(recorded_utc=datetime.now(timezone.utc).isoformat(),
        base_head=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        branch=subprocess.check_output(['git','branch','--show-current'],text=True).strip(),
        fresh_remote_head=subprocess.check_output(['git','rev-parse','origin/silicon-wafer-research'],text=True).strip(),
        all_planned_calculations_complete=False,completion_status='completed and interrupted calculations explicitly separated',
        budget_start_utc='2026-09-20T19:11:55Z',budget_deadline_utc='2026-09-20T23:11:55Z',
        clock_gap_cause_verified=False,sampling=sampling,new_MD=dynamics,
        total_new_MD_ps=sum(r['total_ps'] for r in dynamics),
        total_saved_MCMC_draws=sum(r['saved_draws'] for r in sampling),
        tests=dict(silicon=dict(passed=42,seconds=31.85),material_gate=dict(passed=24,seconds=11.01),
            profile_recheck=dict(passed=1,seconds=.73),unique_related_passed=66,
            full_solver_regression_rerun=False,desktop_smoke_exit=0,compileall_exit=0),
        raw_validation=dict(ensemble_states=len(raw['ensemble_states']),MD_states=len(raw['md_states']),
            independent_triple_states=len(raw['direct_triple_checks']),source_hash_mismatches=raw['source_hash_mismatches'],
            max_ensemble_observation_error=max(r['observation_max_error'] for r in raw['ensemble_states']),
            max_MD_force_energy_error=max(r['force_energy_max_error'] for r in raw['md_states']),
            max_triple_force_error_eV_A=max(r['force_max_error_eV_A'] for r in raw['direct_triple_checks']),
            max_triple_energy_error_eV=max(r['energy_error_eV'] for r in raw['direct_triple_checks'])),
        checked_npz_archives=npz_count,plots_visually_inspected=['thermal_research.png','boundary_diagnostic.png'],
        production_enabled=False,real_wafer_calibrated=False,physical_clock_calibrated=False,
        git_publication='this file records pre-commit verification; final commit and remote equality are reported after push')
    save_json(root/'validation.json',report)
    files=[]
    for path in sorted(root.rglob('*')):
        if not path.is_file() or path.name=='artifact_manifest.json':continue
        if path.suffix not in ('.npz','.json','.csv','.png','.py'):continue
        files.append(dict(path=path.relative_to(root).as_posix(),bytes=path.stat().st_size,
            sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    save_json(root/'artifact_manifest.json',dict(files=files,
        scope='raw/generated artifacts and historical LF source snapshots; prose/logs excluded from byte manifest',
        total_bytes=sum(r['bytes'] for r in files)))
    print('verified',npz_count,'NPZ archives; manifest',len(files),'files;',report['total_new_MD_ps'],'new ps',flush=True)


if __name__=='__main__':main()
