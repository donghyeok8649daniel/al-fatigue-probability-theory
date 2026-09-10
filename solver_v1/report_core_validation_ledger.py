"""Inventory actual completed/failed v22 studies and executed test reports.

No directory name or optimizer success flag substitutes for the force/Morse
gate. Hashes preserve the underlying observations; they are not a physical
material certificate. A fresh output avoids overwriting earlier inventories.
"""
import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

from .run_current_material_core import ROOT
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json


def read_report(path):
    root=ET.parse(path).getroot()
    suites=[root] if root.tag=='testsuite' else list(root.findall('testsuite'))
    if not suites:raise ValueError('actual completed JUnit report required')
    result={key:sum(int(s.get(key,'0')) for s in suites)
            for key in ('tests','failures','errors','skipped')}
    result['elapsed_seconds']=sum(float(s.get('time','0')) for s in suites)
    result['passed']=result['tests']-result['failures']-result['errors']-result['skipped']
    result['sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
    if result['tests']<=0 or result['failures'] or result['errors']:
        raise ValueError('failed/empty regression cannot be certified as passed')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=ROOT/'results/current_material_core_v22')
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--targeted-xml',type=Path,required=True)
    parser.add_argument('--solver-xml',type=Path,required=True)
    parser.add_argument('--app-xml',type=Path,required=True)
    args=parser.parse_args();study=args.root.resolve();out=args.out.resolve()
    out.relative_to(study)
    if out.exists():raise FileExistsError('fresh validation ledger required')
    tests={name:read_report(path) for name,path in [('targeted',args.targeted_xml),
        ('solver_v1',args.solver_xml),('app',args.app_xml)]}
    tick=time.perf_counter()
    smoke=subprocess.run([sys.executable,'-X','utf8','-m','app.desktop_ui','--smoke'],
                         cwd=ROOT,capture_output=True,text=True,encoding='utf8',timeout=120)
    tests['desktop_smoke']=dict(exit_code=smoke.returncode,elapsed_seconds=time.perf_counter()-tick,
                               stdout=smoke.stdout.strip(),stderr=smoke.stderr.strip())
    if smoke.returncode:raise RuntimeError('actual desktop smoke did not pass')
    cases=[];recoveries={};bindings={}
    for path in sorted(study.rglob('metadata.json')):
        meta=json.loads(path.read_bytes());parent=path.parent
        if 'parameter_sha256' not in meta or 'radius_over_L0' not in meta:continue
        parameter=meta['parameter_sha256'];bindings[parameter]=meta['parameter_source']
        recovery=meta.get('checkpoint_recovery')
        if recovery:recoveries.setdefault(recovery['source'],[]).append(parent.relative_to(ROOT).as_posix())
        row=dict(case=parent.relative_to(study).as_posix(),parameter_sha256=parameter,
            radius=meta['radius_over_L0'],ring=meta['ring'],shear_MPa=meta['shear_traction_MPa'],
            force_residual=None,minimum_H=None,optimizer_success=None,error=None,
            material_accepted=False,physical_yield_validated=False,physical_Hz=False)
        if (parent/'summary.json').exists():
            saved=json.loads((parent/'summary.json').read_bytes());mode=saved['final_stability_probe']
            row.update(outcome='completed_fixed_boundary_minimum' if mode['stable_on_tested_fixed_boundary']
                       else 'completed_but_not_a_verified_minimum',
                force_residual=saved['maximum_free_force'],minimum_H=mode['minimum_eigenvalue'],
                optimizer_success=saved['optimizer_success'],elapsed_seconds=saved['elapsed_seconds'])
        elif (parent/'saddle_summary.json').exists():
            saved=json.loads((parent/'saddle_summary.json').read_bytes())
            row.update(outcome='completed_index_one_stationary_point' if saved['force_converged']
                       and saved['negative_eigenvalues']==1 else 'stationary_point_not_verified_index_one',
                force_residual=saved['maximum_force'],minimum_H=saved['minimum_mode_probe']['minimum_eigenvalue'],
                elapsed_seconds=saved['elapsed_seconds'])
        elif (parent/'failure.json').exists():
            saved=json.loads((parent/'failure.json').read_bytes())
            row.update(outcome='failed_attempt_preserved',error=saved['error'],elapsed_seconds=saved['elapsed_seconds'])
        else:row.update(outcome='incomplete_not_certified',elapsed_seconds=None)
        cases.append(row)
    for row in cases:
        path=(study/row['case']).relative_to(ROOT).as_posix()
        row['checkpoint_recovered_by']=';'.join(recoveries.get(path,[]))
    files=[dict(path=p.relative_to(study).as_posix(),bytes=p.stat().st_size,
                sha256=hashlib.sha256(p.read_bytes()).hexdigest())
           for p in sorted(study.rglob('*')) if p.is_file() and out not in p.parents]
    write_csv(out/'case_status.csv',cases);write_csv(out/'artifact_hashes.csv',files)
    save_json(out/'executed_tests.json',tests)
    save_json(out/'status.json',dict(completed=True,recorded_UTC=datetime.now(timezone.utc).isoformat(),
        starting_HEAD='a8e4130934ab8bb67f2a657a489d269d10b20b73',parameter_bindings=bindings,
        actual_atomic_cases=len(cases),failed_attempts=sum(r['outcome']=='failed_attempt_preserved' for r in cases),
        incomplete_cases=sum(r['outcome']=='incomplete_not_certified' for r in cases),
        represented_fixed_boundary_checks_not_material_certification=True,
        production_PDE_changed=False,default_material_changed=False,
        actual_yield_calibrated=False,finite_loop_activation_validated=False,
        production_a_s_mobility_available=False,physical_seconds=False,physical_Hz=False,
        hashes_do_not_certify_external_source_truth=True))


if __name__=='__main__':main()
