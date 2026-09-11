"""Inventory actual v23 results, failed attempts and executed regressions.

No fitting or physical certification. Smoke is really executed here; test
counts are read from completed JUnit reports, never copied from a prior turn.
The final byte manifest excludes its own ledger directory.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

from .report_core_validation_ledger import read_report
from .run_current_material_core import ROOT
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json


START = 'd1ac6695518d3738b8303687d8c6a7b58877203e'
PRESERVED = {
    'results/fcc111_active_interface/tangent_calibration_v20/shape_refinement/old_family/calibration.json':
        '8736cb9d28f1430e3991b1046ef42544cefaa3c9a4f5d743e0d02931329f5c2e',
    'results/current_material_core_v22/wider_probe_validation/research_candidate_snapshot.json':
        '2fd94b89168b173378f4c67d9d6faf52c0363292fc66d8d7e705f214a55de23a',
    '.cache/al-reference/Al99.eam.alloy':
        '60c8a085be79d273324ab421f5b1447578fef55c1acfc6492c0999f15ee8a284',
}


def failure_text(record):
    """Preserve both explicit historical failure schemas; never assume PASS."""
    value = record.get('error', record.get('message'))
    if not isinstance(value, str) or not value:
        raise ValueError('failure record has no explicit error/message text')
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--targeted-xml', type=Path, required=True)
    parser.add_argument('--topology-xml', type=Path, required=True)
    parser.add_argument('--solver-xml', type=Path, required=True)
    parser.add_argument('--app-xml', type=Path, required=True)
    args = parser.parse_args()
    study = ROOT/'results/core_interface_compatibility_v23'
    out = args.out.resolve(); out.relative_to(study)
    if out.exists():
        raise FileExistsError('fresh final ledger required')
    reports = dict(targeted=args.targeted_xml, topology_targeted=args.topology_xml,
                   solver_v1=args.solver_xml, app=args.app_xml)
    tests = {name: dict(read_report(path), report=path.resolve().relative_to(ROOT).as_posix())
             for name, path in reports.items()}
    if any(t['skipped'] for t in tests.values()):
        raise ValueError('unexpected skips must be reported and inspected, not hidden')
    preserved = []
    for name, expected in PRESERVED.items():
        actual = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(f'protected input binding changed: {name}')
        preserved.append(dict(path=name, sha256=actual, unchanged=True))
    for inspection in study.rglob('checkpoint_inspection.json'):
        saved = json.loads(inspection.read_bytes())
        actual = hashlib.sha256((inspection.parent/'checkpoint.csv').read_bytes()).hexdigest()
        if actual != saved['checkpoint_sha256']:
            raise ValueError('recorded checkpoint inspection does not bind the current bytes')
    tick = time.perf_counter()
    smoke = subprocess.run([sys.executable, '-X', 'utf8', '-m', 'app.desktop_ui', '--smoke'],
                           cwd=ROOT, capture_output=True, text=True, encoding='utf8', timeout=120)
    tests['desktop_smoke'] = dict(exit_code=smoke.returncode, elapsed_seconds=time.perf_counter()-tick,
        stdout=smoke.stdout.strip(), stderr=smoke.stderr.strip())
    if smoke.returncode:
        raise RuntimeError('actual smoke did not pass')
    cases, failures = [], []
    for path in sorted(study.rglob('metadata.json')):
        meta = json.loads(path.read_bytes()); parent = path.parent
        if 'radius_over_L0' not in meta or 'parameter_sha256' not in meta:
            continue
        row = dict(case=parent.relative_to(study).as_posix(), parameter_source=meta['parameter_source'],
            parameter_sha256=meta['parameter_sha256'], radius_L0=meta['radius_over_L0'], ring=meta['ring'],
            minimizer=None, outcome='incomplete', force_eV_L0=None, minimum_H_eV_L0sq=None,
            elapsed_seconds=None, error=None, material_accepted=False, physical_Hz=False)
        if (parent/'summary.json').exists():
            saved = json.loads((parent/'summary.json').read_bytes()); stability = saved['final_stability_probe']
            row.update(outcome='verified_fixed_boundary_minimum' if stability['stable_on_tested_fixed_boundary']
                       else 'completed_NOT_verified_minimum', minimizer=saved['minimizer'],
                force_eV_L0=saved['maximum_free_force'], minimum_H_eV_L0sq=stability['minimum_eigenvalue'],
                elapsed_seconds=saved['elapsed_seconds'])
        elif (parent/'failure.json').exists():
            saved = json.loads((parent/'failure.json').read_bytes())
            row.update(outcome='failed_attempt_preserved', error=failure_text(saved),
                       elapsed_seconds=saved.get('elapsed_seconds'))
        cases.append(row)
    for path in sorted(study.rglob('failure.json')):
        saved = json.loads(path.read_bytes())
        failures.append(dict(attempt=path.parent.relative_to(study).as_posix(),
            error_type=saved.get('error_type'), error=failure_text(saved),
            elapsed_seconds=saved.get('elapsed_seconds'), completed=False,
            failure_alone_proves_physical_instability=False))
    artifacts = [dict(path=p.relative_to(study).as_posix(), bytes=p.stat().st_size,
        sha256=hashlib.sha256(p.read_bytes()).hexdigest())
        for p in sorted(study.rglob('*')) if p.is_file() and out not in p.parents]
    write_csv(out/'atomic_case_status.csv', cases)
    write_csv(out/'failed_attempts.csv', failures)
    write_csv(out/'preserved_input_hashes.csv', preserved)
    write_csv(out/'artifact_hashes.csv', artifacts)
    save_json(out/'executed_tests.json', tests)
    save_json(out/'status.json', dict(completed=True, starting_HEAD=START,
        recorded_UTC=datetime.now(timezone.utc).isoformat(),
        recorded_core_cases_including_preflight=len(cases), failure_records=len(failures),
        live_or_unaccounted_atomic_cases=sum(r['outcome'] == 'incomplete' for r in cases),
        material_adopted=False, whole_family_impossibility_proved=False,
        source_Hxx_sign_not_a_candidate_Morse_test=True,
        default_material_changed=False, production_PDE_changed=False,
        actual_yield_validated=False, finite_loop_activation_validated=False,
        production_a_s_mobility_available=False, physical_seconds=False, physical_Hz=False,
        hashes_do_not_certify_physical_source_truth=True))
    print(json.dumps(dict(regressions=tests, recorded_core_cases=len(cases), failures=len(failures)), indent=2))


if __name__ == '__main__':
    main()
