"""Package completed v5 evidence and actual JUnit records; does not run tests."""
from datetime import datetime, timezone
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import zipfile

from .run_atomistic_controls import checksum, save_json


def main():
    root = Path('results/silicon_atomistic_v5')
    cache = Path('.cache/si-atomistic-v5')
    suites, cases = [], {}
    for label, filename in [('related_suite', 'tests.xml'), ('updated_reference', 'reference_final.xml')]:
        tree = ET.parse(cache/filename)
        records = []
        for case in tree.findall('.//testcase'):
            name = case.get('classname')+'::'+case.get('name')
            if case.find('failure') is not None or case.find('error') is not None or case.find('skipped') is not None:
                raise ValueError(f'nonpassing test record: {name}')
            record = dict(node=name, seconds=float(case.get('time', '0')), passed=True)
            records.append(record)
            cases[name] = record
        suite = tree.find('.//testsuite')
        suites.append(dict(label=label, tests=len(records), seconds=float(suite.get('time')), cases=records))
    raw = json.loads((root/'raw_validation.json').read_text())
    controls = json.loads((root/'atomistic_controls.json').read_text())
    dft = json.loads((root/'dft_material_audit.json').read_text())
    assert raw['passed'] and controls['complete'] and dft['complete']
    assert len(cases) == 74
    source_verification = []
    for report in ['atomistic_controls.json', 'dft_material_audit.json']:
        record = json.loads((root/report).read_text())
        for source, expected in record['code_sha256'].items():
            actual_path = Path(source)
            if checksum(actual_path) != expected:
                if source != 'solver_v1/silicon_atomistic_reference.py' or report != 'atomistic_controls.json':
                    raise ValueError(f'unaccounted changed source {source}')
                actual_path = root/'source_snapshots/silicon_atomistic_reference_controls.py'
            if checksum(actual_path) != expected:
                raise ValueError(f'source snapshot mismatch {actual_path}')
            source_verification.append(dict(report=report, runtime_source=source,
                matching_source=str(actual_path).replace('\\', '/'), sha256=expected))
    for report, source in [('reconstruction_polish.json', 'polish_atomistic_reconstruction.py'),
                           ('raw_validation.json', 'validate_atomistic_audit.py')]:
        expected = json.loads((root/report).read_text())['code_sha256']
        path = Path('results/silicon_wafer_feasibility')/source
        if checksum(path) != expected:
            raise ValueError(f'executed source hash changed {source}')
        source_verification.append(dict(report=report, matching_source=str(path).replace('\\', '/'), sha256=expected))
    archives = sorted(root.glob('*.npz'))
    for p in archives:
        with zipfile.ZipFile(p) as z:
            if z.testzip() is not None:
                raise ValueError(f'corrupt result archive {p}')
    validation = dict(recorded_utc=datetime.now(timezone.utc).isoformat(),
        base_head='59e6cdf39f62100767aa13f3a174fcbbda7a82c7', branch='silicon-wafer-research',
        static_audit_complete=True, tests=suites, unique_related_passed=len(cases),
        raw_validation_passed=raw['passed'], checked_npz_archives=len(archives),
        source_verification=source_verification, independent_force_checks=len(raw['independent_evaluations']),
        full_Al_UI_regression_rerun=False, compileall_exit=0,
        material_calibrated=False, kinetic_calibrated=False, production_changed=False,
        new_DFT_calculations=0, new_GAP_calculations=0, new_MD_trajectories=0,
        classical_DFT_structure_evaluations=dft['evaluated_pairs'],
        v4_remaining_sampling_or_boundary_issues_resolved=False,
        figure_visual_inspection='PNG inspected; PDF from same figure',
        repository_status='verify final git status and fresh remote; no commit hash inferred here')
    save_json(root/'validation.json', validation)
    files = sorted(p for p in root.rglob('*') if p.is_file() and p.name != 'artifact_manifest.json')
    files += [Path('solver_v1')/p for p in ['silicon_atomistic_reference.py',
              'test_silicon_atomistic_reference.py', 'SILICON_ATOMISTIC_FOUNDATION_V5.md']]
    files += [Path('results/silicon_wafer_feasibility')/p for p in [
        'run_atomistic_controls.py', 'polish_atomistic_reconstruction.py', 'run_dft_material_audit.py',
        'validate_atomistic_audit.py', 'plot_atomistic_audit.py', 'package_atomistic_audit.py']]
    entries = [dict(path=str(p).replace('\\', '/'), bytes=p.stat().st_size, sha256=checksum(p)) for p in files]
    save_json(root/'artifact_manifest.json', dict(scope='completed v5 result and source bytes; excludes this manifest',
              files=entries, count=len(entries), total_bytes=sum(r['bytes'] for r in entries)))
    print(json.dumps(dict(unique_tests=len(cases), artifacts=len(entries), checked_npz=len(archives))))


if __name__ == '__main__':
    main()
