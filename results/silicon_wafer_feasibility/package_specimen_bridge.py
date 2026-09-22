"""Inventory v6 artifacts and enforce source/output links before Git staging."""
from __future__ import annotations
import json
from pathlib import Path

from .run_atomistic_controls import checksum, save_json


def main():
    root = Path('results/silicon_specimen_v6')
    baseline = json.loads((root/'verified_affine_boundary/specimen_bridge.json').read_text())
    refined = json.loads((root/'verified_internal_boundary/refinement.json').read_text())
    validation = json.loads((root/'validation.json').read_text())
    tests = json.loads((root/'tests.json').read_text())
    assert baseline['complete'] and refined['complete'] and validation['passed'] and tests['passed']
    assert tests['tests'] == 88 and len(validation['actual_recomputed_states']) == 40
    assert (root/'specimen_bridge.png').stat().st_size > 10000
    for folder, data in [(root/'verified_affine_boundary', baseline), (root/'verified_internal_boundary', refined)]:
        for row in data['cases']:
            assert (folder/row['artifact']).is_file()
            if 'stationary_check' in row:
                assert (folder/row['stationary_check']['output']).is_file()
        for path, value in data['code_sha256'].items():
            assert checksum(path) == value
    # Preserve the precise sources of first-run failures and distinguish them
    # from the final verified pipeline. They are historical, not executables.
    original = json.loads((root/'specimen_bridge.json').read_text())
    for name, source in [('solver_v1/silicon_specimen_research.py', 'silicon_specimen_research_initial.py'),
                         ('results/silicon_wafer_feasibility/run_specimen_bridge.py', 'run_specimen_bridge_initial.py')]:
        assert checksum(root/'source_snapshots'/source) == original['code_sha256'][name]
    first_partial = json.loads((root/'internal_boundary/checkpoint.json').read_text())
    assert checksum(root/'source_snapshots/refine_specimen_bridge_initial.py') == first_partial['code_sha256']['results/silicon_wafer_feasibility/refine_specimen_bridge.py']
    sources = [Path(p) for p in sorted(set(baseline['code_sha256']) | set(refined['code_sha256']))]
    sources += [Path('results/silicon_wafer_feasibility/validate_specimen_bridge.py'), Path(__file__),
                Path('solver_v1/test_silicon_specimen_research.py'), Path('solver_v1/test_silicon_specimen_storage.py'),
                Path('solver_v1/SILICON_SPECIMEN_LOADING_V6.md')]
    paths = sorted(set(sources + [p for p in root.rglob('*') if p.is_file() and p.name != 'artifact_manifest.json']), key=lambda p:p.as_posix())
    entries = [dict(path=p.resolve().relative_to(Path.cwd()).as_posix(), bytes=p.stat().st_size, sha256=checksum(p)) for p in paths]
    save_json(root/'artifact_manifest.json', dict(entries=entries, total_bytes=sum(r['bytes'] for r in entries),
        scope='final verified pipeline plus explicitly marked initial-run failures',
        final_static_reference_cases=40, additional_boundary_work_relaxations=8,
        pytest_tests=88, subtests=6, material_calibrated=False, kinetic_calibrated=False))
    print(json.dumps(dict(files=len(entries), bytes=sum(r['bytes'] for r in entries), passed=True)))


if __name__ == '__main__':
    main()
