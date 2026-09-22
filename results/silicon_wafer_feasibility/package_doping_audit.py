"""Check completed v7 evidence and bind exact result/source bytes."""
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'results/silicon_doping_v7'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    test_summary = {}
    for scope in ('solver', 'app'):
        evidence = json.loads((OUT/f'{scope}_tests.json').read_text())
        assert evidence['exit_code'] == 0, scope
        assert not any(r['outcome'] == 'failed' for r in evidence['reports']), scope
        test_summary[scope] = {key: evidence[key] for key in ('passed_calls', 'passed_subtests', 'elapsed_wall_seconds')}
        test_summary[scope]['skipped_reports'] = sum(r['outcome'] == 'skipped' for r in evidence['reports'])
    smoke = json.loads((OUT/'desktop_smoke.json').read_text(encoding='utf-8-sig'))
    assert smoke['exit_code'] == 0
    final = json.loads((OUT/'final_targeted_tests.json').read_text())
    assert final['exit_code'] == 0
    for name, expected in final['source_sha256'].items():
        assert digest(ROOT/name) == expected, name
    audit = json.loads((OUT/'audit.json').read_text())
    for name, expected in audit['source_sha256'].items():
        assert digest(ROOT/name) == expected, name
    replay = json.loads((OUT/'replay_validation.json').read_text())
    assert all(replay['independent_rerun_byte_identical'].values())
    paths = sorted(p for p in OUT.rglob('*') if p.is_file() and p.name != 'artifact_manifest.json')
    paths += [ROOT/'solver_v1'/name for name in ('silicon_doping_research.py', 'test_silicon_doping_research.py', 'SILICON_DOPING_V7.md')]
    paths += [ROOT/'results/silicon_wafer_feasibility'/name for name in ('run_doping_audit.py', 'validate_doping_regression.py', 'package_doping_audit.py')]
    paths += [ROOT/name for name in audit['source_sha256'] if ROOT/name not in paths]
    manifest = dict(scope='v7 code, derivation, sources, results, and cited v6 atom-count records',
                    tests=test_summary, final_targeted_tests=final,
                    desktop_smoke=smoke,
                    files={p.relative_to(ROOT).as_posix(): dict(sha256=digest(p), bytes=p.stat().st_size) for p in sorted(set(paths))})
    (OUT/'artifact_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(dict(files=len(manifest['files']), tests=test_summary, final_targeted_passes=final['passed']), indent=2))


if __name__ == '__main__':
    main()
