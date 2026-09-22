"""Verify and hash v8 evidence, preserving historical v7 result files."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'results/silicon_doping_v8'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    tests = {}
    for name in ['targeted', 'solver', 'app']:
        r = json.loads((OUT/f'{name}_tests.json').read_text())
        assert r['exit_code'] == 0
        assert all(x['outcome'] == 'passed' for x in r['reports'])
        tests[name] = {k:r[k] for k in ['passed_calls', 'elapsed_wall_seconds']}
        tests[name]['passed_subtests'] = sum(x['phase'] == 'call' and x['subtest'] for x in r['reports'])
        for p, sha in r.get('source_sha256', {}).items():
            assert digest(ROOT/p) == sha, p
    smoke = json.loads((OUT/'desktop_smoke.json').read_text())
    assert smoke['exit_code'] == 0
    assert all(json.loads((OUT/'replay_validation.json').read_text()).values())
    inputs = set()
    for p in [OUT/'independent_validation.json', OUT/'revised_continuum/audit.json']:
        audit = json.loads(p.read_text())
        for name, sha in audit['source_sha256'].items():
            assert digest(ROOT/name) == sha, name
            inputs.add(ROOT/name)
    paths = {p for p in OUT.rglob('*') if p.is_file() and p.name != 'artifact_manifest.json'} | inputs
    paths |= {ROOT/'solver_v1'/n for n in ['silicon_doping_research.py', 'silicon_doping_validation.py',
               'test_silicon_doping_research.py', 'test_silicon_doping_validation.py', 'SILICON_DOPING_VALIDATION_V8.md']}
    paths |= {ROOT/'results/silicon_wafer_feasibility'/n for n in ['run_doping_audit.py', 'audit_doping_v8.py',
                                                                              'package_doping_v8.py', 'validate_doping_regression.py']}
    result = dict(base_commit='4d1e3b984e3833f9800c4b09ba1af4b913ea3fc4',
                  scope='v8 mathematical validation; actual doped fracture and clock uncalibrated',
                  tests=tests, desktop_smoke=smoke,
                  files={p.relative_to(ROOT).as_posix():dict(sha256=digest(p), bytes=p.stat().st_size)
                         for p in sorted(paths)})
    (OUT/'artifact_manifest.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(dict(files=len(paths), tests=tests), indent=2))


if __name__ == '__main__':
    main()
