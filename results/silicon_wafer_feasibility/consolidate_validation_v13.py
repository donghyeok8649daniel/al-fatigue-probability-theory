"""Bind already executed disjoint tests; do not count this as another pytest run."""
from __future__ import annotations
import argparse, hashlib, json, math
from datetime import datetime, timezone
from pathlib import Path

DEFAULT_RECORDS = ('resume_implementation_validation.json', 'anharmonic_replay_validation.json',
                   'mpa_qe_replay_guard_validation.json')


def main(args):
    repo = args.repo.resolve(); root = args.results.resolve()
    if root != repo / 'results/silicon_initiation_v13': raise ValueError('v13 results required')
    output = root / 'validation.json'
    if output.exists(): raise ValueError('Preserve existing consolidated validation evidence')
    sources = {}; test_paths = set(); runs = []; count = 0
    for filename in DEFAULT_RECORDS:
        path = root / filename; raw = path.read_bytes(); record = json.loads(raw)
        passed = record['tests_passed']; elapsed = record['elapsed_s']
        if type(passed) is not int or passed <= 0 or not math.isfinite(elapsed) or elapsed < 0:
            raise ValueError('Invalid actual test count or duration')
        entries = record.get('source_files', record.get('source_manifest'))
        if entries is None:
            entries = [dict(path=p, sha256=s) for p,s in record['source_sha256'].items()]
        tests = {entry['path'] for entry in entries if Path(entry['path']).name.startswith('test_')}
        if not tests or tests & test_paths:
            raise ValueError('Missing or overlapping test source; do not double count runs')
        for entry in entries:
            source = repo / entry['path']
            if not source.resolve().is_relative_to(repo): raise ValueError('Source outside repository')
            actual = hashlib.sha256(source.read_bytes()).hexdigest()
            if actual != entry['sha256']: raise ValueError('Source changed after recorded tests: '+entry['path'])
            if entry['path'] in sources and sources[entry['path']] != actual: raise ValueError('Conflicting test-source hashes')
            sources[entry['path']] = actual
        test_paths |= tests; count += passed
        runs.append(dict(path=filename, sha256=hashlib.sha256(raw).hexdigest(), tests_passed=passed,
                         elapsed_seconds=elapsed, test_source_files=sorted(tests), original_record=record))
    result = dict(recorded_UTC=datetime.now(timezone.utc).isoformat(),
        tests_passed=count, disjoint_previously_executed_runs=len(runs), runs=runs,
        tests_rerun_by_consolidation=False, new_test_executions=0,
        source_files=[dict(path=p,sha256=s) for p,s in sorted(sources.items())],
        source_hashes_rechecked=True, old_v12_63_tests_included=False,
        scope='Actual 17+10+11 disjoint test records; implementation validation only. Initial setup failure is preserved in its original record.',
        new_model_calls=0, new_DFT=0, new_MD=0, material_approved=False, physical_clock=None)
    output.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(recorded_tests_passed=count,runs=len(runs),source_files=len(sources),tests_rerun=False)))


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--repo',type=Path,required=True)
    parser.add_argument('--results',type=Path,required=True)
    main(parser.parse_args())
