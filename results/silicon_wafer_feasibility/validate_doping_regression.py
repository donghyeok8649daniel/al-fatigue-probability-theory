"""Persist test outcomes without hostnames or private absolute paths."""
import argparse
import json
from pathlib import Path
import time

import pytest


class Recorder:
    def __init__(self):
        self.reports = []

    def pytest_runtest_logreport(self, report):
        # pytest-subtests reports share a nodeid; keep their context separate.
        self.reports.append(dict(nodeid=report.nodeid, phase=report.when,
                                 outcome=report.outcome, seconds=report.duration,
                                 subtest=type(report).__name__ != 'TestReport'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--scope', choices=['solver', 'app'], required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('choose new test result path')
    recorder = Recorder()
    test_args = ['-q', '-p', 'no:cacheprovider', 'solver_v1' if args.scope == 'solver' else 'app']
    start = time.monotonic()
    result = pytest.main(test_args, plugins=[recorder])
    payload = dict(scope=args.scope, arguments=test_args, exit_code=int(result),
                   elapsed_wall_seconds=time.monotonic()-start,
                   passed_calls=sum(r['phase'] == 'call' and r['outcome'] == 'passed' and not r['subtest'] for r in recorder.reports),
                   passed_subtests=sum(r['phase'] == 'call' and r['outcome'] == 'passed' and r['subtest'] for r in recorder.reports),
                   reports=recorder.reports)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2)+'\n', encoding='utf-8')
    raise SystemExit(result)


if __name__ == '__main__':
    main()
