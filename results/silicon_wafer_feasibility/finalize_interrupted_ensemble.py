"""Record completed chains of a stopped job without claiming a completed run."""
import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path

from .run_local_crack_audit import save_json


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--calculation',type=Path,required=True)
    args=parser.parse_args();out=args.calculation
    if (out/'summary.json').exists() or (out/'interrupted_summary.json').exists():
        parser.error('a final status already exists')
    meta=json.loads((out/'running_summary.json').read_text(encoding='utf-8'))
    records=json.loads((out/'completed_runs.json').read_text(encoding='utf-8'))
    names=meta['configuration']['names'] or [r['name'] for r in meta['reference_summary']['rows']]
    expected={f'{name}_T{round(t)}_seed{seed}' for name in names for t in meta['configuration']['temperatures']
              for seed in meta['configuration']['seeds']}
    saved={r['label'] for r in records}
    if not saved<=expected:raise ValueError('saved labels differ from the declared calculation')
    for row in records:
        if not (out/(row['label']+'.npz')).exists() or not (out/(row['label']+'.json')).exists():
            raise ValueError('incomplete raw files in the completed list')
    save_json(out/'interrupted_summary.json',dict(**meta,records=records,
        recorded_at_utc=datetime.now(timezone.utc).isoformat(),
        computation_completed=False,execution_status='stopped_after_budget_expiry',
        reason='clock check found budget expired after a long unverified wall-clock gap; task-owned remaining workers stopped',
        budget_start_utc='2026-09-20T19:11:55Z',budget_deadline_utc='2026-09-20T23:11:55Z',
        completed_chain_count=len(records),planned_chain_count=len(expected),
        uncompleted_labels=sorted(expected-saved),
        total_completed_draws=sum(r['draws'] for r in records),
        finalizer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
    print('saved complete chains',len(saved),'of',len(expected),flush=True)


if __name__=='__main__':main()
