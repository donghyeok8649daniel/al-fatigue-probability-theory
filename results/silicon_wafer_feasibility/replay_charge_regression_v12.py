"""Verify full v9 synthetic audit outputs after stricter v12 validation."""
import argparse,hashlib,json
from pathlib import Path


def main(args):
    records=[]
    for name in ('charge_rate_convergence.csv','spatial_grid_convergence.csv',
            'absorption_convergence.csv','quenched_mixture_counterexample.csv','dynamics_summary.json'):
        old=(args.baseline/name).read_bytes();new=(args.result/name).read_bytes()
        record=dict(file=name,bytes=len(new),baseline_sha256=hashlib.sha256(old).hexdigest(),
            rerun_sha256=hashlib.sha256(new).hexdigest(),byte_equal=old==new)
        if not record['byte_equal']:raise ValueError('full regression output changed: '+name)
        records.append(record)
    output=args.result/'baseline_replay.json'
    if output.exists():raise ValueError('existing replay must be preserved')
    report=dict(complete=True,byte_equal_files=len(records),files=records,
        new_charge_source_sha256=hashlib.sha256(args.charge_source.read_bytes()).hexdigest(),
        meaning='actual full synthetic audit rerun; validates retained numerical results, not real Si kinetics',
        new_MACE=0,new_DFT=0,new_MD=0,physical_clock=None)
    output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8');print(json.dumps(report,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('baseline','result','charge-source'):p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
