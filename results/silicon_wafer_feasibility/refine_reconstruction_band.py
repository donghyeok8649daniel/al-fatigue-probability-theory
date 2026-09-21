"""Refine a saved band using a smaller numerical FIRE step, preserving history."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

from .run_local_crack_audit import save_json


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--references',type=Path,required=True)
    parser.add_argument('--band',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    if any(args.output.iterdir()):parser.error('use an empty output')
    data=np.load(args.band,allow_pickle=False)
    seed=args.output/'seed.npz'
    np.savez_compressed(seed,**{f'forward_{i:03d}':r for i,r in enumerate(data['positions'])})
    command=[sys.executable,'-u','-m','results.silicon_wafer_feasibility.run_reconstruction_neb',
        '--references',str(args.references),'--path',str(seed),'--output',str(args.output/'band'),
        '--images',str(len(data['positions'])),'--steps','3200','--fire-step','.001',
        '--fire-max-step','.01','--max-image-move','.01','--climb-after','200']
    with (args.output/'run.log').open('w',encoding='utf-8') as stream:
        result=subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT,check=False)
    report=dict(returncode=result.returncode,source=args.band.as_posix(),
        source_sha256=hashlib.sha256(args.band.read_bytes()).hexdigest(),
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    if (args.output/'band/summary.json').exists():
        summary=json.loads((args.output/'band/summary.json').read_text(encoding='utf-8'))
        report.update(band_converged=summary['band_converged'],iterations=summary['iterations'],
            neb_force_max=summary['neb_force_max'],saddle_audit=summary['saddle_audit'])
    save_json(args.output/'summary.json',report)
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':main()
