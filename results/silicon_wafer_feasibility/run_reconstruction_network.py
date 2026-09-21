"""Independent short bands between already audited conditional minima.

This is a finite list of candidate edges, not an exhaustive reaction network.
"""
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
    parser.add_argument('--minima',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--edges',nargs='+',default=['state04:state01','state01:state02'])
    parser.add_argument('--images',type=int,default=21)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    if any(args.output.iterdir()):parser.error('preserve previous results; use an empty output')
    states=np.load(args.minima,allow_pickle=False)
    rows=[]
    for edge in args.edges:
        left,right=edge.split(':')
        seed=args.output/(left+'_'+right+'.npz')
        np.savez_compressed(seed,forward_000=states[left],forward_001=states[right])
        output=args.output/(left+'_'+right)
        command=[sys.executable,'-u','-m','results.silicon_wafer_feasibility.run_reconstruction_neb',
            '--references',str(args.references),'--path',str(seed),'--output',str(output),
            '--images',str(args.images),'--steps','2400','--fire-step','.005',
            '--fire-max-step','.04','--max-image-move','.03','--climb-after','350']
        with (args.output/(left+'_'+right+'.log')).open('w',encoding='utf-8') as stream:
            completed=subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT,check=False)
        row=dict(initial=left,final=right,returncode=completed.returncode)
        if (output/'summary.json').exists():
            summary=json.loads((output/'summary.json').read_text(encoding='utf-8'))
            row.update(band_converged=summary['band_converged'],saddle_audit=summary['saddle_audit'])
        rows.append(row)
        save_json(args.output/'running.json',rows)
        print(left,right,'exit',completed.returncode,'converged',row.get('band_converged'),flush=True)
    save_json(args.output/'summary.json',dict(edges=rows,
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        source_minima_sha256=hashlib.sha256(args.minima.read_bytes()).hexdigest(),
        graph_exhaustive=False,physical_transition_rates=None))


if __name__=='__main__':main()
