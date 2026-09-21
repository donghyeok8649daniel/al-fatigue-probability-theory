"""Quench time snapshots to identify structures, never to count fatigue failures."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import hashlib
import json
from pathlib import Path

import numpy as np

from solver_v1.silicon_crack_research import RelaxedCoordinates,SpatialSW,relax_atoms
from .run_local_crack_audit import save_json
from .run_static_probe import source_parameters


def task(config):
    calculation,references,minima,label=config
    calculation,references,minima=map(Path,(calculation,references,minima))
    ref=np.load(references/'initial.npz',allow_pickle=False)
    source=np.load(calculation/(label+'.npz'),allow_pickle=False)
    known=np.load(minima,allow_pickle=False)
    metadata=json.loads((references/'summary.json').read_text(encoding='utf-8'))
    p,_=source_parameters();model=SpatialSW(p,front_period=metadata['geometry']['front_period_A'])
    coords=RelaxedCoordinates(ref['reference_positions'],ref['fixed'],bond=ref['bond'])
    gap=float(ref['gap_A']);rows=[];quenches=[]
    for i in range(0,len(source['snapshot_steps']),2):
        r,info=relax_atoms(model,coords,initial=source['snapshot_positions'][i],values=[gap],tolerance=2e-6)
        if not info['force_converged']:raise RuntimeError('snapshot quench did not converge')
        distances={name:float(np.linalg.norm(r-known[name])) for name in known.files if name.startswith('state')}
        nearest=min(distances,key=distances.get)
        rows.append(dict(time_ps=float(source['snapshot_steps'][i])*.001,
            energy_from_original_eV=float(np.sum(model.evaluate(r).site_energy-ref['sites'])),
            nearest_saved_minimum=nearest,distance_from_saved_minimum_A=distances[nearest],
            free_gradient_max_eV_A=info['free_gradient_max_eV_A'],
            distance_from_original_A=float(np.linalg.norm(r-ref['positions']))))
        quenches.append(r)
    np.savez_compressed(calculation/(label+'_quenches.npz'),positions=np.asarray(quenches),
        time_ps=[r['time_ps'] for r in rows])
    save_json(calculation/(label+'_structural_analysis.json'),dict(rows=rows,
        source_sha256=hashlib.sha256((calculation/(label+'.npz')).read_bytes()).hexdigest(),
        threshold_for_failure=None,physical_transition_rate_from_counts=None))
    return dict(label=label,rows=rows)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--calculation',type=Path,required=True)
    parser.add_argument('--references',type=Path,required=True)
    parser.add_argument('--minima',type=Path,required=True)
    parser.add_argument('--workers',type=int,default=2)
    args=parser.parse_args()
    summary=json.loads((args.calculation/'summary.json').read_text(encoding='utf-8'))
    configs=[(str(args.calculation),str(args.references),str(args.minima),r['label']) for r in summary['records']]
    rows=[]
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        for future in as_completed([executor.submit(task,c) for c in configs]):
            row=future.result();rows.append(row)
            print(row['label'],'quench energy range',min(r['energy_from_original_eV'] for r in row['rows']),
                max(r['energy_from_original_eV'] for r in row['rows']),flush=True)
    save_json(args.calculation/'structural_analysis.json',dict(records=rows,
        scope='deterministic fixed-q quench of snapshots every 2 ps; closest-state distances retained without a failure threshold',
        failure_probability=None,transition_rate_from_counts=None,
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))


if __name__=='__main__':main()
