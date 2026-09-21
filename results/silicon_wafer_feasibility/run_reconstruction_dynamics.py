"""Actual inertial structural-response probes from explicitly metastable starts.

The sampling box is enlarged for this NONEQUILIBRIUM release. These records do
not estimate canonical memory, a fracture probability or a lifetime from counts.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from solver_v1.silicon_crack_research import RelaxedCoordinates,SpatialSW,relax_atoms
from solver_v1.silicon_conditional_research import SI_MASS_AMU,AMU_EV_PS2_A2,KB_EV_K
from solver_v1.silicon_thermal_dynamics import constrained_verlet
from .run_local_crack_audit import save_json,polish
from .run_static_probe import source_parameters


def task(config):
    references,ensembles,reconstructed,temperature,seed,duration,output=config
    references,ensembles,reconstructed,output=map(Path,(references,ensembles,reconstructed,output))
    start=time.perf_counter()
    ref=np.load(references/'initial.npz',allow_pickle=False)
    branch=np.load(reconstructed,allow_pickle=False)
    meta=json.loads((references/'summary.json').read_text(encoding='utf-8'))
    label=f'initial_T{round(temperature)}_seed{seed}'
    sample=np.load(ensembles/(label+'.npz'),allow_pickle=False)
    source=json.loads((ensembles/(label+'.json')).read_text(encoding='utf-8'))
    coords=RelaxedCoordinates(ref['reference_positions'],ref['fixed'],bond=ref['bond'])
    gap=float(ref['gap_A']);direction=coords.encode(branch['initial']);length=np.linalg.norm(direction)
    direction/=length
    p,_=source_parameters();model=SpatialSW(p,front_period=meta['geometry']['front_period_A'])
    x=sample['snapshot_bath_coordinates_A'][-1]
    # Audit the start's basin BEFORE calling it an original-branch release.
    quench,qinfo=relax_atoms(model,coords,initial=coords.positions(x,[gap]),values=[gap],tolerance=2e-6)
    quench,qcorrection=polish(model,coords,quench,ref['fixed'],values=[gap])
    initial_basin_distance=float(np.linalg.norm(quench-ref['positions']))
    if initial_basin_distance>1e-4:
        raise RuntimeError('selected preparation is not in the original conditional basin')
    mass=SI_MASS_AMU*AMU_EV_PS2_A2
    rng=np.random.default_rng(290000+seed)
    velocity=np.sqrt(KB_EV_K*temperature/mass)*rng.standard_normal(len(x))
    order=[];snapshots=[];snapshot_steps=[]
    class Instrumented:
        def evaluate(self,positions):
            bath=coords.encode(positions)
            index=len(order)
            order.append(float((direction@bath)/length))
            if index%1000==0:
                snapshots.append(positions.copy());snapshot_steps.append(index)
            return model.evaluate(positions)
    # R=2 A leaves room for the verified new minimum at max|x|=1.1544 A.
    run=constrained_verlet(Instrumented(),coords,x,velocity,gap=gap,dt_ps=.001,
        steps=round(duration/.001),halfwidth=2.,stride=2)
    final_minimum,final_info=relax_atoms(model,coords,initial=run['final_positions_A'],values=[gap],tolerance=2e-6)
    final_minimum,final_correction=polish(model,coords,final_minimum,ref['fixed'],values=[gap])
    np.savez_compressed(output/(label+'.npz'),**{k:v for k,v in run.items() if isinstance(v,np.ndarray)},
        initial_bath_coordinates_A=x,initial_bath_velocities_A_ps=velocity,
        reconstruction_fraction_every_step=np.asarray(order),snapshot_positions=np.asarray(snapshots),
        snapshot_steps=snapshot_steps,final_quenched_positions=final_minimum)
    row={k:v for k,v in run.items() if not isinstance(v,np.ndarray)}
    row.update(label=label,temperature_K=temperature,seed=seed,duration_ps=duration,
        original_preparation_halfwidth_A=source['bath_box_halfwidth_A'],dynamics_halfwidth_A=2.,
        initial_basin_distance_A=initial_basin_distance,initial_quench=qinfo,initial_correction=qcorrection,
        final_quench=final_info,final_correction=final_correction,
        final_quenched_energy_from_original_eV=float(np.sum(model.evaluate(final_minimum).site_energy-ref['sites'])),
        final_quenched_distance_from_original_A=float(np.linalg.norm(final_minimum-ref['positions'])),
        final_quenched_distance_from_reconstructed_A=float(np.linalg.norm(final_minimum-branch['initial'])),
        reconstruction_fraction_min=float(min(order)),reconstruction_fraction_max=float(max(order)),
        reconstruction_fraction_final=float(order[-1]),
        elapsed_seconds=time.perf_counter()-start,sample_source_sha256=hashlib.sha256((ensembles/(label+'.npz')).read_bytes()).hexdigest(),
        canonical_stationarity_claimed=False,physical_failure_probability=None,transition_rate_from_counts=None)
    save_json(output/(label+'.json'),row);return row


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--references',type=Path,required=True)
    parser.add_argument('--ensembles',type=Path,required=True)
    parser.add_argument('--reconstructed',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--duration-ps',type=float,default=100.)
    parser.add_argument('--workers',type=int,default=2)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    if any(args.output.iterdir()):parser.error('use an empty output')
    start=time.perf_counter()
    provenance=dict(temperatures_K=[300.,600.],seeds=[211,307],duration_ps=args.duration_ps,
        interpretation='metastable original-basin starts, larger conservative domain; structural reference dynamics only',
        production_t0_seconds=None,physical_failure_probability=None,
        code_sha256={name:hashlib.sha256(Path(name).read_bytes()).hexdigest() for name in (
            'solver_v1/silicon_thermal_dynamics.py','results/silicon_wafer_feasibility/run_reconstruction_dynamics.py')})
    save_json(args.output/'running_summary.json',provenance)
    configs=[(str(args.references),str(args.ensembles),str(args.reconstructed),t,seed,args.duration_ps,str(args.output))
             for t in [300.,600.] for seed in [211,307]]
    rows=[]
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        for future in as_completed([executor.submit(task,config) for config in configs]):
            row=future.result();rows.append(row);save_json(args.output/'completed_runs.json',rows)
            print('release',row['label'],row['final_quenched_energy_from_original_eV'],
                row['reconstruction_fraction_max'],'reflections',row['reflection_count'],flush=True)
    save_json(args.output/'summary.json',dict(**provenance,records=rows,
        total_new_trajectory_ps=sum(row['duration_ps'] for row in rows),elapsed_seconds=time.perf_counter()-start,
        computation_completed=True))


if __name__=='__main__':main()
