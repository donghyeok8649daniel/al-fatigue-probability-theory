"""Run new physical-mass constrained SW motion from saved nonlinear ensembles."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from solver_v1.silicon_crack_research import RelaxedCoordinates,SpatialSW
from solver_v1.silicon_conditional_research import SI_MASS_AMU,AMU_EV_PS2_A2,KB_EV_K
from solver_v1.silicon_thermal_dynamics import constrained_verlet
from .run_local_crack_audit import save_json
from .run_static_probe import source_parameters


def task(configuration):
    calculation,reference_root,label,output,duration,dt=configuration
    start=time.perf_counter()
    calculation,reference_root,output=map(Path,(calculation,reference_root,output))
    record=json.loads((calculation/(label+'.json')).read_text(encoding='utf-8'))
    sample=np.load(calculation/(label+'.npz'),allow_pickle=False)
    reference=np.load(reference_root/(record['reference']+'.npz'),allow_pickle=False)
    meta=json.loads((reference_root/'summary.json').read_text(encoding='utf-8'))
    p,_=source_parameters()
    model=SpatialSW(p,front_period=meta['geometry']['front_period_A'])
    coordinates=RelaxedCoordinates(reference['reference_positions'],reference['fixed'],bond=reference['bond'])
    x=sample['snapshot_bath_coordinates_A'][-1]
    rng=np.random.default_rng(180000+record['seed'])
    mass=SI_MASS_AMU*AMU_EV_PS2_A2
    velocity=np.sqrt(KB_EV_K*record['temperature_K']/mass)*rng.standard_normal(len(x))
    run=constrained_verlet(model,coordinates,x,velocity,gap=record['gap_A'],
        dt_ps=dt,steps=round(duration/dt),halfwidth=record['bath_box_halfwidth_A'],stride=round(.002/dt))
    name=label+f'_dt{dt:.5f}'
    np.savez_compressed(output/(name+'.npz'),**{k:v for k,v in run.items() if isinstance(v,np.ndarray)},
        initial_bath_coordinates_A=x,initial_bath_velocities_A_ps=velocity)
    scalars={k:v for k,v in run.items() if not isinstance(v,np.ndarray)}
    scalars.update(label=name,preparation_label=label,temperature_K=record['temperature_K'],
        reference=record['reference'],seed=record['seed'],duration_ps=duration,
        independent_preparation_seed=record['seed'],velocity_seed=180000+record['seed'],
        source_sample_sha256=hashlib.sha256((calculation/(label+'.npz')).read_bytes()).hexdigest(),
        source_metadata_sha256=hashlib.sha256((calculation/(label+'.json')).read_bytes()).hexdigest(),
        bath_dimension=len(x),elapsed_seconds=time.perf_counter()-start,
        max_energy_error_over_initial_kinetic=run['max_energy_residual_eV']/run['initial_kinetic_energy_eV'],
        fixed_displacement_A=float(np.max(abs(run['final_positions_A'][reference['fixed']]
            -reference['reference_positions'][reference['fixed']]),initial=0)),
        final_gap_error_A=float(run['final_positions_A'][reference['bond'][1],1]
            -run['final_positions_A'][reference['bond'][0],1]-record['gap_A']),
        conditional_force_correlation_is_exact_nonlinear_memory=False)
    save_json(output/(name+'.json'),scalars)
    return scalars


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ensembles',type=Path,required=True)
    parser.add_argument('--references',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--names',nargs='+',default=['initial','saddle','opened_minimum'])
    parser.add_argument('--temperatures',type=float,nargs='+',default=[100.,300.,600.])
    parser.add_argument('--duration-ps',type=float,default=40.)
    parser.add_argument('--dt-ps',type=float,default=.001)
    parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args()
    if (args.duration_ps<=0 or args.dt_ps<=0 or abs(.002/args.dt_ps-round(.002/args.dt_ps))>1e-10):
        parser.error('positive duration and dt dividing 0.002 ps required')
    args.output.mkdir(parents=True,exist_ok=True)
    if any(args.output.iterdir()):
        parser.error('preserve prior runs; choose an empty output')
    start=time.perf_counter()
    meta=json.loads((args.ensembles/'summary.json').read_text(encoding='utf-8'))
    selected=[row for row in meta['records'] if row['reference'] in args.names and row['temperature_K'] in args.temperatures]
    code=['solver_v1/silicon_thermal_dynamics.py','solver_v1/silicon_crack_research.py',
          'results/silicon_wafer_feasibility/run_thermal_dynamics.py']
    provenance=dict(ensemble_source=args.ensembles.as_posix(),reference_source=args.references.as_posix(),
        configuration={k:v for k,v in vars(args).items() if k not in ('ensembles','references','output')},
        code_sha256={name:hashlib.sha256(Path(name).read_bytes()).hexdigest() for name in code},
        source_summary_sha256=hashlib.sha256((args.ensembles/'summary.json').read_bytes()).hexdigest(),
        initialization='last actual Metropolis sample plus independent Maxwell bath velocities',
        thermostat=None,production_t0_seconds=None,constant_drag_identified=False)
    save_json(args.output/'running_summary.json',provenance)
    configs=[(str(args.ensembles),str(args.references),row['label'],str(args.output),args.duration_ps,args.dt_ps) for row in selected]
    rows=[]
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        for future in as_completed([executor.submit(task,config) for config in configs]):
            row=future.result();rows.append(row)
            save_json(args.output/'completed_runs.json',rows)
            print('MD',len(rows),'/',len(configs),row['label'],'drift',row['max_energy_error_over_initial_kinetic'],
                  'reflections',row['reflection_count'],'seconds',row['elapsed_seconds'],flush=True)
    save_json(args.output/'summary.json',dict(**provenance,records=rows,
        total_new_trajectory_ps=sum(row['duration_ps'] for row in rows),
        elapsed_seconds=time.perf_counter()-start,computation_completed=True))


if __name__=='__main__':
    main()
