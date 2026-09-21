"""Fresh Metropolized full-SW conditional ensembles on a common finite box.

No time or fracture probability is inferred from MCMC proposal counts.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np
from scipy.linalg import solve_triangular

from solver_v1.silicon_crack_research import RelaxedCoordinates, SpatialSW
from solver_v1.silicon_thermal_research import (
    GaussianReference, pcn_sample, harmonic_hmc_sample, bounded_force_control, block_statistics,
)
from .run_local_crack_audit import save_json
from .run_static_probe import source_parameters


def task(config):
    start = time.perf_counter()
    source, name, temperature, seed, output, draws, warmup, scale, radius, sampler, hmc_step = config
    source, output = Path(source), Path(output)
    reference_data = np.load(source/(name+'.npz'), allow_pickle=False)
    metadata = json.loads((source/'summary.json').read_text(encoding='utf-8'))
    parameters, parameter_hash = source_parameters()
    model = SpatialSW(parameters, front_period=metadata['geometry']['front_period_A'])
    coordinates = RelaxedCoordinates(reference_data['reference_positions'], reference_data['fixed'],
        bond=reference_data['bond'])
    gaussian = GaussianReference.from_hessian(reference_data['bath_hessian'], temperature, reference_data['mean'])
    gap = float(reference_data['gap_A'])
    free = ~reference_data['fixed']
    response = reference_data['response']
    sites = reference_data['sites']
    rejected_outside = [0]
    evaluations = [0]

    def evaluate(u):
        x = gaussian.transform(u)
        maximum = float(abs(x).max())
        if maximum >= radius:
            rejected_outside[0] += 1
            return np.inf, None, None
        result = model.evaluate(coordinates.positions(x, [gap]))
        gradient, reactions = coordinates.pullback(result.gradient)
        energy = float(np.sum(result.site_energy-sites))
        residual = energy/gaussian.thermal_energy-.5*(u@u)
        control = bounded_force_control(x, gradient, response, temperature, radius)
        evaluations[0] += 1
        grad_phi = (solve_triangular(gaussian.lower,gradient,lower=True,check_finite=False)
                    /np.sqrt(gaussian.thermal_energy)-u) if sampler=='hmc' else None
        return residual, grad_phi, np.array([reactions[0], reactions[0]-control, energy,
            residual, maximum, float(np.mean((x-reference_data['mean'])**2))])

    rng = np.random.default_rng(seed)
    # A fresh Gaussian start, not the same zero vector for every replica.
    initial = rng.standard_normal(coordinates.dimension)
    while not np.isfinite(evaluate(initial)[0]):
        initial = rng.standard_normal(coordinates.dimension)
    if sampler=='hmc':
        result = harmonic_hmc_sample(evaluate,initial,draws=draws,warmup=warmup,rng=rng,step=hmc_step)
    else:
        def pcn_evaluate(u):
            phi,_,observation=evaluate(u)
            return phi,observation
        result = pcn_sample(pcn_evaluate, initial, proposal_scale=scale, draws=draws,
            warmup=warmup, rng=rng)
    label = f'{name}_T{round(temperature)}_seed{seed}'
    snapshot_indices = np.unique(np.r_[np.arange(0, draws, max(1,draws//16)), draws-1])
    snapshots = gaussian.transform(result['samples'][snapshot_indices])
    np.savez_compressed(output/(label+'.npz'), observations=result['observations'],
        potentials=result['potentials'], accepted_at_saved_step=result['accepted_at_saved_step'],
        snapshot_bath_coordinates_A=snapshots, snapshot_indices=snapshot_indices,
        initial_whitened=initial, final_whitened=result['samples'][-1],
        final_positions_A=coordinates.positions(snapshots[-1], [gap]),
        hmc_energy_errors=result.get('hamiltonian_errors',np.empty(0)))
    statistics = block_statistics(result['observations'], blocks=16)
    record = dict(label=label, reference=name, temperature_K=temperature, seed=seed,
        gap_A=gap, draws=draws, warmup=warmup, proposal_scale=scale, sampler=sampler,
        hmc_step=hmc_step if sampler=='hmc' else None,
        bath_box_halfwidth_A=radius, box_center='initial-state orthonormal bath coordinates; independent of q',
        acceptance_fraction=result['acceptance_fraction'],
        warmup_acceptance_fraction=result['warmup_acceptance_fraction'],
        domain_rejected_proposals=rejected_outside[0], force_evaluations=evaluations[0],
        observations=['raw_dU_dq_eV_A','controlled_dU_dq_eV_A','excess_U_eV',
            'anharmonic_residual_over_kBT','max_abs_bath_coordinate_A','bath_mean_square_displacement_A2'],
        statistics={key:np.asarray(value).tolist() for key,value in statistics.items()},
        raw_to_controlled_variance_ratio=float(np.var(result['observations'][:,0])/np.var(result['observations'][:,1])),
        fixed_displacement_A=float(np.max(abs(coordinates.positions(snapshots[-1], [gap])[~free]
            -reference_data['reference_positions'][~free]),initial=0)),
        final_rng_state=result['final_rng_state'],
        source_parameter_sha256=parameter_hash, reference_sha256=hashlib.sha256((source/(name+'.npz')).read_bytes()).hexdigest(),
        elapsed_seconds=time.perf_counter()-start,
        sample_clock_ps=None, finite_T_global_PMF_certified=False, physical_mobility=None)
    save_json(output/(label+'.json'), record)
    return {k:v for k,v in record.items() if k not in ('statistics','final_rng_state')}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--references', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--temperatures', type=float, nargs='+', default=[100.,300.,600.])
    parser.add_argument('--names', nargs='+')
    parser.add_argument('--seeds', type=int, nargs='+', default=[103,211,307,401])
    parser.add_argument('--draws', type=int, default=4096)
    parser.add_argument('--warmup', type=int, default=1024)
    parser.add_argument('--proposal-scale', type=float, default=.6)
    parser.add_argument('--halfwidth', type=float, default=1.)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--sampler',choices=['pcn','hmc'],default='hmc')
    parser.add_argument('--hmc-step',type=float,default=.2)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if any(args.output.iterdir()):
        parser.error('preserve prior runs; output must be empty')
    start = time.perf_counter()
    reference_meta = json.loads((args.references/'summary.json').read_text(encoding='utf-8'))
    names = args.names or [r['name'] for r in reference_meta['rows']]
    code = ['solver_v1/silicon_thermal_research.py','solver_v1/silicon_crack_research.py',
            'results/silicon_wafer_feasibility/run_conditional_ensemble.py']
    metadata = dict(schema='silicon-conditional-ensemble/4',
        head=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        source_references=args.references.as_posix(), reference_summary=reference_meta,
        configuration={k:v for k,v in vars(args).items() if k not in ('references','output')},
        code_sha256={name:hashlib.sha256(Path(name).read_bytes()).hexdigest() for name in code},
        sample_clock_ps=None, production_enabled=False, production_t0_seconds=None)
    save_json(args.output/'running_summary.json', metadata)
    configs = [(str(args.references),name,t,seed,str(args.output),args.draws,args.warmup,
                args.proposal_scale,args.halfwidth,args.sampler,args.hmc_step)
               for name in names for t in args.temperatures for seed in args.seeds]
    records = []
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = [executor.submit(task, config) for config in configs]
        for future in as_completed(futures):
            row = future.result(); records.append(row)
            save_json(args.output/'completed_runs.json', records)
            print('completed',len(records),'/',len(configs),row['label'],
                  'accept',row['acceptance_fraction'],'seconds',row['elapsed_seconds'],flush=True)
    save_json(args.output/'summary.json',dict(**metadata, records=records,
        total_draws=sum(r['draws'] for r in records), elapsed_seconds=time.perf_counter()-start,
        computation_completed=True))


if __name__ == '__main__':
    main()
