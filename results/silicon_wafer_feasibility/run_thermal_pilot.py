"""Fresh joint Gaussian importance pilot, with an explicit common bath box.

This is a sampling/overlap diagnostic, not certified Si thermodynamics.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from solver_v1.silicon_crack_research import RelaxedCoordinates, SpatialSW
from solver_v1.silicon_thermal_research import GaussianReference, exponential_reweight
from .run_conditional_dynamics import reduction
from .run_local_crack_audit import save_json
from .run_static_probe import source_parameters


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--samples', type=int, default=512)
    parser.add_argument('--seed', type=int, default=20260921)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if any(args.output.iterdir()):
        parser.error('preserve prior output; choose an empty folder')
    start = time.perf_counter()
    source = Path('results/silicon_local_crack_v2/front4')
    geometry = np.load(source/'geometry.npz', allow_pickle=False)
    states = np.load(source/'localized/stationary_positions.npz', allow_pickle=False)
    meta = json.loads((source/'geometry.json').read_text(encoding='utf-8'))
    parameters, parameter_hash = source_parameters()
    model = SpatialSW(parameters, front_period=meta['front_period_A'])
    coordinates = RelaxedCoordinates(states['initial'], geometry['fixed'], bond=states['bond'])
    free = np.flatnonzero(~geometry['fixed'])
    rng = np.random.default_rng(args.seed)
    standard = rng.standard_normal((args.samples, coordinates.dimension))
    np.savez_compressed(args.output/'standard_normals.npz', values=standard)
    rows = []
    for name in ('initial', 'saddle', 'opened_minimum'):
        r = states[name]
        sites = model.evaluate(r).site_energy.copy()
        one, _, _ = reduction(model.hessian(r, free_atoms=free), r, geometry['fixed'], states['bond'], (1,))
        mean = coordinates.encode(r)
        gap = float(r[states['bond'][1], 1]-r[states['bond'][0], 1])
        for temperature in (100., 300., 600.):
            reference = GaussianReference.from_hessian(one.bath_hessian, temperature, mean)
            proposals = reference.transform(standard)
            residual, reactions = [], []
            for x, u in zip(proposals, standard):
                evaluation = model.evaluate(coordinates.positions(x, [gap]))
                residual.append(float(np.sum(evaluation.site_energy-sites)-.5*reference.thermal_energy*(u@u)))
                reactions.append(float(coordinates.pullback(evaluation.gradient)[1][0]))
            norms = np.max(abs(proposals), axis=1)
            for radius in (.75, 1., 1.5):
                report = exponential_reweight(residual, temperature, inside=norms < radius, observables=reactions)
                del report['normalized_weights']
                rows.append(dict(state=name, temperature_K=temperature, bath_box_halfwidth_A=radius,
                    gaussian_max_coordinate_A=float(norms.max()), **report))
                print(rows[-1], flush=True)
            np.savez_compressed(args.output/f'{name}_{int(temperature)}K.npz',
                residual_eV=residual, reaction_eV_A=reactions, bath_max_coordinate_A=norms)
            save_json(args.output/'running.json', rows)
    save_json(args.output/'summary.json', dict(rows=rows, seed=args.seed, samples=args.samples,
        source_parameter_sha256=parameter_hash, elapsed_seconds=time.perf_counter()-start,
        code_sha256={name:hashlib.sha256(Path(name).read_bytes()).hexdigest() for name in (
            'solver_v1/silicon_thermal_research.py', 'results/silicon_wafer_feasibility/run_thermal_pilot.py')},
        production_clock_seconds=None, finite_T_PMF_certified=False))


if __name__ == '__main__':
    main()
