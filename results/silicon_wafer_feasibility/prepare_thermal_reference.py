"""Reevaluate conditional SW branches and cache exact Gaussian references."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from solver_v1.silicon_crack_research import RelaxedCoordinates, SpatialSW, relax_atoms
from .run_conditional_dynamics import reduction
from .run_local_crack_audit import save_json, polish
from .run_static_probe import source_parameters


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path('results/silicon_local_crack_v2/front4'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--profile-stride', type=int, default=2)
    parser.add_argument('--profile-stop', type=int, default=24)
    parser.add_argument('--stationary-only', action='store_true')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if any(args.output.iterdir()):
        parser.error('use an empty reference directory')
    start = time.perf_counter()
    source = args.source
    geometry = np.load(source/'geometry.npz', allow_pickle=False)
    states = np.load(source/'localized/stationary_positions.npz', allow_pickle=False)
    profile = np.load(source/'localized/profile_positions.npz', allow_pickle=False)
    meta = json.loads((source/'geometry.json').read_text(encoding='utf-8'))
    parameters, parameter_hash = source_parameters()
    model = SpatialSW(parameters, front_period=meta['front_period_A'])
    fixed, bond = geometry['fixed'], states['bond']
    free = np.flatnonzero(~fixed)
    coordinates = RelaxedCoordinates(states['initial'], fixed, bond=bond)
    reference_sites = model.evaluate(states['initial']).site_energy.copy()
    rows, saved_positions, relaxations = [], {}, {}
    points = {name:states[name] for name in ('initial','saddle','opened_minimum')}
    if not args.stationary_only:
        for index in range(args.profile_stride, min(len(profile['q_A']), args.profile_stop+1), args.profile_stride):
            points[f'profile{index:02d}'] = profile['positions'][index]
        for delta in (-.12, -.06):
            q = profile['q_A'][0]+delta
            position, result = relax_atoms(model, coordinates, values=[q], tolerance=1e-8)
            if not result['force_converged']:
                position, correction = polish(model, coordinates, position, fixed, values=[q])
                result['newton_correction'] = correction
            label = f'lower{round(-100*delta):02d}'
            points[label] = position
            relaxations[label] = result
    for name, r in points.items():
        gap = float(r[bond[1],1]-r[bond[0],1])
        residual = np.max(abs(coordinates.pullback(model.evaluate(r).gradient)[0]))
        if residual > 1e-8:
            r, correction = polish(model, coordinates, r, fixed, values=[gap])
            relaxations.setdefault(name, {})['source_polish'] = correction
        evaluation = model.evaluate(r)
        mean = coordinates.encode(r)
        h = model.hessian(r, free_atoms=free)
        one, b, d = reduction(h, r, fixed, bond, (1,))
        gradient, reaction = coordinates.pullback(evaluation.gradient)
        if np.max(abs(gradient)) > 2e-7:
            raise RuntimeError('source conditional stationary residual too large')
        np.testing.assert_allclose(coordinates.positions(mean, [gap]), r, rtol=0, atol=2e-14)
        np.savez_compressed(args.output/(name+'.npz'), positions=r, fixed=fixed, bond=bond,
            reference_positions=states['initial'], mean=mean, gap_A=gap,
            bath_hessian=one.bath_hessian, response=one.bath_response[:,0],
            cross=one.cross[:,0], bare_curvature=one.bare_curvature,
            relaxed_curvature=one.relaxed_curvature, sites=evaluation.site_energy,
            full_hessian=h.toarray() if name == 'initial' else np.empty((0,0)))
        row = dict(name=name, gap_A=gap, bath_dimension=len(mean),
            energy_difference_eV=float(np.sum(evaluation.site_energy-reference_sites)),
            logdet_bath=one.logdet_bath, static_reaction_eV_A=float(reaction[0]),
            force_residual_eV_A=float(np.max(abs(gradient))),
            relaxed_curvature_eV_A2=float(one.relaxed_curvature[0,0]),
            bath_mean_max_coordinate_A=float(abs(mean).max()))
        rows.append(row); saved_positions[name] = r
        save_json(args.output/'running.json', rows)
        print(row, flush=True)
    np.savez_compressed(args.output/'branch_positions.npz', **saved_positions,
        reference=states['initial'], fixed=fixed, bond=bond)
    save_json(args.output/'summary.json', dict(source=source.as_posix(), geometry=meta,
        source_parameter_sha256=parameter_hash, rows=rows, relaxations=relaxations,
        input_sha256={p.as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in (
            source/'geometry.npz', source/'geometry.json', source/'localized/stationary_positions.npz',
            source/'localized/profile_positions.npz')},
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        elapsed_seconds=time.perf_counter()-start))


if __name__ == '__main__':
    main()
