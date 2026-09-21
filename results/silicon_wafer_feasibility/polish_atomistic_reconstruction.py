"""Force-based correction and conditional Hessian audit of Tersoff endpoints.

Keeps the first L-BFGS results (including their failed force tolerances) intact.
The two finite-difference steps test local curvature, not material accuracy.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy.linalg import eigh, solve

from solver_v1.silicon_atomistic_reference import LammpsSilicon
from solver_v1.silicon_crack_research import RelaxedCoordinates
from .run_atomistic_controls import sources, checksum, save_json


def hessian(force, x, step):
    columns = []
    for j in range(len(x)):
        d = np.zeros_like(x)
        d[j] = step
        columns.append((force(x+d)-force(x-d))/(2*step))
    h = np.column_stack(columns)
    skew = float(np.max(abs(h-h.T)))
    return (h+h.T)/2, skew


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--calculation', type=Path, default=Path('results/silicon_atomistic_v5'))
    args = parser.parse_args()
    output = args.calculation/'reconstruction_polish.json'
    if output.exists():
        raise FileExistsError(output)
    started = time.perf_counter()
    specs, _ = sources(args.calculation)
    style, potential, source_hash = specs['tersoff_1989']
    source = args.calculation/'tersoff_1989_reconstruction.npz'
    with np.load(source) as data:
        states = {k:data[k].copy() for k in data.files}
    reference, fixed, bond = states['state00_initial'], states['fixed'], states['bond']
    origin, cell = states['origin'], states['cell']
    coordinates = RelaxedCoordinates(reference, fixed, bond=bond)
    gap = float(reference[bond[1], 1]-reference[bond[0], 1])
    rows, saved = [], {}
    with LammpsSilicon(style, potential, sha256=source_hash) as engine:
        def evaluate(x):
            return engine.evaluate(coordinates.positions(x, [gap])-origin, cell,
                                   periodic=(False, False, True))

        def force(x):
            return coordinates.pullback(evaluate(x).gradient)[0]

        for name in ['state00', 'state04', 'state01', 'state02']:
            x = coordinates.encode(states[name+'_relaxed'])
            residuals = [float(np.max(abs(force(x))))]
            hc, skew = hessian(force, x, 1e-4)
            coarse_eigen = eigh(hc, subset_by_index=(0, 4), eigvals_only=True)
            if coarse_eigen[0] <= 0:
                raise RuntimeError('endpoint has a nonpositive conditional Hessian')
            for _ in range(3):
                if residuals[-1] < 1e-9:
                    break
                x -= solve(hc, force(x), assume_a='pos')
                residuals.append(float(np.max(abs(force(x)))))
            if residuals[-1] >= 1e-8:
                raise RuntimeError(f'force correction failed: {name}, {residuals}')
            # Both steps at the SAME corrected endpoint; the first pre-correction
            # Hessian is recorded only for the Newton correction.
            h1, skew1 = hessian(force, x, 1e-4)
            h2, skew2 = hessian(force, x, 5e-5)
            values, vectors = eigh(h2, subset_by_index=(0, 4))
            coarse_final = eigh(h1, subset_by_index=(0, 4), eigvals_only=True)
            current = evaluate(x)
            g, reaction = coordinates.pullback(current.gradient)
            row = dict(name=name, residual_history_eV_A=residuals,
                       correction_hessian_skew_eV_A2=skew,
                       steps_A=[1e-4, 5e-5], skew_eV_A2=[skew1, skew2],
                       lowest_eigenvalues_eV_A2=values.tolist(),
                       coarse_eigenvalues_eV_A2=coarse_final.tolist(),
                       eigenpair_residuals_eV_A2=np.linalg.norm(h2 @ vectors-vectors*values, axis=0).tolist(),
                       hessian_step_max_difference_eV_A2=float(np.max(abs(h1-h2))),
                       energy_eV=current.energy, reaction_eV_A=reaction.tolist(),
                       final_force_max_eV_A=float(np.max(abs(g))))
            rows.append(row)
            saved[name] = coordinates.positions(x, [gap])
            saved[name+'_sites'] = current.site_energy
            print(json.dumps(row), flush=True)
    for row in rows:
        row['energy_vs_state00_eV'] = float(np.sum(saved[row['name']+'_sites']-saved['state00_sites']))
    saved.update(fixed=fixed, bond=bond, cell=cell, origin=origin)
    np.savez_compressed(args.calculation/'reconstruction_polished.npz', **saved)
    save_json(output, dict(rows=rows, source_sha256=checksum(source),
              code_sha256=checksum(__file__), elapsed_seconds=time.perf_counter()-started,
              scope='fixed-q bath minima of ordinary Tersoff, same SW grips; not full-q stability',
              complete=True, material_calibrated=False))


if __name__ == '__main__':
    main()
