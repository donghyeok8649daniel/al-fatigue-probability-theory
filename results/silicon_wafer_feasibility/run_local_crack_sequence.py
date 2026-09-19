"""Verify a complete four-repeat front advance through local stationary hops.

Each step frees ALL unconstrained atoms, including previously opened bonds.
Nothing prevents healing. These are static connected states, not trajectories
at a physical temperature, frequency or speed.
"""
import argparse
import json
from pathlib import Path
import time

import numpy as np

from solver_v1.silicon_crack_research import SpatialSW, diamond_crack_strip
from .run_local_crack_audit import energy_difference, profile, save_json
from .run_static_probe import source_parameters


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--calculation', type=Path, default=Path('results/silicon_local_crack_v2'))
    args = parser.parse_args()
    started = time.perf_counter()
    root = args.calculation
    metadata = json.loads((root/'running_summary.json').read_text(encoding='utf-8'))
    p, _ = source_parameters()
    strip = diamond_crack_strip(metadata['lattice_A'], nx=6, ny=4, nz=4)
    model = SpatialSW(p, front_period=strip.front_period)
    first = json.loads((root/'front4/localized/summary.json').read_text(encoding='utf-8'))
    data = np.load(root/'front4/localized/stationary_positions.npz')
    initial, current = data['initial'], data['opened_minimum']
    selected = first['states']['saddle']['frontier']['bond_index']
    middle = strip.positions[strip.crossing_bonds].mean(axis=1)
    front = np.flatnonzero(np.isclose(middle[:, 0], middle[selected, 0]))
    front = front[np.argsort((middle[front, 2]-middle[selected, 2]) % strip.front_period)]
    if len(front) != 4 or front[0] != selected:
        raise RuntimeError('sequence requires the audited four-bond front')
    output = root/'sequence'
    output.mkdir()
    rows = [dict(step=1, bond_index=int(selected),
        saddle_energy_from_initial_eV=first['forward_barrier_eV'],
        endpoint_energy_from_initial_eV=energy_difference(model, current, initial),
        forward_barrier_eV=first['forward_barrier_eV'],
        reverse_barrier_eV=first['reverse_barrier_eV'])]
    for step, index in enumerate(front[1:], 2):
        summary, states = profile(model, strip, current, int(index), coherent=False,
                                  intervals=32, output=output/f'step{step}')
        current = states['opened_minimum']
        row = dict(step=step, bond_index=int(index),
            saddle_energy_from_initial_eV=energy_difference(model, states['saddle'], initial),
            endpoint_energy_from_initial_eV=energy_difference(model, current, initial),
            forward_barrier_eV=summary['forward_barrier_eV'],
            reverse_barrier_eV=summary['reverse_barrier_eV'])
        rows.append(row)
        save_json(output/'running_summary.json', rows)
        print(json.dumps(row), flush=True)
    coherent = np.load(root/'front4/coherent/stationary_positions.npz')['opened_minimum']
    max_error = float(np.max(abs(current-coherent)))
    if max_error > 2e-5:
        raise RuntimeError('sequential advance did not connect the coherent control endpoint')
    summary = {'steps':rows, 'complete_front_endpoint_max_difference_A':max_error,
        'maximum_static_path_energy_above_initial_eV':max(row['saddle_energy_from_initial_eV'] for row in rows),
        'scope':'one verified sequential path in the finite cell; no global minimum-path proof',
        'physical_time_seconds':None, 'elapsed_seconds':time.perf_counter()-started, 'completed':True}
    save_json(output/'summary.json', summary)


if __name__ == '__main__':
    main()
