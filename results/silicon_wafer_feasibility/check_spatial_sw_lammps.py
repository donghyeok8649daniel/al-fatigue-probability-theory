"""Independent optional LAMMPS run-0 energy/force comparison; no MD is run.

Requires a separately installed LAMMPS Python library with the MANYBODY package.
The production UI/solver does not import it. Uses the hash-bound original Si.sw.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from solver_v1.silicon_crack_research import SpatialSW, diamond_crack_strip
from .run_static_probe import ROOT, source_parameters


def lammps_evaluation(positions, period, potential):
    from lammps import lammps
    lmp = lammps(cmdargs=['-log', 'none', '-screen', 'none'])
    commands = []

    def command(line):
        commands.append(line)
        lmp.command(line)

    try:
        command('units metal')
        command('atom_style atomic')
        command('atom_modify map array')
        command('boundary f f p')
        lo, hi = positions.min(axis=0)-12., positions.max(axis=0)+12.
        command(f'region box block {lo[0]:.17g} {hi[0]:.17g} {lo[1]:.17g} {hi[1]:.17g} 0 {period:.17g} units box')
        command('create_box 1 box')
        current = positions.copy()
        current[:, 2] %= period
        n = len(positions)
        made = lmp.create_atoms(n, list(range(1, n+1)), [1]*n, current.ravel().tolist())
        if made != n:
            raise RuntimeError('LAMMPS did not create every supplied atom')
        command('mass 1 28.0855')  # Static energy/force only; not an overdamped clock.
        command('pair_style sw')
        command(f'pair_coeff * * "{potential.resolve().as_posix()}" Si')
        command('neighbor 0.6 bin')
        command('neigh_modify every 1 delay 0 check yes')
        command('thermo_modify norm no')
        command('run 0 post no')
        energy = float(lmp.get_thermo('pe'))
        forces = np.ctypeslib.as_array(lmp.gather_atoms('f', 1, 3), shape=(3*n,)).reshape(n, 3).copy()
        library = Path(lmp.lib._name)
        metadata = {'lammps_version':int(lmp.version()), 'library_filename':library.name,
                    'library_sha256':hashlib.sha256(library.read_bytes()).hexdigest()}
        metadata['commands'] = [line.replace(potential.resolve().as_posix(), 'source_Si.sw') for line in commands]
        return energy, forces, metadata
    finally:
        lmp.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--calculation', type=Path, default=Path('results/silicon_local_crack_v2'))
    parser.add_argument('--output', type=Path, default=Path('results/silicon_local_crack_v2/lammps_validation.json'))
    args = parser.parse_args()
    started = time.perf_counter()
    params, checksum = source_parameters()
    summary = json.loads((args.calculation/'running_summary.json').read_text(encoding='utf-8'))
    strip = diamond_crack_strip(summary['lattice_A'], nx=6, ny=4, nz=4)
    model = SpatialSW(params, front_period=strip.front_period)
    rng = np.random.default_rng(19284)
    states = {'perfect_strip':strip.positions, 'opened_seed':strip.seed(3.5, tip=-3., transition=5.),
              'nonuniform_perturbation':strip.positions+rng.normal(0, .09, strip.positions.shape)}
    stationary = args.calculation/'front4/localized/stationary_positions.npz'
    if stationary.exists():
        with np.load(stationary) as data:
            states['localized_saddle'] = data['saddle']
            states['localized_opened_minimum'] = data['opened_minimum']
    rows = []
    for label, positions in states.items():
        actual = model.evaluate(positions)
        independent = model.evaluate(positions, representation='direct')
        e, f, metadata = lammps_evaluation(positions, strip.front_period, ROOT/'source_Si.sw')
        row = {'state':label, 'atoms':len(positions), 'moment_energy_eV':actual.energy,
            'lammps_energy_eV':e, 'energy_abs_difference_eV':abs(e-actual.energy),
            'force_max_abs_difference_eV_A':float(np.max(abs(f+actual.gradient))),
            'direct_moment_energy_difference_eV':abs(independent.energy-actual.energy),
            'direct_moment_force_max_difference_eV_A':float(np.max(abs(independent.gradient-actual.gradient)))}
        if row['energy_abs_difference_eV'] > 1e-8 or row['force_max_abs_difference_eV_A'] > 1e-10:
            raise RuntimeError(f'independent LAMMPS disagreement: {row}')
        rows.append(row)
        print(json.dumps(row), flush=True)
    result = {'scope':'independent run-0 energy and forces, no MD or kinetic inference',
        'source_sha256':checksum, 'library':metadata, 'cases':rows,
        'elapsed_seconds':time.perf_counter()-started, 'passed':True}
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')


if __name__ == '__main__':
    main()
