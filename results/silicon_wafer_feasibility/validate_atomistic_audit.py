"""Recompute saved Si/DFT metrics and independently check actual input cells."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import time

import numpy as np

from solver_v1.silicon_environment_research import environment_jet
from .run_static_probe import source_parameters
from .run_atomistic_controls import save_json, checksum, sources
from .run_dft_material_audit import read_frames, reference_fields


def independent_periodic_sw(atoms, parameters):
    """ASE image enumeration + explicit per-center triple jets, including self images."""
    from ase.neighborlist import neighbor_list
    i, j, vectors = neighbor_list('ijD', atoms, parameters['sigma']*parameters['a'])
    gradient = np.zeros_like(atoms.positions)
    sites = []
    for center in range(len(atoms)):
        selected = i == center
        jet = environment_jet(vectors[selected], parameters, representation='direct')
        edge_gradient = jet.gradient.reshape(-1, 3)
        np.add.at(gradient, j[selected], edge_gradient)
        gradient[center] -= edge_gradient.sum(axis=0)
        sites.append(jet.value)
    return float(np.sum(sites)), -gradient


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, default=Path('.cache/si-atomistic-v5/Si_PRX_GAP.zip'))
    parser.add_argument('--calculation', type=Path, default=Path('results/silicon_atomistic_v5'))
    args = parser.parse_args()
    started = time.perf_counter()
    frames, _ = read_frames(args.archive)
    records = list(csv.DictReader((args.calculation/'dft_frame_metrics.csv').open(encoding='utf-8')))
    summary = json.loads((args.calculation/'dft_material_audit.json').read_text())
    specs, params = sources(args.calculation)
    comparisons, maximum_metric_error = [], 0.
    # Materialize each compressed array once. Accessing NpzFile[key] in the
    # per-frame loop repeatedly decompresses the entire 171815-atom field.
    with np.load(args.calculation/'dft_predictions.npz') as archive:
        data = {key:archive[key] for key in archive.files}
        offsets = data['offsets']
        assert offsets.tolist() == np.r_[0, np.cumsum([len(a) for a in frames])].tolist()
        for frame, atoms in enumerate(frames):
            ed, fd, _ = reference_fields(atoms)
            np.testing.assert_array_equal(data['reference_force'][offsets[frame]:offsets[frame+1]], fd)
            assert data['reference_energy'][frame] == ed
        for row in records:
            frame, label = int(row['frame']), row['model']
            sl = slice(offsets[frame], offsets[frame+1])
            error = data[label+'_force'][sl]-data['reference_force'][sl]
            expected = np.sqrt(np.sum(error*error)/error.size)
            maximum_metric_error = max(maximum_metric_error, abs(expected-float(row['component_rmse_eV_A'])))
            assert abs(expected-float(row['component_rmse_eV_A'])) < 1e-12
            assert data[label+'_energy'][frame] == float(row['model_energy_eV'])
        for group in summary['group_results']:
            subset = [int(r['frame']) for r in records if r['model'] == group['model']
                      and r['config_type'] == group['config_type'] and r['declared_xc'] == group['declared_xc']]
            errors = np.vstack([data[group['model']+'_force'][offsets[i]:offsets[i+1]]
                                - data['reference_force'][offsets[i]:offsets[i+1]] for i in subset])
            assert abs(np.sqrt(np.mean(errors**2))-group['component_rmse_eV_A']) < 1e-12
        # Actual source structures; each independent calculator sees the ORIGINAL
        # cell, including a negative determinant. No LAMMPS transform is reused.
        selected = {1, int(summary['energy_anchor']['frame']),
                    int(summary['cell_basis_changes'][0]['frame'])}
        for kind in ['surface_111', 'surface_111_pandey', 'crack_111_1-10',
                     'crack_110_1-10', 'screw_disloc']:
            selected.add(next(i for i, a in enumerate(frames) if a.info['config_type'] == kind))
        from ase.calculators.tersoff import Tersoff
        for i in sorted(selected):
            atoms = frames[i].copy()
            sl = slice(offsets[i], offsets[i+1])
            for label in specs:
                if label == 'original_sw_1985':
                    energy, force = independent_periodic_sw(atoms, params)
                else:
                    atoms.calc = Tersoff.from_lammps(specs[label][1])
                    energy, force = atoms.get_potential_energy(), atoms.get_forces()
                row = dict(frame=i, config_type=atoms.info['config_type'], model=label,
                           original_cell_determinant_A3=float(np.linalg.det(atoms.cell.array)),
                           energy_error_eV=abs(energy-data[label+'_energy'][i]),
                           force_max_error_eV_A=float(np.max(abs(force-data[label+'_force'][sl]))))
                comparisons.append(row)
                assert row['energy_error_eV'] < 1e-8 and row['force_max_error_eV_A'] < 1e-9, row
                print(json.dumps(row), flush=True)
    polish = json.loads((args.calculation/'reconstruction_polish.json').read_text())
    assert all(r['final_force_max_eV_A'] < 1e-8 and r['lowest_eigenvalues_eV_A2'][0] > 0
               for r in polish['rows'])
    report = dict(passed=True, checked_reference_frames=len(frames), checked_model_frame_pairs=len(records),
                  checked_groups=len(summary['group_results']), max_rmse_recomputation_error=maximum_metric_error,
                  independent_evaluations=comparisons, corrected_conditional_minima=len(polish['rows']),
                  code_sha256=checksum(__file__), elapsed_seconds=time.perf_counter()-started,
                  scope='numerical/source/metric validation only; not material approval')
    save_json(args.calculation/'raw_validation.json', report)


if __name__ == '__main__':
    main()
