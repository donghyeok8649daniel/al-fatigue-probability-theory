"""Fixed-cell B/P concentration sweep in a neutral ML rigid-cleavage model.

New potential evaluations, no DFT, MD, fitted dopant model, or Si probability.
All relaxed bulk structures and rigid opening configurations are preserved.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from solver_v1.silicon_concentration_research import (
    chemical_concentration_cm3, constrained_energy_derivative,
    deterministic_substitution_order, rigid_opening_geometry)
from results.silicon_wafer_feasibility.mace_force_only_v9 import force_calculator
from results.silicon_wafer_feasibility.run_mace_boron_v9 import MODEL_SHA256


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def csv_write(path, rows):
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def main(args):
    import torch
    from ase import Atoms
    from ase.io import write
    from ase.lattice.cubic import Diamond
    from ase.optimize import LBFGS, FIRE
    from mace.calculators import MACECalculator

    if hashlib.sha256(args.model.read_bytes()).hexdigest() != MODEL_SHA256:
        raise ValueError('model identity mismatch')
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError('fresh output required')
    args.output.mkdir(parents=True, exist_ok=True)
    torch.set_num_threads(2); torch.set_num_interop_threads(1)
    class CountedMACE(MACECalculator):
        def __init__(self, **kwargs):
            self.calls = 0
            super().__init__(**kwargs)
        def calculate(self, *a, **kw):
            self.calls += 1
            return super().calculate(*a, **kw)
    model = CountedMACE(model_paths=str(args.model), device='cpu', default_dtype='float64')
    calc = force_calculator(model)
    lattice = 5.470992313871248
    base = Diamond('Si', directions=[[1,-1,0],[1,1,-2],[1,1,1]],
                   size=(3,2,3), latticeconstant=lattice, pbc=True)
    # Put the periodic cut in the WIDE (111) shuffle interplane gap.
    z = np.unique(np.round(base.positions[:, 2], 10))
    gaps = np.diff(np.r_[z, z[0]+base.cell[2,2]])
    index = int(np.argmax(gaps))
    cut = z[index]+gaps[index]/2
    base.positions[:, 2] = (base.positions[:, 2]-cut) % base.cell[2,2]
    # ASE's default wrap eps may leave slightly negative coordinates.
    base.positions[:] %= np.diag(base.cell)
    volume, area = base.get_volume(), base.cell[0,0]*base.cell[1,1]
    openings = np.array([0., .03, .1, .2, .35, .5, .7, .9, 1.1, 1.4, 1.8, 2.3, 3., 4., 6., 8., 10.])
    plans = [('Si', 0, 'pure', [])]
    for species in args.species:
        for arrangement in args.arrangements:
            order = deterministic_substitution_order(base.positions, base.cell.array, arrangement=arrangement)
            for count in args.counts:
                if count > len(order): raise ValueError('not enough selected substitution sites')
                plans.append((species, count, arrangement, order[:count].tolist()))
    manifest = dict(model_sha256=MODEL_SHA256, lattice_A=lattice, atoms=len(base),
        cell_A=base.cell.array.tolist(), reference_material_volume_A3=volume, area_A2=area,
        original_wide_gap_A=float(gaps[index]), opening_grid_A=openings.tolist(),
        plans=[dict(species=s, count=n, arrangement=a, substitution_indices=i) for s,n,a,i in plans],
        scope='neutral ML; fixed material volume; static internally relaxed bulk then rigid (111) separation',
        carrier_density_cm3=None, temperature_K=None,
        temperature_note='0 K potential-energy calculation, no thermal sampling or PMF',
        physical_mobility=None, specimen_failure_probability=None, production_enabled=False)
    dump(args.output/'protocol.json', manifest)
    write(args.output/'perfect_reference.extxyz', base)
    rows, results = [], []
    for species, count, arrangement, indices in plans:
        started = time.perf_counter(); standard_start = model.calls; force_start = calc.evaluations
        tag = f'{species}_{count:02d}_{arrangement}'
        folder = args.output/tag; folder.mkdir()
        atoms = base.copy()
        for i in indices: atoms[i].symbol = species
        atoms.calc = calc
        write(folder/'initial.extxyz', atoms)
        attempts = []
        for method, steps in [('LBFGS', 220), ('FIRE', 300)]:
            if method == 'LBFGS': opt = LBFGS(atoms, logfile=str(folder/'lbfgs.log'), maxstep=.08)
            else: opt = FIRE(atoms, logfile=str(folder/'fire.log'), maxstep=.05, dt=.025)
            converged = opt.run(fmax=2e-4, steps=steps)
            attempts.append(dict(method=method, steps=opt.nsteps, converged=bool(converged),
                force_max_eV_A=float(np.linalg.norm(atoms.get_forces(), axis=1).max())))
            if converged: break
        atoms.positions[:] %= np.diag(atoms.cell)
        atoms.calc = model
        energy0 = float(atoms.get_potential_energy())
        force0 = atoms.get_forces().copy(); stress0 = atoms.get_stress(voigt=False).copy()
        write(folder/'bulk.extxyz', atoms)
        bulk_r, bulk_c, numbers = atoms.positions.copy(), atoms.cell.array.copy(), atoms.numbers.copy()
        bulk = dict(species=species, dopant_count=count, arrangement=arrangement,
            indices=indices, chemical_concentration_cm3=chemical_concentration_cm3(count, reference_material_volume_A3=volume),
            volume_A3=volume, energy_eV=energy0, stress_GPa=(stress0*160.2176634).tolist(),
            force_max_eV_A=float(np.linalg.norm(force0, axis=1).max()), converged=bool(converged),
            attempts=attempts, full_Hessian_checked=False)
        dump(folder/'bulk.json', bulk)
        if not converged: raise RuntimeError('bulk force convergence failed: '+tag)
        case_rows = []
        for q in openings:
            r, c, dr, dc = rigid_opening_geometry(bulk_r, bulk_c, q)
            opened = Atoms(numbers=numbers, positions=r, cell=c, pbc=True, calculator=model)
            energy = float(opened.get_potential_energy())
            force = opened.get_forces().copy(); stress = opened.get_stress(voigt=False).copy()
            derivative = constrained_energy_derivative(r, force, c, stress, position_tangent=dr, cell_tangent=dc)
            row = dict(case=tag, species=species, count=count, arrangement=arrangement,
                chemical_concentration_cm3=bulk['chemical_concentration_cm3'], opening_A=float(q),
                energy_eV=energy, delta_energy_eV=energy-energy0, area_A2=area,
                separation_work_J_m2=(energy-energy0)/area*16.02176634,
                conjugate_force_eV_A=derivative, rigid_traction_GPa=derivative/area*160.2176634,
                affine_stress_zz_GPa=float(stress[2,2]*160.2176634),
                nonaffine_correction_GPa=float((derivative/area-stress[2,2])*160.2176634),
                total_force_norm_eV_A=float(np.linalg.norm(force.sum(axis=0))))
            case_rows.append(row); rows.append(row)
            np.savez_compressed(folder/f'opening_{q:.3f}.npz', positions=r,cell=c,numbers=numbers,
                energy=energy,forces=force,stress=stress,position_tangent=dr,cell_tangent=dc)
        csv_write(folder/'curve.csv', case_rows)
        top = max(case_rows, key=lambda r:r['rigid_traction_GPa'])
        result = dict(**bulk, case=tag, sampled_peak_traction_GPa=top['rigid_traction_GPa'],
            sampled_peak_opening_A=top['opening_A'], rigid_separation_work_J_m2=case_rows[-1]['separation_work_J_m2'],
            plateau_8_to_10_eV=case_rows[-1]['energy_eV']-case_rows[-2]['energy_eV'],
            standard_calculator_calls=model.calls-standard_start,
            force_only_calculator_calls=calc.evaluations-force_start,
            elapsed_s=time.perf_counter()-started, status='rigid constrained path; not relaxed cleavage or strength')
        dump(folder/'result.json', result); results.append(result)
        csv_write(args.output/'all_curves.csv', rows)
        dump(args.output/'progress.json', dict(completed=len(results),planned=len(plans),cases=results))
        print('CASE', tag, 'C',bulk['chemical_concentration_cm3'],'W',result['rigid_separation_work_J_m2'],
              'peak',result['sampled_peak_traction_GPa'],'seconds',result['elapsed_s'],flush=True)
    dump(args.output/'summary.json',dict(cases=results,complete=True,model_sha256=MODEL_SHA256,
        standard_calculator_calls=model.calls,force_only_calculator_calls=calc.evaluations,
        new_DFT_runs=0,new_MD_runs=0,production_enabled=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--model', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--species', nargs='+', default=['B','P'])
    parser.add_argument('--counts', type=int, nargs='+', default=[1,2,4,8])
    parser.add_argument('--arrangements', nargs='+', default=['interface_cluster','interface_spread','bulk_spread'])
    main(parser.parse_args())
