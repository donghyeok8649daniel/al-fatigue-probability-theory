"""Actual neutral MACE relaxation of a finite, initially uncracked Si prism.

Numerical force tolerance, trial step and specimen size are explicit. A load
drop is not automatically classified as crack initiation. All states survive.
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
from solver_v1.silicon_initiation_research import (
    intact_prism, prescribed_grips, grip_observables, connectivity_diagnostics)
from results.silicon_wafer_feasibility.mace_force_only_v9 import force_calculator, force_only
from results.silicon_wafer_feasibility.run_mace_boron_v9 import MODEL_SHA256


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def main(args):
    import torch
    from ase.constraints import FixAtoms
    from ase.io import write, read
    from ase.optimize import LBFGS, FIRE, LBFGSLineSearch
    from mace.calculators import MACECalculator

    started = time.perf_counter()
    if hashlib.sha256(args.model.read_bytes()).hexdigest() != MODEL_SHA256:
        raise ValueError('unexpected MACE weights')
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError('fresh output required; preserve failed trials')
    if args.strains[0] != 0 or np.any(np.diff(args.strains) <= 0):
        raise ValueError('increasing strain schedule starting at zero required')
    args.output.mkdir(parents=True, exist_ok=True)
    torch.set_num_threads(2); torch.set_num_interop_threads(1)
    model = MACECalculator(model_paths=str(args.model), device='cpu', default_dtype='float64')
    calc = force_calculator(model)
    atoms, geometry = intact_prism(repeats=args.repeats, grip_width_A=args.grip_width)
    lower, upper, free = (geometry[x] for x in ('lower', 'upper', 'free'))
    fixed = lower | upper
    reference = atoms.positions.copy()
    write(args.output/'uncut_lattice_prism.extxyz', atoms)
    rng = np.random.default_rng(args.seed)
    atoms.positions[free] += rng.normal(0., args.perturbation, (np.count_nonzero(free), 3))
    initial_record=None
    if args.initial is not None:
        if args.initial.suffix=='.npz':
            with np.load(args.initial) as data:
                prior_numbers=data['numbers'].copy();prior_positions=data['positions'].copy()
        else:
            prior=read(args.initial);prior_numbers=prior.numbers;prior_positions=prior.positions
        # ASE extxyz writes eight decimal places by default. Check that scale,
        # then restore prescribed grips exactly; never propagate text rounding.
        if (not np.array_equal(prior_numbers,atoms.numbers)
                or not np.allclose(prior_positions[fixed],reference[fixed],atol=1e-7,rtol=0)):
            raise ValueError('restart is not the same intact-prism geometry and zero-load grips')
        rounding=float(np.max(abs(prior_positions[fixed]-reference[fixed])))
        atoms.positions[:]=prior_positions;atoms.positions[fixed]=reference[fixed]
        initial_record=dict(file_name=args.initial.name,sha256=hashlib.sha256(args.initial.read_bytes()).hexdigest(),
            fixed_coordinate_roundtrip_error_A=rounding,fixed_coordinates_restored_exactly=True,
            meaning='continuation of the preserved zero-load relaxation, not a new initial crack')
    atoms.set_constraint(FixAtoms(mask=fixed)); atoms.calc = calc
    np.savez_compressed(args.output/'geometry.npz', reference=reference, **geometry)
    initial_graph = connectivity_diagnostics(reference, lower=lower, upper=upper)
    if not all(row['components'] == 1 and row['grip_connected'] for row in initial_graph):
        raise RuntimeError('initial reference is not intact and connected')
    protocol = dict(repeats=args.repeats, atoms=len(atoms), free_atoms=int(free.sum()),
        seed=args.seed, perturbation_A=args.perturbation, grip_width_A=args.grip_width,
        restart=initial_record,initial_optimizer=args.optimizer,
        gauge_length_A=geometry['gauge_length_A'], area_A2=geometry['area_A2'],
        nominal_lengths_A=geometry['nominal_lengths_A'].tolist(),
        nominal_area_convention='uncut periodic-cell rectangle before external surface creation',
        initial_graph=initial_graph, strains=args.strains, fmax_eV_A=args.fmax,
        model_sha256=MODEL_SHA256, temperature_K=None,
        scope='0 K bare finite prism; rigid grips; no initial internal cut, void or crack seed',
        initial_surface='unpassivated ideal (001) and (1-10) facets, followed by local relaxation',
        oxide_included=False, measured_wafer_strength=False,
        initiation_probability=None, first_crack_event_certified=False, physical_clock=None)
    dump(args.output/'protocol.json', protocol)
    if args.validate_initial:
        probe=atoms.copy();probe.set_constraint();probe.calc=model
        standard_energy=float(probe.get_potential_energy());standard_force=probe.get_forces()
        energy,force=force_only(model,probe)
        checks=dict(energy_difference_eV=energy-standard_energy,
                    force_max_difference_eV_A=float(np.max(abs(force-standard_force))),differences=[])
        for name,direction in [('free_random',rng.normal(size=reference.shape)),('grip_extension',np.zeros_like(reference))]:
            if name=='free_random':
                direction[fixed]=0;direction/=np.linalg.norm(direction)
            else:direction[lower,2]=-.5;direction[upper,2]=.5
            exact=-float(np.sum(force*direction))
            for h in (1e-3,5e-4):
                plus=probe.copy();minus=probe.copy()
                plus.positions[:]+=h*direction;minus.positions[:]-=h*direction
                ep,_=force_only(model,plus);em,_=force_only(model,minus)
                got=(ep-em)/(2*h)
                checks['differences'].append(dict(direction=name,step_A=h,analytic=exact,finite_difference=got,
                    absolute_error_eV_A=abs(got-exact)))
        dump(args.output/'initial_force_validation.json',checks)
        if (abs(checks['energy_difference_eV'])>1e-8 or checks['force_max_difference_eV_A']>1e-8
                or max(x['absolute_error_eV_A'] for x in checks['differences'])>5e-5):
            raise RuntimeError('initial geometry force/energy validation failed')
    rows, previous_strain, previous_positions = [], 0., atoms.positions.copy()
    write(args.output/'initial.extxyz', atoms)
    for index, strain in enumerate(args.strains):
        folder = args.output/f'state_{index:03d}'; folder.mkdir()
        state_started, calls_before = time.perf_counter(), calc.evaluations
        extension = strain*geometry['gauge_length_A']
        target = prescribed_grips(reference, lower=lower, upper=upper, extension_A=extension)
        current = previous_positions.copy()
        # Affine predictor only for free atoms; the actual endpoints are rigid.
        zlo, zhi = reference[lower, 2].mean(), reference[upper, 2].mean()
        fraction = np.clip((reference[:, 2]-zlo)/(zhi-zlo), 0, 1)-.5
        current[:, 2] += fraction*(strain-previous_strain)*geometry['gauge_length_A']
        current[fixed] = target[fixed]
        atoms.set_positions(current, apply_constraint=False)
        write(folder/'before_relax.extxyz', atoms)
        attempts = []
        interrupted = None
        try:
            for method in (args.optimizer, 'FIRE'):
                if method=='LBFGS':opt=LBFGS(atoms,logfile=str(folder/'lbfgs.log'),maxstep=.06)
                elif method=='LBFGSLineSearch':opt=LBFGSLineSearch(atoms,logfile=str(folder/'lbfgs_linesearch.log'),maxstep=.08)
                else:opt=FIRE(atoms,logfile=str(folder/'fire.log'),maxstep=.04,dt=.02)
                def checkpoint():
                    if opt.nsteps % 20 == 0:
                        write(folder/'checkpoint.extxyz', atoms)
                        np.savez_compressed(folder/'checkpoint.npz',positions=atoms.positions,
                                            numbers=atoms.numbers,cell=atoms.cell.array)
                        dump(args.output/'running.json', dict(state=index, strain=strain,
                            method=method, step=opt.nsteps, calls=calc.evaluations,
                            elapsed_s=time.perf_counter()-started))
                    if time.perf_counter()-started >= args.max_seconds:
                        raise TimeoutError('explicit run wall-time budget reached')
                opt.attach(checkpoint, interval=1)
                converged = bool(opt.run(fmax=args.fmax, steps=args.max_steps))
                attempts.append(dict(method=method, steps=opt.nsteps, converged=converged))
                if converged: break
        except TimeoutError as exc:
            interrupted = str(exc); converged = False
        energy = float(atoms.get_potential_energy())
        force = atoms.get_forces(apply_constraint=False)
        observables = grip_observables(atoms.positions, force, lower=lower, upper=upper,
                                      area_A2=geometry['area_A2'])
        if not np.allclose(atoms.positions[fixed], target[fixed], rtol=0, atol=1e-12):
            raise RuntimeError('grip constraints drifted')
        graph = connectivity_diagnostics(atoms.positions, lower=lower, upper=upper)
        row = dict(state=index, strain=float(strain), extension_A=float(extension), energy_eV=energy,
            converged=converged, force_calls=calc.evaluations-calls_before,
            elapsed_s=time.perf_counter()-state_started, **observables)
        np.savez_compressed(folder/'raw.npz', positions=atoms.positions, forces=force,
                            energy=energy, numbers=atoms.numbers, cell=atoms.cell.array)
        write(folder/'relaxed.extxyz', atoms)
        dump(folder/'result.json', dict(**row, attempts=attempts, connectivity=graph,
            interruption=interrupted, crack_initiation_certified=False))
        rows.append(row)
        with (args.output/'curve.csv').open('w', newline='', encoding='utf-8') as stream:
            scalar = [{k:v for k,v in row.items() if not isinstance(v, list)} for row in rows]
            writer = csv.DictWriter(stream, fieldnames=list(scalar[0])); writer.writeheader(); writer.writerows(scalar)
        dump(args.output/'progress.json', dict(completed_states=len(rows), planned_states=len(args.strains),
            last_converged=converged, elapsed_s=time.perf_counter()-started, force_calls=calc.evaluations,
            interruption=interrupted, protocol=protocol))
        print('STATE', index, 'strain', strain, 'stress_GPa', observables['nominal_stress_GPa'],
              'fmax', observables['free_force_max_eV_A'], 'converged', converged,
              'elapsed_s', time.perf_counter()-started, flush=True)
        if not converged:
            raise RuntimeError('state did not converge; raw state and completed prefix preserved')
        previous_strain, previous_positions = strain, atoms.positions.copy()
    dump(args.output/'summary.json', dict(completed_states=len(rows), complete=True,
        force_calls=calc.evaluations, elapsed_s=time.perf_counter()-started,
        actual_DFT_runs=0, actual_MD_runs=0, first_crack_event_certified=False,
        status='completed force-relaxed load schedule; no automatic initiation classification'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--model', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--repeats', type=int, nargs=3, default=[3,3,10])
    parser.add_argument('--strains', type=float, nargs='+', default=[0,.02,.04,.06,.08,.10,.12,.14,.16,.18,.20])
    parser.add_argument('--seed', type=int, default=1101)
    parser.add_argument('--perturbation', type=float, default=.01)
    parser.add_argument('--grip-width', type=float, default=6.)
    parser.add_argument('--fmax', type=float, default=.001)
    parser.add_argument('--max-steps', type=int, default=400)
    parser.add_argument('--max-seconds', type=float, default=7200)
    parser.add_argument('--initial',type=Path)
    parser.add_argument('--optimizer',choices=['LBFGS','LBFGSLineSearch'],default='LBFGSLineSearch')
    parser.add_argument('--validate-initial',action='store_true')
    main(parser.parse_args())
