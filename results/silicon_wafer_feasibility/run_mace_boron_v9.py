"""Independent neutral MACE-MP-0b3 research on source B/Si configurations.

This is new ML-potential minimization, NOT new DFT or a charged-dopant model.
The published comparison uses same-composition neutral 3B configurations only.
"""
from __future__ import annotations
import argparse
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import time
import numpy as np

MODEL_SHA256 = '2f2be696351ac9e94fbe01cdfb6f017679acdbd2db7645209ef55fec9826b012'
MODEL_URL = 'https://github.com/ACEsuit/mace-foundations/releases/download/mace_mp_0b3/mace-mp-0b3-medium.model'


def write_json(path, obj):
    path.write_text(json.dumps(obj, indent=2)+'\n', encoding='utf-8')


def main(args):
    import torch
    from ase import Atoms
    from ase.build import bulk
    from ase.io import write
    from ase.optimize import LBFGS, FIRE
    from mace.calculators import MACECalculator
    from scipy.optimize import minimize_scalar
    if hashlib.sha256(args.model.read_bytes()).hexdigest() != MODEL_SHA256:
        raise ValueError('unexpected model file')
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError('use a fresh output folder; do not overwrite research evidence')
    args.output.mkdir(parents=True, exist_ok=True)
    torch.set_num_threads(2)
    torch.set_num_interop_threads(1)
    calc = MACECalculator(model_paths=str(args.model), device='cpu', default_dtype='float64')
    count = [0]
    original = calc.calculate
    def counted(*a, **kw):
        count[0] += 1
        return original(*a, **kw)
    calc.calculate = counted
    started = time.perf_counter()
    source = json.loads(args.structures.read_text(encoding='utf-8'))
    def source_atoms(key):
        s = source[key]
        return Atoms(symbols=s['symbols'], scaled_positions=s['fractional_positions'],
                     cell=s['cell_A'], pbc=True, calculator=calc)
    def energy(atoms):
        atoms.calc = calc
        return float(atoms.get_potential_energy())
    def state(atoms):
        atoms.calc = calc
        e = energy(atoms); forces = atoms.get_forces(); stress = atoms.get_stress()
        return dict(energy_eV=e, max_force_eV_A=float(np.linalg.norm(forces, axis=1).max()),
                    total_force_eV_A=forces.sum(axis=0).tolist(), stress_eV_A3=stress.tolist(),
                    atoms=len(atoms), volume_A3=atoms.get_volume(), evaluations=count[0])

    # Smooth non-equilibrium source geometry probes both Si and B derivatives.
    atoms = source_atoms('3boron_int_cluster_confi_1')
    initial = state(atoms)
    forces = atoms.get_forces().copy(); stress = atoms.get_stress().copy()
    checks = []
    for index, axis in [(0, 0), (17, 1), (64, 2), (66, 0)]:
        for step in [1e-3, 3e-4, 1e-4]:
            plus, minus = atoms.copy(), atoms.copy()
            plus.positions[index, axis] += step; minus.positions[index, axis] -= step
            fd = -(energy(plus)-energy(minus))/(2*step)
            checks.append(dict(kind='force', atom=index, axis=axis, step=step,
                               autodiff=float(forces[index, axis]), finite_difference=fd,
                               absolute_error=abs(fd-forces[index, axis])))
    for axis in range(3):
        for step in [1e-3, 3e-4, 1e-4]:
            plus, minus = atoms.copy(), atoms.copy()
            deformation = np.eye(3); deformation[axis, axis] += step
            plus.set_cell(atoms.cell.array@deformation.T, scale_atoms=True)
            deformation[axis, axis] -= 2*step
            minus.set_cell(atoms.cell.array@deformation.T, scale_atoms=True)
            fd = (energy(plus)-energy(minus))/(2*step*atoms.get_volume())
            checks.append(dict(kind='stress', axis=axis, step=step,
                               autodiff=float(stress[axis]), finite_difference=fd,
                               absolute_error=abs(fd-stress[axis])))
    translated = atoms.copy(); translated.positions += [.123, -.231, .057]
    permuted = atoms[np.arange(len(atoms))[::-1]]
    checks.append(dict(kind='translation', energy_difference_eV=energy(translated)-initial['energy_eV']))
    checks.append(dict(kind='permutation', energy_difference_eV=energy(permuted)-initial['energy_eV']))
    write_json(args.output/'derivative_checks.json', checks)
    # Pure Si reference volume is determined independently; source relaxations
    # below keep the archive's 10.862 A cube to match its stated protocol.
    bulk_curve = []
    def bulk_energy(lattice):
        b = bulk('Si', 'diamond', a=lattice, cubic=True)
        e = energy(b)/len(b)
        bulk_curve.append(dict(a_A=float(lattice), energy_per_atom_eV=e))
        return e
    optimum = minimize_scalar(bulk_energy, bounds=(5.2, 5.7), method='bounded', options={'xatol':1e-7})
    bulk_eq = bulk('Si', 'diamond', a=optimum.x, cubic=True); bulk_eq.calc = calc
    reference = dict(a_A=float(optimum.x), optimizer_success=bool(optimum.success),
                     state=state(bulk_eq), curve=bulk_curve)
    write_json(args.output/'bulk_reference.json', reference)
    print('BULK', json.dumps(reference['state']), 'a', optimum.x, flush=True)

    results = []
    for key in sorted(source):
        case_start = time.perf_counter(); before = count[0]
        atoms = source_atoms(key)
        record = dict(structure=key, model='MACE-MP-0b3-medium', imposed_charge='none; neutral ML model',
                      source_geometry='Durham B-2025 input', cell_constraint='fixed source cube',
                      initial=state(atoms), attempts=[])
        write(args.output/(key+'_initial.extxyz'), atoms)
        dynamics = LBFGS(atoms, logfile=str(args.output/(key+'_lbfgs.log')), maxstep=.10)
        ok = dynamics.run(fmax=args.fmax, steps=args.steps)
        record['attempts'].append(dict(optimizer='LBFGS', steps=dynamics.nsteps, converged=bool(ok), state=state(atoms)))
        write(args.output/(key+'_lbfgs.extxyz'), atoms)
        if not ok:
            dynamics = FIRE(atoms, logfile=str(args.output/(key+'_fire.log')), maxstep=.08, dt=.03)
            ok = dynamics.run(fmax=args.fmax, steps=args.steps)
            record['attempts'].append(dict(optimizer='FIRE', steps=dynamics.nsteps, converged=bool(ok), state=state(atoms)))
            write(args.output/(key+'_fire.extxyz'), atoms)
        final = state(atoms)
        record.update(final=final, converged=bool(ok and final['max_force_eV_A'] <= args.fmax),
                      evaluations=count[0]-before, elapsed_reported_s=time.perf_counter()-case_start)
        write(args.output/(key+'_final.extxyz'), atoms)
        write_json(args.output/(key+'_result.json'), record)
        results.append(record)
        print('CASE', key, json.dumps(final), 'converged', record['converged'], flush=True)
    tri = [next(r for r in results if r['structure'] == f'3boron_int_cluster_confi_{i}') for i in [1,2,3]]
    reference_deltas = [0., 1.42, 4.51]
    comparisons = []
    for r, published in zip(tri, reference_deltas):
        delta = r['final']['energy_eV']-tri[0]['final']['energy_eV']
        comparisons.append(dict(structure=r['structure'], MACE_relative_energy_eV=delta,
            published_LDA_relative_formation_energy_eV=published, difference_eV=delta-published,
            all_three_force_converged=all(s['converged'] for s in tri),
            status='cross-functional transfer check; not common-level DFT validation'))
    report = dict(model_name='MACE-MP-0b3-medium', model_sha256=MODEL_SHA256, model_URL=MODEL_URL,
        license='MIT', training='MPTrj PBE+U; transfer to B defect LDA targets is not guaranteed',
        dtype='float64', device='CPU', torch_threads=2,
        versions={p:version(p) for p in ['mace-torch','torch','e3nn','ase','numpy','scipy','matscipy']},
        force_tolerance_eV_A=args.fmax, max_steps_per_optimizer=args.steps,
        cases=results, comparisons=comparisons, total_calculator_evaluations=count[0],
        elapsed_reported_s=time.perf_counter()-started,
        status='research reference only; no charge control, doping fracture or clock calibration',
        new_DFT_runs=0, new_MD_runs=0)
    write_json(args.output/'run_summary.json', report)
    print('DONE', json.dumps(comparisons), 'evaluations', count[0], flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--model', type=Path, required=True)
    p.add_argument('--structures', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--fmax', type=float, default=1e-4)
    p.add_argument('--steps', type=int, default=600)
    main(p.parse_args())
