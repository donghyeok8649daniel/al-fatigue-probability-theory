"""Run independent numerical controls, bulk elasticity and same-grip Si audit.

Requires optional ASE 3.26.0 and serial LAMMPS. No parameter fit, new DFT,
physical time calibration, or production solver change is performed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time
from types import SimpleNamespace
from urllib.request import urlopen

import numpy as np
from scipy.optimize import minimize_scalar, root

from solver_v1.silicon_atomistic_reference import LammpsSilicon
from solver_v1.silicon_crack_research import SpatialSW, RelaxedCoordinates, relax_atoms
from .run_static_probe import source_parameters, ROOT, EV_A3_TO_GPA


TERSOFF_URL = ('https://raw.githubusercontent.com/lammps/lammps/'
               'stable_22Jul2025_update4/potentials/SiC.tersoff')
TERSOFF_SHA256 = '157a155a3fc263d43b30dbdf4ccc28f0d633ce7112c68de42fc09b5ad070ce64'
V4 = Path('results/silicon_thermal_v4/references/multibasin')


def checksum(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def sources(output):
    p, h = source_parameters()
    target = Path(output)/'sources'/'SiC.tersoff'
    target.parent.mkdir(parents=True, exist_ok=True)
    if not target.exists():
        with urlopen(TERSOFF_URL, timeout=30) as r:
            raw = r.read()
        if hashlib.sha256(raw).hexdigest() != TERSOFF_SHA256:
            raise ValueError('upstream Tersoff source changed')
        target.write_bytes(raw)
    if checksum(target) != TERSOFF_SHA256:
        raise ValueError('cached Tersoff source differs')
    return {'original_sw_1985': ('sw', ROOT/'source_Si.sw', h),
            'tersoff_1989': ('tersoff', target, TERSOFF_SHA256)}, p


def strip_box(positions, period):
    low, high = positions.min(axis=0)-12., positions.max(axis=0)+12.
    low[2], high[2] = 0., period
    return low, np.diag(high-low)


def controls(engine, style, parameters):
    from ase.build import bulk
    from ase.calculators.tersoff import Tersoff
    rng = np.random.default_rng(281441)
    atoms = bulk('Si', a=5.43)*(2, 2, 2)
    r = atoms.positions+rng.normal(0., .035, atoms.positions.shape)
    cell = atoms.cell.array @ np.array([[1.01, .03, -.01], [0, .99, .02], [0, 0, 1.02]])
    first = engine.evaluate(r, cell)
    direction = rng.normal(size=r.shape)
    direction /= np.linalg.norm(direction)
    gradient = float(np.sum(first.gradient*direction))
    force_checks = []
    for step in [1e-4, 5e-5, 2.5e-5]:
        fd = (engine.evaluate(r+step*direction, cell).energy
              - engine.evaluate(r-step*direction, cell).energy)/(2*step)
        force_checks.append(dict(step_A=step, abs_error_eV_A=abs(fd-gradient)))
    rotation, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    if np.linalg.det(rotation) < 0:
        rotation[:, 0] *= -1
    rotated = engine.evaluate(r @ rotation, cell @ rotation)
    translated = engine.evaluate(r+np.array([2.31, -.41, 1.1]), cell)
    repeated = atoms.copy()
    repeated.set_cell(cell)
    repeated.set_positions(r)
    repeated = repeated.repeat((2, 1, 1))
    large = engine.evaluate(repeated.positions, repeated.cell.array)
    deformation = np.array([[.2, .3, -.1], [.3, -.4, .2], [-.1, .2, .1]])
    affine_exact = float(np.linalg.det(cell)*np.sum(first.stress*deformation))
    h = 1e-5
    fp, fm = np.eye(3)+h*deformation, np.eye(3)-h*deformation
    affine_fd = (engine.evaluate(r @ fp.T, cell @ fp.T).energy
                 - engine.evaluate(r @ fm.T, cell @ fm.T).energy)/(2*h)
    row = dict(force_finite_differences=force_checks,
        rotation_energy_error_eV=abs(rotated.energy-first.energy),
        rotation_force_error_eV_A=float(np.max(abs(rotated.gradient-first.gradient @ rotation))),
        rotation_stress_error_eV_A3=float(np.max(abs(rotated.stress-rotation.T @ first.stress @ rotation))),
        translation_energy_error_eV=abs(translated.energy-first.energy),
        translation_force_error_eV_A=float(np.max(abs(translated.gradient-first.gradient))),
        supercell_energy_per_copy_error_eV=abs(large.energy/2-first.energy),
        supercell_force_error_eV_A=float(np.max(abs(large.gradient-np.tile(first.gradient, (2, 1))))),
        affine_stress_derivative_error_eV=abs(affine_fd-affine_exact),
        net_force_eV_A=np.sum(-first.gradient, axis=0).tolist())
    if style == 'tersoff':
        independent = atoms.copy()
        independent.set_cell(cell)
        independent.set_positions(r)
        independent.calc = Tersoff.from_lammps(engine.potential)
        row['independent_ase_energy_error_eV'] = abs(independent.get_potential_energy()-first.energy)
        row['independent_ase_force_error_eV_A'] = float(np.max(abs(independent.get_forces()+first.gradient)))
    else:
        with np.load(V4/'branch_positions.npz') as d:
            r = d['state01'].copy()
        geometry = json.loads((V4/'summary.json').read_text())['geometry']
        period = geometry['front_period_A']
        origin, box = strip_box(r, period)
        direct = SpatialSW(parameters, front_period=period).evaluate(r, representation='direct')
        checked = engine.evaluate(r-origin, box, periodic=(False, False, True))
        row['independent_sw_triples_energy_error_eV'] = abs(direct.energy-checked.energy)
        row['independent_sw_triples_force_error_eV_A'] = float(np.max(abs(direct.gradient-checked.gradient)))
    errors = [value for key, value in row.items() if 'error' in key]
    errors += [item['abs_error_eV_A'] for item in force_checks]
    if max(errors) > 2e-6:
        raise AssertionError(f'atomistic numerical control failed: {row}')
    row['passed'] = True
    return row


def bulk_properties(engine):
    from ase.build import bulk

    def e_at(a):
        atoms = bulk('Si', a=a)
        return engine.evaluate(atoms.positions, atoms.cell.array).energy/len(atoms)

    optimum = minimize_scalar(e_at, bounds=(5.2, 5.7), method='bounded',
                              options={'xatol':1e-10})
    if not optimum.success:
        raise RuntimeError('bulk scalar relaxation failed')
    a = float(optimum.x)
    atoms = bulk('Si', a=a)
    p, c = atoms.positions, atoms.cell.array
    e0, v0 = engine.evaluate(p, c).energy, float(np.linalg.det(c))
    raw = []

    def strain_energy(f, relax):
        current, box = p @ f.T, c @ f.T
        if not relax:
            return engine.evaluate(current, box).energy

        def objective(x):
            trial = current.copy()
            trial[1] += x
            val = engine.evaluate(trial, box)
            return val.energy-e0, val.gradient[1]
        # Resolve the three internal force equations directly: near zero strain
        # extensive-energy line searches hit roundoff before force convergence.
        fit = root(lambda x:objective(x)[1], np.zeros(3), method='hybr',
                   options={'xtol':1e-9, 'eps':1e-10})
        residual = float(np.max(abs(objective(fit.x)[1])))
        if residual > 1e-7:
            raise RuntimeError(f'bulk internal relaxation residual {residual}')
        step = 1e-4
        hessian = np.column_stack([(objective(fit.x+step*v)[1]
                                  -objective(fit.x-step*v)[1])/(2*step)
                                  for v in np.eye(3)])
        if np.linalg.eigvalsh((hessian+hessian.T)/2).min() <= 0:
            raise RuntimeError('bulk internal stationary point is not a local minimum')
        return float(objective(fit.x)[0]+e0)

    for step in [1e-3, 5e-4, 2.5e-4]:
        for relaxed in [False, True]:
            curvatures = {}
            for name, mode in [('axial', np.diag([1., 0., 0.])), ('hydro', np.eye(3)),
                                ('shear', np.array([[0., 1., 0.], [0., 0., 0.], [0., 0., 0.]]))]:
                ep = strain_energy(np.eye(3)+step*mode, relaxed)
                em = strain_energy(np.eye(3)-step*mode, relaxed)
                curvature = (ep+em-2*e0)/(step*step*v0)*EV_A3_TO_GPA
                curvatures[name] = curvature
                raw.append(dict(step=step, relaxed=relaxed, mode=name, plus_eV=ep,
                                minus_eV=em, zero_eV=e0, energy_curvature_GPa=curvature))
            raw.append(dict(step=step, relaxed=relaxed, C11_GPa=curvatures['axial'],
                            C12_GPa=(curvatures['hydro']/3-curvatures['axial'])/2,
                            C44_GPa=curvatures['shear']))
    return dict(lattice_A=a, energy_per_atom_eV=e0/2, volume_per_atom_A3=v0/2,
                stress_eV_A3=engine.evaluate(p, c).stress.tolist(), rows=raw,
                state='0 K perfect pure-Si diamond; clamped and internally relaxed basis')


def reconstruction_audit(engine, output):
    data = np.load(V4/'branch_positions.npz')
    reference, fixed, bond = data['reference'], data['fixed'], data['bond']
    geometry = json.loads((V4/'summary.json').read_text())['geometry']
    origin, cell = strip_box(reference, geometry['front_period_A'])
    model = SimpleNamespace(evaluate=lambda r:engine.evaluate(r-origin, cell,
                                          periodic=(False, False, True)))
    coordinates = RelaxedCoordinates(reference, fixed, bond=bond)
    gap = float(reference[bond[1], 1]-reference[bond[0], 1])
    rows, saved, energies = [], {}, []
    base = model.evaluate(data['state00'])
    for name in ['state00', 'state04', 'state01', 'state02']:
        initial = data[name]
        before = model.evaluate(initial)
        g, reaction = coordinates.pullback(before.gradient)
        final, info = relax_atoms(model, coordinates, initial=initial, values=[gap],
                                 tolerance=2e-6, maxiter=3000)
        after = model.evaluate(final)
        row = dict(name=name, same_geometry_delta_eV=float(np.sum(before.site_energy-base.site_energy)),
                   same_geometry_bath_force_max_eV_A=float(np.max(abs(g))),
                   same_geometry_reaction_eV_A=reaction.tolist(), relaxation=info,
                   relaxed_vs_initial_state00_eV=float(np.sum(after.site_energy-base.site_energy)),
                   rms_displacement_A=float(np.sqrt(np.mean((final-initial)**2))),
                   fixed_atom_max_change_A=float(np.max(abs(final[fixed]-initial[fixed]))),
                   gap_error_A=float(final[bond[1], 1]-final[bond[0], 1]-gap))
        rows.append(row)
        saved[name+'_initial'] = initial
        saved[name+'_relaxed'] = final
        energies.append(after.site_energy.copy())
        print(json.dumps(dict(event='reconstruction', **row)), flush=True)
    for i, row in enumerate(rows):
        row['relaxed_delta_vs_relaxed_state00_eV'] = float(np.sum(energies[i]-energies[0]))
    saved.update(fixed=fixed, bond=bond, cell=cell, origin=origin)
    np.savez_compressed(output, **saved)
    return dict(gap_A=gap, rows=rows, input_sha256=checksum(V4/'branch_positions.npz'),
                scope='same fixed grips and gap as SW; not each potential equilibrium loading',
                hessian_stability_certified=False, transition_barriers_certified=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('results/silicon_atomistic_v5'))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    target = args.output/'atomistic_controls.json'
    if target.exists():
        raise FileExistsError('refusing to replace an existing completed audit')
    started = time.perf_counter()
    specs, parameters = sources(args.output)
    result = dict(scope='static atomistic reference audit; not material/kinetic approval', models={})
    for name, (style, potential, source_hash) in specs.items():
        with LammpsSilicon(style, potential, sha256=source_hash) as engine:
            metadata = engine.metadata()
            checked = controls(engine, style, parameters)
            print(json.dumps(dict(event='controls', model=name, **checked)), flush=True)
            properties = bulk_properties(engine)
            print(json.dumps(dict(event='bulk', model=name, **properties)), flush=True)
            states = reconstruction_audit(engine, args.output/(name+'_reconstruction.npz'))
            result['models'][name] = dict(metadata=metadata, controls=checked,
                                         bulk=properties, reconstruction=states)
            save_json(args.output/'atomistic_controls_checkpoint.json', result)
    result.update(elapsed_seconds=time.perf_counter()-started,
                  code_sha256={str(p).replace('\\', '/'):checksum(p) for p in
                      [Path(__file__), Path('solver_v1/silicon_atomistic_reference.py')]},
                  complete=True)
    # Store repository-relative source keys, never a machine-specific script path.
    result['code_sha256'] = {str(Path(k).resolve().relative_to(Path.cwd())).replace('\\', '/'):v
                            for k, v in result['code_sha256'].items()}
    save_json(target, result)
    print(json.dumps(dict(event='complete', elapsed_seconds=result['elapsed_seconds'])), flush=True)


if __name__ == '__main__':
    main()
