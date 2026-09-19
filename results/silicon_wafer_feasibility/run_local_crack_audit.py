"""Actual spatial Si crack calculations; never a replay of a fitted barrier.

Run: python -m results.silicon_wafer_feasibility.run_local_crack_audit
The finite fixed-grip original-SW strip is a research control, not a wafer fit.
No potential parameter is altered to make fracture easier. No bond is disabled.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.sparse.linalg import eigsh, spsolve

from solver_v1.silicon_crack_research import (
    RelaxedCoordinates, SpatialSW, diamond_crack_strip, relax_atoms,
)
from .run_static_probe import PARAMETER_URL, pair, source_parameters


def save_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def atomic_hessian(model, coordinates, positions, fixed):
    free = np.flatnonzero(~fixed)
    dofs = (free[:, None]*3+np.arange(3)).ravel()
    h = model.hessian(positions, free_atoms=free)
    b = coordinates.tangent_matrix()[dofs]
    return (b.T@h@b).tocsr(), h, b


def polish(model, coordinates, positions, fixed, *, values=(), tolerance=1e-9):
    """Newton root correction near an already located minimum or saddle.

    This is NOT minimization: negative curvature is retained at a saddle.
    Reject large corrections and record residuals, not just a solver flag.
    """
    x = coordinates.encode(positions)
    history = []
    for iteration in range(6):
        positions = coordinates.positions(x, values)
        evaluation = model.evaluate(positions)
        g, _ = coordinates.pullback(evaluation.gradient)
        residual = float(np.max(abs(g), initial=0))
        history.append(residual)
        if residual <= tolerance:
            break
        h, _, _ = atomic_hessian(model, coordinates, positions, fixed)
        step = spsolve(h, -g)
        if not np.all(np.isfinite(step)) or np.max(abs(step), initial=0) > .2:
            raise RuntimeError('Newton correction left the local stationary branch')
        for scale in (1., .5, .25, .125):
            proposed = coordinates.positions(x+scale*step, values)
            gg, _ = coordinates.pullback(model.evaluate(proposed).gradient)
            if np.linalg.norm(gg) < np.linalg.norm(g):
                x += scale*step
                break
        else:
            raise RuntimeError('stationary correction failed to reduce the force residual')
    if history[-1] > tolerance:
        raise RuntimeError(f'stationary residual {history[-1]:.4g} exceeds {tolerance}')
    return positions, {'residual_history_eV_A':history, 'tolerance_eV_A':tolerance,
                       'converged':True, 'newton_steps':len(history)-1}


def spectrum(h, count):
    k = min(count, h.shape[0]-1)
    values, vectors = eigsh(h, k=k, which='SA', tol=2e-10,
        v0=np.sin(np.arange(h.shape[0])+.731))
    residual = np.linalg.norm(h@vectors-vectors*values, axis=0)
    if values[-1] <= 0 or max(residual) > 2e-7:
        raise RuntimeError('lowest spectrum does not resolve the complete negative index')
    return {'lowest_eigenvalues_eV_A2':values.tolist(),
            'eigenpair_residuals_eV_A2':residual.tolist(),
            'negative_index':int(sum(values < -1e-8)),
            'near_zero_eigenvalues':int(sum(abs(values) <= 1e-8))}, vectors[:, 0]


def energy_difference(model, positions, reference):
    return float(np.sum(model.evaluate(positions).site_energy-model.evaluate(reference).site_energy))


def frontier(strip, model, positions, bond_index):
    bonds = strip.crossing_bonds
    middle = strip.positions[bonds].mean(axis=1)
    gaps = strip.bond_gaps(positions)
    selected = np.isclose(middle[:, 0], middle[bond_index, 0], atol=1e-7)
    order = np.argsort(middle[selected, 2])
    return {'bond_index':int(bond_index), 'front_x_A':float(middle[bond_index, 0]),
        'front_z_A':middle[selected, 2][order].tolist(),
        'front_normal_gaps_A':gaps[selected][order].tolist(),
        'original_bonds_with_normal_gap_beyond_SW_cutoff':int(sum(gaps > model.cutoff)),
        'interpretation':'reversible geometric gap observable; not an irreversible damage flag'}


def select_tip_bond(strip, model, positions):
    bonds = strip.crossing_bonds
    middle = strip.positions[bonds].mean(axis=1)
    gaps = strip.bond_gaps(positions)
    eligible = (~np.any(strip.fixed[bonds], axis=1)) & (gaps < model.cutoff) & (middle[:, 0] > -1.)
    if not np.any(eligible):
        raise RuntimeError('no intact mobile bond in front of the seeded crack')
    ids = np.flatnonzero(eligible)
    return int(ids[np.argmin(middle[ids, 0]+.01*abs(middle[ids, 2]-.5*strip.front_period))])


def profile(model, strip, initial, bond_index, *, coherent, intervals, output):
    bond = strip.crossing_bonds[bond_index]
    groups = strip.front_groups if coherent else None
    constrained = RelaxedCoordinates(initial, strip.fixed, groups=groups, bond=bond)
    unconstrained = RelaxedCoordinates(initial, strip.fixed, groups=groups)
    start_gap = float(strip.bond_gaps(initial)[bond_index])
    cache, rows, configurations = {}, [], []
    q_grid = np.linspace(start_gap, 4.65, intervals+1)
    previous = initial.copy()

    def evaluate(q, start=None):
        nonlocal previous
        key = float(q)
        if key in cache:
            return cache[key]
        if start is None:
            start = cache[min(cache, key=lambda x: abs(x-key))][1] if cache else previous
        current, info = relax_atoms(model, constrained, initial=start, values=[q],
                                    tolerance=2e-6, maxiter=2000)
        if not info['force_converged']:
            current, correction = polish(model, constrained, current, strip.fixed, values=[q])
            info['correction_after_failed_force_check'] = correction
        evaluation = model.evaluate(current)
        g, reaction = constrained.pullback(evaluation.gradient)
        record = dict(q_A=key, relative_energy_eV=energy_difference(model, current, initial),
            reaction_eV_A=float(reaction[0]), free_gradient_max_eV_A=float(max(abs(g))),
            relaxation=info)
        cache[key] = record, current
        return record, current

    for q in q_grid:
        record, previous = evaluate(q, previous)
        rows.append(record)
        configurations.append(previous.copy())
    output.mkdir(parents=True)
    save_json(output/'profile.json', rows)
    np.savez_compressed(output/'profile_positions.npz', q_A=q_grid,
                        positions=np.asarray(configurations), bond=bond)
    with (output/'profile.csv').open('w', newline='', encoding='utf-8') as stream:
        fields = ['q_A', 'relative_energy_eV', 'reaction_eV_A', 'free_gradient_max_eV_A']
        writer = csv.DictWriter(stream, fields)
        writer.writeheader()
        writer.writerows({key:row[key] for key in fields} for row in rows)
    saddle_brackets = [(a['q_A'], b['q_A']) for a, b in zip(rows, rows[1:])
        if a['reaction_eV_A'] > 1e-5 and b['reaction_eV_A'] < -1e-5]
    if not saddle_brackets:
        raise RuntimeError('no resolved positive-to-negative reaction root')
    lo, hi = saddle_brackets[0]
    saddle_q = brentq(lambda q: evaluate(q)[0]['reaction_eV_A'], lo, hi, xtol=2e-8)
    minimum_brackets = [(a['q_A'], b['q_A']) for a, b in zip(rows, rows[1:])
        if a['q_A'] > saddle_q and a['reaction_eV_A'] < -1e-5 and b['reaction_eV_A'] > 1e-5]
    if not minimum_brackets:
        raise RuntimeError('no first local minimum after the opening saddle')
    lo, hi = minimum_brackets[0]
    minimum_q = brentq(lambda q: evaluate(q)[0]['reaction_eV_A'], lo, hi, xtol=2e-8)
    states, details = {}, {}
    for name, q in [('saddle', saddle_q), ('opened_minimum', minimum_q)]:
        current = evaluate(q)[1]
        current, correction = polish(model, unconstrained, current, strip.fixed)
        h_restricted, h_full, tangent = atomic_hessian(model, unconstrained, current, strip.fixed)
        local_spectrum, unstable = spectrum(h_restricted, 6)
        full_spectrum = spectrum(h_full, max(6, int(round(strip.front_period/strip.basis_period))+2))[0]
        if local_spectrum['negative_index'] != (1 if name == 'saddle' else 0):
            raise RuntimeError('stationary point has the wrong restricted Morse index')
        if not coherent and full_spectrum['negative_index'] != local_spectrum['negative_index']:
            raise RuntimeError('independent atom coordinates lost a free direction')
        details[name] = {'q_A':float(strip.bond_gaps(current)[bond_index]),
            'relative_energy_eV':energy_difference(model, current, initial),
            'stationary_correction':correction, 'restricted_spectrum':local_spectrum,
            'full_atom_spectrum':full_spectrum, 'frontier':frontier(strip, model, current, bond_index)}
        states[name] = current
        if name == 'saddle':
            unstable_mode = np.asarray(unconstrained.tangent_matrix()@unstable).reshape(-1, 3)
    # Verify both downhill connections from the actual saddle, rather than
    # identifying a maximum of an arbitrary sampled interpolation as a barrier.
    connections = []
    for sign in (-1, 1):
        current, info = relax_atoms(model, unconstrained,
            initial=states['saddle']+sign*.06*unstable_mode, tolerance=2e-6)
        current, correction = polish(model, unconstrained, current, strip.fixed)
        h, _, _ = atomic_hessian(model, unconstrained, current, strip.fixed)
        spec, _ = spectrum(h, 4)
        distances = {'initial':float(np.max(abs(current-initial))),
            'opened_minimum':float(np.max(abs(current-states['opened_minimum'])))}
        target = min(distances, key=distances.get)
        if distances[target] > 2e-5 or spec['negative_index']:
            raise RuntimeError('saddle downhill branch did not reach the expected stable minimum')
        connections.append(dict(sign=sign, target=target, max_position_errors_A=distances,
            relaxation=info, stationary_correction=correction, spectrum=spec))
    if {x['target'] for x in connections} != {'initial', 'opened_minimum'}:
        raise RuntimeError('the two saddle branches did not connect different minima')
    # Independent relaxed-energy derivative refinement below the saddle.
    qcheck = start_gap+.4*(saddle_q-start_gap)
    center, center_positions = evaluate(qcheck)
    derivative_checks = []
    for step in (.004, .002, .001):
        a = evaluate(qcheck-step, center_positions)[0]
        b = evaluate(qcheck+step, center_positions)[0]
        derivative = (b['relative_energy_eV']-a['relative_energy_eV'])/(2*step)
        derivative_checks.append({'step_A':step, 'difference_eV_A':derivative,
            'reaction_eV_A':center['reaction_eV_A'],
            'absolute_error_eV_A':abs(derivative-center['reaction_eV_A'])})
    if derivative_checks[-1]['absolute_error_eV_A'] > 1e-4:
        raise RuntimeError('relaxed energy derivative does not match the constraint reaction')
    np.savez_compressed(output/'stationary_positions.npz', initial=initial,
        saddle=states['saddle'], opened_minimum=states['opened_minimum'],
        unstable_mode=unstable_mode, bond=bond)
    summary = {'coherent_front_control':coherent, 'intervals':intervals,
        'normalization':'total eV of the explicit supercell; no Ac, area or fitted multiplicity',
        'states':details, 'connections':connections, 'derivative_checks':derivative_checks,
        'forward_barrier_eV':details['saddle']['relative_energy_eV'],
        'reverse_barrier_eV':details['saddle']['relative_energy_eV']-details['opened_minimum']['relative_energy_eV'],
        'physical_activation_free_energy_eV':None, 'kinetic_closure_certified':False}
    save_json(output/'summary.json', summary)
    return summary, states


def run_case(parameters, lattice, *, nx, ny, nz, opening, intervals, output):
    started = time.perf_counter()
    strip = diamond_crack_strip(lattice, nx=nx, ny=ny, nz=nz, grip_width=4.)
    model = SpatialSW(parameters, front_period=strip.front_period)
    seed = strip.seed(opening, tip=-3., transition=5.)
    coherent_coordinates = RelaxedCoordinates(seed, strip.fixed, groups=strip.front_groups)
    initial, initial_info = relax_atoms(model, coherent_coordinates, tolerance=2e-6)
    initial, correction = polish(model, coherent_coordinates, initial, strip.fixed)
    full_coordinates = RelaxedCoordinates(seed, strip.fixed)
    h, _, _ = atomic_hessian(model, full_coordinates, initial, strip.fixed)
    stable, _ = spectrum(h, 5)
    if stable['negative_index'] or stable['near_zero_eigenvalues']:
        raise RuntimeError('initial crack is not a resolved full-atom local minimum')
    bond_index = select_tip_bond(strip, model, initial)
    output.mkdir(parents=True)
    np.savez_compressed(output/'geometry.npz', reference=strip.positions, fixed=strip.fixed,
        front_groups=strip.front_groups, frame=strip.frame, crossing_bonds=strip.crossing_bonds,
        seed=seed, initial=initial)
    metadata = {'nx':nx, 'ny':ny, 'nz':nz, 'atoms':len(initial),
        'free_atoms':int(sum(~strip.fixed)), 'front_period_A':strip.front_period,
        'width_A':strip.width, 'height_A':strip.height, 'grip_width_A':4.,
        'imposed_opening_A':opening, 'seed_tip_A':-3., 'seed_transition_A':5.,
        'cut':'111 shuffle', 'propagation_cubic':[1, 1, -2], 'front_cubic':[1, -1, 0],
        'initial_relaxation':initial_info, 'initial_correction':correction,
        'initial_full_spectrum':stable, 'tip_bond_index':bond_index,
        'initial_frontier':frontier(strip, model, initial, bond_index)}
    save_json(output/'geometry.json', metadata)
    summaries = {}
    for coherent, label in [(False, 'localized'), (True, 'coherent')]:
        summary, _ = profile(model, strip, initial, bond_index, coherent=coherent,
                             intervals=intervals, output=output/label)
        summaries[label] = summary
        print(json.dumps({'case':output.name, 'path':label,
            'barrier_eV':summary['forward_barrier_eV'],
            'reverse_barrier_eV':summary['reverse_barrier_eV'],
            'full_saddle_index':summary['states']['saddle']['full_atom_spectrum']['negative_index'],
            'seconds':time.perf_counter()-started}), flush=True)
    result = {'geometry':metadata, 'paths':summaries, 'elapsed_seconds':time.perf_counter()-started}
    save_json(output/'summary.json', result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('results/silicon_local_crack_v2'))
    parser.add_argument('--front-repeats', type=int, nargs='+', default=[4, 8, 12])
    parser.add_argument('--nx', type=int, default=6)
    parser.add_argument('--ny', type=int, default=4)
    parser.add_argument('--opening', type=float, default=3.5)
    parser.add_argument('--intervals', type=int, default=32)
    args = parser.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        raise FileExistsError('use a fresh output folder; preserve previous calculations')
    args.output.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    parameters, checksum = source_parameters()
    bond = minimize_scalar(lambda r: float(pair(r, parameters)), bounds=(2.2, 2.5),
                           method='bounded', options={'xatol':1e-14}).x
    lattice = 4*bond/np.sqrt(3)
    paths = [Path(__file__), Path('solver_v1/silicon_crack_research.py'),
             Path('solver_v1/silicon_environment_research.py')]
    summary = {'scope':'0 K finite fixed-grip original-SW spatial crack; not wafer strength or fatigue calibration',
        'source_url':PARAMETER_URL, 'source_sha256':checksum, 'parameters':parameters,
        'lattice_A':float(lattice), 'bond_A':float(bond),
        'code_sha256':{str(p.relative_to(Path.cwd()) if p.is_absolute() else p).replace('\\', '/'):
            hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'production_parameters_changed':False, 'physical_time_seconds':None,
        'wafer_strength_GPa':None, 'physical_activation_free_energy_eV':None, 'cases':[]}
    save_json(args.output/'running_summary.json', summary)
    for nz in args.front_repeats:
        result = run_case(parameters, lattice, nx=args.nx, ny=args.ny, nz=nz,
            opening=args.opening, intervals=args.intervals, output=args.output/f'front{nz}')
        summary['cases'].append(result)
        save_json(args.output/'running_summary.json', summary)
    summary['elapsed_seconds'] = time.perf_counter()-started
    summary['completed'] = True
    save_json(args.output/'summary.json', summary)
    print(json.dumps({'completed':True, 'elapsed_seconds':summary['elapsed_seconds']}), flush=True)


if __name__ == '__main__':
    main()
