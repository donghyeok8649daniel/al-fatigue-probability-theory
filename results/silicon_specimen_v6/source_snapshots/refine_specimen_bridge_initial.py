"""Audit optical basis, force convergence, boundary work and free-atom modes.

The initial 18 v6 calculations are preserved, including three failed force
checks. New corrected-boundary calculations are kept in a separate directory.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy.sparse.linalg import LinearOperator, cg, eigsh

from solver_v1.silicon_atomistic_reference import LammpsSilicon
from solver_v1.silicon_crystal_boundary import internal_strain_response, with_internal_relaxation
from solver_v1.silicon_specimen_research import AnisotropicModeI, atomic_crack_boundary
from .run_specimen_bridge import V5, load_case, observables
from .run_atomistic_controls import checksum, save_json
from .run_static_probe import EV_A3_TO_GPA


def optical_checks(engine, lattice, target_c44):
    response, info = internal_strain_response(engine, lattice, step=1e-4)
    response2, fine = internal_strain_response(engine, lattice, step=5e-5)
    cell = lattice*.5*np.array([[0., 1., 1.], [1., 0., 1.], [1., 1., 0.]])
    original = lattice*np.array([[0., 0., 0.], [.25, .25, .25]])
    zero = engine.evaluate(original, cell).energy
    rows = []
    for step in [1e-3, 5e-4, 2.5e-4]:
        energies, forces, affine = [], [], []
        for sign in [-1., 1.]:
            strain = np.zeros(6)
            strain[-1] = sign*step
            transform = np.eye(3)
            transform[0, 1] = transform[1, 0] = sign*step/2
            rr, cc = original@transform.T, cell@transform.T
            affine.append(engine.evaluate(rr, cc).energy)
            shift = response2@strain
            rr += np.array([-.5, .5])[:, None]*shift
            evaluated = engine.evaluate(rr, cc)
            energies.append(evaluated.energy)
            forces.append(float(np.max(abs(evaluated.gradient))))
        conversion = EV_A3_TO_GPA/(step**2*abs(np.linalg.det(cell)))
        rows.append(dict(engineering_shear=step, force_after_linear_correction_eV_A=max(forces),
                         C44_affine_GPa=(sum(affine)-2*zero)*conversion,
                         C44_optical_GPa=(sum(energies)-2*zero)*conversion))
    error = abs(rows[-1]['C44_optical_GPa']-target_c44)
    difference = float(np.max(abs(response-response2)))
    if difference > 1e-6 or error > 2e-3:
        raise AssertionError('internal strain response does not reproduce relaxed bulk elasticity')
    info.update(fine=fine, response_step_max_difference_A=difference,
                independent_bulk_shear_checks=rows, v5_C44_error_GPa=error)
    return response2, info


def operator(engine, state, step):
    free = ~state['fixed']
    initial = state['positions'].copy()
    x = initial[free].ravel().copy()

    def evaluate(value):
        r = initial.copy()
        r[free] = value.reshape(-1, 3)
        return engine.evaluate(r-state['origin'], state['cell'], periodic=(False, False, True))

    def force(value):
        return evaluate(value).gradient[free].ravel()

    def hessian_at(value):
        def action(direction):
            size = np.linalg.norm(direction)
            if size == 0:
                return np.zeros_like(direction)
            scale = step/size
            return (force(value+scale*direction)-force(value-scale*direction))/(2*scale)
        return LinearOperator((len(value), len(value)), matvec=action, dtype=float)
    return x, evaluate, force, hessian_at


def polish_and_modes(engine, path, output, *, modes=False):
    with np.load(path) as data:
        state = {k:data[k].copy() for k in data.files}
    x, evaluate, force, hessian_at = operator(engine, state, 1e-4)
    residuals = [float(np.max(abs(force(x))))]
    iterations = []
    for _ in range(4):
        if residuals[-1] < 2e-8:
            break
        step, status = cg(hessian_at(x), -force(x), rtol=2e-7, atol=1e-12, maxiter=1500)
        iterations.append(int(status))
        if status or not np.all(np.isfinite(step)) or np.max(abs(step)) > .02:
            break
        for scale in [1., .5, .25, .125]:
            next_x = x+scale*step
            next_residual = float(np.max(abs(force(next_x))))
            if next_residual < residuals[-1]:
                x = next_x
                residuals.append(next_residual)
                break
        else:
            break
    final = evaluate(x)
    state['positions'][~state['fixed']] = x.reshape(-1, 3)
    state.update(gradient=final.gradient, site_energy=final.site_energy)
    np.savez_compressed(output, **state)
    record = dict(input=path.name, output=output.name, residuals_eV_A=residuals,
                   cg_status=iterations, force_converged=bool(residuals[-1] < 2e-8),
                   max_position_change_A=None, full_q_stability_certified=False)
    with np.load(path) as source:
        record['max_position_change_A'] = float(np.max(abs(state['positions']-source['positions'])))
    if modes:
        spectra = []
        for delta in [1e-4, 5e-5]:
            _, _, _, hessian_at = operator(engine, state, delta)
            h = hessian_at(x)
            vals, vecs = eigsh(h, k=6, which='SA', tol=2e-7, maxiter=3000,
                              v0=np.sin(np.arange(len(x))+.273))
            rr = np.linalg.norm(h@vecs-vecs*vals, axis=0)
            spectra.append(dict(step_A=delta, eigenvalues_eV_A2=vals.tolist(), eigen_residuals_eV_A2=rr.tolist()))
        record['spectra'] = spectra
        record['eigenvalue_step_max_difference_eV_A2'] = max(abs(np.array(spectra[0]['eigenvalues_eV_A2'])-spectra[1]['eigenvalues_eV_A2']))
        record['full_q_stability_certified'] = bool(record['force_converged'] and
            min(spectra[-1]['eigenvalues_eV_A2']) > 1e-5 and max(spectra[-1]['eigen_residuals_eV_A2']) < 2e-5
            and record['eigenvalue_step_max_difference_eV_A2'] < 2e-5)
    return state, record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path('results/silicon_specimen_v6'))
    parser.add_argument('--output', type=Path, default=Path('results/silicon_specimen_v6/internal_boundary'))
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('refusing to overwrite an existing refinement')
    args.output.mkdir(parents=True)
    started = time.perf_counter()
    old = json.loads((args.source/'specimen_bridge.json').read_text())
    result = dict(scope='new optically relaxed boundary statics plus full free-atom force/Hessian/work checks',
                  source_sha256=checksum(args.source/'specimen_bridge.json'), models={}, cases=[],
                  boundary_work=[], original_failed_force_corrections=[],
                  code_sha256={p:checksum(p) for p in ['solver_v1/silicon_crystal_boundary.py',
                      'solver_v1/silicon_specimen_research.py', 'results/silicon_wafer_feasibility/refine_specimen_bridge.py',
                      'results/silicon_wafer_feasibility/run_specimen_bridge.py']})
    potentials = dict(original_sw_1985=Path('results/silicon_wafer_feasibility/source_Si.sw'), tersoff_1989=V5/'sources/SiC.tersoff')
    for label, info in old['models'].items():
        md = info['engine']
        elastic = AnisotropicModeI(*[info['bulk_input'][k] for k in ['C11_GPa', 'C12_GPa', 'C44_GPa']])
        lattice, kg = info['lattice_A'], info['griffith_K_MPa_sqrt_m']
        cutoff = 3.77118 if md['style'] == 'sw' else 3.2
        with LammpsSilicon(md['style'], potentials[label], sha256=md['potential_sha256']) as engine:
            response, optical = optical_checks(engine, lattice, info['bulk_input']['C44_GPa'])
            result['models'][label] = dict(internal_strain=optical, metadata=engine.metadata())
            print(json.dumps(dict(event='optical_check', model=label, C44_error=optical['v5_C44_error_GPa'])), flush=True)
            for case in old['cases']:
                if case['model'] == label and not case['relaxation']['force_converged']:
                    _, corrected = polish_and_modes(engine, args.source/case['artifact'],
                                                     args.output/('corrected_initial_'+case['artifact']))
                    result['original_failed_force_corrections'].append(corrected)
                    print(json.dumps(dict(event='old_force_polish', model=label, **corrected)), flush=True)
            for radius, front in [(28., 4), (40., 4), (56., 4), (28., 8)]:
                base = atomic_crack_boundary(lattice, radius_A=radius, front_repeats=front, grip_width_A=8.)
                boundary = with_internal_relaxation(base, lattice, response)
                previous, previous_k = None, 0.
                factors = [.8, 1., 1.2] if front == 4 else [1., 1.2]
                for factor in factors:
                    k = factor*kg
                    if previous is not None:
                        previous = previous+(boundary.displaced(elastic, k)-boundary.displaced(elastic, previous_k))
                    stem = f'{label}_R{radius:g}_front{front}_K{factor:g}'
                    path = args.output/(stem+'.npz')
                    positions, row = load_case(engine, elastic, boundary, k, previous, path)
                    row.update(model=label, factor_of_rigid_griffith=factor, front_repeats=front, artifact=path.name,
                               boundary_basis='linear optical response from same potential',
                               front=observables(boundary, positions, cutoff))
                    # R=28 at 1 and 1.2 KG: test all free Cartesian directions,
                    # including front-varying modes (not a coherent constraint).
                    needs_modes = radius == 28. and factor >= 1.
                    if needs_modes or not row['relaxation']['force_converged']:
                        state, checked = polish_and_modes(engine, path, args.output/(stem+'_polished.npz'), modes=needs_modes)
                        row['stationary_check'] = checked
                        positions = state['positions']
                        row['front_after_polish'] = observables(boundary, positions, cutoff)
                    result['cases'].append(row)
                    previous, previous_k = positions, k
                    save_json(args.output/(stem+'.json'), row)
                    save_json(args.output/'checkpoint.json', result)
                    print(json.dumps(dict(event='internal_boundary', model=label, radius=radius, front=front,
                        factor=factor, force_max=row['relaxation']['free_gradient_max_eV_A'],
                        full_stability=row.get('stationary_check', {}).get('full_q_stability_certified'))), flush=True)
                    if radius == 28. and front == 4 and factor == 1.:
                        unit = boundary.displaced(elastic, 1.)-boundary.reference
                        exact = float(np.sum(state['gradient'][boundary.fixed]*unit[boundary.fixed]))
                        for relative_step in [1e-3, 5e-4]:
                            h = kg*relative_step
                            energies, work_rows, changes = [], [], []
                            for sign in [-1, 1]:
                                wp = args.output/(stem+f'_work{relative_step:g}_{sign}.npz')
                                shifted, work = load_case(engine, elastic, boundary, k+sign*h, positions, wp, tolerance=2e-7)
                                # A force-root correction handles extensive energy roundoff.
                                ws, check = polish_and_modes(engine, wp, wp.with_name(wp.stem+'_polished.npz'))
                                energies.append(ws['site_energy'])
                                work_rows.append(check)
                                changes.append(float(np.max(abs(ws['positions']-positions))))
                            numerical = float(np.sum(energies[1]-energies[0])/(2*h))
                            result['boundary_work'].append(dict(model=label, relative_step=relative_step,
                                analytic_eV_per_K=exact, numerical_eV_per_K=numerical,
                                absolute_error=abs(numerical-exact), relative_error=abs(numerical-exact)/max(1., abs(exact)),
                                max_branch_displacement_A=max(changes), force_checks=work_rows))
            save_json(args.output/'checkpoint.json', result)
    result.update(complete=True, elapsed_seconds=time.perf_counter()-started, material_calibrated=False,
                  kinetic_calibrated=False, physical_clock_available=False)
    save_json(args.output/'refinement.json', result)
    print(json.dumps(dict(event='complete', cases=len(result['cases']), elapsed_seconds=result['elapsed_seconds'])), flush=True)


if __name__ == '__main__':
    main()
