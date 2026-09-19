"""Derivative, conditional stability and boundary work audit of actual Si states."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy.sparse.linalg import spsolve

from solver_v1.silicon_crack_research import RelaxedCoordinates, SpatialSW, diamond_crack_strip, relax_atoms
from .run_local_crack_audit import atomic_hessian, energy_difference, polish, save_json, spectrum
from .run_static_probe import source_parameters


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--calculation', type=Path, default=Path('results/silicon_local_crack_v2'))
    args = parser.parse_args()
    started = time.perf_counter()
    metadata = json.loads((args.calculation/'running_summary.json').read_text(encoding='utf-8'))
    p, _ = source_parameters()
    strip = diamond_crack_strip(metadata['lattice_A'], nx=6, ny=4, nz=4)
    model = SpatialSW(p, front_period=strip.front_period)
    data = np.load(args.calculation/'front4/localized/stationary_positions.npz')
    bond = data['bond']
    fixed, free = strip.fixed, np.flatnonzero(~strip.fixed)
    dofs = (free[:, None]*3+np.arange(3)).ravel()
    results = {}
    for name in ['initial', 'saddle', 'opened_minimum']:
        r = data[name]
        q = float(r[bond[1], 1]-r[bond[0], 1])
        coordinates = RelaxedCoordinates(r, fixed, bond=bond)
        hzz, h, tangent = atomic_hessian(model, coordinates, r, fixed)
        spec, _ = spectrum(hzz, 4)
        if spec['negative_index'] or spec['near_zero_eigenvalues']:
            raise RuntimeError('omitted Cartesian coordinates are not conditionally stable')
        zeros = np.zeros(coordinates.dimension)
        lift = (coordinates.positions(zeros, [q+1])-coordinates.positions(zeros, [q])).ravel()[dofs]
        cross = tangent.T@h@lift
        schur = float(lift@h@lift-cross@spsolve(hzz, cross))
        # The full analytic Hessian is checked against forces, not its own
        # energy expression differentiated with the same Jet class.
        direction = np.sin(np.arange(len(dofs))+.431)
        direction /= np.linalg.norm(direction)
        displacement = np.zeros_like(r)
        displacement[free] = direction.reshape(-1, 3)
        hv_checks = []
        for step in [.004, .002, .001]:
            plus = model.evaluate(r+step*displacement).gradient[free].ravel()
            minus = model.evaluate(r-step*displacement).gradient[free].ravel()
            derivative = (plus-minus)/(2*step)
            hv_checks.append(dict(step_A=step, max_error_eV_A2=float(max(abs(derivative-h@direction)))))
        curvature_checks = []
        for step in [.004, .002, .001]:
            reactions = []
            for sign in [-1, 1]:
                current, info = relax_atoms(model, coordinates, initial=r, values=[q+sign*step], tolerance=2e-6)
                current, correction = polish(model, coordinates, current, fixed, values=[q+sign*step])
                _, reaction = coordinates.pullback(model.evaluate(current).gradient)
                reactions.append(float(reaction[0]))
            curvature = (reactions[1]-reactions[0])/(2*step)
            curvature_checks.append(dict(step_A=step, difference_eV_A2=curvature,
                                        absolute_error_eV_A2=abs(curvature-schur)))
        if hv_checks[-1]['max_error_eV_A2'] > 1e-5 or curvature_checks[-1]['absolute_error_eV_A2'] > 1e-3:
            raise RuntimeError('spatial/relaxed Hessian refinement failed')
        results[name] = {'conditional_spectrum':spec, 'schur_curvature_eV_A2':schur,
            'force_derivative_checks':hv_checks, 'reaction_derivative_checks':curvature_checks}
        print(name, 'conditional lambda_min', spec['lowest_eigenvalues_eV_A2'][0], 'Schur', schur, flush=True)
    # Boundary loading work is conjugate to the actual imposed displacement
    # field. It includes the side grips: do not label it uniform remote stress.
    initial = data['initial']
    load_direction = strip.seed(4.5, tip=-3., transition=5.)-strip.seed(3.5, tip=-3., transition=5.)
    evaluation = model.evaluate(initial)
    reaction = float(np.sum(evaluation.gradient[fixed]*load_direction[fixed]))
    work_checks = []
    for step in [.004, .002, .001]:
        energies = []
        for sign in [-1, 1]:
            reference = initial.copy()
            reference[fixed] += sign*step*load_direction[fixed]
            coordinates = RelaxedCoordinates(reference, fixed, groups=strip.front_groups)
            current, info = relax_atoms(model, coordinates, initial=reference, tolerance=2e-6)
            current, correction = polish(model, coordinates, current, fixed)
            energies.append(energy_difference(model, current, initial))
        derivative = (energies[1]-energies[0])/(2*step)
        work_checks.append(dict(step_A=step, energy_derivative_eV_A=derivative,
                               absolute_error_eV_A=abs(derivative-reaction)))
    if work_checks[-1]['absolute_error_eV_A'] > 1e-4:
        raise RuntimeError('boundary virtual-work check failed')
    result = {'states':results, 'boundary_reaction_eV_A':reaction,
        'boundary_work_refinement':work_checks,
        'boundary_work_interpretation':'conjugate to prescribed nonuniform fixed-grip displacements, not uniform nominal stress',
        'total_force_balance_eV_A':evaluation.gradient.sum(axis=0).tolist(),
        'finite_temperature_PMF_calibrated':False, 'mobility_calibrated':False,
        'elapsed_seconds':time.perf_counter()-started, 'passed':True}
    save_json(args.calculation/'reduction_validation.json', result)


if __name__ == '__main__':
    main()
