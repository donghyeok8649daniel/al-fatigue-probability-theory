"""v20 independent static tests of completed, positive-LJ candidates.

No optimizer, no production PDE, no physical-time or yield inference. Target
jets are the 48 observations declared before the v20 controls/refits.
"""
import argparse
import hashlib
from pathlib import Path
import time

import numpy as np

from .interface_tangent_calibration import declared_tangent_validation
from .validate_coordination_screening import load_candidate
from .coordination_screening import CoordinationScreenedInterface, CoordinationScreenedBulk
from .run_vector_material_calibration import source_and_targets
from .vector_material_calibration import UNITS, LENGTH_M
from .run_vector_registry_audit import save_json, point_row
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_tensor_traction_audit import declared_tensor_scenarios
from .interface_static_scenarios import resolved_tensor_tractions
from .vector_registry_audit import stationary_state
from .run_even_environment_validation import resolved_opening_extrema, curvature_refinement
from .tail_constrained_material import declared_wavepoints


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--models', type=Path, nargs='+', required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError('fresh independent validation directory required')
    began = time.perf_counter()
    source, _, states = source_and_targets()
    observations = declared_tangent_validation(source)
    q_states = list(dict.fromkeys(o.state for o in observations))
    models = [('source_Al99_target_only', source, None)]
    bindings = []
    for directory in args.models:
        model, shape, c, law, definition = load_candidate(directory)
        if definition['source_sha256'] != source.reference.sha256:
            raise ValueError('source identity mismatch')
        if definition.get('new_validation_used_for_selection', False):
            raise ValueError('validation data were used for model selection')
        name = directory.parent.name+'/'+directory.name
        if name in [item[0] for item in models]:
            raise ValueError('unique model names required')
        models.append((name, model, (shape, c, law)))
        bindings.append(dict(model=name, shape=shape, coefficients=c, law=law,
            calibration_sha256=hashlib.sha256((directory/'calibration.json').read_bytes()).hexdigest()))
    save_json(args.out/'definition.json', dict(source_sha256=source.reference.sha256,
        model_bindings=bindings, new_validation_states=q_states,
        target_units='eV/cell, eV/L0, eV/L0^2; physical traction MPa',
        L0_m=LENGTH_M, density_and_reciprocal_tolerances=[2e-11, 2e-13],
        independent_radii=[12., 16.], independent_stretches=[.98, .995, 1., 1.007, 1.02],
        wavepoint_step=.08, derivative_steps=[4e-5, 2e-5, 1e-5],
        opening_extrema_bracket_samples=[65, 129], material_adoption_automatic=False,
        validation_used_for_selection=False, production_PDE_changed=False))
    tables = {name: [] for name in ('heldout_observations', 'heldout_summary', 'stationary_states',
        'static_tensor_states', 'static_returns', 'interface_curves', 'opening_extrema',
        'derivative_refinement', 'per_site_direct', 'finite_q_refinement', 'nonlinear_bloch_check',
        'term_budgets')}
    cases, normal, slip, _ = declared_tensor_scenarios()
    for name, model, parameters in models:
        zero = stationary_state(model, states['perfect'], expected_index=0)
        if not zero['valid']:
            raise ArithmeticError(f'{name}: pristine stable equilibrium not verified')
        for label, index in (('perfect', 0), ('fault', 0), ('saddle', 1)):
            result = stationary_state(model, states[label], expected_index=index)
            tables['stationary_states'].append(dict(
                **point_row(name, label, result['q'], result['evaluation'], UNITS),
                verified=result['valid'], force_residual=result['force_residual']))
        for case, tensor in cases.items():
            q = zero['q']; valid = True
            for fraction in (0., .5, 1., .5, 0., -.5, -1., -.5, 0.):
                traction = resolved_tensor_tractions(fraction*tensor, normal, slip)
                force = UNITS.traction_mpa_to_force(traction)
                result = stationary_state(model, q, force=force, expected_index=0)
                q = result['q']; delta = (q-zero['q'])*LENGTH_M/1e-10
                linear = np.linalg.solve(zero['evaluation'].hessian, force)*LENGTH_M/1e-10
                tables['static_tensor_states'].append(dict(model=name, case=case, fraction=fraction,
                    normal_MPa=traction[0], shear1_MPa=traction[1], shear2_MPa=traction[2],
                    delta_a_angstrom=delta[0], delta_x_angstrom=delta[1], delta_y_angstrom=delta[2],
                    linear_delta_a_angstrom=linear[0], linear_delta_x_angstrom=linear[1],
                    linear_delta_y_angstrom=linear[2], local_root_verified=result['valid'],
                    minimum_H=result['eigenvalues'][0], force_residual=result['force_residual'],
                    dynamic_hold=False, actual_yield_inferred=False))
                if not result['valid']:
                    valid = False
                    break
            tables['static_returns'].append(dict(model=name, case=case, path_verified=valid,
                return_error_L0=float(np.linalg.norm(q-zero['q'])) if valid else None,
                dynamic_hold_performed=False, residual_plasticity_validated=False))
        for path in ('opening', 'direct110', 'shockley112'):
            fractions = (np.r_[np.linspace(1., 2., 41), 2.5, 3., 5., 20., 40.]
                         if path == 'opening' else np.linspace(0., 1., 51))
            for fraction in fractions:
                q = ((fraction*model.h, 0., 0.) if path == 'opening' else
                     (model.h, fraction, 0.) if path == 'direct110' else
                     (model.h, .5*fraction, np.sqrt(3)*fraction/6))
                value = model.evaluate(q)
                tables['interface_curves'].append(dict(model=name, path=path, fraction=fraction,
                    energy_J_m2=float(UNITS.energy_to_surface(value.energy)),
                    traction_MPa=float(UNITS.force_to_traction_mpa(value.gradient[0])),
                    Haa=value.hessian[0, 0], Hxx=value.hessian[1, 1]))
        for samples in (65, 129):
            tables['opening_extrema'].extend(dict(model=name, **row)
                for row in resolved_opening_extrema(model, samples))
        print(name, 'signed MPa states, curves and refined normal extrema completed', flush=True)
        if parameters is None:
            continue
        shape, c, law = parameters
        for number, q in enumerate([(model.h, 0., 0.)]+q_states):
            value = model.evaluate(q)
            for component, jet in value.components.items():
                tables['term_budgets'].append(dict(model=name, state=number, component=component,
                    energy_eV_cell=jet[0], normal_force_eV_L0=jet[1],
                    Haa_eV_L0sq=jet[4], Hxx_eV_L0sq=jet[7]))
        for o in observations:
            value = model.evaluate(o.state)
            jet = np.r_[value.energy, value.gradient, value.hessian[0],
                        value.hessian[1, 1:], value.hessian[2, 2]]
            prediction = float(jet@o.jet_weights)
            tables['heldout_observations'].append(dict(model=name, observable=o.name,
                a_L0=o.state[0], x_L0=o.state[1], y_L0=o.state[2],
                target=o.target, prediction=prediction, scale=o.scale, units=o.units,
                normalized_residual=(prediction-o.target)/o.scale, fitted=False))
        for field in ('energy', 'force', 'Haa', 'Hxx'):
            rows = [r for r in tables['heldout_observations'] if r['model'] == name
                    and r['observable'].endswith('_'+field)]
            residual = np.array([r['normalized_residual'] for r in rows])
            tables['heldout_summary'].append(dict(model=name, field=field, count=len(rows),
                normalized_rms=float(np.sqrt(np.mean(residual**2))),
                maximum_absolute_normalized_error=float(max(abs(residual))),
                within_declared_scale=int(np.count_nonzero(abs(residual) <= 1))))
        fine = CoordinationScreenedInterface(shape, c, law=law, tolerance=2e-13)
        for number, q in enumerate(q_states[::3]):
            q = np.array(q); value = model.evaluate(q); tighter = fine.evaluate(q)
            for step in (4e-5, 2e-5, 1e-5):
                grad, H = [], []
                for e in np.eye(3):
                    hi, lo = model.evaluate(q+step*e), model.evaluate(q-step*e)
                    grad.append((hi.energy-lo.energy)/(2*step))
                    H.append((hi.gradient-lo.gradient)/(2*step))
                tables['derivative_refinement'].append(dict(model=name, state=number, step=step,
                    gradient_max_error=float(np.max(abs(grad-value.gradient))),
                    hessian_max_error=float(np.max(abs(np.asarray(H).T-value.hessian))),
                    reciprocal_H_change=float(np.max(abs(tighter.hessian-value.hessian))),
                    reciprocal_energy_change=abs(tighter.energy-value.energy)))
            exact, _ = model.screened_jet(q)
            for radius in (6, 10, 16):
                direct = model.direct_screened_energy(q, radius=radius, layers=radius)
                tables['per_site_direct'].append(dict(model=name, state=number, radius=radius,
                    reciprocal_unit_channel=exact[0], direct_unit_channel=direct,
                    absolute_error=abs(direct-exact[0])))
        points = declared_wavepoints(.08)
        for stretch in (.98, .995, 1., 1.007, 1.02):
            previous = None
            for radius in (12., 16.):
                bulk = CoordinationScreenedBulk(shape, radius=radius, stretch=stretch, law=law)
                current = []
                for j, q in enumerate(points):
                    columns, tails = bulk.evaluate(q)
                    H = np.einsum('c,cij->ij', c, columns); current.append(H)
                    bound = float(tails@abs(c)); minimum = float(np.linalg.eigvalsh(H)[0])
                    tables['finite_q_refinement'].append(dict(model=name, stretch=stretch, radius=radius,
                        qx=q[0], qy=q[1], qz=q[2], minimum_H=minimum, tail_bound=bound,
                        margin=minimum-bound, radius_change=None if previous is None else
                        float(np.linalg.norm(H-previous[j], 2))))
                previous = current
            print(name, 'independent dilation', stretch, 'completed', flush=True)
        bulk = CoordinationScreenedBulk(shape, radius=12., stretch=1.007, law=law)
        for axis in range(3):
            tables['nonlinear_bloch_check'].extend(dict(model=name, axis=axis, **row)
                for row in curvature_refinement(bulk, c, alpha=shape[3], polarization=np.eye(3)[axis]))
        save_json(args.out/'progress.json', dict(completed=False, last_model=name))
    for filename, rows in tables.items():
        write_csv(args.out/(filename+'.csv'), rows)
    save_json(args.out/'completion.json', dict(completed=True,
        actual_tensor_states=len(tables['static_tensor_states']),
        all_tensor_roots_verified=all(r['local_root_verified'] for r in tables['static_tensor_states']),
        max_final_FD_H_error=max(r['hessian_max_error'] for r in tables['derivative_refinement'] if r['step'] == 1e-5),
        min_sampled_spectral_margin=min(r['margin'] for r in tables['finite_q_refinement']),
        material_accepted=False, kinetics_calibrated=False, physical_PDE_Hz=False,
        actual_yield_validated=False, dynamic_hold_performed=False, whole_zone_proof=False,
        elapsed_seconds=time.perf_counter()-began))


if __name__ == '__main__':
    main()
