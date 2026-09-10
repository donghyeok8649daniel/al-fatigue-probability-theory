"""Independent static validation of v19 fits; no refit and no production gate."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from .coordination_screening import CoordinationScreenedInterface, CoordinationScreenedBulk
from .run_vector_material_calibration import source_and_targets
from .run_vector_registry_audit import save_json, point_row
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_tensor_traction_audit import declared_tensor_scenarios
from .interface_static_scenarios import resolved_tensor_tractions
from .vector_material_calibration import IDEAL_H, UNITS, LENGTH_M
from .vector_registry_audit import stationary_state
from .run_even_environment_validation import resolved_opening_extrema, curvature_refinement
from .tail_constrained_material import declared_wavepoints


def load_candidate(directory):
    directory = Path(directory)
    data = json.loads((directory/'calibration.json').read_bytes())
    definition = json.loads((directory.parent/'definition.json').read_bytes())
    if not data['completed'] or data['best'] is None or not data['best']['strictly_positive_LJ']:
        raise ValueError('completed positive-LJ candidate required, not a zero-pair closure')
    shape = np.asarray(data.get('shape', data['best'].get('shape')), float)
    law = data.get('screening_law', definition.get('screening_law', 'rational'))
    c = np.asarray(data['best']['coefficients'], float)
    model = CoordinationScreenedInterface(shape, c, law=law)
    return model, shape, c, law, definition


def independent_states():
    """Additional validation declared before joint-shape completion, not fit."""
    states = [(r*IDEAL_H, 0., 0.) for r in (1.009, 1.079, 1.209, 1.469, 1.649, 1.929)]
    states += [(1.061*IDEAL_H, .287, -.057), (1.121*IDEAL_H, .427, .117),
               (1.031*IDEAL_H, .371, .201), (.987*IDEAL_H, .097, -.033)]
    return states


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--models', type=Path, nargs='+', required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError('preserve prior validation; fresh output directory required')
    began = time.perf_counter()
    source, _, states = source_and_targets()
    validation_q = independent_states()
    save_json(args.out/'definition.json', dict(independent_states=validation_q,
        cutoff_radii=[12., 16.], stretches=[.98, .995, 1., 1.007, 1.02],
        wavepoint_step=.08, bracket_samples=[65, 129],
        numerical_derivative_steps=[4e-5, 2e-5, 1e-5],
        source_sha256=source.reference.sha256, production_changed=False))
    models = [('source_Al99', source, None)]
    bindings = []
    for directory in args.models:
        model, shape, c, law, definition = load_candidate(directory)
        if definition['source_sha256'] != source.reference.sha256:
            raise ValueError('source provenance mismatch')
        old_observations = definition['observations']
        old_states = {tuple(o['state']) for o in old_observations if o['state'] is not None}
        if any(q in old_states for q in validation_q):
            raise ValueError('new independent state overlaps prior data')
        name = directory.parent.name+'/'+directory.name
        if name in [row[0] for row in models]:
            raise ValueError('unique result model names required')
        models.append((name, model, (shape, c, law)))
        bindings.append(dict(model=name, calibration_sha256=hashlib.sha256(
            (directory/'calibration.json').read_bytes()).hexdigest(), shape=shape,
            coefficients=c, law=law))
    save_json(args.out/'model_bindings.json', bindings)
    derivative_rows, direct_rows, spectral_rows, bloch_rows = [], [], [], []
    root_rows, static_rows, extrema_rows, curves, budgets, heldout = [], [], [], [], [], []
    cases, n, m, _ = declared_tensor_scenarios()
    roots_by_model = {}
    for name, model, parameters in models:
        root_data = {}
        for label, index in (('perfect', 0), ('fault', 0), ('saddle', 1)):
            root = stationary_state(model, states[label], expected_index=index)
            root_data[label] = dict(q=root['q'], valid=root['valid'], morse_index=root['morse_index'])
            root_rows.append(dict(model=name, label=label, valid=root['valid'],
                **{k: v for k, v in point_row(name, label, root['q'], root['evaluation'], UNITS).items()
                   if k not in ('model', 'state')}, force_residual=root['force_residual']))
        roots_by_model[name] = root_data
        zero = stationary_state(model, states['perfect'], expected_index=0)
        if not zero['valid']:
            raise ArithmeticError('pristine mechanical equilibrium failed')
        for case, tensor in cases.items():
            q = zero['q']
            for fraction in (0., .5, 1., .5, 0., -.5, -1., -.5, 0.):
                traction = resolved_tensor_tractions(fraction*tensor, n, m)
                out = stationary_state(model, q, force=UNITS.traction_mpa_to_force(traction), expected_index=0)
                q = out['q']; delta = q-zero['q']
                static_rows.append(dict(model=name, case=case, load_fraction=fraction,
                    normal_MPa=traction[0], shear1_MPa=traction[1], shear2_MPa=traction[2],
                    delta_a_angstrom=delta[0]*LENGTH_M/1e-10,
                    delta_x_angstrom=delta[1]*LENGTH_M/1e-10,
                    delta_y_angstrom=delta[2]*LENGTH_M/1e-10,
                    local_root_verified=out['valid'], force_residual=out['force_residual'],
                    minimum_H=out['eigenvalues'][0], dynamic_hold=False))
                if not out['valid']:
                    break
        for samples in (65, 129):
            extrema_rows.extend(dict(model=name, **row) for row in resolved_opening_extrema(model, samples))
        for path in ('opening', 'direct110', 'shockley112'):
            fractions = np.r_[np.linspace(1., 2., 31), 2.5, 3., 5., 20., 40.] if path == 'opening' else np.linspace(0., 1., 41)
            for f in fractions:
                q = ((f*model.h, 0., 0.) if path == 'opening' else
                     (model.h, f, 0.) if path == 'direct110' else
                     (model.h, .5*f, np.sqrt(3)*f/6))
                value = model.evaluate(q)
                curves.append(dict(model=name, path=path, fraction=f,
                    energy_J_m2=float(UNITS.energy_to_surface(value.energy)),
                    normal_MPa=float(UNITS.force_to_traction_mpa(value.gradient[0])),
                    Haa=value.hessian[0, 0], Hxx=value.hessian[1, 1]))
        print(name, 'stationary/tensor/curve/extrema checks completed', flush=True)
        if parameters is None:
            continue
        shape, c, law = parameters
        fine = CoordinationScreenedInterface(shape, c, law=law, tolerance=2e-13)
        for i, q in enumerate(validation_q):
            result, reference = model.evaluate(q), source.evaluate(q)
            for field, predicted, target, scale in (
                ('energy', result.energy, reference.energy, max(.1*abs(reference.energy), .01/float(UNITS.energy_to_surface(1.)))),
                ('normal_force', result.gradient[0], reference.gradient[0], float(UNITS.traction_mpa_to_force(250.))),
                ('Haa', result.hessian[0, 0], reference.hessian[0, 0], max(.1*abs(reference.hessian[0, 0]), .1)),
                ('Hxx', result.hessian[1, 1], reference.hessian[1, 1], max(.1*abs(reference.hessian[1, 1]), .1))):
                heldout.append(dict(model=name, state=i, field=field, target=target, prediction=predicted,
                                    scale=scale, normalized_residual=(predicted-target)/scale))
            for term, jet in result.components.items():
                budgets.append(dict(model=name, state=i, component=term, energy=jet[0],
                                    normal_force=jet[1], Haa=jet[4], Hxx=jet[7]))
        for i, q in enumerate(validation_q[::2]):
            q = np.asarray(q); result = model.evaluate(q)
            finer = fine.evaluate(q)
            for step in (4e-5, 2e-5, 1e-5):
                grad, H = [], []
                for e in np.eye(3):
                    hi, lo = model.evaluate(q+step*e), model.evaluate(q-step*e)
                    grad.append((hi.energy-lo.energy)/(2*step)); H.append((hi.gradient-lo.gradient)/(2*step))
                derivative_rows.append(dict(model=name, state=i, step=step,
                    gradient_max_error=float(np.max(abs(np.asarray(grad)-result.gradient))),
                    hessian_max_error=float(np.max(abs(np.asarray(H).T-result.hessian))),
                    reciprocal_H_change=float(np.max(abs(finer.hessian-result.hessian))),
                    reciprocal_energy_change=abs(finer.energy-result.energy)))
            exact, diag = model.screened_jet(q)
            for radius in (6, 10, 16):
                direct = model.direct_screened_energy(q, radius=radius, layers=radius)
                direct_rows.append(dict(model=name, state=i, radius=radius, layers=radius,
                    reciprocal_unit_channel=exact[0], direct_unit_channel=direct,
                    absolute_error=abs(direct-exact[0]), reciprocal_layers=diag['vector']['layers'],
                    reciprocal_shells=diag['vector']['shells'], empirical_tail_indicator=diag['vector']['last_layer_jet']))
        points = declared_wavepoints(.08)
        for stretch in (.98, .995, 1., 1.007, 1.02):
            previous = None
            for radius in (12., 16.):
                bulk = CoordinationScreenedBulk(shape, radius=radius, stretch=stretch, law=law)
                current = []
                for j, q in enumerate(points):
                    cols, tail = bulk.evaluate(q)
                    H = np.einsum('c,cij->ij', c, cols); current.append(H)
                    bound = float(tail@abs(c)); minimum = float(np.linalg.eigvalsh(H)[0])
                    spectral_rows.append(dict(model=name, stretch=stretch, radius=radius,
                        qx=q[0], qy=q[1], qz=q[2], minimum_H=minimum, tail_bound=bound,
                        margin=minimum-bound, radius_change=None if previous is None else
                        float(np.linalg.norm(H-previous[j], 2))))
                previous = current
            print(name, 'independent dilation', stretch, 'finished', flush=True)
        bulk = CoordinationScreenedBulk(shape, radius=12., stretch=1.007, law=law)
        for axis in range(3):
            bloch_rows.extend(dict(model=name, axis=axis, **r) for r in curvature_refinement(
                bulk, c, alpha=shape[3], polarization=np.eye(3)[axis]))
        print(name, 'analytic/direct/dilated Bloch checks completed', flush=True)
        save_json(args.out/'progress.json', dict(completed=False, last_model=name))
    tables = dict(stationary_states=root_rows, static_tensor_states=static_rows,
        opening_extrema=extrema_rows, interface_curves=curves, independent_observations=heldout,
        derivative_refinement=derivative_rows, per_site_direct=direct_rows,
        finite_q_refinement=spectral_rows, displaced_site_harmonic=bloch_rows, term_budgets=budgets)
    for filename, rows in tables.items():
        write_csv(args.out/(filename+'.csv'), rows)
    save_json(args.out/'completion.json', dict(completed=True, roots=roots_by_model,
        analytic_jet_error=max(r['hessian_max_error'] for r in derivative_rows if r['step'] == 1e-5),
        per_site_direct_fine_error=max(r['absolute_error'] for r in direct_rows if r['radius'] == 16),
        minimum_sampled_spectral_margin=min(r['margin'] for r in spectral_rows),
        whole_zone_proof=False, material_accepted=False, actual_yield_validated=False,
        physical_PDE_seconds=False, physical_PDE_Hz=False,
        elapsed_seconds=time.perf_counter()-began))


if __name__ == '__main__':
    main()
