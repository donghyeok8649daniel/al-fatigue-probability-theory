"""Actual stationary topology, not a sign test at another model's saddle.

Reuses the unchanged full-vector interface and branch/root routines. This is
static validation of rejected research candidates, not a PDE, dynamic hold,
finite-core activation, experimental yield or a physical-time prediction.
"""
import argparse
from pathlib import Path
import time

import numpy as np

from .run_current_material_core import ROOT, load_current_material
from .run_vector_material_calibration import source_and_targets
from .run_even_environment_validation import resolved_opening_extrema
from .run_vector_registry_audit import point_row, save_json
from .run_low_stress_cyclic_diagnostic import write_csv
from .vector_material_calibration import UNITS
from .vector_registry_audit import stationary_state, saddle_downhill_endpoints
from .vector_interface_reference import MishinVectorInterfaceReference


def registry_basis_reduced(model):
    """Source geometry is angstrom; candidate geometry already uses L0."""
    if isinstance(model, MishinVectorInterfaceReference):
        geometry, scale = model.reference.geometry, model.length
    else:
        geometry, scale = model.geometry, 1.
    return np.column_stack([geometry.a1, geometry.a2])/scale


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError('fresh actual-topology output required')
    began = time.perf_counter()
    source, _, states = source_and_targets()
    models = [('source_Al99_target_only', source)]
    bindings = []
    for label, path in (
        ('v22_wider', ROOT/'results/current_material_core_v22/wider_probe_validation/research_candidate_snapshot.json'),
        ('v23_radial_joint', ROOT/'results/core_interface_compatibility_v23/radial_full_validation/joint_LS/research_candidate_snapshot.json')):
        model, _, binding = load_current_material(path)
        models.append((label, model)); bindings.append(dict(model=label, **binding))
    cases = [(0., 0., 0.), (50., 4., -4.), (0., 0., 0.), (-50., -4., 4.), (0., 0., 0.)]
    save_json(args.out/'definition.json', dict(source_sha256=source.reference.sha256,
        bindings=bindings, source_stationary_seeds=states, static_traction_MPa=cases,
        actual_candidates_held_fixed=True, source_not_a_production_potential=True,
        source_coordinate_H_sign_is_not_candidate_stationary_topology=True,
        normal_extrema_bracket_samples=[33, 65], derivative_steps=[4e-5, 2e-5, 1e-5],
        all_results_diagnostic=True, new_blind_target_data=False,
        no_new_fitting=True, physical_time=False,
        dynamic_hold=False, material_accepted=False, actual_yield_validated=False))
    tables = {k: [] for k in ('stationary_points', 'saddle_connectivity', 'topology_summary',
                             'opening_extrema', 'static_traction_states', 'source_jet_FD')}
    # Check the source's nontrivial curvatures too; a numerical source jet is
    # not an experimentally exact Al datum or evidence against analytic sums.
    q_open = np.array([1.392126670481773, 0., 0.])
    for label, q in (('source_saddle', states['saddle']), ('inspected_opening', q_open)):
        q = np.asarray(q); exact = source.evaluate(q)
        for step in (4e-5, 2e-5, 1e-5):
            numerical = np.column_stack([(source.evaluate(q+step*e).gradient-source.evaluate(q-step*e).gradient)/(2*step)
                                         for e in np.eye(3)])
            tables['source_jet_FD'].append(dict(state=label, step_L0=step,
                Haa=exact.hessian[0, 0], Hxx=exact.hessian[1, 1],
                maximum_H_error=float(np.max(abs(numerical-exact.hessian))),
                physical_reference='0K Al99 matched rigid geometry, not measured curvature'))
    for name, model in models:
        roots = {}
        for label, index in (('perfect', 0), ('fault', 0), ('saddle', 1)):
            root = stationary_state(model, states[label], expected_index=index)
            roots[label] = root
            tables['stationary_points'].append(dict(
                **point_row(name, label, root['q'], root['evaluation'], UNITS),
                expected_index=index, actual_index=root['morse_index'], verified=root['valid'],
                stationarity_residual=root['force_residual'],
                distance_from_source_seed_L0=float(np.linalg.norm(root['q']-states[label]))))
        connectivity = False
        connectivity_error = None
        if roots['saddle']['valid']:
            try:
                endpoints = saddle_downhill_endpoints(model, roots['saddle'])
                for i, end in enumerate(endpoints):
                    tables['saddle_connectivity'].append(dict(
                        **point_row(name, f'downhill_{i}', end['q'], end['evaluation'], UNITS),
                        stationarity_residual=end['force_residual'], verified=end['valid'],
                        energy_decrease_eV_cell=end['energy_decrease']))
                # Distinct basins, modulo one triangular lattice cell. An
                # energy decrease alone could return to the SAME minimum.
                difference = endpoints[1]['q']-endpoints[0]['q']
                lattice = registry_basis_reduced(model)
                registry_fraction = np.linalg.solve(lattice, difference[1:])
                distinct = (abs(difference[0]) > 1e-7 or
                            np.linalg.norm(registry_fraction-np.rint(registry_fraction)) > 1e-7)
                connectivity = bool(all(e['valid'] for e in endpoints) and distinct)
            except (ArithmeticError, ValueError, RuntimeError) as error:
                connectivity_error = str(error)
        valid = all(r['valid'] for r in roots.values())
        e0 = roots['perfect']['evaluation'].energy
        summary = dict(model=name, local_roots_verified=valid,
            distinct_downhill_minima_verified=connectivity, connectivity_error=connectivity_error,
            fault_energy_J_m2=float(UNITS.energy_to_surface(roots['fault']['evaluation'].energy-e0)) if valid else None,
            saddle_energy_J_m2=float(UNITS.energy_to_surface(roots['saddle']['evaluation'].energy-e0)) if valid else None,
            reverse_barrier_J_m2=float(UNITS.energy_to_surface(roots['saddle']['evaluation'].energy-roots['fault']['evaluation'].energy)) if valid else None,
            rigid_opening_40h_J_m2=float(UNITS.energy_to_surface(model.evaluate((40*model.h, 0., 0.)).energy-e0)),
            opening_40h_is_finite_sample_not_infinite_limit=True,
            saddle_at_source_Hxx=model.evaluate(states['saddle']).hessian[1, 1],
            candidate_own_saddle_Hxx=roots['saddle']['evaluation'].hessian[1, 1],
            material_accepted=False)
        tables['topology_summary'].append(summary)
        for samples in (33, 65):
            tables['opening_extrema'].extend(dict(model=name, fixed_registry_not_coupled_spinodal=True, **r)
                for r in resolved_opening_extrema(model, samples))
        if roots['perfect']['valid']:
            q = roots['perfect']['q'].copy(); q0 = q.copy()
            for step, traction in enumerate(cases):
                result = stationary_state(model, q, force=UNITS.traction_mpa_to_force(traction), expected_index=0)
                q = result['q']
                tables['static_traction_states'].append(dict(model=name, step=step,
                    normal_MPa=traction[0], shear_x_MPa=traction[1], shear_y_MPa=traction[2],
                    delta_a_L0=q[0]-q0[0], delta_x_L0=q[1]-q0[1], delta_y_L0=q[2]-q0[2],
                    residual=result['force_residual'], verified=result['valid'],
                    minimum_H=result['eigenvalues'][0], static_return=step in (2, 4),
                    dynamic_hold=False, residual_plasticity_claimed=False))
                if not result['valid']:
                    break
        for filename, rows in tables.items():
            if rows:
                write_csv(args.out/(filename+'.csv'), rows)
        save_json(args.out/'progress.json', dict(completed=False, last_model=name,
            elapsed_seconds=time.perf_counter()-began))
        print(name, summary, flush=True)
    save_json(args.out/'completion.json', dict(completed=True,
        local_roots_verified=all(r['local_roots_verified'] for r in tables['topology_summary']),
        matched_saddle_connectivity_verified=all(r['distinct_downhill_minima_verified'] for r in tables['topology_summary']),
        actual_static_traction_states=len(tables['static_traction_states']),
        all_static_roots_verified=all(r['verified'] for r in tables['static_traction_states']),
        elapsed_seconds=time.perf_counter()-began, material_accepted=False,
        physical_time=False, actual_yield_validated=False, production_changed=False))


if __name__ == '__main__':
    main()
