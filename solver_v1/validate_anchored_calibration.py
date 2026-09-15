"""Compare a completed training-selected shape with its parent on excluded states."""
import argparse
import json
from pathlib import Path
import time

import numpy as np

from .coordination_screening import CoordinationScreenedInterface, CoordinationScreenedBulk
from .periodic_plane_covariance import SourceEAMBlochHessian
from .run_finite_q_compatibility_v31 import points
from .run_joint_shape_v31 import decode_shape
from .run_low_frequency_forcing_v29 import sha
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_source_core_reference import load_source_material
from .run_vector_registry_audit import save_json
from .vector_interface_reference import MishinVectorInterfaceReference
from .vector_material_calibration import LENGTH_M


def run(candidate, baseline, out):
    candidate, baseline, out = map(Path, (candidate, baseline, out))
    if out.exists():
        raise FileExistsError('fresh independent validation directory required')
    started = time.perf_counter()
    plan = json.loads((candidate/'validation_plan.json').read_bytes())
    observations = json.loads((candidate/'observations.json').read_bytes())
    if sha(candidate/'observations.json') != sha(baseline.parent/'observations.json'):
        raise ValueError('baseline and candidate must use identical training targets and scales')
    fitted_states = {tuple(o['state']) for o in observations if o.get('state') is not None}
    if fitted_states & {tuple(q) for q in plan['interface_states']}:
        raise ValueError('declared new validation must not overlap the training states')
    if {tuple(q) for _, q in points()} & {tuple(q) for q in plan['finite_q']}:
        raise ValueError('declared new wavevectors must not overlap the old eight points')
    source, _, binding = load_source_material()
    factor = LENGTH_M/1e-10
    ref_interface = MishinVectorInterfaceReference(source.source, factor)
    ref_bulk = SourceEAMBlochHessian(source.source)
    state_targets = [ref_interface.evaluate(q) for q in plan['interface_states']]
    all_points = list(points())+[('new_excluded', tuple(q)) for q in plan['finite_q']]
    q_targets = [ref_bulk.evaluate(np.asarray(q)*2*np.pi/source.source.geometry.lattice_constant)*factor**2
                 for _, q in all_points]
    interface_rows, qrows, summaries = [], [], []
    for label, path in (('baseline', baseline), ('candidate', candidate/'joint_summary.json')):
        report = json.loads(path.read_bytes())
        profile = report['final_profile']
        if not report['completed'] or not profile['strictly_positive_LJ']:
            raise ValueError('completed positive-LJ training-selected candidates required')
        shape = decode_shape(report['final_shape_coordinates'], 1., True)
        c = np.asarray(profile['coefficients'])
        model = CoordinationScreenedInterface(shape, c, law='power')
        tight = CoordinationScreenedInterface(shape, c, law='power', tolerance=2e-13)
        rows = []
        for index, (q, reference) in enumerate(zip(plan['interface_states'], state_targets)):
            value, refined = model.evaluate(q), tight.evaluate(q)
            row = dict(model=label, state_index=index, a=q[0], x=q[1], y=q[2],
                energy_eV=value.energy, source_energy_eV=reference.energy,
                gradient=value.gradient.tolist(), source_gradient=reference.gradient.tolist(),
                hessian=value.hessian.tolist(), source_hessian=reference.hessian.tolist(),
                minimum_H_eigenvalue=float(np.linalg.eigvalsh(value.hessian).min()),
                source_minimum_H_eigenvalue=float(np.linalg.eigvalsh(reference.hessian).min()),
                energy_absolute_error_eV=abs(value.energy-reference.energy),
                gradient_error_norm=float(np.linalg.norm(value.gradient-reference.gradient)),
                relative_H_error=float(np.linalg.norm(value.hessian-reference.hessian)/np.linalg.norm(reference.hessian)),
                analytic_tolerance_H_change=float(np.max(abs(value.hessian-refined.hessian))))
            # Independently differentiate energy/gradient at two off-axis states.
            if index in (6, 9):
                q = np.asarray(q)
                differences = []
                for step in (2e-5, 1e-5):
                    high = [tight.evaluate(q+step*e) for e in np.eye(3)]
                    low = [tight.evaluate(q-step*e) for e in np.eye(3)]
                    g = np.array([(a.energy-b.energy)/(2*step) for a, b in zip(high, low)])
                    H = np.column_stack([(a.gradient-b.gradient)/(2*step) for a, b in zip(high, low)])
                    differences.append((g, H))
                g, H = [(4*differences[1][k]-differences[0][k])/3 for k in (0, 1)]
                row.update(gradient_fd_error=float(np.max(abs(g-refined.gradient))),
                           hessian_fd_error=float(np.max(abs(H-refined.hessian))))
            else:
                row.update(gradient_fd_error=None, hessian_fd_error=None)
            rows.append(row)
        interface_rows.extend(rows)
        coarse = CoordinationScreenedBulk(shape, radius=12., law='power')
        fine = CoordinationScreenedBulk(shape, radius=16., law='power')
        matrix_rows = []
        for (role, q), reference in zip(all_points, q_targets):
            high, tail = fine.evaluate(q)
            low, _ = coarse.evaluate(q)
            H = np.einsum('c,cij->ij', c, high)
            radius_difference = np.einsum('c,cij->ij', c, high-low)
            matrix_rows.append(dict(model=label, role=role, q=list(q),
                hessian=H.tolist(), source_hessian=reference.tolist(),
                relative_H_error=float(np.linalg.norm(H-reference)/np.linalg.norm(reference)),
                minimum_eigenvalue=float(np.linalg.eigvalsh(H).min()),
                tail_bound=float(tail@abs(c)),
                radius_change=float(np.linalg.norm(radius_difference)),
                radius_change_operator_norm=float(np.linalg.norm(radius_difference, 2))))
        qrows.extend(matrix_rows)
        summaries.append(dict(model=label, source_report_sha256=sha(path), shape=shape,
            coefficients=c, training_eta=profile['minimax_normalized_error'],
            maximum_new_interface_H_error=max(r['relative_H_error'] for r in rows),
            maximum_analytic_tolerance_H_change=max(r['analytic_tolerance_H_change'] for r in rows),
            maximum_gradient_fd_error=max(r['gradient_fd_error'] for r in rows if r['gradient_fd_error'] is not None),
            maximum_hessian_fd_error=max(r['hessian_fd_error'] for r in rows if r['hessian_fd_error'] is not None),
            maximum_previous_excluded_q_error=max(r['relative_H_error'] for r in matrix_rows if r['role']=='excluded'),
            maximum_new_excluded_q_error=max(r['relative_H_error'] for r in matrix_rows if r['role']=='new_excluded'),
            minimum_tested_q_eigenvalue=min(r['minimum_eigenvalue'] for r in matrix_rows),
            maximum_q_tail_bound=max(r['tail_bound'] for r in matrix_rows),
            maximum_q_radius_change=max(r['radius_change'] for r in matrix_rows),
            maximum_q_radius_change_operator_norm=max(r['radius_change_operator_norm'] for r in matrix_rows)))
        print('validated', label, summaries[-1]['maximum_new_interface_H_error'], flush=True)
    out.mkdir(parents=True)
    write_csv(out/'interface.csv', interface_rows)
    write_csv(out/'finite_q.csv', qrows)
    save_json(out/'interface.json', interface_rows)
    save_json(out/'finite_q.json', qrows)
    # The source loader also serves a straight-row study. Its row-energy and
    # elastic-tensor unit labels do not describe these interface/Bloch jets.
    provenance = {key: binding[key] for key in ('parameter_source', 'parameter_sha256',
        'source_url', 'source_nominal_lattice_angstrom', 'source_reference_only', 'kernel')}
    summary = dict(completed=True, source=provenance,
        units=dict(length_scale_m=LENGTH_M, interface_coordinates='physical distances / L0',
            interface_energy='eV per primitive interface cell',
            interface_gradient='derivative in reduced coordinates; eV per interface cell',
            interface_hessian='second derivative in reduced coordinates; eV per interface cell',
            bulk_hessian='harmonic stiffness in reduced displacements; eV per atom',
            wavevector='cubic reciprocal coordinates in units 2*pi/a_lat',
            reference_state='static ideal FCC; no finite-temperature or surface relaxation fit'),
        validation_plan_sha256=sha(candidate/'validation_plan.json'),
        observations_sha256=sha(candidate/'observations.json'),
        norms=dict(relative_H_error='Frobenius', tail_bound='operator 2-norm',
            radius_change='Frobenius', analytic_tolerance_H_change='maximum absolute entry',
            gradient_fd_error='maximum absolute component', hessian_fd_error='maximum absolute entry'),
        selection='candidate selected using training eta before computing excluded validation',
        comparisons=summaries, elapsed_seconds=time.perf_counter()-started,
        full_history_blindness_certified=False, all_q_stability_proved=False,
        material_accepted=False, physical_time_calibrated=False, production_changed=False)
    save_json(out/'summary.json', summary)
    print((out/'summary.json').read_text(encoding='utf-8'), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(__doc__)
    for name in ('candidate', 'baseline', 'out'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    run(args.candidate, args.baseline, args.out)
