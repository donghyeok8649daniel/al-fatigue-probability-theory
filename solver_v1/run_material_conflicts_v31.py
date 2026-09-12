"""Identify competing measured targets at the completed radial shape.

Every profile is a diagnostic ablation of existing targets, never a promoted
calibration. No signs, exact anchors, source values or material terms change.
The LP dual gives a fixed-shape lower bound only.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from .core_interface_compatibility import minimax_compatibility
from .coordination_screening import CoordinationScreenedBulk
from .interface_tangent_calibration import TangentCalibrationProblem, impose_tangents, NONNEGATIVE
from .periodic_plane_covariance import SourceEAMBlochHessian
from .run_current_material_core import ROOT
from .run_finite_q_compatibility_v31 import points
from .run_source_core_reference import load_source_material
from .run_vector_registry_audit import save_json
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_low_frequency_forcing_v29 import sha
from .vector_material_calibration import LENGTH_M, MaterialObservation


def run(study, out):
    study, out = Path(study), Path(out)
    if out.exists():
        raise FileExistsError('fresh conflict report required')
    report = json.loads((study / 'summary.json').read_bytes())
    if not report['completed'] or not report['best']:
        raise ValueError('completed same-family search required')
    shape = np.asarray(report['best']['shape'])
    problem = TangentCalibrationProblem(ROOT / 'results/fcc111_active_interface/coordination_screening_v19/joint_refinement')
    obs = list(impose_tangents(problem.observations, problem.source_tangents))
    base = problem.matrix(tuple(shape))
    interface = [i for i, o in enumerate(obs) if o.role != 'exact' and o.units == 'eV/L0^2']
    source, _, binding = load_source_material()
    operator = SourceEAMBlochHessian(source.source)
    bulk = CoordinationScreenedBulk(shape, radius=16., law='power')
    extra, group = [], {'K': [], 'L': []}
    for k, (role, q) in enumerate(points()):
        if role != 'fit':
            continue
        H = operator.evaluate(np.asarray(q)*2*np.pi/source.source.geometry.lattice_constant)*(LENGTH_M/1e-10)**2
        columns, _ = bulk.evaluate(q)
        key = 'K' if q[2] == 0 else 'L'
        for i, j in ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)):
            weight = 1 if i == j else np.sqrt(2)
            group[key].append(len(obs))
            obs.append(MaterialObservation(f'v31_{key}{k}_H{i}{j}', float(weight*H[i, j]),
                                           float(.05*np.linalg.norm(H)), 'eV/L0^2', 'fit'))
            extra.append(weight*columns[:, i, j])
    matrix = np.vstack([base, extra])
    qrows = group['K'] + group['L']
    profiles = [('K_only', group['K']), ('L_only', group['L']), ('finite_q_only', qrows),
                ('interface_only', interface), ('interface_plus_K', interface+group['K']),
                ('interface_plus_L', interface+group['L']), ('joint', interface+qrows)]
    out.mkdir(parents=True)
    reports, dual_rows, residual_rows = [], [], []
    for label, rows in profiles:
        fit = minimax_compatibility(matrix, obs, rows, nonnegative=NONNEGATIVE)
        if not fit['completed']:
            raise ArithmeticError('fixed-shape LP did not complete')
        reports.append(dict(profile=label, **fit))
        prediction = np.asarray(fit['predictions'])
        for i, o in enumerate(obs):
            residual_rows.append(dict(profile=label, target=o.name, prediction=prediction[i],
                reference=o.target, scale=o.scale, residual=(prediction[i]-o.target)/o.scale,
                selected_in_profile=i in rows, exact=o.role == 'exact'))
        for j, i in enumerate(rows):
            dual_rows.append(dict(profile=label, target=obs[i].name,
                upper_error_dual=fit['upper_error_duals'][j], lower_error_dual=fit['lower_error_duals'][j],
                absolute_dual_weight=abs(fit['upper_error_duals'][j])+abs(fit['lower_error_duals'][j])))
        print(label, fit['minimax_normalized_error'], flush=True)
    joint = reports[-1]
    if abs(joint['minimax_normalized_error']-report['best']['value']) > 2e-5:
        raise ArithmeticError('joint objective replay failed')
    save_json(out / 'profiles.json', reports)
    write_csv(out / 'dual_target_weights.csv', dual_rows)
    write_csv(out / 'profile_residuals.csv', residual_rows)
    save_json(out / 'summary.json', dict(completed=True, search_sha256=sha(study/'summary.json'),
        source=binding, shape=shape, profiles=[dict(profile=r['profile'], eta=r['minimax_normalized_error'],
            dual_lower_bound=r['dual_lower_bound'], duality_gap=r['duality_gap'],
            positive_LJ=r['strictly_positive_LJ']) for r in reports],
        interpretation='fixed-shape diagnostic ablations, not fits eligible for adoption',
        full_radial_family_impossibility_proved=False, material_accepted=False, production_changed=False))


if __name__ == '__main__':
    p = argparse.ArgumentParser(__doc__)
    p.add_argument('--study', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    run(a.study, a.out)
