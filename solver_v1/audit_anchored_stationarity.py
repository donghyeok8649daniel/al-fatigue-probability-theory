"""Check local first-order stationarity from a saved analytic design/Jacobian."""
import argparse
import json
from pathlib import Path

import numpy as np

from .anchored_shape_calibration import ExactAnchorCoordinates, tangent_descent
from .interface_tangent_calibration import NONNEGATIVE
from .run_low_frequency_forcing_v29 import sha
from .run_vector_registry_audit import save_json


def run(study, out):
    study, out = Path(study), Path(out)
    if out.exists():
        raise FileExistsError('fresh local stationarity report required')
    report = json.loads((study/'joint_summary.json').read_bytes())
    if not report['completed']:
        raise ValueError('completed shape study required')
    data = json.loads((study/'final_matrix.json').read_bytes())
    definition = json.loads((study/'anchored_definition.json').read_bytes())
    profile = report['final_profile']
    M, dM = np.asarray(data['matrix']), np.asarray(data['derivatives'])
    D = np.asarray(definition['coefficient_scales'])
    chart = ExactAnchorCoordinates(M, data['target'], data['scales'], profile['exact_rows'], D)
    point = np.r_[report['final_shape_coordinates'],
                  (np.asarray(profile['coefficients'])/D)[chart.free], profile['minimax_normalized_error']]
    value = chart.evaluate(M, dM, point[6:-1], point[-1], profile['tested_rows'], NONNEGATIVE)
    lower = np.r_[report['lower'], np.full(len(chart.free), -np.inf), 0.]
    upper = np.r_[report['upper'], np.full(len(chart.free)+1, np.inf)]
    # The second audit removes only this run's local trust box. It searches
    # infinitesimal directions at the SAME point, not finite infeasible shapes.
    # The original isolated-neighbour energy-decay condition remains included.
    x = point[:6]
    decay = 2*np.exp(x[4])+min(x[5], 0.)*np.exp(x[0])
    if decay <= 0:
        raise ArithmeticError('isolated-neighbour energy does not decay at the audited point')
    decay_gradient = np.zeros(len(point))
    decay_gradient[0] = min(x[5], 0.)*np.exp(x[0])
    decay_gradient[4] = 2*np.exp(x[4])
    decay_gradient[5] = np.exp(x[0]) if x[5] < 0 else 0.
    constraints = np.r_[value['constraints'], decay]
    jacobian = np.vstack([value['jacobian'], decay_gradient])
    global_lower = np.r_[np.log([.7, 2., 2., 1e-5, 2.]), -1., lower[6:]]
    global_upper = np.r_[np.log([8., 24., 24., 1e6, 24.]), 1., upper[6:]]
    audits = []
    for label, lo, hi in (('executed_box', lower, upper),
                          ('original_domain_tangent', global_lower, global_upper)):
        for tol in (2e-6, 1e-5):
            audits.append(dict(bounds=label, **tangent_descent(
                constraints, jacobian, point, lo, hi, tol)))
    out.mkdir(parents=True)
    result = dict(completed=True, study_sha256=sha(study/'joint_summary.json'),
        matrix_sha256=sha(study/'final_matrix.json'), audits=audits,
        exact_residual=value['exact_residual'], pivot_condition=value['pivot_condition'],
        numerical_coordinate_scales=D, free_coefficients=chart.free,
        isolated_neighbour_decay_exponent=float(decay),
        criterion='active feasible-cone LP; finite matrix derivatives and declared search bounds',
        global_shape_optimum_certified=False, material_accepted=False, physical_time_calibrated=False)
    save_json(out/'summary.json', result)
    print(json.dumps(dict(minimum_linearized_deta=[r['minimum_linearized_deta'] for r in audits],
                         exact_residual=result['exact_residual']), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--study', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    run(args.study, args.out)
