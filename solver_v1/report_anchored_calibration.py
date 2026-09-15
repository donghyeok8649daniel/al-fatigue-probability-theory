"""Select completed material candidates by training error before validation."""
import argparse
import json
from pathlib import Path

import numpy as np

from .run_low_frequency_forcing_v29 import sha
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json


def run(fits, baseline, out):
    fits, baseline, out = list(map(Path, fits)), Path(baseline), Path(out)
    if not fits or out.exists():
        raise ValueError('completed fits and a fresh training selection directory required')
    observations = json.loads((fits[0]/'observations.json').read_bytes())
    if any(sha(path/'observations.json') != sha(fits[0]/'observations.json')
           for path in [baseline.parent, *fits]):
        raise ValueError('candidate targets, scales and row order must be identical')
    reports = [('baseline', baseline)]+[(path.name, path/'joint_summary.json') for path in fits]
    if len({label for label, _ in reports}) != len(reports):
        raise ValueError('distinct model labels required')
    rows, residual_rows = [], []
    for label, path in reports:
        report = json.loads(path.read_bytes())
        profile = report['final_profile']
        if not report['completed'] or not profile['completed'] or not profile['strictly_positive_LJ']:
            raise ValueError('completed certified positive-LJ profiles required')
        predictions = np.asarray(profile['predictions'])
        if predictions.shape != (len(observations),):
            raise ValueError('matching training prediction rows required')
        x = np.asarray(report['final_shape_coordinates'])
        at_lower = np.flatnonzero(abs(x-np.asarray(report['lower'])) < 2e-6).tolist() if 'lower' in report else []
        at_upper = np.flatnonzero(abs(x-np.asarray(report['upper'])) < 2e-6).tolist() if 'upper' in report else []
        signs = 0
        for i, (obs, prediction) in enumerate(zip(observations, predictions)):
            diagonal = any(obs['name'].endswith('_'+key) for key in ('Haa', 'Hxx', 'Hyy', 'H00', 'H11', 'H22'))
            opposite = bool(diagonal and obs['role'] != 'exact' and prediction*obs['target'] < 0
                            and min(abs(prediction), abs(obs['target'])) > obs['scale'])
            signs += int(opposite)
            residual_rows.append(dict(model=label, index=i, name=obs['name'], role=obs['role'],
                units=obs['units'], state=obs['state'], target=obs['target'], scale=obs['scale'],
                prediction=float(prediction), residual=float((prediction-obs['target'])/obs['scale']),
                diagonal_curvature_opposite_beyond_discrepancy_scale=opposite))
        rows.append(dict(model=label, report=path.as_posix(), report_sha256=sha(path),
            training_eta=profile['minimax_normalized_error'], positive_LJ=True,
            optimizer_success=report['optimizer_success'], message=report['message'],
            iterations=report['iterations'], shape=np.r_[np.exp(x[:5]), x[5]].tolist(), coefficients=profile['coefficients'],
            exact_residual=profile['equality_residual'], lp_duality_gap=profile['duality_gap'],
            diagonal_curvature_opposite_beyond_discrepancy_scale=signs,
            active_lower_coordinates=at_lower, active_upper_coordinates=at_upper))
    selected = min(rows, key=lambda row: row['training_eta'])
    baseline_eta = rows[0]['training_eta']
    out.mkdir(parents=True)
    write_csv(out/'training_residuals.csv', residual_rows)
    summary = dict(completed=True, selection_rule='minimum certified positive-LJ training eta; baseline retained',
        observations_sha256=sha(fits[0]/'observations.json'), comparisons=rows,
        selected_model=selected['model'], selected_report=selected['report'],
        relative_training_improvement=1-selected['training_eta']/baseline_eta,
        excluded_validation_used_for_selection=False, material_accepted=False,
        physical_time_calibrated=False, production_changed=False)
    save_json(out/'training_selection.json', summary)
    print(json.dumps(dict(selected_model=selected['model'], selected_report=selected['report'],
                          relative_training_improvement=summary['relative_training_improvement']), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--fits', type=Path, nargs='+', required=True)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    run(args.fits, args.baseline, args.out)
