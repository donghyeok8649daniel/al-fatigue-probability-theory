"""Joint fixed-range STF4 capacity diagnostic; NOT a selected material law.

Rebuild all columns from the analytic energy and independently replay each
saved single-column prediction before combining the five existing ranges.
No inspected/excluded state is represented as a new blind validation point.
"""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from .coordination_screening import CoordinationScreenedBulk
from .rank_four_environment import RankFourEnvironment, rank_four_bulk_column
from .core_interface_compatibility import minimax_compatibility
from .interface_tangent_calibration import TangentCalibrationProblem, NONNEGATIVE
from .run_current_material_core import ROOT
from .run_finite_q_compatibility_v31 import points
from .run_radial_channels_v32 import conditional_svd
from .run_low_frequency_forcing_v29 import sha
from .run_vector_registry_audit import save_json
from .vector_material_calibration import MaterialObservation
from .yield_elastic_metric import cubic_to_mode_matrix


def run(studies, out):
    studies, out = Path(studies), Path(out)
    if out.exists():
        raise FileExistsError('fresh capacity output required')
    names = ('fixed_shape', 'decay_3.0', 'decay_5.0', 'decay_8.0', 'decay_11.0')
    definition = json.loads((studies/names[0]/'definition.json').read_bytes())
    raw_obs = json.loads((studies/names[0]/'observations.json').read_bytes())
    obs = [MaterialObservation(**o) for o in raw_obs]
    shape = definition['shape']
    decays = [shape[0], 3., 5., 8., 11.]
    inputs = {}
    for name in names:
        saved_definition = json.loads((studies/name/'definition.json').read_bytes())
        saved_summary = json.loads((studies/name/'summary.json').read_bytes())
        if saved_definition['shape'] != shape or not saved_summary['completed']:
            raise ValueError('same parent shape and completed studies required')
        if json.loads((studies/name/'observations.json').read_bytes()) != raw_obs:
            raise ValueError('observations or units changed between studies')
        inputs[name] = sha(studies/name/'profiles.json')
    out.mkdir(parents=True)
    save_json(out/'definition.json', dict(parent_profile_hashes=inputs, shape=shape,
        decays=decays, purpose='joint capacity only, not minimal model selection',
        nonnegative_control=True, signed_control=True,
        excluded_q_used_in_loss=False, new_blind_validation=False,
        production_changed=False))
    started = time.perf_counter()
    problem = TangentCalibrationProblem(ROOT/'results/fcc111_active_interface/coordination_screening_v19/joint_refinement')
    base = problem.matrix(tuple(shape))
    bulk = CoordinationScreenedBulk(shape, radius=16., law='power')
    components = ((0,0),(1,1),(2,2),(0,1),(0,2),(1,2))
    weights = np.array([1.,1.,1.,np.sqrt(2),np.sqrt(2),np.sqrt(2)])
    qbase = []
    for _, q in points():
        columns, _ = bulk.evaluate(q)
        qbase.append(weights[:,None]*np.array([columns[:,i,j] for i,j in components]))
    matrix = np.vstack([base, *qbase])
    extra_columns, replay, tails = [], [], []
    rows = json.loads((studies/names[0]/'profiles.json').read_bytes())[0]['tested_rows']
    for name, decay in zip(names, decays):
        model = RankFourEnvironment(decay)
        column = np.zeros(len(obs))
        column[:5] = model.bulk_column()
        column[2:5] = np.linalg.solve(cubic_to_mode_matrix(), column[2:5])
        for i, o in enumerate(problem.raw_observations):
            if o.bulk_index is None:
                column[i] = np.asarray(o.jet_weights) @ model.interface_jet(o.state)[0]
        for j, (_, q) in enumerate(points()):
            H, bound = rank_four_bulk_column(model, bulk, q)
            column[len(base)+6*j:len(base)+6*(j+1)] = [w*H[i,k] for w,(i,k) in zip(weights,components)]
            tails.append(dict(decay=decay, q=q, coefficient_H_tail=bound))
        for profile in json.loads((studies/name/'profiles.json').read_bytes()):
            design = matrix if len(profile['coefficients']) == 10 else np.column_stack([matrix,column])
            error = float(np.max(abs(design @ profile['coefficients'] - profile['predictions'])))
            if error > 2e-8:
                raise ArithmeticError(f'independent single-column replay failed: {name}, {error}')
            replay.append(dict(study=name, profile=profile['label'], maximum_error=error))
        extra_columns.append(column)
    design = np.column_stack([matrix, *extra_columns])
    profiles = []
    for signed in (False, True):
        signs = NONNEGATIVE + (() if signed else tuple(range(10,15)))
        fit = minimax_compatibility(design, obs, rows, nonnegative=signs)
        if not fit['completed']:
            raise ArithmeticError('joint capacity LP failed')
        excluded = []
        for j, (role, _) in enumerate(points()):
            if role == 'excluded':
                ix = np.arange(len(base)+6*j, len(base)+6*(j+1))
                target = np.array([obs[i].target for i in ix])
                excluded.append(float(np.linalg.norm(fit['predictions'][ix]-target)/np.linalg.norm(target)))
        profiles.append(dict(label='signed_capacity_DIAGNOSTIC' if signed else 'nonnegative_capacity_DIAGNOSTIC',
            **fit, identifiability=conditional_svd(design,obs,rows), excluded_q_errors=excluded))
    save_json(out/'profiles.json', profiles)
    save_json(out/'column_replay.json', replay)
    save_json(out/'unit_coefficient_tails.json', tails)
    save_json(out/'summary.json', dict(completed=True, elapsed_seconds=time.perf_counter()-started,
        maximum_single_column_replay_error=max(r['maximum_error'] for r in replay),
        profiles=[dict(label=p['label'], eta=p['minimax_normalized_error'],
            excluded_q_error=max(p['excluded_q_errors']), positive_LJ=p['strictly_positive_LJ'],
            duality_gap=p['duality_gap'], conditional_rank=p['identifiability']['conditional_rank']) for p in profiles],
        interpretation='joint fixed-shape capacity; no all-q stability or material acceptance',
        material_accepted=False, physical_time_calibrated=False, actual_yield_validated=False,
        production_changed=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--studies', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    run(args.studies,args.out)
