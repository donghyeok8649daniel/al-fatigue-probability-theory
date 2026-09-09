"""Deterministic one-extra-range nested fit; no strength or kinetic target."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize

from .range_resolved_material import RangeObservationCache, radial_identifiability
from .run_low_stress_cyclic_diagnostic import FIT, write_csv
from .run_vector_registry_audit import save_json
from .run_vector_material_calibration import source_and_targets
from .vector_material_calibration import fit_coefficients
from .yield_elastic_metric import cubic_metric_problem, fit_exact_bulk


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/fcc111_active_interface/range_core_v12/material'


def execute(out, local_evaluations):
    started = time.perf_counter()
    source, observations, states = source_and_targets()
    cache = RangeObservationCache(observations)
    previous = json.loads((ROOT/'results/fcc111_active_interface/yield_bridge_v11/material_metric_refined/calibration.json').read_text())['best']
    profiles = []; optimizers = []
    def evaluate(decays, stage, metric=None):
        matrix, obs = cubic_metric_problem(cache.matrix(decays), observations)
        for mode in ([metric] if metric else ['cubic_5pct', 'bulk_exact']):
            fit = fit_coefficients(matrix, obs) if mode == 'cubic_5pct' else fit_exact_bulk(matrix, obs)
            row = dict(scalar_decay=float(decays[0]), odd_decay=float(decays[1]),
                       quadrupole_decay=float(decays[2]), stage=stage, metric=mode, **fit)
            if fit['coefficients'] is not None:
                row['cubic_GPa'] = matrix[2:5]@fit['coefficients']
                hold = [i for i, o in enumerate(obs) if o.role == 'heldout']
                res = (matrix@fit['coefficients']-np.array([o.target for o in obs]))/np.array([o.scale for o in obs])
                row['heldout_normalized_rms'] = float(np.sqrt(np.mean(res[hold]**2)))
            profiles.append(row)
            print(f"{len(profiles)} {stage} {mode} {decays}: loss={row.get('squared_loss')} feasible={row['admissible']}", flush=True)
        save_json(out/'progress.json', dict(completed=False, profiles=profiles,
                                           elapsed_seconds=time.perf_counter()-started))
        return row
    # Predeclared one-dimensional ablations around both actual v11 starts, plus
    # a small coarse3D grid. These are fixed starts, never randomized samples.
    starts = []
    for old in previous.values():
        d, k = old['scalar_decay'], old['angular_decay']
        starts.extend((d, k, q) for q in (2., 3.2, k, 6.4, 9., 12.))
    starts.extend((d, o, e) for d in (1.3, 2.7, 5.3) for o in (3.2, 6.4, 10.)
                  for e in (2., 5., 10.))
    starts = list(dict.fromkeys(starts))
    for point in starts:
        evaluate(point, 'predeclared_nested_grid')
    for metric in ('cubic_5pct', 'bulk_exact'):
        feasible = [r for r in profiles if r['metric'] == metric and r['admissible'] and r['strictly_positive_LJ_resolved']]
        if not feasible:
            optimizers.append(dict(metric=metric, no_feasible_start=True)); continue
        start = min(feasible, key=lambda r: r['squared_loss'])
        p = [start[k] for k in ('scalar_decay', 'odd_decay', 'quadrupole_decay')]
        if local_evaluations:
            def objective(log_parameters):
                row = evaluate(np.exp(log_parameters), 'local_log_Powell', metric)
                return row.get('squared_loss', np.inf)
            run = minimize(objective, np.log(p), method='Powell',
                           bounds=np.log([[.7, 8.], [2., 12.], [2., 12.]]),
                           options=dict(maxfev=local_evaluations, ftol=1e-6, xtol=1e-4))
            optimizers.append(dict(metric=metric, start=p, success=bool(run.success),
                                   message=str(run.message), nfev=run.nfev,
                                   returned_objective=float(run.fun)))
    bests = {}; residuals = []; identities = {}
    for metric in ('cubic_5pct', 'bulk_exact'):
        rows = [r for r in profiles if r['metric'] == metric and r['admissible'] and r['strictly_positive_LJ_resolved']]
        if not rows:
            continue
        best = min(rows, key=lambda r: r['squared_loss']); bests[metric] = best
        parameters = [best[k] for k in ('scalar_decay', 'odd_decay', 'quadrupole_decay')]
        matrix, obs = cubic_metric_problem(cache.matrix(parameters), observations)
        for value, o in zip(matrix@best['coefficients'], obs):
            residuals.append(dict(candidate=metric, observable=o.name, target=o.target,
                                  prediction=value, scale=o.scale, units=o.units, role=o.role,
                                  normalized_residual=(value-o.target)/o.scale))
        identities[metric] = radial_identifiability(parameters, best['coefficients'], observations)
    save_json(out/'calibration.json', dict(completed=True, best=bests, profiles=profiles,
        predeclared_starts=starts, local_optimizations=optimizers,
        elapsed_seconds=time.perf_counter()-started, material_accepted=False,
        extra_shape_parameters=1, physical_yield_validated=False, production_changed=False))
    save_json(out/'target_definition.json', dict(source_sha256=source.reference.sha256,
        historical_parameter_sha256=hashlib.sha256(FIT.read_bytes()).hexdigest(),
        observations=[asdict(o) for o in observations], states=states, temperature_K=0,
        force_and_cohesion_exact=True, no_heldout_in_loss=True,
        strengths_kinetics_source_lengths_in_loss=False))
    write_csv(out/'residuals.csv', residuals)
    save_json(out/'identifiability.json', identities)
    save_json(out/'progress.json', dict(completed=True, evaluated_profiles=len(profiles)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUT)
    parser.add_argument('--local-evaluations', type=int, default=180)
    args = parser.parse_args()
    execute(args.output, args.local_evaluations)


if __name__ == '__main__':
    main()
