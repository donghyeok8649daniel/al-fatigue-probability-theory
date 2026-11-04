"""Retain boundary cells using time-node-specific positive liftings."""
from __future__ import annotations

import json
from pathlib import Path
import time

import numpy as np

from coarse_memory import aggregation, cdf_bracket, exact_memory, markov_audit


def main():
    started = time.perf_counter()
    root = Path(__file__).resolve().parent
    data = np.load(root / 'SOURCE_OPERATOR_PAIR.npz')
    steps, valid, paths = data['steps'], data['valid'], data['fine_paths']
    initial = [data['base'], data['mixed']]
    well = np.tile(np.floor(data['s'] + 0.5).astype(int), data['a'].size)
    cases = []
    checks = []

    def check(name, condition, value=None):
        checks.append({'name': name, 'passed': bool(condition), 'value': value})
        if not condition:
            raise AssertionError(name)

    for bands in [1, 3, 9]:
        normal = np.arange(data['a'].size) * bands // data['a'].size
        coordinates = np.c_[np.repeat(normal, data['s'].size), well]
        _, labels = np.unique(coordinates, axis=0, return_inverse=True)
        labels = np.r_[labels, labels.max() + 1]
        liftings = []
        empty_groups = []
        for node_valid in valid:
            weights = node_valid.astype(float)
            empty = []
            for group in np.unique(labels):
                cells = labels == group
                if weights[cells].sum() == 0:
                    empty.append(int(group))
                    weights[cells] = 1.0
            p, r = aggregation(labels, weights)
            liftings.append(r)
            empty_groups.append(empty)
        check(f'right_inverses_{bands}', max(float(np.max(abs(p @ r - np.eye(p.shape[0])))) for r in liftings) < 1e-14)
        group_cases = []
        for case, (p0, truth) in enumerate(zip(initial, paths)):
            audit = markov_audit(steps, truth, p, liftings, valid)
            lower, upper = cdf_bracket(steps, p, labels, p @ p0, valid)
            memory, _, _ = exact_memory(steps, p, liftings, p0)
            memory_error = float(np.max(abs(memory - truth @ p.T)))
            check(f'moving_memory_{bands}_{case}', memory_error < 2e-12, memory_error)
            check(f'moving_bounds_{bands}_{case}', lower - 2e-13 <= truth[-1, -1] <= upper + 2e-13)
            check(f'moving_tv_{bands}_{case}', all(row['cdf_error'] <= row['posterior_tv_bound'] + 2e-13 for row in audit))
            for n, empty in enumerate(empty_groups):
                check_mass = p @ truth[n]
                if empty:
                    check(f'empty_mass_{bands}_{case}_{n}', float(np.max(abs(check_mass[empty]))) == 0.0)
            group_cases.append({'lower': lower, 'upper': upper, 'memory_max_error': memory_error, 'final_audit': audit[-1]})
        cases.append({'normal_bands': bands, 'macro_states': p.shape[0], 'cases': group_cases,
                      'empty_groups_by_time': empty_groups})

    # Regression of the constant-lifting memory formula after the extension.
    labels = np.r_[well - well.min(), well.max() - well.min() + 1]
    p, r = aggregation(labels, np.r_[np.logical_and.reduce(valid[:, :-1]).astype(float), 1.0])
    constant_memory, _, _ = exact_memory(steps, p, r, initial[1])
    check('constant_lifting_regression', float(np.max(abs(constant_memory - paths[1] @ p.T))) < 2e-12)
    check('nested_bracket_widths', all(
        cases[j + 1]['cases'][i]['upper'] - cases[j + 1]['cases'][i]['lower']
        <= cases[j]['cases'][i]['upper'] - cases[j]['cases'][i]['lower'] + 2e-13
        for j in range(2) for i in range(2)))
    receipt = {'status': 'PASS_MOVING_LIFTING_BOUNDARY_RETENTION', 'dateKST': '2026-11-04',
               'wall_seconds': time.perf_counter() - started, 'source_rebuild_or_old_suite': False,
               'input': 'SOURCE_OPERATOR_PAIR.npz', 'levels': cases,
               'checks': checks, 'passed_checks': len(checks),
               'empty_group_convention': 'Uniform dummy lifting inside its own block with zero actual mass; bounds [0,1]; original masks keep incoming mass zero',
               'independent_review': 'NOT_PERFORMED', 'outward_roundoff_certificate': False}
    (root / 'MOVING_LIFTING_RESULT.json').write_text(json.dumps(receipt, indent=2), encoding='utf-8')
    print(json.dumps({'status': receipt['status'], 'checks': len(checks),
                      'mixed_final': [r['cases'][1]['final_audit']['cdf_markov'] for r in cases]}))


if __name__ == '__main__':
    main()
