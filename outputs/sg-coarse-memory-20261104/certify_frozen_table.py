"""Exact-rational Bellman envelopes for captured binary64 SG tables."""
from __future__ import annotations

from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import time

import numpy as np


def decimal_bound(value, rounding):
    with localcontext() as context:
        context.prec = 60
        context.rounding = rounding
        return str(Decimal(value.numerator) / Decimal(value.denominator))


def encode(value):
    return {'numerator': str(value.numerator), 'denominator': str(value.denominator),
            'decimal_lower': decimal_bound(value, ROUND_FLOOR),
            'decimal_upper': decimal_bound(value, ROUND_CEILING)}


def apply(columns, vector):
    out = [F(0)] * len(vector)
    for i, entries in enumerate(columns):
        for j, coefficient in entries:
            out[j] += coefficient * vector[i]
    return out


def main():
    # These integers are generated from this bounded, local binary64 fixture.
    # Exact products can exceed Python's default decimal serialization limit.
    sys.set_int_max_str_digits(50000)
    started = time.perf_counter()
    root = Path(__file__).resolve().parent
    data = np.load(root / 'SOURCE_OPERATOR_PAIR.npz')
    steps, valid = data['steps'], data['valid']
    size = steps.shape[1]
    captured = []
    stochastic = []
    perturbations = []
    norms = []
    for step in steps:
        original_columns = []
        corrected_columns = []
        maximum_change = F(0)
        maximum_norm = F(0)
        for i in range(size):
            entries = {int(j): F.from_float(float(step[j, i])) for j in np.flatnonzero(step[:, i])}
            if any(value < 0 for value in entries.values()):
                raise ValueError('Captured table contains a negative coefficient')
            original_columns.append(list(entries.items()))
            column_sum = sum(entries.values(), F(0))
            maximum_norm = max(maximum_norm, column_sum)
            anchor = max(entries, key=entries.get)
            adjustment = F(1) - column_sum
            entries[anchor] += adjustment
            if entries[anchor] < 0:
                raise ValueError('Column residual construction lost positivity')
            corrected_columns.append(list(entries.items()))
            maximum_change = max(maximum_change, abs(adjustment))
        captured.append(original_columns)
        stochastic.append(corrected_columns)
        perturbations.append(maximum_change)
        norms.append(maximum_norm)

    boundary = json.loads((root / 'BOUNDARY_MIXTURE_RESULT.json').read_text())['boundary_cell_index']
    base = [F.from_float(float(x)) for x in data['base']]
    mass = sum(base, F(0))
    mixed = [F(3, 5) * x for x in base]
    mixed[boundary] += F(2, 5) * sum(base[:-1], F(0))
    mixed[-1] = base[-1]
    well = np.tile(np.floor(data['s'] + 0.5).astype(int), data['a'].size)
    coarse_labels = np.r_[well - well.min(), well.max() - well.min() + 1]
    for group in np.unique(coarse_labels):
        assert sum((base[i] for i in np.flatnonzero(coarse_labels == group)), F(0)) == sum((mixed[i] for i in np.flatnonzero(coarse_labels == group)), F(0))

    normal = np.repeat(np.arange(data['a'].size), data['s'].size)
    _, labels = np.unique(np.c_[normal, well], axis=0, return_inverse=True)
    labels = np.r_[labels, labels.max() + 1]
    groups = int(labels.max() + 1)
    lower = [F(0)] * groups
    upper = [F(0)] * groups
    lower[-1] = upper[-1] = F(1)
    for n in range(len(steps) - 1, -1, -1):
        low_cells = [sum((value * lower[labels[j]] for j, value in column), F(0)) for column in stochastic[n]]
        high_cells = [sum((value * upper[labels[j]] for j, value in column), F(0)) for column in stochastic[n]]
        next_low, next_high = [], []
        for group in range(groups):
            cells = np.flatnonzero((labels == group) & valid[n])
            next_low.append(min((low_cells[i] for i in cells), default=F(0)))
            next_high.append(max((high_cells[i] for i in cells), default=F(1)))
        lower, upper = next_low, next_high

    cases = []
    captured_cdfs = []
    for initial in (base, mixed):
        original = initial.copy()
        corrected = initial.copy()
        perturbation_bound = F(0)
        for n in range(len(steps)):
            perturbation_bound += perturbations[n] * sum(original, F(0))
            original = apply(captured[n], original)
            corrected = apply(stochastic[n], corrected)
        aggregate = [sum((initial[i] for i in np.flatnonzero(labels == g)), F(0)) for g in range(groups)]
        low = sum((x * y for x, y in zip(lower, aggregate)), F(0))
        high = sum((x * y for x, y in zip(upper, aggregate)), F(0))
        assert low <= corrected[-1] <= high
        assert abs(original[-1] - corrected[-1]) <= perturbation_bound
        assert sum(corrected, F(0)) == sum(initial, F(0))
        captured_cdfs.append(original[-1])
        cases.append({'corrected_table_absorption': encode(corrected[-1]),
                      'captured_table_exact_product_absorption': encode(original[-1]),
                      'corrected_lower': encode(low), 'corrected_upper': encode(high),
                      'captured_product_lower': encode(max(F(0), low - perturbation_bound)),
                      'captured_product_upper': encode(high + perturbation_bound),
                      'table_perturbation_bound': encode(perturbation_bound)})

    cdf_gap = captured_cdfs[1] - captured_cdfs[0]
    assert cdf_gap > F(39, 100)
    receipt = {
        'status': 'EXACT_RATIONAL_FROZEN_TABLE_CDF_ENVELOPES_VERIFIED',
        'dateKST': '2026-11-04', 'wall_seconds': time.perf_counter() - started,
        'target': 'Exact-real products of captured binary64 transition entries and raw initial masses',
        'column_construction': 'Change only the largest positive entry by the exact column residual; preserve every column support and record perturbation',
        'initial_measure_renormalized': False,
        'exact_initial_raw_mass': encode(mass), 'equal_initial_macro_laws_exact': True,
        'maximum_column_l1_perturbation': encode(max(perturbations)),
        'macro_states': groups, 'fine_states': size,
        'cases': cases, 'exact_cdf_gap': encode(cdf_gap),
        'checks': {'exact_equal_macro': True, 'positive_stochastic_reference': True,
                   'both_exact_cdf_enclosures': True, 'both_perturbation_bounds': True,
                   'raw_mass_preserved': True, 'cdf_gap_gt_0p39': True},
        'original_energy_table_generation_roundoff': 'OUTSIDE_THIS_CERTIFICATE',
        'floating_solver_execution_roundoff': 'OUTSIDE_THIS_CERTIFICATE',
        'continuum_event_and_physical_calibration': 'OPEN',
        'independent_review': 'NOT_PERFORMED',
    }
    (root / 'EXACT_TABLE_CERTIFICATE.json').write_text(json.dumps(receipt, indent=2), encoding='utf-8')
    print(json.dumps({'status': receipt['status'], 'seconds': receipt['wall_seconds'],
                      'mixed_lower': cases[1]['captured_product_lower']['decimal_lower'],
                      'mixed_upper': cases[1]['captured_product_upper']['decimal_upper'],
                      'maximum_column_change': receipt['maximum_column_l1_perturbation']['decimal_upper']}), flush=True)


if __name__ == '__main__':
    main()
