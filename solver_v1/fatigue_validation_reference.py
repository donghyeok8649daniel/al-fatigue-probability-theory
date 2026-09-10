"""Condition/endpoint-aware source fatigue validation, never an S-N fit."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np


REFERENCE = Path(__file__).with_name('data')/'aluminum_fatigue_validation_v18.json'


def load_fatigue_reference():
    data = json.loads(REFERENCE.read_bytes())
    rows = data['records']
    if (len({row['id'] for row in rows}) != len(rows)
            or data['used_in_kinetic_or_energy_fit']):
        raise ValueError('unique independent validation references required')
    for row in rows:
        if row['control'] not in ('axial_stress', 'total_axial_strain'):
            raise ValueError('unknown source control mode')
        if (row['ratio'] != -1 or row['frequency_Hz'] <= 0 or row['full_range'] <= 0
                or row['failure_cycles_min'] <= 0
                or row['failure_cycles_max'] < row['failure_cycles_min']):
            raise ValueError('invalid source range/time/cycle record')
        expected_unit = 'MPa' if row['control'] == 'axial_stress' else '1'
        if row['full_range_unit'] != expected_unit:
            raise ValueError('source control mode and units disagree')
    return data


def normalized_fatigue_records():
    """Retain pooled source ranges, not invented per-specimen sample pairs."""
    data = load_fatigue_reference()
    rows = []
    for row in data['records']:
        rate = 2*row['frequency_Hz']*row['full_range']
        # Include floating arithmetic but do not reinterpret printed rounding
        # as an experimental error bar or enlarge it to force agreement.
        roundoff = np.finfo(float).eps*max(rate, row['reported_rate'])*8
        if abs(rate-row['reported_rate']) > row['rate_rounding_halfwidth']+roundoff:
            raise ValueError('source rates inconsistent with declared triangular range convention')
        rows.append(dict(row, amplitude=row['full_range']/2,
            minimum=-row['full_range']/2, maximum=row['full_range']/2,
            derived_triangular_rate=rate, rate_rounding_residual=rate-row['reported_rate'],
            laboratory_failure_time_min_seconds=row['failure_cycles_min']/row['frequency_Hz'],
            laboratory_failure_time_max_seconds=row['failure_cycles_max']/row['frequency_Hz'],
            endpoint='final_specimen_fracture', initiation_cycles=None,
            local_absorbed_probability=None, source_doi=data['doi']))
    return rows


def fatigue_comparison_gate(*, prediction_endpoint, physical_clock_calibrated,
                            specimen_mapping_validated, control_protocol_matched,
                            microstructure_matched, probability_resolution_certified):
    """Refuse to turn a local/noncalibrated run into a specimen lifetime fit.

    A line-level kinetic calibration alone cannot satisfy the PDE-clock flag.
    This checks declared prerequisites, not the truth of external evidence.
    """
    reasons = []
    if prediction_endpoint != 'final_specimen_fracture':
        reasons.append('local opening / AE growth / final fracture are distinct endpoints')
    for passed, reason in ((physical_clock_calibrated, 'matched PDE physical clock unavailable'),
            (specimen_mapping_validated, 'specimen-scale mapping unvalidated'),
            (control_protocol_matched, 'loading control/waveform/frequency not matched'),
            (microstructure_matched, 'grain/orientation/defect state not matched'),
            (probability_resolution_certified, 'probability resolution uncertified')):
        if not passed:
            reasons.append(reason)
    return dict(comparable=not reasons, reasons=reasons,
                lifetime_error_available=False,  # needs an actual matched prediction
                empirical_lifetime_fit_performed=False)
