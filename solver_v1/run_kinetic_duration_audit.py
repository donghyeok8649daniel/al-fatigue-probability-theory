"""Audit duration scaling and heating using the completed v31 MD campaign."""
import argparse
import csv
import json
from pathlib import Path
import time

import numpy as np

from .kinetic_duration_audit import (KB_EV_K, cycle_windows, sine_weights,
    conditional_precision_budget, thermal_work_summary)
from .run_low_frequency_forcing_v29 import sha
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json
from .work_phase_audit import duration_sensitivity

PARTS = (1, 2, 4, 5, 10, 19, 20, 38, 76)
THERMO_COLUMNS = ['temperature_K', 'potential_eV', 'kinetic_eV', 'total_eV', 'pressure_bar']


def _bound_records(study, previous, work_audit):
    """Validate a complete signed/null campaign before creating any output."""
    protocol = json.loads((study/'protocol.json').read_bytes())
    cases = json.loads((study/'cases.json').read_bytes())
    prior = json.loads((previous/'summary.json').read_bytes())
    audit = json.loads((work_audit/'summary.json').read_bytes())
    if (not prior['completed'] or prior['protocol'] != protocol
            or not audit['completed'] or audit['new_MD_run']
            or len(set(protocol['seeds'])) != len(protocol['seeds'])):
        raise ValueError('completed matching previous analyses required')
    expected = {(seed, axis, sign) for seed in protocol['seeds']
                for axis, sign in ((None, 0), (0, -1), (0, 1), (1, -1), (1, 1))}
    names = [case['name'] for case in cases]
    keys = [(case['seed'], case['axis'], case['sign']) for case in cases]
    hashes = {row['case']: row['trajectory_sha256'] for row in prior['checks']}
    if (len(set(names)) != len(names) or len(set(keys)) != len(keys)
            or set(keys) != expected or set(names) != set(hashes)
            or len(prior['checks']) != len(names) or audit['raw_hashes'] != hashes):
        raise ValueError('complete unique signed/null source binding required')
    records = []
    for case in cases:
        name = case['name']
        if Path(name).name != name or name in ('.', '..'):
            raise ValueError('simple source case names required')
        folder = study/name
        metadata = json.loads((folder/'summary.json').read_bytes())
        restart = study/f'seed{case["seed"]}_init'/'equilibrated.restart'
        force = None if case['axis'] is None else (
            case['sign']*protocol['force_scale_eV_A'][case['axis']])
        drive = None if force is None else dict(axis=case['axis'], force_eV_A=force,
                                               frequency_per_ps=protocol['frequency_per_ps'])
        if (not metadata['completed'] or metadata['ensemble'] != 'nve'
                or metadata['thermostat_damping_ps'] is not None
                or any(metadata[k] != protocol[k] for k in ('dt_ps', 'frame_ps', 'duration_ps'))
                or metadata['conjugate_drive'] != drive
                or (force is not None and (case['force_eV_A'] != force
                                          or not metadata['internal_step_work']))
                or metadata['restart_sha256'] != case['restart_sha256']
                or sha(restart) != case['restart_sha256']
                or sha(folder/'plane_coordinates.npz') != hashes[name]
                or metadata['thermo_columns'] != THERMO_COLUMNS
                or metadata['atom_count'] != metadata['repeats']*metadata['atoms_per_plane']):
            raise ValueError('completed unchanged NVE trajectory and restart binding required')
        records.append((case, metadata, folder/'plane_coordinates.npz'))
    return protocol, records, hashes, audit


def run(study, previous, work_audit, out, *, parts=PARTS, thermal_window_ps=100.):
    study, previous, work_audit, out = map(Path, (study, previous, work_audit, out))
    if out.exists():
        raise FileExistsError('fresh report directory required')
    if (not parts or len(set(parts)) != len(parts) or 1 not in parts
            or any(not isinstance(n, int) or n <= 0 for n in parts)):
        raise ValueError('unique positive integer partitions including the full record required')
    started = time.perf_counter()
    protocol, records, hashes, old_work = _bound_records(study, previous, work_audit)
    frequency = protocol['frequency_per_ps']
    duration = protocol['duration_ps']-protocol['exclude_ps']
    null_blocks, scaling, thermal_blocks, thermal, provenance = [], [], [], [], []
    null_temperatures, atom_counts, wall_rates = {}, set(), []
    zero_sum = 0.
    for case, metadata, raw in records:
        name, seed = case['name'], case['seed']
        atom_counts.add(metadata['atom_count'])
        wall_rates.append(metadata['elapsed_seconds']/metadata['duration_ps'])
        provenance.append(dict(case=name, seed=seed, axis=case['axis'], sign=case['sign'],
            trajectory_sha256=hashes[name], restart_sha256=case['restart_sha256'],
            metadata_sha256=sha(raw.parent/'summary.json'), atom_count=metadata['atom_count'],
            atoms_per_plane=metadata['atoms_per_plane'], ensemble=metadata['ensemble']))
        with np.load(raw, allow_pickle=False) as data:
            t, q, h = data['time_seconds']*1e12, data['coordinates_m']*1e10, data['thermo']
            w = None if case['axis'] is None else data['internal_step_work_eV']
            expected_t = np.arange(round(protocol['duration_ps']/protocol['frame_ps'])+1)*protocol['frame_ps']
            if (t.shape != expected_t.shape or not np.allclose(t, expected_t, rtol=0, atol=1e-8)
                    or q.shape != (len(t), metadata['repeats'], 3) or not np.all(np.isfinite(q))):
                raise ValueError('complete finite trajectory with protocol sampling required')
            first, last = cycle_windows(t, frequency, duration, protocol['exclude_ps'])[0]
            checks = thermal_work_summary(t[first:last+1], h[first:last+1],
                                          None if w is None else w[first:last+1])
            blocks = cycle_windows(t, frequency, thermal_window_ps, protocol['exclude_ps'])
            means = []
            for block, (i, j) in enumerate(blocks):
                values = thermal_work_summary(t[i:j+1], h[i:j+1], None if w is None else w[i:j+1])
                means.append(values['mean_temperature_K'])
                thermal_blocks.append(dict(case=name, seed=seed, axis=case['axis'], sign=case['sign'],
                    block=block, start_ps=float(t[i]), end_ps=float(t[j]), **values))
            thermal.append(dict(case=name, seed=seed, axis=case['axis'], sign=case['sign'],
                analyzed_duration_ps=duration, first_window_mean_K=means[0], last_window_mean_K=means[-1],
                last_minus_first_window_K=means[-1]-means[0], **checks))
            if case['axis'] is not None:
                continue
            null_temperatures[seed] = checks['mean_temperature_K']
            zero_sum = max(zero_sum, float(np.max(abs(q[first:last+1].sum(axis=1)))))
            for axis in (0, 1):
                force = protocol['force_scale_eV_A'][axis]
                for count in parts:
                    width = duration/count
                    loss_blocks = []
                    for block, (i, j) in enumerate(cycle_windows(t, frequency, width, protocol['exclude_ps'])):
                        losses = sine_weights(t[i:j+1], frequency, force) @ q[i:j+1, :, axis]
                        loss_blocks.append(losses)
                        for plane, loss in enumerate(losses):
                            null_blocks.append(dict(seed=seed, axis=axis, parts=count, block=block,
                                plane=plane, start_ps=float(t[i]), end_ps=float(t[j]),
                                duration_ps=width, signed_loss_A2_eV=float(loss)))
                    losses = np.asarray(loss_blocks)
                    rms = float(np.sqrt(np.mean(losses**2)))
                    # An RMS over equivalent planes is descriptive, not SE/sqrt(planes).
                    scaling.append(dict(seed=seed, axis=axis, parts=count, duration_ps=width,
                        temporal_blocks=len(losses), planes=losses.shape[1],
                        rms_loss_A2_eV=rms, max_abs_loss_A2_eV=float(np.max(abs(losses))),
                        rms_times_sqrt_duration=rms*np.sqrt(width),
                        max_abs_mean_over_planes_A2_eV=float(np.max(abs(losses.mean(axis=1)))),
                        independent_sample_count=None, confidence_interval=False))
    with (previous/'fdt.csv').open(newline='', encoding='utf8') as stream:
        fdt = [row for row in csv.DictReader(stream) if int(row['parts']) == 1]
    expected_fdt = {(seed, axis) for seed in protocol['seeds'] for axis in (0, 1)}
    if {(int(r['seed']), int(r['axis'])) for r in fdt} != expected_fdt:
        raise ValueError('full-record FDT predictions for every seed and axis required')
    precision = []
    for row in fdt:
        seed, axis = int(row['seed']), int(row['axis'])
        loss = -float(row['imag_A2_eV'])
        full = next(r for r in scaling if (r['seed'], r['axis'], r['parts']) == (seed, axis, 1))
        budget = conditional_precision_budget(protocol['force_scale_eV_A'][axis], frequency,
            loss, null_temperatures[seed], duration)
        precision.append(dict(seed=seed, axis=axis, nw=float(row['nw']),
            half_band_cycles_ps=float(row['half_band_cycles_ps']), independent_fdt_loss_A2_eV=loss,
            null_temperature_K=null_temperatures[seed], observed_full_rms_A2_eV=full['rms_loss_A2_eV'],
            observed_rms_to_conditional_std=full['rms_loss_A2_eV']/budget['predicted_loss_std_A2_eV'],
            **budget))
    # Keep the prior MAX-envelope criterion; do not replace it by a standard deviation.
    heating = []
    if len(atom_counts) != 1:
        raise ValueError('one verified crystal size required for shared heating sensitivity')
    atoms = next(iter(atom_counts))
    heat_capacity = (3*atoms-3)*KB_EV_K  # Harmonic fixed-COM assumption, not a measured Al Cv.
    if sorted(row['axis'] for row in old_work['planning']) != [0, 1]:
        raise ValueError('one prior envelope plan per axis required')
    for row in old_work['planning']:
        axis = row['axis']
        duration_required = row['conditional_record_duration_ps']
        predicted = min(r['independent_fdt_loss_A2_eV'] for r in precision if r['axis'] == axis)
        null = max(r['max_abs_loss_A2_eV'] for r in scaling if r['axis'] == axis and r['parts'] == 1)
        expected = duration_sensitivity(duration, null, predicted, row['desired_ratio'])
        if not np.allclose([row['independent_fdt_loss_A2_eV'], row['observed_imaginary_null_A2_eV'],
                            duration_required], [predicted, null, expected], rtol=1e-8, atol=0):
            raise ValueError('previous planning must match source FDT and full-record null')
        budget = conditional_precision_budget(protocol['force_scale_eV_A'][axis], frequency,
            row['independent_fdt_loss_A2_eV'], max(null_temperatures.values()), duration_required)
        heating.append(dict(axis=axis, force_eV_A=protocol['force_scale_eV_A'][axis],
            predicted_loss_A2_eV=row['independent_fdt_loss_A2_eV'],
            observed_imaginary_null_A2_eV=row['observed_imaginary_null_A2_eV'],
            desired_envelope_ratio=row['desired_ratio'],
            conditional_record_duration_ps=duration_required,
            expected_work_eV=budget['expected_work_eV'],
            assumed_heat_capacity_eV_K=heat_capacity,
            conditional_temperature_rise_K=budget['expected_work_eV']/heat_capacity,
            estimated_single_record_wall_hours=row['estimated_single_record_wall_hours'],
            heat_capacity_measured=False, stationarity_certified=False,
            production_clock_calibrated=False))
    durations = {r['axis']: r['conditional_record_duration_ps'] for r in heating}
    campaign_ps = len(protocol['seeds'])*(max(durations.values())+2*sum(durations.values()))
    ranges = []
    for seed in protocol['seeds']:
        for axis in (0, 1):
            chosen = [r for r in scaling if r['seed'] == seed and r['axis'] == axis]
            products = [r['rms_times_sqrt_duration'] for r in chosen]
            ranges.append(dict(seed=seed, axis=axis, min_rms_sqrt_duration=min(products),
                max_rms_sqrt_duration=max(products), max_to_min=max(products)/min(products),
                duration_min_ps=min(r['duration_ps'] for r in chosen), duration_max_ps=duration))
    summary = dict(completed=True, new_MD_run=False, records_reanalyzed=len(records),
        protocol=protocol, partitions=list(parts), thermal_window_ps=thermal_window_ps,
        raw_hashes=hashes, source_records=provenance,
        analysis_input_hashes={label: sha(path) for label, path in (
            ('protocol', study/'protocol.json'), ('cases', study/'cases.json'),
            ('previous_summary', previous/'summary.json'), ('previous_fdt', previous/'fdt.csv'),
            ('work_audit_summary', work_audit/'summary.json'))},
        max_null_plane_sum_A=zero_sum, noise_scaling_ranges=ranges,
        full_record_null=[r for r in scaling if r['parts'] == 1], thermal=thermal,
        heating_sensitivity=heating,
        conditional_signed_null_campaign_duration_ns=campaign_ps/1000,
        conditional_sum_record_wall_hours=campaign_ps*float(np.median(wall_rates))/3600,
        campaign_scheduled=False, distinct_initializations=len(protocol['seeds']),
        statistical_independence_certified=False, stationarity_certified=False,
        confidence_interval=False, heat_capacity_measured=False, production_clock_calibrated=False,
        conclusion='Duration and heating sensitivity only; no new resolved drag or production clock.',
        elapsed_seconds=time.perf_counter()-started)
    out.mkdir(parents=True)
    for name, rows in (('null_blocks', null_blocks), ('noise_scaling', scaling),
                       ('thermal_blocks', thermal_blocks), ('thermal', thermal),
                       ('precision', precision), ('heating_sensitivity', heating)):
        write_csv(out/(name+'.csv'), rows)
    save_json(out/'summary.json', summary)
    print(json.dumps(dict(noise_scaling_ranges=ranges, heating_sensitivity=heating,
        conditional_signed_null_campaign_duration_ns=campaign_ps/1000,
        elapsed_seconds=summary['elapsed_seconds'], production_clock_calibrated=False), indent=2))
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(__doc__)
    for name in ('study', 'previous', 'work-audit', 'out'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    run(args.study, args.previous, args.work_audit, args.out)
