"""v28: execute a constrained-scope mobility calibration on completed records.

Uses v25 finite-band spectra and v26/v27 forced responses, without rerunning
or modifying MD/production. No random sampling, fitted mass, or physical clock.
All low bands are retained, including failed/unresolved tests.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from .impedance_mobility import (
    MOBILITY_SI_PER_A2_EV_PS, scalar_impedance, inverse_disk_bounds,
    scalar_low_band_drag, fit_constant_drag,
)
from .run_vector_registry_audit import save_json
from .run_low_stress_cyclic_diagnostic import write_csv

ROOT = Path(__file__).resolve().parents[1]


def read_csv(path):
    with path.open(encoding='utf-8', newline='') as stream:
        return list(csv.DictReader(stream))


def run(output):
    output = Path(output)
    if output.exists(): raise FileExistsError('fresh results directory required')
    inputs = []

    def load_json(path):
        raw = path.read_bytes()
        inputs.append(dict(path=path.relative_to(ROOT).as_posix(), sha256=hashlib.sha256(raw).hexdigest()))
        return json.loads(raw)

    def load_csv(path):
        inputs.append(dict(path=path.relative_to(ROOT).as_posix(), sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        return read_csv(path)

    r27 = ROOT/'results/phase_resolution_v27'
    summary = load_json(r27/'summary.json')
    if not summary['completed']: raise ValueError('completed v27 required')
    protocol = summary['protocol']
    comparisons = load_csv(r27/'phase_comparison.csv')
    static = {int(r['axis']): float(r['static_covariance_chi_A2_eV']) for r in comparisons}
    pairs = load_csv(r27/'paired_responses.csv')
    nulls = load_csv(r27/'null_windows.csv')
    v26 = load_json(ROOT/'results/collective_forcing_v26/summary.json')
    if not v26['completed'] or v26['protocol']['source_projection_sha256'] != protocol['source_sha256']:
        raise ValueError('matched completed v26 required')
    for r in v26['paired']:
        # v27 supplies empirical null windows only for 180/360 ps. Legacy
        # 90 ps quarter fits remain in v26, not assigned an invented floor.
        if r['parts'] not in (1,2): continue
        pairs.append(dict(r, frequency_per_ps=.05, dt_ps=.0025,
            sign_disagreement_A2_eV=r['signed_response_disagreement_A2_eV']))

    # Exact inverse-map sensitivity, including controls which fail positivity.
    driven = []
    for r in pairs:
        axis, parts = int(r['axis']), int(r['parts'])
        f, fraction = float(r['frequency_per_ps']), float(r['fraction'])
        chi = complex(float(r['real_A2_eV']), float(r['imag_A2_eV']))
        F = protocol['force_scale_eV_A'][axis]*fraction
        if f == .05:
            F = next(c['force_eV_A'] for c in v26['protocol']['cases']
                if c['axis']==axis and c['fraction']==fraction and c['sign']==1)
        null = [x for x in nulls if int(x['axis']) == axis
            and float(x['frequency_per_ps']) == f and float(x['duration_ps']) == 360/parts]
        if not null: raise ValueError('matched unforced null windows missing')
        radius = max(np.hypot(float(x['real_A']), float(x['imag_A'])) for x in null)/F
        # These empirical differences mix numerical, nonlinear and thermal error.
        # Add them conservatively; never pretend independent quadrature errors.
        full = [x for x in pairs if int(x['axis'])==axis and float(x['frequency_per_ps'])==f
            and int(x['parts'])==1]
        as_complex = lambda x: complex(float(x['real_A2_eV']),float(x['imag_A2_eV']))
        base_response = next(x for x in full if float(x['dt_ps'])==.0025 and float(x['fraction'])==1)
        half_response = next(x for x in full if float(x['dt_ps'])==.0025 and float(x['fraction'])==.5)
        fine = [x for x in full if float(x['dt_ps'])==.00125 and float(x['fraction'])==1]
        dt_change = abs(as_complex(fine[0])-as_complex(base_response)) if fine else 0.
        amp_change = abs(as_complex(half_response)-as_complex(base_response))
        same = [x for x in pairs if int(x['axis']) == axis
            and float(x['frequency_per_ps']) == f and float(x['fraction']) == fraction
            and float(x['dt_ps']) == float(r['dt_ps'])]
        center = next(x for x in same if int(x['parts']) == 1)
        central_chi = complex(float(center['real_A2_eV']), float(center['imag_A2_eV']))
        block_change = max(abs(complex(float(x['real_A2_eV']),float(x['imag_A2_eV']))-central_chi)
            for x in same)
        sign_change = float(r['sign_disagreement_A2_eV'])
        for label, error in [('null_only', radius),
                ('null_plus_observed_controls', radius+dt_change+amp_change+block_change+sign_change)]:
            driven.append(dict(axis=axis, frequency_per_ps=f, dt_ps=float(r['dt_ps']),
                fraction=fraction, parts=parts, block=int(r['block']), envelope=label,
                real_chi_A2_eV=chi.real, imag_chi_A2_eV=chi.imag,
                radius_A2_eV=error, null_radius_A2_eV=radius,
                dt_difference_A2_eV=dt_change, amplitude_difference_A2_eV=amp_change,
                block_difference_A2_eV=block_change, sign_difference_A2_eV=sign_change,
                matched_dt_control_available=bool(fine),matched_amplitude_control_available=True,
                **scalar_impedance(chi, f, static_chi=static[axis]),
                **{k:v for k,v in inverse_disk_bounds(chi,error,f).items()
                   if k != 'production_clock_calibrated'}))

    # Three pre-existing broad low bands fit one scalar drag, not a matrix entry.
    # Lower bands are post-fit validation, not new/blind experimental data.
    bands = []
    base = ROOT/'results/modal_kinetic_calibration_v25'
    for dt, prefix in [(dt,p) for dt in (.0025,.005) for p in (1000,2000,5000)]:
        directory = f'long_record_{"dt5_" if dt == .005 else ""}prefix_{prefix}'
        data = load_json(base/directory/'summary.json')
        if (not data['completed'] or data['source_dt_ps'] != dt
                or data['source_repeats'] != 6 or data['source_ensemble'] != 'nve'):
            raise ValueError('matched completed N6 low-frequency controls required')
        if dt == .0025 and data['source_projection_sha256'] != protocol['source_sha256']:
            raise ValueError('v25/v27 source hash mismatch')
        for r in data['records']:
            if not r['valid_finite_band']: continue
            C = np.asarray(r['covariance_m2']); K = np.asarray(r['band']['integral_proxy_m2_seconds'])
            for axis in (0,1):
                gamma = scalar_low_band_drag(C[axis,axis], K[axis,axis], r['temperature_K'])
                bands.append(dict(axis=axis, dt_ps=dt, prefix_ps=prefix,
                    parts=r['parts'],block=r['block'],nw=r['nw'],
                    lower_per_ps=r['lower_cycle_THz'],upper_per_ps=r['upper_cycle_THz'],
                    bins=r['band']['bins'], temperature_K=r['temperature_K'],
                    drag_J_s_m2=gamma,scalar_mobility_m2_J_s=1/gamma,
                    coupled_matrix_diagonal_mobility_m2_J_s=r['response']['mobility_m2_per_J_second'][axis][axis],
                    fitted_band=r['lower_cycle_THz'] >= .02,
                    zero_frequency_limit_certified=False))

    fits, residuals, fit_controls, predictions = [], [], [], []
    for axis in (0,1):
        rows = [r for r in bands if r['axis'] == axis]
        training = [r for r in rows if r['prefix_ps']==5000 and r['dt_ps']==.0025
            and r['nw']==3 and r['parts']==1 and r['fitted_band']]
        scales = {}
        for r in training:
            control = [x['drag_J_s_m2'] for x in rows if x['lower_per_ps']==r['lower_per_ps']]
            scales[r['lower_per_ps']] = max(control)-min(control)
        y = [r['drag_J_s_m2'] for r in training]
        s = [scales[r['lower_per_ps']] for r in training]
        fit = fit_constant_drag(y,s)
        if not fit['admissible']: raise ValueError('nonpositive low-band fit; not repaired')
        gamma = fit['drag']
        fits.append(dict(axis=axis,**fit,mobility_m2_per_J_s=1/gamma,
            interpretation='finite-band scalar low-frequency hypothesis, not certified M(0)'))
        for r in rows:
            scale = scales.get(r['lower_per_ps'])
            residuals.append(dict(**r, fitted_drag_J_s_m2=gamma,
                prediction_relative_error=(gamma-r['drag_J_s_m2'])/r['drag_J_s_m2'],
                normalized_residual=(gamma-r['drag_J_s_m2'])/scale if scale else None,
                sensitivity_scale_J_s_m2=scale))
        keys = sorted(set((r['dt_ps'],r['prefix_ps'],r['nw'],r['parts'],r['block']) for r in rows))
        for dt,prefix,nw,parts,block in keys:
            rr=[r for r in rows if (r['dt_ps'],r['prefix_ps'],r['nw'],r['parts'],r['block'])
                == (dt,prefix,nw,parts,block) and r['fitted_band']]
            if len(rr) != 3: continue
            cf=fit_constant_drag([r['drag_J_s_m2'] for r in rr],[scales[r['lower_per_ps']] for r in rr])
            fit_controls.append(dict(axis=axis,dt_ps=dt,prefix_ps=prefix,nw=nw,parts=parts,block=block,
                drag_J_s_m2=cf['drag'],mobility_m2_J_s=1/cf['drag'] if cf['admissible'] else None,
                objective=cf['objective']))
        for r in driven:
            if r['axis']!=axis or r['envelope']!='null_only' or r['parts']!=1: continue
            omega=2*np.pi*r['frequency_per_ps']
            predicted=1/(1/static[axis]+1j*omega*gamma*MOBILITY_SI_PER_A2_EV_PS)
            measured=complex(r['real_chi_A2_eV'],r['imag_chi_A2_eV'])
            predictions.append(dict(axis=axis,frequency_per_ps=r['frequency_per_ps'],
                dt_ps=r['dt_ps'],fraction=r['fraction'],
                predicted_real_A2_eV=predicted.real,predicted_imag_A2_eV=predicted.imag,
                predicted_phase_deg=np.degrees(np.angle(predicted)),
                measured_real_A2_eV=measured.real,measured_imag_A2_eV=measured.imag,
                complex_error_A2_eV=abs(predicted-measured),
                observed_null_radius_A2_eV=r['null_radius_A2_eV'],
                error_to_observed_null=abs(predicted-measured)/r['null_radius_A2_eV']))

    output.mkdir(parents=True)
    for name, records in [('forced_impedance',driven),('low_band_residuals',residuals),
            ('fit_control_sensitivity',fit_controls),('overdamped_predictions',predictions)]:
        write_csv(output/(name+'.csv'), records)
    save_json(output/'calibration_candidate.json',dict(completed=True,
        status='research_only_finite_band_candidate', fits=fits,
        fitted_band_edges_per_ps=[[.02,.04],[.04,.08],[.08,.16]],
        residual_scale='full observed control span per band; not statistical uncertainty',
        source_inputs=inputs,conjugate_coordinate='periodic gap of 144-atom plane means',
        finite_frequency_scalar_inverse_is_not_matrix_mobility=True,
        zero_frequency_limit_certified=False, material_kinetics_calibrated=False,
        local_coordinate_transfer_validated=False,production_clock_calibrated=False,
        t0_seconds=None))
    print(json.dumps(dict(fits=[dict(axis=r['axis'],M=r['mobility_m2_per_J_s'],objective=r['objective'])
        for r in fits],forced_rows=len(driven),band_rows=len(bands),control_fits=len(fit_controls)),indent=2))


if __name__ == '__main__':
    p=argparse.ArgumentParser(__doc__)
    p.add_argument('--out',type=Path,required=True)
    run(p.parse_args().out)
