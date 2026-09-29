"""Audit stored Si Hessians and replay older SW force histories. No new atom MD."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh, solve
from scipy.signal import correlate

from solver_v1.silicon_charge_dynamics import sg_chain
from solver_v1.silicon_quantum_diagnostic import quantum_harmonic_difference
from solver_v1.silicon_thermal_validity import KB, oscillator_statistics, validate_fixed_cartesian


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def main(root, output):
    if output.exists():
        raise ValueError('use a fresh output directory')
    output.mkdir(parents=True)
    inputs = {}

    def source(relative):
        path = root / relative
        inputs[relative] = sha(path)
        return path

    geometry = np.load(source('results/silicon_initiation_v11/intact_prism_360_linesearch_v2/geometry.npz'))
    free = geometry['free']
    reference = geometry['reference']
    # An explicitly labelled internal extension observable, NOT a crack coordinate.
    lower = free & (reference[:, 2] < np.median(reference[free, 2]))
    upper = free & ~lower
    if not lower.any() or not upper.any():
        raise ValueError('both free-atom halves required')
    observation = np.zeros_like(reference)
    observation[upper, 2] = 1 / upper.sum()
    observation[lower, 2] = -1 / lower.sum()
    published = json.loads(source('results/silicon_initiation_v13/local_harmonic/summary.json').read_text())
    rows, spectra, arrays, records = [], {}, {}, []
    model_hashes = set()
    for name in ('loading8', 'return8', 'loading10'):
        prefix = 'results/silicon_initiation_v13/dense_' + name
        protocol = json.loads(source(prefix + '/protocol.json').read_text())
        summary = json.loads(source(prefix + '/summary.json').read_text())
        if not summary['complete'] or protocol['ensemble'] != 'displacement':
            raise ValueError('complete fixed-grip source required')
        model_hashes.add(protocol['model_sha256'])
        d = np.load(source(prefix + '/raw_hessian.npz'))
        if np.any(d['numbers'] != 14):
            raise ValueError('equal Si reference masses only')
        h, checks = validate_fixed_cartesian(d['hessian'], d['basis'], free)
        lam, vectors = eigh(h)
        spectrum_numpy = np.linalg.eigvalsh(h)
        checks['independent_spectrum_max_abs'] = float(np.max(abs(lam - spectrum_numpy)))
        if checks['independent_spectrum_max_abs'] > 1e-10:
            raise ValueError('independent eigensolver mismatch')
        inv = solve(h, np.eye(len(h)), assume_a='pos')
        checks['inverse_identity_error'] = float(np.max(abs(h @ inv - np.eye(len(h)))))
        lift = d['basis'] @ vectors
        c = d['basis'].T @ observation.ravel()
        weights = (vectors.T @ c) ** 2
        checks['observation_inverse_identity_error'] = float(abs(weights @ (1 / lam) - c @ inv @ c))
        shift = -inv @ d['gradient']
        checks['quadratic_center_shift_max_atom_A'] = float(np.linalg.norm((d['basis'] @ shift).reshape(-1, 3), axis=1).max())
        spectra[name] = lam
        records.append(dict(state=name, dimension=len(h), free_atoms=int(free.sum()),
                            fixed_atoms=int((~free).sum()), stored_energy_eV=float(d['energy']),
                            maximum_gradient_eV_A=float(np.max(abs(d['gradient']))), checks=checks))
        for temperature in (50., 150., 300., 600.):
            stats = oscillator_statistics(lam, temperature)
            vc, vq = stats['classical_variance_A2'], stats['quantum_variance_A2']
            atom_c = ((lift * lift) @ vc).reshape(-1, 3).sum(axis=1)
            atom_q = ((lift * lift) @ vq).reshape(-1, 3).sum(axis=1)
            freq = stats['omega_rad_ps'] / (2 * np.pi)
            ratio = stats['variance_ratio']
            row = dict(state=name, temperature_K=temperature,
                       frequency_min_THz=float(freq.min()), frequency_max_THz=float(freq.max()),
                       period_min_ps=float(1 / freq.max()), period_max_ps=float(1 / freq.min()),
                       modes_hbar_omega_above_kBT=int(np.sum(stats['quantum_energy_eV'] > KB * temperature)),
                       mode_variance_ratio_min=float(ratio.min()), mode_variance_ratio_max=float(ratio.max()),
                       modes_variance_correction_above_10pct=int(np.sum(ratio > 1.1)),
                       classical_free_atom_RMS_A=float(np.sqrt(atom_c[free].mean())),
                       quantum_free_atom_RMS_A=float(np.sqrt(atom_q[free].mean())),
                       classical_maximum_atom_RMS_A=float(np.sqrt(atom_c.max())),
                       quantum_maximum_atom_RMS_A=float(np.sqrt(atom_q.max())),
                       fixed_atom_maximum_variance_A2=float(max(atom_c[~free].max(), atom_q[~free].max())),
                       classical_extension_variance_A2=float(weights @ vc),
                       quantum_extension_variance_A2=float(weights @ vq),
                       classical_potential_energy_mean_eV=float(.5 * lam @ vc),
                       equipartition_identity_error_eV=float(abs(.5 * lam @ vc - len(h) * KB * temperature / 2)),
                       classical_vibrational_free_energy_eV=float(stats['classical_free_energy_eV'].sum()),
                       quantum_vibrational_free_energy_eV=float(stats['quantum_free_energy_eV'].sum()))
            if name in ('loading8', 'return8'):
                old = next(v for v in published['temperature_diagnostics'] if v['temperature_K'] == temperature) if temperature in (50., 150., 300.) else None
                if old is not None:
                    error = abs(row['classical_free_atom_RMS_A'] - old[name + '_Gaussian_free_atom_RMS_A'])
                    if error > 1e-10:
                        raise ValueError('previous classical result did not reproduce')
                    row['previous_classical_RMS_replay_error_A'] = error
            rows.append(row)
            if temperature == 300:
                arrays[name + '_curvatures_eV_A2'] = lam
                arrays[name + '_frequencies_THz'] = freq
                arrays[name + '_classical_atom_variance_A2'] = atom_c
                arrays[name + '_quantum_atom_variance_A2'] = atom_q
                time_ps = np.linspace(0, 2, 4001)
                # Analytic, undamped harmonic correlation, not sampled real dynamics.
                normalized = np.cos(np.outer(time_ps, stats['omega_rad_ps'])) @ (weights * vc) / (weights @ vc)
                arrays[name + '_extension_normalized_harmonic_correlation'] = normalized
                arrays['harmonic_time_ps'] = time_ps
                row['undamped_harmonic_extension_correlation_min_0_2ps'] = float(normalized.min())
                row['uncoupled_diagonal_RMS_A'] = float(np.sqrt(KB * temperature * np.sum(1 / np.diag(h)) / free.sum()))
        print(name, 'harmonic audit complete', flush=True)
    if len(model_hashes) != 1:
        raise ValueError('mixed force fields')

    differences = []
    for temperature in (50., 150., 300., 600.):
        a, b = [next(v for v in rows if v['state'] == name and v['temperature_K'] == temperature) for name in ('loading8', 'return8')]
        dq = b['quantum_vibrational_free_energy_eV'] - a['quantum_vibrational_free_energy_eV']
        dc = b['classical_vibrational_free_energy_eV'] - a['classical_vibrational_free_energy_eV']
        replay = quantum_harmonic_difference(spectra['return8'], spectra['loading8'], temperature)
        if abs(dq - replay) > 1e-10:
            raise ValueError('quantum free-energy independent formula mismatch')
        differences.append(dict(temperature_K=temperature, quantum_vibrational_difference_eV=dq,
                                classical_vibrational_difference_eV=dc, quantum_minus_classical_difference_eV=dq - dc,
                                formula_replay_error_eV=abs(dq - replay), is_transition_barrier=False))

    # Actual SG generator identity, on a declared synthetic coupled-coordinate energy.
    grid = np.linspace(-1, 1, 101)
    f = .12 * grid**2 + .04 * grid**4
    mobility = 1 + .2 * (grid[:-1] + grid[1:]) / 2
    kt = KB * 300
    base = sg_chain(f, thermal_energy=kt, spacing=grid[1]-grid[0], face_mobility=mobility)
    normalization = []
    for n in (2, 216, 360):
        scaled = sg_chain(f/n, thermal_energy=kt/n, spacing=grid[1]-grid[0], face_mobility=n*mobility)
        err = float(np.max(abs((base-scaled).toarray())))
        bad = sg_chain(f/n, thermal_energy=kt, spacing=grid[1]-grid[0], face_mobility=n*mobility)
        gibbs = np.exp(-f/kt); gibbs /= gibbs.sum()
        wrong = np.exp(-f/(n*kt)); wrong /= wrong.sum()
        normalization.append(dict(N=n, maximum_generator_error=err,
                                  relative_generator_error=err / abs(base.toarray()).max(),
                                  gibbs_stationarity_error=float(np.max(abs(base@gibbs))),
                                  wrong_thermal_generator_difference=float(np.max(abs((base-bad).toarray()))),
                                  wrong_thermal_equilibrium_variance_ratio=float((wrong@grid**2)/(gibbs@grid**2))))
        if err / abs(base.toarray()).max() > 1e-12:
            raise ValueError('energy normalization changed generator')

    # Replay all four initial/300 K histories in the previous SW control group.
    prefix = 'results/silicon_thermal_v4/thermal_dynamics/'
    old_analysis = json.loads(source(prefix+'analysis.json').read_text())
    old = next(v for v in old_analysis['rows'] if v['reference']=='initial' and v['temperature_K']==300)
    meta = json.loads(source(prefix+'summary.json').read_text())
    runs = sorted([v for v in meta['records'] if v['reference']=='initial' and v['temperature_K']==300], key=lambda v:v['seed'])
    replay_rows, integrals = [], []
    center = old['mean_from_independent_preparation_eV_A']
    for r in runs:
        d = np.load(source(prefix+r['label']+'.npz'))
        dt = float(d['time_ps'][1]-d['time_ps'][0])
        centered = d['observations'][:,0] - center
        n = len(centered); maxlag = min(round(10/dt), n//4)
        corr = correlate(centered, centered, mode='full', method='fft')[n-1:n+maxlag] / np.arange(n,n-maxlag-1,-1) / kt
        direct_error = max(abs(corr[k] - np.dot(centered[:n-k], centered[k:])/(n-k)/kt) for k in (0,1,10,100,maxlag))
        integral = np.r_[0., np.cumsum((corr[1:]+corr[:-1])*dt/2)]
        kinetic = float(2*d['kinetic_energy_eV'].mean()/(r['bath_dimension']*KB))
        replay_rows.append(dict(label=r['label'], original_duration_ps=r['duration_ps'],
                                reflections=r['reflection_count'], kinetic_temperature_K=kinetic,
                                max_direct_correlation_error=direct_error, new_MD=False))
        integrals.append(integral)
    integrals = np.asarray(integrals)
    windows = []
    for w in old['finite_window_integrals']:
        vals = integrals[:,round(w['lag_ps']/dt)]
        error = float(np.max(abs(vals - w['replica_values'])))
        if error > 1e-10:
            raise ValueError('old SW correlation did not reproduce')
        windows.append(dict(lag_ps=w['lag_ps'], mean=float(vals.mean()), minimum=float(vals.min()), maximum=float(vals.max()), replay_error=error))
    result = dict(complete=True, scope='stored MP-0b3 fixed-grip Hessians; separate historical SW dynamics replay',
                  states=records, temperatures=rows, vibrational_differences=differences,
                  energy_normalization_synthetic_checks=normalization,
                  historical_SW=dict(records=replay_rows, finite_integrals_eV_ps_A2=windows,
                                     center_from_stored_independent_preparation=center,
                                     exact_nonlinear_memory=False, transfers_to_current_model=False),
                  mass_amu=28.0855, model_sha256=next(iter(model_hashes)),
                  extension_observable='mean z of upper free reference half minus mean z of lower free half; not crack coordinate',
                  thresholds_10pct_and_hbaromega_kBT='diagnostic bins, not acceptance criteria',
                  excluded_force_ensemble='649-dimensional grip coordinate has no independently specified loading-device inertia',
                  new_potential_evaluations=0, new_MD=0, new_DFT=0,
                  measured_current_thermal_distribution=False, measured_current_time_correlation=False,
                  material_approved=False, overdamped_closure_validated=False, physical_clock=None,
                  limitations=['local harmonic approximation', 'finite-temperature anharmonicity and basin mixing unresolved',
                               'reference mass and empirical potential quantum sensitivity only',
                               'harmonic correlation has no damping and is not measured dynamics',
                               'historical SW geometry/potential differ from current intact MP-0b3 specimen'])
    save(output/'summary.json', result)
    np.savez_compressed(output/'spectra_and_covariance.npz', **arrays)
    fields = sorted(set().union(*(r.keys() for r in rows)))
    with (output/'temperature_diagnostics.csv').open('w', newline='', encoding='utf-8') as stream:
        writer=csv.DictWriter(stream, fieldnames=fields); writer.writeheader(); writer.writerows(rows)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4), layout='constrained')
    for name in spectra:
        current = [r for r in rows if r['state']==name]
        t = [r['temperature_K'] for r in current]
        axes[0].plot(t, [r['classical_free_atom_RMS_A'] for r in current], '--', label=name+' classical')
        axes[0].plot(t, [r['quantum_free_atom_RMS_A'] for r in current], '-', label=name+' quantum')
        r = next(r for r in current if r['temperature_K']==300)
        axes[1].plot(arrays[name+'_frequencies_THz'], oscillator_statistics(spectra[name],300)['variance_ratio'], label=name)
    axes[0].set(xlabel='Assumed temperature (K)', ylabel='Harmonic free-atom RMS (Angstrom)')
    axes[1].set(xlabel='Mode frequency (THz)', ylabel='Quantum / classical variance at 300 K')
    axes[2].plot([w['lag_ps'] for w in windows], [w['mean'] for w in windows], 'o-', color='black', label='4-history mean')
    axes[2].fill_between([w['lag_ps'] for w in windows], [w['minimum'] for w in windows], [w['maximum'] for w in windows], alpha=.25, label='Observed min-max (not CI)')
    axes[2].axhline(0,color='gray',ls='--');axes[2].set(xscale='log',xlabel='Correlation integration cutoff (ps)',ylabel='SW finite integral (eV ps / Angstrom^2)')
    for ax in axes: ax.legend(fontsize=7)
    fig.suptitle('Current MP-0b3 harmonic sensitivity | Historical SW replay\nNo new MD; no physical mobility or first-crack probability',fontsize=11)
    fig.savefig(output/'thermal_validity.png',dpi=170);plt.close(fig)
    for relative in ('solver_v1/silicon_thermal_validity.py','solver_v1/silicon_quantum_diagnostic.py',
                     'solver_v1/silicon_conditional_research.py','solver_v1/silicon_charge_dynamics.py',
                     'solver_v1/probability_pde_2d.py','results/silicon_wafer_feasibility/audit_thermal_validity_v14.py'):
        source(relative)
    save(output/'input_manifest.json', dict(sha256=inputs))
    save(output/'artifact_manifest.json', dict(sha256={p.name:sha(p) for p in output.iterdir() if p.is_file()}))
    print(json.dumps({'complete':True,'states':len(records),'temperature_rows':len(rows),'historical_MD_replayed':len(runs),'new_MD':0}))


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path('.'))
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();main(a.root,a.output)
