"""Analyze fresh v3 trajectories and candidate retained coordinates, without MD."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import cho_factor, cho_solve

from solver_v1.silicon_crack_research import SpatialSW
from .run_conditional_dynamics import reduction
from .run_local_crack_audit import save_json
from .run_static_probe import source_parameters


def explained_force_fraction(bath_hessian, force_coupling, observables):
    """Conditional Gaussian variance removed by independent linear observables.

    kBT cancels from this ratio. This is a static coordinate diagnostic, not a
    fit, a dynamical timescale certificate or an argument for a large state.
    Singular/duplicate observable combinations are rejected, not regularized.
    """
    coupling = np.asarray(force_coupling, float).reshape(-1)
    c = np.atleast_2d(np.asarray(observables, float))
    factor = cho_factor(bath_hessian, lower=True)
    variance = float(coupling@cho_solve(factor, coupling))
    covariance = c@cho_solve(factor, c.T)
    cross = c@cho_solve(factor, coupling)
    explained = float(cross@cho_solve(cho_factor(covariance, lower=True), cross))
    if variance <= 0 or explained < -1e-12 or explained > variance+1e-10:
        raise ValueError('invalid positive Gaussian variance decomposition')
    return dict(observable_count=len(c), explained_force_variance_fraction=explained/variance,
                remaining_force_variance_fraction=(variance-explained)/variance)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--calculation', type=Path, default=Path('results/silicon_conditional_v3'))
    args = parser.parse_args()
    out = args.calculation
    summary = json.loads((out/'summary.json').read_text(encoding='utf-8'))
    for name, expected_hash in summary['code_hashes'].items():
        if hashlib.sha256(Path(name).read_bytes()).hexdigest() != expected_hash:
            raise RuntimeError('research source changed since calculation: '+name)
    for name, expected_hash in summary['input_hashes'].items():
        if hashlib.sha256(Path(name.replace('\\', '/')).read_bytes()).hexdigest() != expected_hash:
            raise RuntimeError('research input changed since calculation: '+name)
    source = Path(next(k for k in summary['input_hashes'] if k.endswith('geometry.npz')).replace('\\', '/')).parent
    geometry = np.load(source/'geometry.npz', allow_pickle=False)
    states = np.load(source/'localized/stationary_positions.npz', allow_pickle=False)
    r, fixed, bond = states['initial'], geometry['fixed'], states['bond']
    free = np.flatnonzero(~fixed)
    lookup = {int(atom):i for i, atom in enumerate(free)}
    p, _ = source_parameters()
    model = SpatialSW(p, front_period=summary['geometry']['front_period_A'])
    one, b, d = reduction(model.hessian(r, free_atoms=free), r, fixed, bond, (1,))

    def observable(pair, component):
        row = np.zeros(3*len(free))
        row[3*lookup[int(pair[0])]+component] = -1.
        row[3*lookup[int(pair[1])]+component] = 1.
        return np.asarray(row@b).reshape(-1)

    crossing = geometry['crossing_bonds']
    eligible = crossing[np.all(~fixed[crossing], axis=1)]
    eligible = eligible[np.any(eligible != bond, axis=1)]
    centers = geometry['reference'][eligible].mean(axis=1)
    selected_x = geometry['reference'][bond, 0].mean()
    along_front = eligible[np.isclose(centers[:, 0], selected_x, atol=1e-8)]
    ahead_x = centers[centers[:, 0] > selected_x+1e-8, 0].min()
    next_front = eligible[np.isclose(centers[:, 0], ahead_x, atol=1e-8)]
    observations = {
        'local_slip':np.array([observable(bond, 0)]),
        'other_openings_on_same_front':np.array([observable(pair, 1) for pair in along_front]),
        'next_front_openings':np.array([observable(pair, 1) for pair in next_front]),
        'all_other_mobile_shuffle_openings':np.array([observable(pair, 1) for pair in eligible]),
    }
    observations['local_slip_and_same_front'] = np.vstack([
        observations['local_slip'], observations['other_openings_on_same_front']])
    coordinates = {name:explained_force_fraction(one.bath_hessian, one.cross, c)
                   for name, c in observations.items()}
    # Independent block-Schur check for adding the same pair's slip coordinate.
    two, _, _ = reduction(model.hessian(r, free_atoms=free), r, fixed, bond, (1, 0))
    expected = (two.relaxed_curvature[0, 0]-one.relaxed_curvature[0, 0])/(one.bare_curvature[0, 0]-one.relaxed_curvature[0, 0])
    np.testing.assert_allclose(coordinates['local_slip']['explained_force_variance_fraction'], expected, atol=1e-12)
    thermal = list(csv.DictReader((out/'thermal_mode_slices.csv').open(encoding='utf-8')))
    worst_slices = {str(t):max(abs(float(row['residual_over_kBT'])) for row in thermal
        if float(row['temperature_K']) == t) for t in (100, 300, 600)}
    harmonic = np.load(out/'harmonic_release.npz', allow_pickle=False)
    time, expected_response = harmonic['time_ps'], harmonic['response_per_A']

    def odd(amplitude, dt):
        plus = np.load(out/f'gap_{amplitude:+.4f}_dt_{dt:.5f}.npz', allow_pickle=False)['observations']
        minus = np.load(out/f'gap_{-amplitude:+.4f}_dt_{dt:.5f}.npz', allow_pickle=False)['observations']
        return (plus-minus)/(2*amplitude)

    extrapolated = {str(a):(4*odd(a, .00025)-odd(a, .0005))/3 for a in (.004, .002)}
    # This cancels leading O(dt^2) numerics; it is not another MD trajectory.
    richardson_errors = {a:float(np.max(abs(value-expected_response))) for a, value in extrapolated.items()}
    direct = odd(.002, .00025)
    k = summary['stationary']['initial']['relaxed_curvature_a_eV_A2']
    mq = summary['stationary']['initial']['q_mass_eV_ps2_A2']
    path_mass = summary['stationary']['initial']['relaxed_path_mass_eV_ps2_A2']
    adiabatic = np.cos(np.sqrt(k/mq)*time)
    relaxed_path = np.cos(np.sqrt(k/path_mass)*time)
    error_orders = []
    for a in (.004, -.004, .002, -.002):
        records = sorted([r for r in summary['md_runs'] if r['amplitude_A'] == a], key=lambda r:-r['dt_ps'])
        ratios = [records[j]['max_energy_residual_eV']/records[j+1]['max_energy_residual_eV'] for j in (0, 1)]
        error_orders.append(dict(amplitude_A=a, energy_error_ratios=ratios,
            orders_base2=np.log2(ratios).tolist()))
    negative = np.flatnonzero(direct[:, 0] < 0)
    analysis = dict(input_summary_sha256=hashlib.sha256((out/'summary.json').read_bytes()).hexdigest(),
        analysis_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        coordinate_variance_diagnostics=coordinates, mode_slice_max_abs_error_kBT=worst_slices,
        dt_energy_orders=error_orders, richardson_estimated_odd_errors=richardson_errors,
        amplitude_error_ratio_after_time_extrapolation=richardson_errors['0.004']/richardson_errors['0.002'],
        direct_fine_normalized_gap_min=float(direct[:, 0].min()),
        first_sampled_negative_response_ps=float(time[negative[0]]) if len(negative) else None,
        max_static_bath_bare_mass_response_error=float(np.max(abs(direct[:, 0]-adiabatic))),
        max_static_bath_path_mass_response_error=float(np.max(abs(direct[:, 0]-relaxed_path))),
        thermal_joint_distribution_sampled=False, full_PMF_certified=False,
        constant_drag_identified=False, production_t0_seconds=None)
    save_json(out/'analysis.json', analysis)
    print(json.dumps(analysis, indent=2), flush=True)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':10, 'svg.hashsalt':'silicon-conditional-v3'})
    rows = list(csv.DictReader((out/'gaussian_profile.csv').open(encoding='utf-8')))
    gaps = np.array([float(row['gap_A']) for row in rows])
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    ax = axes[0, 0]
    ax.plot(gaps, [float(row['energy_difference_eV']) for row in rows], 'o-', ms=3, label='0 K energy branch')
    for t, color in [(300, '#c45d20'), (600, '#25886c')]:
        values = [float(row[f'gaussian_{t}K_eV']) if row[f'gaussian_{t}K_eV'] else np.nan for row in rows]
        ax.plot(gaps, values, 'o--', ms=3, label=f'Classical Gaussian candidate, {t} K', color=color)
    ax.set(xlabel='Local normal gap [Angstrom]', ylabel='Difference from initial geometry [eV]',
           title='A. Conditional Gaussian correction (not full PMF)')
    ax.legend(fontsize=8)
    ax = axes[0, 1]
    spectral = np.load(out/'memory_spectra.npz', allow_pickle=False)
    for name, label in [('initial','Intact minimum'), ('opened_minimum','Opened minimum')]:
        kernel = spectral[name+'_kernel_eV_A2']
        ax.plot(time, kernel/kernel[0], label=label)
    ax.axhline(0, color='.65', lw=.6)
    ax.set(xlim=(0, .5), xlabel='Newtonian reference time [ps]', ylabel='Memory / initial memory',
           title='B. Finite harmonic bath retains memory')
    ax.legend(fontsize=8)
    ax = axes[1, 0]
    ax.plot(time, expected_response[:, 0], color='#222222', lw=2, label='Full harmonic Cartesian response')
    ax.plot(time[::8], direct[::8, 0], '.', ms=3, color='#c45d20', label='New nonlinear atom dynamics, 0.25 fs')
    ax.plot(time, adiabatic, '--', color='#64799c', alpha=.8, label='Instantaneous bath, bare mass')
    ax.plot(time, relaxed_path, ':', color='#25886c', label='Instantaneous bath, relaxed-path mass')
    ax.set(xlim=(0, .2), xlabel='Newtonian reference time [ps]', ylabel='Odd gap response / initial perturbation',
           title='C. Exact bath dynamics versus static elimination')
    ax.legend(fontsize=7, loc='lower right')
    ax = axes[1, 1]
    for a, marker in [(.004, 'o'), (.002, 's')]:
        controls = [row for row in summary['md_controls'] if row['amplitude_A'] == a]
        ax.loglog([row['dt_ps']*1000 for row in controls],
                  [row['max_odd_minus_harmonic'] for row in controls], marker+'-', label=f'+/- {a:g} Angstrom')
    ax.set(xlabel='Integration step [fs]', ylabel='Maximum normalized response discrepancy',
           title='D. Timestep refinement of actual atom trajectories')
    ax.legend(fontsize=8)
    for ax in axes.ravel():
        ax.grid(alpha=.2)
    fig.suptitle('Pure Si / original SW / fixed grips: conditional measure and dynamics audit', fontsize=14)
    fig.savefig(out/'conditional_dynamics.png', dpi=180)
    fig.savefig(out/'conditional_dynamics.svg', metadata={'Date':None})
    plt.close(fig)


if __name__ == '__main__':
    main()
