"""Fresh Si Gaussian-measure and conservative atom-dynamics audit (research).

Reads the saved v2 coordinates, reevaluates all energies/Hessians, and runs new
velocity-Verlet trajectories. It does NOT replay MD, count fracture trajectories,
fit damping, certify a finite-T PMF, or activate the production Si/Hz gates.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np
from scipy.constants import hbar, electron_volt

from solver_v1.silicon_crack_research import RelaxedCoordinates, SpatialSW
from solver_v1.silicon_conditional_research import (
    AMU_EV_PS2_A2, KB_EV_K, SI_MASS_AMU, conditional_harmonic,
    harmonic_free_energy_difference, harmonic_memory, harmonic_release, velocity_verlet,
)
from .run_local_crack_audit import save_json
from .run_static_probe import source_parameters


TEMPERATURES_K = (100., 300., 600.)  # Diagnostic probes, not fitted wafer conditions.


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def reduction(hessian, positions, fixed, bond, components):
    coordinates = RelaxedCoordinates(positions, fixed, bond=bond, components=components)
    dofs = (np.flatnonzero(~fixed)[:, None]*3+np.arange(3)).ravel()
    b = coordinates.tangent_matrix()[dofs]
    d = np.zeros((len(dofs), len(components)))
    for k, (a, z, _) in enumerate(coordinates.constraints):
        d[a, k], d[z, k] = -.5, .5
    return conditional_harmonic(hessian, b, d), b, d


def write_csv(path, rows):
    with Path(path).open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path('results/silicon_local_crack_v2'))
    parser.add_argument('--output', type=Path, default=Path('results/silicon_conditional_v3'))
    parser.add_argument('--profile-stride', type=int, default=2)
    parser.add_argument('--duration-ps', type=float, default=1.)
    args = parser.parse_args()
    if args.profile_stride < 1 or not np.isfinite(args.duration_ps) or args.duration_ps <= 0:
        parser.error('positive profile stride and duration required')
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        parser.error('output must be empty; preserve previous research runs')
    start = time.perf_counter()
    p, parameter_hash = source_parameters()
    geometry_path = args.source/'front4/geometry.npz'
    geometry = np.load(geometry_path, allow_pickle=False)
    meta = json.loads((args.source/'front4/geometry.json').read_text(encoding='utf-8'))
    states_path = args.source/'front4/localized/stationary_positions.npz'
    states = np.load(states_path, allow_pickle=False)
    profile_path = args.source/'front4/localized/profile_positions.npz'
    profile = np.load(profile_path, allow_pickle=False)
    fixed, bond = geometry['fixed'], states['bond']
    free = np.flatnonzero(~fixed)
    model = SpatialSW(p, front_period=meta['front_period_A'])
    reference = states['initial']
    reference_sites = model.evaluate(reference).site_energy.copy()
    provenance = dict(schema='silicon-conditional-dynamics/3',
        head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        branch=subprocess.check_output(['git', 'branch', '--show-current'], text=True).strip(),
        units=dict(energy='total eV', position='Angstrom', newtonian_time='ps', mass='amu'),
        reference_mass_amu=SI_MASS_AMU, mass_eV_ps2_A2=SI_MASS_AMU*AMU_EV_PS2_A2,
        kB_eV_K=KB_EV_K, temperatures_K=TEMPERATURES_K, duration_ps=args.duration_ps,
        atom_count=len(reference), free_atom_count=len(free), geometry=meta,
        source_parameter_sha256=parameter_hash,
        input_hashes={str(path):digest(path) for path in (geometry_path, states_path, profile_path)},
        code_hashes={name:digest(name) for name in (
            'solver_v1/silicon_crack_research.py', 'solver_v1/silicon_conditional_research.py',
            'results/silicon_wafer_feasibility/run_conditional_dynamics.py')},
        production_t0_seconds=None, physical_mobility=None, finite_T_PMF_certified=False,
        wafer_strength_calibrated=False, silicon_production_enabled=False,
        temperature_interpretation='classical Gaussian diagnostic; no finite-T ensemble sampled',
        dynamics_interpretation='new zero-temperature conservative small-displacement releases; no thermostat or damping fit')
    save_json(output/'running_summary.json', provenance)
    times = np.arange(round(args.duration_ps/.001)+1)*.001
    reports, modes, kernels, spectra, thermal_rows = {}, {}, {}, {}, []
    initial_h = initial_reduction = initial_d = None
    for name in ('initial', 'saddle', 'opened_minimum'):
        r = states[name]
        evaluation = model.evaluate(r)
        h = model.hessian(r, free_atoms=free)
        one, basis, d = reduction(h, r, fixed, bond, (1,))
        two, _, _ = reduction(h, r, fixed, bond, (1, 0))
        memory = harmonic_memory(one, times)
        values, weights = memory['eigenvalues_eV_A2'], memory['weights_eV_A2'][:, 0, 0]
        mass_q = float(memory['q_mass_eV_ps2_A2'][0, 0])
        k = float(one.relaxed_curvature[0, 0])
        rate = np.sqrt(abs(k)/mass_q)
        marginal_error = one.logdet_bath-two.logdet_bath-np.log(2*two.relaxed_curvature[1, 1])
        if abs(marginal_error) > 1e-9:
            raise RuntimeError('one/two-coordinate Gaussian marginalization mismatch')
        reports[name] = dict(gap_A=float(r[bond[1], 1]-r[bond[0], 1]),
            energy_difference_eV=float(np.sum(evaluation.site_energy-reference_sites)),
            full_force_residual_eV_A=float(np.max(abs(evaluation.gradient[free]))),
            bath_dimension=len(values), bath_min_eigenvalue_eV_A2=float(values[0]),
            logdet_a=one.logdet_bath, logdet_as=two.logdet_bath,
            logdet_marginalization_error=float(marginal_error),
            relaxed_curvature_a_eV_A2=k, relaxed_curvature_as_eV_A2=two.relaxed_curvature.tolist(),
            q_mass_eV_ps2_A2=mass_q,
            relaxed_path_mass_eV_ps2_A2=float(SI_MASS_AMU*AMU_EV_PS2_A2*np.sum(one.relaxed_lift**2)),
            memory_zero_eV_A2=float(weights.sum()),
            memory_below_abs_curvature_rate_fraction=float(weights[memory['angular_frequency_ps'] < rate].sum()/weights.sum()),
            abs_curvature_angular_rate_ps=float(rate),
            slowest_bath_angular_frequency_ps=float(memory['angular_frequency_ps'][0]),
            max_hbar_omega_over_kBT_300K=float((hbar/electron_volt/1e-12)*memory['angular_frequency_ps'][-1]/(KB_EV_K*300)))
        if name == 'initial':
            initial_h, initial_reduction, initial_d = h, one, d
        kernels[name] = memory['kernel_eV_A2'][:, 0, 0]
        spectra[name+'_lambda_eV_A2'] = values
        spectra[name+'_omega_ps'] = memory['angular_frequency_ps']
        spectra[name+'_weights_eV_A2'] = weights
        # Slice probes are necessary checks on a Gaussian approximation, never
        # a replacement for a joint conditional finite-temperature integral.
        selected = sorted(set([0, 1, 2, int(np.argmax(weights))]))
        for mode in selected:
            displacement = np.zeros_like(r)
            displacement[free] = np.asarray(basis@memory['eigenvectors'][:, mode]).reshape(-1, 3)
            for temperature in TEMPERATURES_K:
                sigma = np.sqrt(KB_EV_K*temperature/values[mode])
                for multiple in (-2., -1., -.5, .5, 1., 2.):
                    delta = multiple*sigma
                    exact = float(np.sum(model.evaluate(r+delta*displacement).site_energy-evaluation.site_energy))
                    gaussian = .5*values[mode]*delta**2
                    thermal_rows.append(dict(state=name, mode=mode, temperature_K=temperature,
                        sigma_coordinate_A=float(sigma), sigma_multiple=multiple,
                        max_atom_displacement_A=float(np.max(np.linalg.norm(delta*displacement, axis=1))),
                        exact_energy_difference_eV=exact, gaussian_energy_eV=float(gaussian),
                        residual_over_kBT=float((exact-gaussian)/(KB_EV_K*temperature))))
        print('state', name, 'logdet', one.logdet_bath, 'Gamma0/K', weights.sum()/abs(k), flush=True)
    for report in reports.values():
        report['gaussian_differences_at_static_geometries_eV'] = {
            str(int(t)):float(harmonic_free_energy_difference(report['energy_difference_eV'],
                report['logdet_a'], reports['initial']['logdet_a'], t)) for t in TEMPERATURES_K}
    save_json(output/'stationary_gaussian.json', reports)
    write_csv(output/'thermal_mode_slices.csv', thermal_rows)
    np.savez_compressed(output/'memory_spectra.npz', time_ps=times, **spectra,
                        **{name+'_kernel_eV_A2':v for name, v in kernels.items()})
    profile_rows = []
    for index in sorted(set(range(0, len(profile['q_A']), args.profile_stride)) | {len(profile['q_A'])-1}):
        r = profile['positions'][index]
        e = model.evaluate(r)
        row = dict(source_index=index, gap_A=float(profile['q_A'][index]),
                   energy_difference_eV=float(np.sum(e.site_energy-reference_sites)))
        try:
            reduced, _, _ = reduction(model.hessian(r, free_atoms=free), r, fixed, bond, (1,))
        except np.linalg.LinAlgError:
            row.update(conditional_stable=False, logdet_bath=None,
                **{f'gaussian_{int(t)}K_eV':None for t in TEMPERATURES_K})
        else:
            row.update(conditional_stable=True, logdet_bath=reduced.logdet_bath,
                **{f'gaussian_{int(t)}K_eV':float(harmonic_free_energy_difference(row['energy_difference_eV'],
                    reduced.logdet_bath, reports['initial']['logdet_a'], t)) for t in TEMPERATURES_K})
        profile_rows.append(row)
        print('profile', index, 'conditional stable', row['conditional_stable'], flush=True)
    write_csv(output/'gaussian_profile.csv', profile_rows)
    observation = lambda r: r[bond[1], [1, 0]]-r[bond[0], [1, 0]]
    q0 = observation(reference)
    observe_free = np.zeros((2, 3*len(free)))
    for row, component in enumerate((1, 0)):
        observe_free[row, 3*np.searchsorted(free, bond[0])+component] = -1.
        observe_free[row, 3*np.searchsorted(free, bond[1])+component] = 1.
    expected = harmonic_release(initial_h, initial_reduction.relaxed_lift[:, 0], observe_free, times)
    np.savez_compressed(output/'harmonic_release.npz', time_ps=times, response_per_A=expected)
    runs, histories = [], {}
    displacement = np.zeros_like(reference)
    displacement[free] = initial_reduction.relaxed_lift[:, 0].reshape(-1, 3)
    for amplitude in (0., .004, -.004, .002, -.002):
        for dt in ((.00025,) if amplitude == 0 else (.001, .0005, .00025)):
            label = f'gap_{amplitude:+.4f}_dt_{dt:.5f}'
            r = reference+amplitude*displacement
            energy = float(np.sum(model.evaluate(r).site_energy-reference_sites))
            run = velocity_verlet(model, r, fixed, dt_ps=dt, steps=round(args.duration_ps/dt),
                stride=round(.001/dt), observation=observation,
                progress=lambda step: print('MD', label, step, flush=True))
            np.testing.assert_allclose(run['time_ps'], times, atol=1e-14, rtol=0)
            np.savez_compressed(output/(label+'.npz'), **{k:v for k,v in run.items() if isinstance(v, np.ndarray)})
            record = {k:v for k,v in run.items() if not isinstance(v, np.ndarray)}
            record.update(label=label, amplitude_A=amplitude, initial_excess_energy_eV=energy,
                preparation='linear conditional minimum tangent; zero initial velocity',
                max_energy_residual_over_excess=run['max_energy_residual_eV']/energy if energy > 0 else None)
            if record['max_fixed_displacement_A'] != 0:
                raise RuntimeError('fixed grips moved')
            if amplitude and record['max_energy_residual_over_excess'] > .02:
                raise RuntimeError('conservative energy error exceeds declared 2% control limit')
            runs.append(record)
            histories[(amplitude, dt)] = run['observations']
            save_json(output/'md_runs.json', runs)
            print('completed', label, 'energy relative', record['max_energy_residual_over_excess'], flush=True)
    controls = []
    baseline = histories[(0., .00025)]
    for amplitude in (.004, .002):
        for dt in (.001, .0005, .00025):
            plus, minus = histories[(amplitude, dt)], histories[(-amplitude, dt)]
            odd = (plus-minus)/(2*amplitude)
            even = (plus+minus-2*baseline)/(2*amplitude)
            controls.append(dict(amplitude_A=amplitude, dt_ps=dt,
                max_odd_minus_harmonic=float(np.max(abs(odd-expected))),
                rms_odd_minus_harmonic=float(np.sqrt(np.mean((odd-expected)**2))),
                max_even_over_amplitude=float(np.max(abs(even)))))
    save_json(output/'md_controls.json', controls)
    save_json(output/'summary.json', dict(**provenance, stationary=reports, md_runs=runs,
        md_controls=controls, computed_profile_points=len(profile_rows),
        stable_profile_points=sum(row['conditional_stable'] for row in profile_rows),
        total_new_trajectory_ps=len(runs)*args.duration_ps,
        baseline_max_coordinate_drift_A=float(np.max(abs(baseline-q0))),
        elapsed_seconds=time.perf_counter()-start, computation_completed=True,
        scientific_gate='Gaussian candidate and finite conservative memory only; no validated finite-T PMF or Markov mobility'))
    print('complete', time.perf_counter()-start, flush=True)


if __name__ == '__main__':
    main()
