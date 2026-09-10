"""Execute line-level physical-time checks without granting an a/s clock."""
from __future__ import annotations

import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from .dislocation_line_kinetics import PinnedLineKinetics
from .fatigue_validation_reference import normalized_fatigue_records, fatigue_comparison_gate
from .fcc111_geometry import fcc111_geometry_from_b
from .finite_source_reference import line_energy_coefficient, solve_pinned_graph, pinned_source_branch
from .nonlocal_interface_elasticity import cubic_elastic_tensor, rotate_elastic_tensor
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json
from .vector_material_calibration import LENGTH_M

ROOT = Path(__file__).resolve().parents[1]
MATERIAL = ROOT/'results/fcc111_active_interface/normal_response_v17/nested_range/calibration.json'


def run(kinetics_file, output):
    if output.exists():
        raise FileExistsError('fresh research result directory required')
    started = time.perf_counter()
    kinetic = json.loads(kinetics_file.read_bytes())
    material = json.loads(MATERIAL.read_bytes())
    constants = material['best']['cubic_GPa']
    tensor = rotate_elastic_tensor(cubic_elastic_tensor(*(np.array(constants)*1e9)),
        fcc111_geometry_from_b(1.).plane_basis_in_stacked_cubic_axes())
    lines = {n: line_energy_coefficient(tensor, [LENGTH_M, 0., 0.], samples=n) for n in (32, 64)}
    # Explicit edge-oriented line tangent. The experiment pools edge/mixed
    # characters: this orientation is a hypothesis, not measured pin geometry.
    oriented = {n: replace(line, coefficients=line.coefficients*np.exp(1j*line.frequencies*np.pi/2))
                for n, line in lines.items()}
    drag = LENGTH_M/kinetic['velocity_per_stress_m_per_Pa_s']
    frequency_rows, spatial_rows, static_rows, time_rows, histories = [], [], [], [], []
    for span_um in (.2, 1., 5.):
        span = span_um*1e-6
        # Declared sensitivity geometry. Both factors remain fixed during variation.
        logarithm = np.log(span/(2*LENGTH_M))
        tension = float(oriented[64].stiffness(0)*logarithm)
        model = PinnedLineKinetics(drag, tension, LENGTH_M, span)
        critical = pinned_source_branch(oriented[64], np.pi/2, span_m=span,
            outer_log_ratio=logarithm)['critical_shear_outer_only_Pa']
        for stress_mpa in (.01, .1, 1., 4.):
            stress = stress_mpa*1e6
            nonlinear = None
            if stress < critical:
                nonlinear = solve_pinned_graph(oriented[64], span_m=span,
                    outer_log_ratio=logarithm, shear_Pa=stress, segments=256)
                if not nonlinear['converged']:
                    raise RuntimeError('independent nonlinear static branch failed')
            prediction = float(model.static_mean_bow_m(stress))
            measured = None if nonlinear is None else nonlinear['swept_area_m2']/span
            static_rows.append(dict(span_um=span_um, shear_MPa=stress_mpa,
                linear_mean_bow_m=prediction, nonlinear_mean_bow_m=measured,
                linearization_relative_error=None if measured is None else (prediction-measured)/measured,
                predicted_maximum_slope=float(model.static_maximum_slope(stress)),
                outer_only_fold_MPa=critical/1e6,
                experimental_yield_MPa=None,
                status='subcritical bow, no yield claim' if nonlinear else 'outside certified graph branch'))
        for frequency in (.1, 1., 25., 100., 1/(2*np.pi*model.slowest_seconds)):
            exact = model.frequency_response(frequency, tail_tolerance=1e-10)
            G = complex(exact['transfer'])
            frequency_rows.append(dict(span_um=span_um, frequency_Hz=frequency,
                slowest_line_time_seconds=model.slowest_seconds, drag_Pa_s=drag,
                line_stiffness_J_m=tension, omega_tau=2*np.pi*frequency*model.slowest_seconds,
                amplitude_ratio=abs(G), phase_degrees=float(np.angle(G, deg=True)),
                reciprocal_angular_stiffness_change_J_m=float(abs(oriented[32].stiffness(0)*logarithm-tension)),
                maximum_series_mode=exact['maximum_mode'], absolute_series_tail=exact['absolute_tail_bound'],
                time_scope='measured line-velocity slope + hypothetical pins + fixed 0K elastic comparator; NOT production a/s'))
            for segments in (16, 32, 64, 128):
                result = model.discrete_harmonic(frequency, .01e6, segments=segments)
                spatial_rows.append(dict(span_um=span_um, frequency_Hz=frequency, segments=segments,
                    absolute_transfer_error=abs(result['transfer']-G),
                    phase_error_degrees=float(np.angle(result['transfer']/G, deg=True)),
                    work_per_cycle_J=result['work_per_cycle_J'],
                    dissipation_per_cycle_J=result['dissipation_per_cycle_J'],
                    work_dissipation_relative_residual=result['energy_balance_residual_J']/result['dissipation_per_cycle_J']))
        if span_um != 1.:
            continue
        period = 2*np.pi*model.slowest_seconds
        frequency = 1/period
        exact = complex(model.frequency_response(frequency, tail_tolerance=1e-11)['transfer'])
        equilibrium = float(model.static_mean_bow_m(.01e6))
        for steps in (64, 128, 256, 512):
            load_time = np.linspace(0, 8*period, 8*steps+1)
            load = .01e6*np.cos(2*np.pi*load_time/period)
            hold_dt = period/steps
            hold = np.arange(1, int(np.ceil(20*model.slowest_seconds/hold_dt))+1)*hold_dt
            times = np.r_[load_time, load_time[-1]+hold]
            shear = np.r_[load, np.zeros(len(hold))]
            result = model.integrate_backward_euler(times, shear, segments=128)
            chosen = (times >= 6*period) & (times < 8*period)
            X = np.column_stack([np.cos(2*np.pi*times[chosen]/period),
                                 np.sin(2*np.pi*times[chosen]/period), np.ones(chosen.sum())])
            coefficients = np.linalg.lstsq(X, result['mean_bow_m'][chosen], rcond=None)[0]
            observed = (coefficients[0]-1j*coefficients[1])/equilibrium
            time_rows.append(dict(steps_per_cycle=steps, segments=128, dt_seconds=hold_dt,
                absolute_transfer_error=abs(observed-exact), amplitude_ratio=abs(observed),
                phase_degrees=float(np.angle(observed, deg=True)),
                final_hold_seconds=float(hold[-1]), final_mean_bow_m=result['mean_bow_m'][-1],
                final_bow_over_static_mean=result['mean_bow_m'][-1]/equilibrium,
                permanent_slip_generated=False, opening_probability=None))
            if steps == 512:
                for t, stress, bow in zip(times, shear, result['mean_bow_m']):
                    histories.append(dict(line_time_seconds=t, applied_shear_Pa=stress,
                        mean_bow_m=bow, phase='cyclic' if t <= 8*period else 'zero_stress_hold'))
    for name, rows in [('frequency_response', frequency_rows), ('spatial_refinement', spatial_rows),
                       ('stress_scenarios', static_rows), ('time_refinement', time_rows),
                       ('cyclic_unload_hold', histories), ('fatigue_reference', normalized_fatigue_records())]:
        write_csv(output/f'{name}.csv', rows)
    gate = fatigue_comparison_gate(prediction_endpoint='local_atomic_opening_absorption',
        physical_clock_calibrated=False, specimen_mapping_validated=False,
        control_protocol_matched=False, microstructure_matched=False,
        probability_resolution_certified=False)
    save_json(output/'scope.json', dict(completed=True, elapsed_seconds=time.perf_counter()-started,
        material_file=str(MATERIAL.relative_to(ROOT)), material_sha256=hashlib.sha256(MATERIAL.read_bytes()).hexdigest(),
        material_accepted=material['material_accepted'], cubic_GPa=constants,
        material_temperature='0 K static comparator; NOT matched to measured 23 C elasticity',
        mobility_source_doi=kinetic['source_doi'], kinetic_file_sha256=hashlib.sha256(kinetics_file.read_bytes()).hexdigest(),
        velocity_slope_m_per_Pa_s=kinetic['velocity_per_stress_m_per_Pa_s'],
        conditional_drag_Pa_s=drag, burgers_m=LENGTH_M,
        burgers_convention='fixed model 0K b; not a new measured 23 C Burgers vector',
        line_character='explicit edge reference; measured kinetic data pool edge/mixed',
        source_spans_are_measured=False, core_radius_over_b=2, outer_radius_over_span=1,
        geometry_or_drag_fitted_to_yield_or_fatigue=False,
        source_coordinate_seconds_available=True, matched_production_seconds_available=False,
        production_M_a_phys=None, production_M_s_phys=None, production_t0_seconds=None,
        fatigue_comparison=gate, actual_yield_calibrated=False, fatigue_calibrated=False,
        production_changed=False))
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), layout='constrained')
    for span in (.2, 1., 5.):
        rows = [r for r in frequency_rows if r['span_um'] == span]
        rows.sort(key=lambda r: r['frequency_Hz'])
        axes[0].loglog([r['frequency_Hz'] for r in rows],
                      [-r['phase_degrees'] for r in rows], 'o-', label=f'L={span} um (hypothetical)')
    axes[0].set(xlabel='Line-reference forcing [Hz]', ylabel='Lag magnitude [degrees]')
    axes[0].legend(fontsize=7)
    axes[1].plot(np.array([r['line_time_seconds'] for r in histories])*1e6,
                 np.array([r['mean_bow_m'] for r in histories])*1e9)
    axes[1].set(xlabel='Line-reference time [microseconds]', ylabel='Mean bow [nm]')
    for ax in axes:
        ax.grid(alpha=.2)
    fig.suptitle('Source-informed line reference: reversible bowing, NOT calibrated fatigue / a-s clock')
    fig.savefig(output/'line_response.png', dpi=150)
    fig.savefig(output/'line_response.svg')
    plt.close(fig)
    return time_rows


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--kinetics', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.kinetics, args.out), default=lambda x: x.tolist()), flush=True)
