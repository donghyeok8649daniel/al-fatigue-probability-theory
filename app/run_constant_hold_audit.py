"""Reproducible constant-stress control; model time, not Al creep calibration."""
from dataclasses import asdict
import json
from pathlib import Path
import time
import numpy as np
from .solver_adapter import UIAnalysisConfig, run_ui_analysis


def main():
    output = Path('results/constant_hold_audit')
    output.mkdir(parents=True, exist_ok=True)
    rows = []
    cases = [('reference',41,60,.0005,1000.,25.,10.),
             ('s_refined',41,90,.0005,1000.,25.,10.),
             ('a_refined',61,60,.0005,1000.,25.,10.),
             ('dt_refined',41,60,.00025,1000.,25.,10.),
             ('same_duration',41,60,.0005,1000.,50.,20.),
             ('zero_stress',41,60,.0005,0.,25.,10.)]
    for name, na, ns, dt, mean, frequency, cycles in cases:
        config = UIAnalysisConfig(stress_mean_mpa=mean, stress_amplitude_mpa=0,
            grid_n_a=na, grid_n_s=ns, max_dt=dt, integration_method='implicit',
            model_frequency=frequency, cycles=cycles, steps_per_cycle=20 if frequency==25 else 10)
        start = time.perf_counter()
        r = run_ui_analysis(config)
        fields = ['model_time','strain','normal_strain','intrawell_strain','plastic_strain',
                  'cumulative_absorbed_mass','first_passage_flux','mass_balance_residual',
                  'accumulated_net_registry_transfer','cumulative_gross_registry_activity',
                  'cumulative_negative_mass_correction','flux_consistency_residual']
        history = {key: np.asarray(r[key]).tolist() for key in fields}
        row = dict(name=name, config=asdict(config), elapsed_seconds=time.perf_counter()-start, history=history)
        rows.append(row)
        (output/'runs.json').write_text(json.dumps(rows, indent=2), encoding='utf-8')
        print(name, {key: history[key][-1] for key in fields}, flush=True)
    summarize(rows, output)


def summarize(rows, output):
    summary = []
    reference = rows[0]['history']
    for row in rows:
        h = row['history']
        assert np.all(np.diff(h['cumulative_absorbed_mass']) >= 0)
        if row['name'] == 'same_duration':
            for field in reference:
                np.testing.assert_array_equal(h[field], reference[field])
        mass_error = float(np.max(np.abs(h['mass_balance_residual'])))
        probability_change = float(np.max(np.abs(np.asarray(h['cumulative_absorbed_mass'])-reference['cumulative_absorbed_mass'])))
        plastic_change = float(np.max(np.abs(np.asarray(h['plastic_strain'])-reference['plastic_strain'])))
        intervals = []
        times = np.asarray(h['model_time'])
        for left, right in zip(np.arange(0, .4, .1), np.arange(.1, .5, .1)):
            i, j = np.argmin(abs(times-left)), np.argmin(abs(times-right))
            intervals.append(dict(start=float(times[i]), end=float(times[j]),
                absorbed_increment=h['cumulative_absorbed_mass'][j]-h['cumulative_absorbed_mass'][i],
                strain_increment=h['strain'][j]-h['strain'][i],
                plastic_increment=h['plastic_strain'][j]-h['plastic_strain'][i]))
        summary.append(dict(name=row['name'], mass_error=mass_error,
            identical_duration_arrays_equal=True if row['name']=='same_duration' else None,
            probability_difference_from_reference=probability_change,
            plastic_difference_from_reference=plastic_change, intervals=intervals,
            first_strain=h['strain'][0], final_strain=h['strain'][-1],
            final_probability=h['cumulative_absorbed_mass'][-1],
            final_plastic=h['plastic_strain'][-1],
            interpretation='No material/kinetic/creep certification; compare only matching loads.'))
    (output/'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')


if __name__ == '__main__': main()
