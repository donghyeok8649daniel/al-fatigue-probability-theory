"""Actual local-PDE stress/unit/resolution audit; no parameter fitting."""
from dataclasses import replace, asdict
import json
from pathlib import Path
import time
import numpy as np
from .solver_adapter import UIAnalysisConfig, run_ui_analysis, physical_load_conversion


def main():
    output = Path('results/ui_load_audit')
    output.mkdir(parents=True, exist_ok=True)
    rows = []
    for name, mean, amplitude in [('default', 50., 100.), ('compression', -500., 400.),
                                  ('extreme_mechanism', 900., 2000.)]:
        base = UIAnalysisConfig(stress_mean_mpa=mean, stress_amplitude_mpa=amplitude)
        histories = []
        for quality, na, ns, dt, method in [('preview',21,31,.001,'explicit'),
            ('aligned',41,60,.0005,'implicit'), ('refined',61,90,.00025,'implicit')]:
            config = replace(base, grid_n_a=na, grid_n_s=ns, max_dt=dt, integration_method=method)
            start = time.perf_counter()
            result = run_ui_analysis(config)
            maximum = lambda key: float(np.nanmax(np.abs(result[key])))
            final = lambda key: float(np.asarray(result[key])[-1])
            conversion = physical_load_conversion(config)
            row = dict(case=name, grid=quality, config=asdict(config),
                seconds_runtime=time.perf_counter()-start,
                kappa=conversion['relaxed_axial_kappa'], force_max=conversion['force_max'],
                sigma_over_E_max=config.stress_max_mpa/config.young_mpa,
                local_probability=final('cumulative_absorbed_mass'),
                local_floor_lower_bound=maximum('mass_balance_residual'),
                max_plastic=maximum('plastic_strain'), final_plastic=final('plastic_strain'),
                net_transfer=final('accumulated_net_registry_transfer'),
                gross_grid_traffic=final('cumulative_gross_registry_activity'),
                decomposition_residual=maximum('strain_decomposition_residual'),
                flux_residual=maximum('flux_consistency_residual'),
                registry_residual=maximum('registry_moment_balance_residual'),
                negative_repair=maximum('cumulative_negative_mass_correction'))
            assert np.array_equal(result['initiation_probability'], result['cumulative_absorbed_mass'])
            assert np.all(np.diff(result['cumulative_absorbed_mass']) >= 0)
            np.testing.assert_allclose(result['survival'], 1-np.asarray(result['cumulative_absorbed_mass']), atol=0, rtol=0)
            histories.append(result)
            rows.append(row)
            print(json.dumps(row), flush=True)
            (output/'runs.json').write_text(json.dumps(rows, indent=2), encoding='utf-8')
        p_change = float(np.max(np.abs(np.asarray(histories[-1]['cumulative_absorbed_mass'])-histories[-2]['cumulative_absorbed_mass'])))
        plastic_change = float(np.max(np.abs(np.asarray(histories[-1]['plastic_strain'])-histories[-2]['plastic_strain'])))
        rows[-1]['spatial_dt_probability_difference'] = p_change
        rows[-1]['spatial_dt_plastic_difference'] = plastic_change
        rows[-1]['interpretation'] = 'Combined grid/dt comparison only; no material/kinetic certification or residual-hold claim.'
        (output/'runs.json').write_text(json.dumps(rows, indent=2), encoding='utf-8')


if __name__ == '__main__': main()
