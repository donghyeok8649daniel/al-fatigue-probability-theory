"""Check v3 raw trajectories, numerical controls and immutable source hashes."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import __version__ as scipy_version

from solver_v1.silicon_conditional_research import AMU_EV_PS2_A2
from solver_v1.silicon_crack_research import SpatialSW
from .run_conditional_dynamics import reduction
from .run_local_crack_audit import save_json
from .run_static_probe import source_parameters


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--calculation', type=Path, default=Path('results/silicon_conditional_v3'))
    args = parser.parse_args()
    out = args.calculation
    summary = json.loads((out/'summary.json').read_text(encoding='utf-8'))
    analysis = json.loads((out/'analysis.json').read_text(encoding='utf-8'))
    for name, digest in {**summary['code_hashes'], **summary['input_hashes']}.items():
        assert hashlib.sha256(Path(name.replace('\\', '/')).read_bytes()).hexdigest() == digest, name
    assert hashlib.sha256((out/'summary.json').read_bytes()).hexdigest() == analysis['input_summary_sha256']
    assert hashlib.sha256(Path('results/silicon_wafer_feasibility/analyze_conditional_dynamics.py').read_bytes()).hexdigest() == analysis['analysis_code_sha256']
    source = Path(next(k for k in summary['input_hashes'] if k.endswith('geometry.npz')).replace('\\', '/')).parent
    assert json.loads((source/'geometry.json').read_text(encoding='utf-8')) == summary['geometry']
    geometry = np.load(source/'geometry.npz', allow_pickle=False)
    states = np.load(source/'localized/stationary_positions.npz', allow_pickle=False)
    fixed, reference, bond = geometry['fixed'], states['initial'], states['bond']
    free = np.flatnonzero(~fixed)
    params, parameter_hash = source_parameters()
    assert parameter_hash == summary['source_parameter_sha256']
    model = SpatialSW(params, front_period=summary['geometry']['front_period_A'])
    reduced, _, _ = reduction(model.hessian(reference, free_atoms=free), reference, fixed, bond, (1,))
    displacement = np.zeros_like(reference)
    displacement[free] = reduced.relaxed_lift[:, 0].reshape(-1, 3)
    energy_checks, force_checks = [], []
    assert len(summary['md_runs']) == 13
    assert len({(r['amplitude_A'], r['dt_ps']) for r in summary['md_runs']}) == 13
    for record in summary['md_runs']:
        data = np.load(out/(record['label']+'.npz'), allow_pickle=False)
        for key in data.files:
            assert np.isfinite(data[key]).all(), (record['label'], key)
        assert len(data['time_ps']) == record['steps']//record['stride']+1
        np.testing.assert_allclose(data['time_ps'][-1], summary['duration_ps'], atol=1e-14)
        start = reference+record['amplitude_A']*displacement
        final = data['final_positions_A']
        np.testing.assert_array_equal(final[fixed], reference[fixed])
        np.testing.assert_array_equal(data['final_velocities_A_ps'][fixed], 0.)
        exact_potential = float(np.sum(model.evaluate(final).site_energy-model.evaluate(start).site_energy))
        exact_kinetic = float(.5*record['mass_amu']*AMU_EV_PS2_A2*np.sum(data['final_velocities_A_ps'][free]**2))
        error = max(abs(exact_potential-data['potential_difference_eV'][-1]),
                    abs(exact_kinetic-data['kinetic_energy_eV'][-1]))
        assert error < 2e-10
        np.testing.assert_allclose(data['energy_residual_eV'],
            data['potential_difference_eV']+data['kinetic_energy_eV'], atol=2e-15)
        np.testing.assert_allclose(data['observations'][-1],
            final[bond[1], [1, 0]]-final[bond[0], [1, 0]], atol=2e-15)
        assert record['max_energy_residual_eV']+1e-20 >= np.max(abs(data['energy_residual_eV']))
        energy_checks.append(dict(label=record['label'], final_energy_recompute_error_eV=error))
    # The explicit pair/triplet reference is independent of the moment force
    # evaluation used for integration. This is a new endpoint check, not LAMMPS MD.
    for label, position in [('initial', reference), ('fine_release_endpoint', final)]:
        direct = SpatialSW(params, front_period=summary['geometry']['front_period_A']).evaluate(position, representation='direct')
        moments = model.evaluate(position)
        force_error = float(np.max(abs(direct.gradient-moments.gradient)))
        energy_error = float(abs(np.sum(direct.site_energy-moments.site_energy)))
        assert force_error < 1e-10 and energy_error < 1e-9
        force_checks.append(dict(state=label, force_error_eV_A=force_error, total_energy_error_eV=energy_error))
    for row in analysis['dt_energy_orders']:
        assert all(3.8 < ratio < 4.2 for ratio in row['energy_error_ratios'])
    assert 3.8 < analysis['amplitude_error_ratio_after_time_extrapolation'] < 4.2
    assert analysis['richardson_estimated_odd_errors']['0.002'] < 1e-5
    assert summary['baseline_max_coordinate_drift_A'] < 1e-10
    assert summary['physical_mobility'] is summary['production_t0_seconds'] is None
    assert not any(summary[key] for key in ('finite_T_PMF_certified','wafer_strength_calibrated','silicon_production_enabled'))
    report = dict(numerical_consistency_passed=True, actual_MD_runs_checked=len(energy_checks),
        total_new_trajectory_ps=summary['total_new_trajectory_ps'], final_energy_recomputation=energy_checks,
        direct_triple_force_checks=force_checks,
        geometry_metadata_sha256=hashlib.sha256((source/'geometry.json').read_bytes()).hexdigest(),
        validator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        numpy_version=np.__version__, scipy_version=scipy_version,
        material_calibrated=False, finite_T_PMF_certified=False, kinetic_calibration_approved=False)
    save_json(out/'numerical_validation.json', report)
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
