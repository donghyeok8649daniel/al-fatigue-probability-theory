"""Execute fixed-shape LP/dual audits on existing source interface jets.

No material refit, new energy, production parameter or kinetic inference.
All previously examined held-out states used below are retrospective diagnostic
data, NOT blind validation. Output must be a fresh directory.
"""
import argparse
import hashlib
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np

from .core_interface_compatibility import minimax_compatibility
from .current_core_coefficient_basis import NAMES
from .interface_tangent_calibration import TangentCalibrationProblem, impose_tangents, NONNEGATIVE
from .run_current_material_core import ROOT, DEFAULT_MATERIAL, load_current_material
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--models', type=Path, nargs='+', default=[DEFAULT_MATERIAL,
        ROOT/'results/current_material_core_v22/wider_probe_validation/research_candidate_snapshot.json'])
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError('fresh compatibility output directory required')
    began = time.perf_counter()
    problem = TangentCalibrationProblem(ROOT/'results/fcc111_active_interface/coordination_screening_v19/joint_refinement')
    observations = impose_tangents(problem.observations, problem.source_tangents)
    lookup = {o.name: i for i, o in enumerate(observations)}
    critical = [lookup['saddle_Hxx'], lookup['v19_new_4_Haa']]
    groups = dict(saddle_Hxx=[critical[0]], opening_Haa=[critical[1]], critical_pair=critical,
        fit_curvatures=[i for i, o in enumerate(observations) if o.role == 'fit' and o.units == 'eV/L0^2'],
        inspected_curvatures=[i for i, o in enumerate(observations) if o.role != 'exact' and o.units == 'eV/L0^2'],
        old_fit_all=[i for i, o in enumerate(observations) if o.role == 'fit'])
    exact7 = [i for i, o in enumerate(observations) if o.role == 'exact']
    exact5 = [i for i in exact7 if o_name(observations, i) not in ('perfect_Haa', 'perfect_Hxx')]
    save_json(args.out/'definition.json', dict(starting_HEAD='d1ac6695518d3738b8303687d8c6a7b58877203e',
        observations=[asdict(o) for o in observations], groups=groups, coefficient_names=NAMES,
        source_sha256=problem.source.reference.sha256, exact7=exact7, exact5_control=exact5,
        models=[p.resolve().relative_to(ROOT).as_posix() for p in args.models],
        nonlinear_shapes_optimized=False, scales_are_discrepancies_not_uncertainties=True,
        all_inspected_rows_are_retrospective=True, spectral_constraints_included=False,
        negative_result_is_only_fixed_shape=True, energy_changed=False, production_changed=False))
    summary, raw, dual_rows, tradeoffs = [], [], [], []
    for number, path in enumerate(args.models):
        tick = time.perf_counter()
        model, _, binding = load_current_material(path)
        M = problem.matrix(tuple(model.screened_shape))
        identity = 'original' if path.resolve() == DEFAULT_MATERIAL else f'candidate_{number}'
        write_csv(args.out/f'{identity}_design.csv', [dict(observable=o.name, role=o.role,
            target=o.target, scale=o.scale, units=o.units, **dict(zip(NAMES, row)))
            for o, row in zip(observations, M)])
        for anchor, eq in (('bulk5_and_tangents', exact7), ('bulk5_only_CONTROL', exact5)):
            for sign_name, signs in (('declared_signs', NONNEGATIVE), ('all_signs_relaxed_CONTROL', ())):
                for group, indices in groups.items():
                    result = minimax_compatibility(M, observations, indices,
                                                   nonnegative=signs, exact_rows=eq)
                    case = dict(model=identity, model_sha256=binding['parameter_sha256'],
                                anchor=anchor, signs=sign_name, group=group)
                    raw.append(dict(**case, shape=model.screened_shape, **result))
                    summary.append(dict(**case, completed=result['completed'],
                        eta=result.get('minimax_normalized_error'), lower_bound=result.get('dual_lower_bound'),
                        primal=result.get('primal_violation'), equality=result.get('equality_residual'),
                        stationarity=result.get('stationarity_residual'), dual_gap=result.get('duality_gap'),
                        target_box_feasible=result.get('target_box_feasible'),
                        positive_LJ=result.get('strictly_positive_LJ'), material_accepted=False))
                    if result['completed']:
                        for i, high, low in zip(indices, result['upper_error_duals'], result['lower_error_duals']):
                            dual_rows.append(dict(**case, constraint='target_interval',
                                observable=observations[i].name, upper_error_dual=high, lower_error_dual=low,
                                equality_dual=None, bound_dual=None))
                        for i, multiplier in zip(eq, result['equality_duals']):
                            dual_rows.append(dict(**case, constraint='exact_anchor', observable=observations[i].name,
                                upper_error_dual=None, lower_error_dual=None, equality_dual=multiplier, bound_dual=None))
                        for label, multiplier in zip((*NAMES, 'eta'), result['lower_bound_duals']):
                            dual_rows.append(dict(**case, constraint='nonnegative_bound', observable=label,
                                upper_error_dual=None, lower_error_dual=None, equality_dual=None, bound_dual=multiplier))
                        r = (result['predictions']-np.array([o.target for o in observations]))/np.array([o.scale for o in observations])
                        for i in critical:
                            tradeoffs.append(dict(**case, observable=observations[i].name,
                                target=observations[i].target, baseline=float(M[i]@model.coefficients),
                                prediction=result['predictions'][i], normalized_error=r[i],
                                selected_for_this_LP=i in indices, units=observations[i].units))
                    print(f'{identity}/{anchor}/{sign_name}/{group}: eta={summary[-1]["eta"]}', flush=True)
        save_json(args.out/'checkpoint.json', dict(completed=False, profiles=raw))
        print(f'{identity} elapsed {time.perf_counter()-tick:.2f}s', flush=True)
    write_csv(args.out/'compatibility_summary.csv', summary)
    write_csv(args.out/'dual_certificate.csv', dual_rows)
    write_csv(args.out/'critical_tradeoffs.csv', tradeoffs)
    save_json(args.out/'profiles.json', raw)
    save_json(args.out/'completion.json', dict(completed=True, actual_profiles=len(raw),
        unsuccessful_LPs=sum(not r['completed'] for r in raw), elapsed_seconds=time.perf_counter()-began,
        nonlinear_family_impossibility_proved=False, material_accepted=False,
        physical_yield_calibrated=False, physical_Hz=False, production_changed=False,
        files_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in args.out.glob('*.csv')}))


def o_name(observations, index):
    return observations[index].name


if __name__ == '__main__':
    main()
