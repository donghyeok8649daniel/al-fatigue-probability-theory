"""Actually execute v20 fixed-shape tangent controls; never overwrite v19."""
import argparse
from pathlib import Path
import time

import numpy as np

from .interface_tangent_calibration import TangentCalibrationProblem, tangent_summary
from .run_vector_registry_audit import save_json
from .run_low_stress_cyclic_diagnostic import write_csv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent-joint', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError('new output directory required; preserve prior experiments')
    began = time.perf_counter()
    problem = TangentCalibrationProblem(args.parent_joint)
    save_json(args.out/'definition.json', problem.definition | dict(stage='fixed_shape_controls',
        interpolation_fractions=[0., .25, .5, .75, 1.], actual_optimization=True,
        actual_radial_optimization=False, fixed_shape_control=True, joint_refinement=False,
        inherited_optimizer_metadata_is_historical=True))
    summaries, all_tables = [], []
    for family, parent in problem.parents.items():
        save_json(args.out/family/'definition.json', problem.definition | dict(
            stage='fixed_shape_controls', fixed_family=family, actual_radial_optimization=False,
            fixed_shape_control=True, joint_refinement=False,
            inherited_optimizer_metadata_is_historical=True))
        shape = parent['shape']
        replay = problem.matrix(tuple(shape))@np.array(parent['coefficients'])
        source_haa = problem.source_tangents['perfect_Haa']
        source_hxx = problem.source_tangents['perfect_Hxx']
        parent_haa = float(replay[problem.tangent_rows['perfect_Haa']])
        cases = [('baseline', {}), ('normal_exact', {'perfect_Haa': source_haa}),
                 ('registry_exact', {'perfect_Hxx': source_hxx}),
                 ('both_exact', dict(problem.source_tangents))]
        for fraction in (0., .25, .5, .75):
            cases.append((f'normal_scan_{fraction:g}', dict(perfect_Haa=(1-fraction)*parent_haa+fraction*source_haa,
                                                          perfect_Hxx=source_hxx)))
        for case, imposed in cases:
            started = time.perf_counter()
            fit = problem.fit(shape, imposed)
            table = fit.pop('residual_table', [])
            fit['elapsed_seconds'] = time.perf_counter()-started
            save_json(args.out/family/case/'calibration.json', dict(completed=fit['completed'], best=fit,
                screening_law='power', shape=shape, material_accepted=False))
            if table:
                write_csv(args.out/family/case/'residuals.csv', table)
                all_tables.extend(dict(family=family, case=case, **row) for row in table)
            row = tangent_summary(case, family, fit, problem)
            summaries.append(row)
            write_csv(args.out/'summary.csv', summaries)
            print(family, case, row, flush=True)
    write_csv(args.out/'all_residuals.csv', all_tables)
    save_json(args.out/'completion.json', dict(completed=True, executed_profiles=len(summaries),
        verified_profiles=sum(r['completed'] for r in summaries), summary=summaries,
        elapsed_seconds=time.perf_counter()-began, material_accepted=False, physical_kinetics=False))


if __name__ == '__main__':
    main()
