"""Predeclared 0 -> +/-5 -> +/-20 -> +/-50 -> 0 MPa static core branches.

No stress is chosen to manufacture a yield value. Each signed branch starts
from the same verified zero-load core; failures are preserved, not continued
as if stable. Zero-load return is STATIC, not a kinetic hold or fatigue run.
"""
import argparse
import json
from pathlib import Path
import time

from .run_current_material_core import argument_parser, run, ROOT
from .run_source_core_reference import load_source_material, build_source_core
from .reference_eam_targets import DEFAULT_CACHE
from .run_vector_registry_audit import save_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--radius', type=float, required=True)
    parser.add_argument('--ring', type=int, required=True)
    parser.add_argument('--source-reference', action='store_true')
    parser.add_argument('--material', type=Path)
    parser.add_argument('--max-seconds-per-case', type=float, default=1800.)
    args = parser.parse_args()
    if args.source_reference and args.material:
        raise ValueError('source and candidate material inputs are separate')
    if args.out.exists():
        raise FileExistsError('new protocol directory required')
    if not json.loads((args.baseline/'summary.json').read_bytes())['final_stability_probe']['stable_on_tested_fixed_boundary']:
        raise ValueError('the zero-load baseline is not a verified fixed-boundary minimum')
    if json.loads((args.baseline/'metadata.json').read_bytes())['shear_traction_MPa'] != 0:
        raise ValueError('zero-load baseline required')
    began = time.perf_counter()
    records = []
    schedule = [5.,20.,50.,0.]
    save_json(args.out/'protocol.json', dict(schedule_absolute_MPa=schedule, signs=[1,-1],
        source_reference_only=args.source_reference, zero_return_is_static=True,
        physical_hold=False, fatigue_cycles=False, actual_yield_calibrated=False,
        baseline=args.baseline.resolve().relative_to(ROOT).as_posix()))
    for sign, label in ((1,'positive'),(-1,'negative')):
        previous = args.baseline
        for index, magnitude in enumerate(schedule):
            name = f'{label}_{index}_{magnitude:g}MPa'
            destination = args.out/name
            options = argument_parser().parse_args(['--out',str(destination),'--radius',str(args.radius),
                '--ring',str(args.ring),'--from-case',str(previous),'--shear-mpa',str(sign*magnitude),
                '--iterations','160','--max-seconds',str(args.max_seconds_per_case),
                '--minimizer','newton_cg','--stability-method','assembled'])
            extra = {}
            if args.material:
                options.material = args.material
            if args.source_reference:
                options.material = DEFAULT_CACHE
                extra.update(material_loader=load_source_material, core_builder=build_source_core)
            try:
                run(options, **extra)
                result = json.loads((destination/'summary.json').read_bytes())
                stable = result['final_stability_probe']['stable_on_tested_fixed_boundary']
                records.append(dict(case=name, executed=True, stable_fixed_boundary=stable))
                if not stable:
                    break
                previous = destination
            except Exception as exc:
                records.append(dict(case=name, executed=False, stable_fixed_boundary=False,
                    error_type=type(exc).__name__, error=str(exc)))
                break
            finally:
                save_json(args.out/'progress.json', dict(completed=False, records=records,
                    elapsed_seconds=time.perf_counter()-began))
    save_json(args.out/'completion.json', dict(completed=True, records=records,
        all_eight_cases_stable=len(records)==8 and all(r['stable_fixed_boundary'] for r in records),
        actual_yield_calibrated=False, infinite_domain_motion_validated=False,
        elapsed_seconds=time.perf_counter()-began))


if __name__ == '__main__':
    main()
