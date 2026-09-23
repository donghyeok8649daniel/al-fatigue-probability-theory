"""Record completed low-cost validation from its actual stored output.

This does not run tests, certify incomplete atom calculations, or publish raw
local tracebacks. A test implementation hash pins what the log describes.
"""
import argparse,hashlib,json,re,sys
from pathlib import Path


def main(args):
    root=Path(__file__).resolve().parents[2]
    text=args.test_log.read_text(encoding='utf-8-sig')
    matches=re.findall(r'(\d+) passed in ([0-9.]+)s',text)
    if len(matches)!=1 or int(matches[0][0])!=args.expected_tests or 'failed' in text.lower():
        raise ValueError(f'expected actual completed {args.expected_tests}-test run')
    files=[
        'solver_v1/silicon_initiation_probability.py','solver_v1/silicon_initiation_research.py',
        'solver_v1/silicon_oxide_research.py','solver_v1/silicon_charge_dynamics.py',
        'solver_v1/test_silicon_initiation_probability.py','solver_v1/test_silicon_initiation_research.py',
        'solver_v1/test_silicon_charge_dynamics.py']
    record=dict(tests_passed=int(matches[0][0]),test_elapsed_seconds=float(matches[0][1]),
        command='python -m pytest solver_v1/test_silicon_initiation_research.py solver_v1/test_silicon_initiation_probability.py solver_v1/test_silicon_charge_dynamics.py -q',
        log_summary=matches[0][0]+' passed in '+matches[0][1]+'s',
        source_files=[dict(path=p,sha256=hashlib.sha256((root/p).read_bytes()).hexdigest()) for p in files],
        scope='geometry, work/force consistency, composite stress, first passage, frozen-mixture, time-unit-invariant conservation and charge-dynamics regressions',
        new_atomic_calculations=0,material_approved=False,kinetic_clock_calibrated=False)
    args.output.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(record,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--test-log',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--expected-tests',type=int,default=41)
    main(parser.parse_args())
