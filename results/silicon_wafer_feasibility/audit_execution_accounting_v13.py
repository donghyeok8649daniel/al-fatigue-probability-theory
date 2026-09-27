"""Count reported new work once, keeping reuse and incomplete stages separate.

This is a reporting audit of completed, independently checked result records.
It makes no potential calls and does not infer scientific completion from time.
"""
import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def main(args):
    root = args.results.resolve()
    sources = {}
    def read(relative):
        path = root / relative
        raw = path.read_bytes()
        sources[relative] = hashlib.sha256(raw).hexdigest()
        return json.loads(raw)
    rows = []
    def add(stage, forward, ad=0, reused_rows=0, reused_frames=0, note='', elapsed=None):
        counts = (forward, ad, reused_rows, reused_frames)
        if any(type(v) is not int or v < 0 for v in counts):
            raise ValueError('Non-integer or negative counter: ' + stage)
        rows.append(dict(stage=stage, status='completed', new_forward_evaluations=forward,
            new_autograd_rows=ad, reused_hessian_rows=reused_rows,
            reused_prediction_frames=reused_frames, note=note,
            reported_numerical_elapsed_s=elapsed, elapsed_is_not_cpu_time=True))
    for state in ('loading8', 'return8', 'loading10', 'return_zero'):
        name = 'dense_' + state
        if not (root / name / 'summary.json').exists():
            if args.final:
                raise ValueError('Final accounting requires the zero-Hessian outcome')
            rows.append(dict(stage=name, status='pending_not_counted',
                new_forward_evaluations=None, new_autograd_rows=None))
            continue
        item = read(name + '/summary.json')
        if not item['complete'] or item['completed_rows'] != item['dimension'] or item['independent_checked_directions'] != 6:
            raise ValueError('Partial Hessian needs explicit manual accounting: ' + name)
        reused = item['reused_vjp_rows']
        checks = item['resume_validation_vjp_rows']
        if item['autograd_vjp_rows'] != item['dimension'] - reused + checks:
            raise ValueError('Hessian new/reused row mismatch: ' + name)
        if state == 'return_zero':
            if not read('dense_return_zero_replay/summary.json')['complete']:
                raise ValueError('Zero-Hessian replay incomplete')
        add(name, item['new_graph_evaluations'], item['autograd_vjp_rows'], reused,
            note=f'Resume validation {checks} rows are included in new AD, not added twice.',
            elapsed=item['elapsed_s'])
    old = read('dense_force100/reuse_provenance.json')
    previous = read('dense_force100/summary.json')
    if old['new_atomic_evaluations'] != 0 or old['new_hessian_rows'] != 0:
        raise ValueError('Old 100 MPa result was not pure reuse')
    add('dense_force100', 0, reused_rows=previous['dimension'],
        note='Entire v12 result reused. Old summary counters and elapsed excluded.')
    grip = read('grip_augmented/summary.json')
    if not grip['complete']:
        raise ValueError('Grip augmentation incomplete')
    add('grip_augmented', grip['new_graph_evaluations'], grip['new_grip_autograd_rows'],
        grip['reused_fixed_grip_rows'], note='Rows reused from current dense runs; no extra AD.',
        elapsed=grip['elapsed_s'])
    mpa = read('mpa_comparison/summary.json')
    full = read('mpa_qe_full/summary.json')
    if not mpa['complete'] or not full['complete'] or full['baseline_calls'] != 0:
        raise ValueError('Model comparison incomplete or baseline reevaluated')
    if full['new_model_calls'] != full['new_prediction_frames'] + len(full['standard_controls']):
        raise ValueError('Full-QE control accounting mismatch')
    if full['frames'] != full['new_prediction_frames'] + full['reused_model_frames']:
        raise ValueError('Full-QE frame accounting mismatch')
    add('mpa_comparison', mpa['new_model_calls'],
        note='199 predictions plus 9 standard-calculator controls.', elapsed=mpa['elapsed_s'])
    add('mpa_qe_full', full['new_model_calls'], reused_frames=full['reused_model_frames'],
        note='Only complement predictions and new controls count again.', elapsed=full['elapsed_s'])
    zero = read('return_zero_force/summary.json')
    if not zero['complete'] or zero['force_only_calls'] != sum(x['force_calls'] for x in zero['states']):
        raise ValueError('Zero-force call accounting mismatch')
    add('return_zero_force', zero['force_only_calls'], elapsed=zero['elapsed_s'])
    probes = read('anharmonic_probes/summary.json')
    if not probes['complete'] or probes['new_model_calls'] != probes['actual_points'] + len(probes['completed_source_replays']):
        raise ValueError('Nonlinear probe accounting mismatch')
    add('anharmonic_probes', probes['new_model_calls'],
        note='Source replays included in forward count.', elapsed=probes['elapsed_s'])
    for folder in sorted(root.glob('large_oxide_321_1200*')):
        name = folder.name
        if (folder / 'summary.json').exists():
            oxide = read(name + '/summary.json')
            if not oxide['complete']:
                raise ValueError('Large-oxide partial summary requires explicit accounting')
            add(name, oxide['force_only_calls'] + len(oxide['standard_controls']),
                note='New predictions and standard controls; previous allocation failure preserved.',
                elapsed=oxide['elapsed_s'])
        elif (folder / 'partial.json').exists():
            item = read(name + '/partial.json')
            if item != dict(reason='less than 12 GB available before model allocation', completed=0):
                raise ValueError('Partial model work needs explicit accounting: ' + name)
            add(name, 0, note='Failed before model allocation; no new potential evaluation.')
            rows[-1]['status'] = 'failed_before_allocation'
        else:
            raise ValueError('Large-oxide result is still being written: ' + name)
    counted = [r for r in rows if r['new_forward_evaluations'] is not None]
    result = dict(recorded_utc=datetime.now(timezone.utc).isoformat(),
        complete_accounting=all(r['status'] != 'pending_not_counted' for r in rows),
        requested_final=args.final, stages=rows,
        completed_or_failed_stage_new_forward=sum(r['new_forward_evaluations'] for r in counted),
        completed_or_failed_stage_new_autograd_rows=sum(r['new_autograd_rows'] for r in counted),
        actual_new_DFT=0, actual_new_MD=0,
        source_sha256=sources,
        scope='Reported executed calls, not independent specimens, physical time, or model adoption. Pending work excluded.',
        postprocessing_new_potential_calls=0,
        new_model_calls_in_this_audit=0, new_pytest_runs_in_this_audit=0)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print(json.dumps({k: v for k, v in result.items() if k not in ('stages', 'source_sha256')}, ensure_ascii=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--results', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--final', action='store_true')
    main(parser.parse_args())
