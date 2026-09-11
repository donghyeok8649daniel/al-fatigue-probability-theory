'''Predeclared real-MD record/cutoff comparisons; not a production clock fit.'''
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from .public_aluminum_kinetics import correlation_audit
from .zero_frequency_kinetics import correlation_integral_audit
from .run_vector_registry_audit import save_json


def study_plan(frames):
    '''Plan without inspecting the measured signal; independent cutoff axis.'''
    if not isinstance(frames, int) or frames < 1024:
        raise ValueError('at least 1024 complete frames required')
    plans = []
    count = 1024
    while count <= frames:
        for lag in (126, 254, 510, 1022):
            if count // 4 >= 2 * (lag + 1):
                plans.append((f'prefix_{count}_lag_{lag}', 0, count, lag, 4))
        count *= 2
    width = frames // 4
    if width >= 1024:
        for block in range(4):
            plans.append((f'quarter_{block}_lag_126', block * width,
                          (block + 1) * width, 126, 4))
    return plans


def audit(source, output):
    source, output = Path(source), Path(output)
    if output.exists():
        raise FileExistsError('fresh refinement output required')
    metadata = json.loads((source / 'source_metadata.json').read_bytes())
    file = source / 'plane_coordinates.npz'
    with np.load(file, allow_pickle=False) as raw:
        time, q = raw['time_seconds'], raw['coordinates_m']
    step = metadata['frame_seconds']
    if (q.shape != (metadata['frames'], metadata['plane_classes'], 3)
            or time.shape != (len(q),) or not np.isfinite(step) or step <= 0
            or np.any(~np.isfinite(time)) or np.any(~np.isfinite(q))
            or not np.allclose(np.diff(time), step, rtol=1e-10, atol=0)):
        raise ValueError('complete finite source record with declared sampling required')
    plan = study_plan(len(q))
    output.mkdir(parents=True)
    save_json(output / 'protocol.json', dict(studies=plan, source_metadata=metadata,
        projection_sha256=hashlib.sha256(file.read_bytes()).hexdigest(),
        protocol='fixed powers-of-two records; independent cutoff and quarter checks',
        blocks_are_confidence_intervals=False, positivity_alone_certifies_mobility=False,
        matched_production_coordinate=False))
    summaries = []
    for name, start, stop, lag, blocks in plan:
        result = correlation_integral_audit(q[start:stop], frame_seconds=step,
                                           max_lag=lag, blocks=blocks)
        save_json(output / f'{name}.json', result)
        late = slice(len(result['cutoff_seconds']) // 2, None)
        summary = dict(study=name, first_frame=start, frames=stop-start,
            record_seconds=float(time[stop-1]-time[start]),
            cutoff_seconds=float(result['cutoff_seconds'][-1]),
            eigenvalues_seconds=result['integral_eigenvalues_seconds'][-1],
            empirical_floor_seconds=float(result['empirical_floor_seconds'][-1]),
            half_cutoff_change_seconds=float(result['half_cutoff_change_seconds'][-1]),
            late_positive_fraction=float(np.mean(result['positive_integral_above_floor'][late])),
            covariance_m2=result['covariance_m2'])
        summaries.append(summary)
        print(json.dumps(summary, default=lambda x: x.tolist()), flush=True)
    short = correlation_audit(q, frame_seconds=step, max_lag=31, blocks=8)
    save_json(output / 'short_lag_markov_check.json', short)
    decision = dict(completed=True, summaries=summaries,
        resolved_short_lag_negative_modes=bool(np.any(short['negative_beyond_empirical_floor'])),
        source_time_is_physical=True, source_is_300K_Al99_NVT_reference=True,
        matched_production_coordinate=False, thermostat_independence_tested=False,
        M_a_phys=None, M_s_phys=None, t0_seconds=None, production_clock_calibrated=False,
        scope='Actual source-plane memory audit; neither a yield nor fatigue prediction')
    save_json(output / 'decision.json', decision)
    return decision


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    audit(args.source, args.out)
