"""Replay v31 native-step work against phase, with raw-source hash binding.

This performs no new MD and does not write production calibration files.
"""
import argparse
import csv
import json
from pathlib import Path
import time
import numpy as np
from .collective_forcing import harmonic_response
from .work_phase_audit import work_phase, duration_sensitivity
from .run_low_frequency_forcing_v29 import sha
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json


def run(study, previous, out):
    study, previous, out = map(Path, (study, previous, out))
    if out.exists():
        raise FileExistsError('fresh report directory required')
    started = time.perf_counter()
    protocol = json.loads((study/'protocol.json').read_bytes())
    cases = json.loads((study/'cases.json').read_bytes())
    old = json.loads((previous/'summary.json').read_bytes())
    if not old['completed'] or old['protocol'] != protocol:
        raise ValueError('completed matching prior audit required')
    hashes = {row['case']: row['trajectory_sha256'] for row in old['checks']}
    frequency = protocol['frequency_per_ps']
    rows, nulls, elapsed = [], [], []
    for case in cases:
        source = study/case['name']
        metadata = json.loads((source/'summary.json').read_bytes())
        raw = source/'plane_coordinates.npz'
        if (not metadata['completed'] or sha(raw) != hashes[case['name']]
                or metadata['restart_sha256'] != case['restart_sha256']):
            raise ValueError('completed unchanged trajectory binding required')
        elapsed.append(metadata['elapsed_seconds']/metadata['duration_ps'])
        with np.load(raw, allow_pickle=False) as data:
            t = data['time_seconds']*1e12
            q = data['coordinates_m']*1e10
            w = data['internal_step_work_eV'] if case['axis'] is not None else None
            temp = data['thermo'][:, 0]
            for parts in (1, 2, 4):
                width = (protocol['duration_ps']-protocol['exclude_ps'])/parts
                for block in range(parts):
                    lo = protocol['exclude_ps']+block*width
                    hi = lo+width
                    i, j = (int(round(x/protocol['frame_ps'])) for x in (lo, hi))
                    ts = t[i:j+1]
                    if not np.allclose([ts[0],ts[-1]], [lo,hi], atol=1e-8,rtol=0):
                        raise ValueError('missing block endpoints')
                    if case['axis'] is None:
                        for axis in (0,1):
                            for plane in range(q.shape[1]):
                                fit = harmonic_response(ts, q[i:j+1,plane,axis],frequency,
                                                        protocol['force_scale_eV_A'][axis])
                                nulls.append(dict(seed=case['seed'], axis=axis, parts=parts,
                                    block=block, plane=plane, duration_ps=width,
                                    signed_loss_A2_eV=-fit['imag_A2_eV']))
                        continue
                    axis, force = case['axis'], case['force_eV_A']
                    if metadata['conjugate_drive'] != dict(axis=axis,force_eV_A=force,
                                                          frequency_per_ps=frequency):
                        raise ValueError('conjugate drive mismatch')
                    base = dict(case=case['name'], seed=case['seed'], axis=axis,
                                sign=case['sign'], parts=parts, block=block)
                    for stride in (1,2,4):
                        qs = q[i:j+1:stride,0,axis]
                        fit = harmonic_response(ts[::stride], qs, frequency, force)
                        result = work_phase(ts[::stride], qs, w[i:j+1:stride],force,frequency)
                        rows.append(dict(**base, save_stride=stride,
                            mean_temperature_K=float(np.mean(temp[i:j+1])),
                            loss_lstsq_A2_eV=-fit['imag_A2_eV'], **result))
    pairs=[]
    for seed in protocol['seeds']:
        for axis in (0,1):
            for parts in (1,2,4):
                floor=max(abs(r['signed_loss_A2_eV']) for r in nulls
                          if r['seed']==seed and r['axis']==axis and r['parts']==parts)
                for block in range(parts):
                    sample=[r for r in rows if r['seed']==seed and r['axis']==axis
                            and r['parts']==parts and r['block']==block and r['save_stride']==1]
                    native=float(np.mean([r['loss_native_A2_eV'] for r in sample]))
                    coordinate=float(np.mean([r['loss_coordinate_A2_eV'] for r in sample]))
                    pairs.append(dict(seed=seed,axis=axis,parts=parts,block=block,
                        loss_native_A2_eV=native,loss_coordinate_A2_eV=coordinate,
                        observed_imaginary_null_A2_eV=floor,
                        native_loss_to_null=native/floor,
                        discrepancy_to_null=abs(native-coordinate)/floor,
                        confidence_interval=False,production_clock_calibrated=False))
    with (previous/'fdt.csv').open(newline='',encoding='utf8') as stream:
        fdt=list(csv.DictReader(stream))
    planning=[]
    for axis in (0,1):
        predicted=min(-float(r['imag_A2_eV']) for r in fdt
                      if int(r['axis'])==axis and int(r['parts'])==1)
        floor=max(r['observed_imaginary_null_A2_eV'] for r in pairs
                  if r['axis']==axis and r['parts']==1)
        required=duration_sensitivity(protocol['duration_ps']-protocol['exclude_ps'],floor,predicted)
        planning.append(dict(axis=axis,independent_fdt_loss_A2_eV=predicted,
            observed_imaginary_null_A2_eV=floor,desired_ratio=3.,
            conditional_record_duration_ps=required,
            estimated_single_record_wall_hours=required*float(np.median(elapsed))/3600,
            scaling_assumption='stationary noise proportional to inverse sqrt duration',
            guaranteed_precision=False,production_clock_calibrated=False))
    out.mkdir(parents=True)
    for name, values in [('work_phase',rows),('null',nulls),('signed_pairs',pairs),('planning',planning)]:
        write_csv(out/(name+'.csv'),values)
    base=[r for r in rows if r['save_stride']==1]
    summary=dict(completed=True,new_MD_run=False,records_reanalyzed=len(cases),
        raw_hashes=hashes,frequency_cycles_per_ps=frequency,
        max_abs_native_coordinate_loss_difference_A2_eV=max(abs(r['loss_discrepancy_A2_eV']) for r in base),
        max_pair_discrepancy_to_null=max(r['discrepancy_to_null'] for r in pairs),
        full_record_pairs=[r for r in pairs if r['parts']==1],planning=planning,
        elapsed_seconds=time.perf_counter()-started,production_clock_calibrated=False,
        conclusion='Native work is an estimator cross-check, not extra independent kinetic data.')
    save_json(out/'summary.json',summary)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(__doc__)
    for name in ('study','previous','out'):parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args()
    run(args.study,args.previous,args.out)
