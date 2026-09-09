"""Unit/protocol audit of actual public Al wire relaxation measurements.

This is validation data, not fatigue fitting or a cell-mobility calibration.
Source specimen area converts measured load to stress. It is neither atomic
interface area nor the statistical correlation area, and never enters W/PDE.
Downloaded Mathematica notebooks are NOT executed.
"""
import argparse
import hashlib
import io
import itertools
import json
from pathlib import Path
import re
import zipfile

import numpy as np

from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json


EXPECTED_COLUMNS = ('Time[s]', 'Motordispl[µm]', 'Force[mN]', 'EngStrain',
                    'EngStress[MPa]', 'Shearstrain', 'Shearstress [MPa]')


def parse_relaxation(raw):
    lines = [s.strip() for s in raw.decode('cp1252').splitlines() if s.strip()]
    columns = tuple(lines[0].split('\t'))
    if columns != EXPECTED_COLUMNS:
        raise ValueError(f'unrecognized source column units: {columns}')
    values = np.loadtxt(io.StringIO('\n'.join(lines[1:])), delimiter='\t', ndmin=2)
    if values.shape[1] < 7 or len(values) < 2 or np.any(~np.isfinite(values)):
        raise ValueError('finite rows with the seven explicitly named measurements required')
    # Some original files contain two unnamed extra columns and repeated DAQ
    # timestamps. Retain all numbers and expose these facts; do not resample,
    # average duplicates or invent labels for the unnamed columns. We compute
    # no time derivative at repeated timestamps.
    if np.any(np.diff(values[:, 0]) < 0) or values[-1, 0] == values[0, 0]:
        raise ValueError('physical measurement time cannot decrease or remain constant')
    return values


def maximum_fcc_schmid(direction):
    """Geometric maximum over the 12 {111}<110> systems, not a fitted chi.

    Measured orientations are rounded in the paper. Agreement with this
    maximum does not prove which system was active or its source geometry.
    """
    e = np.asarray(direction, float)
    if e.shape != (3,) or not np.all(np.isfinite(e)) or np.linalg.norm(e) == 0:
        raise ValueError('finite nonzero loading direction required')
    e = e/np.linalg.norm(e)
    normals = np.array([(1, a, b) for a in (-1, 1) for b in (-1, 1)])
    directions = sorted(set(itertools.permutations((1, 1, 0))) |
                        set(itertools.permutations((1, -1, 0))))
    candidates = []
    systems = set()
    for n in normals:
        for m in directions:
            m = np.asarray(m)
            if n@m != 0:
                continue
            first = m[np.flatnonzero(m)[0]]
            canonical = tuple((m*first).tolist())
            if (tuple(n), canonical) in systems:
                continue
            systems.add((tuple(n), canonical))
            nn = n/np.linalg.norm(n); mm = m/np.linalg.norm(m)
            candidates.append((abs(float(e@nn*(e@mm))), n.tolist(), m.tolist(), float((e@nn)**2)))
    if len(systems) != 12:
        raise ArithmeticError('incorrect FCC slip-family enumeration')
    factor, normal, slip, normal_factor = max(candidates, key=lambda row: row[0])
    return dict(schmid_factor=factor, slip_plane_normal=normal,
                slip_direction=slip, normal_traction_factor=normal_factor)


def summarize_trace(values, *, diameter_um, orientation):
    if not np.isfinite(diameter_um) or diameter_um <= 0:
        raise ValueError('positive measured specimen diameter required')
    q = np.asarray(values)
    cross_section_um2 = np.pi/4*diameter_um**2
    from_load = q[:, 2]*1000/cross_section_um2
    factor = maximum_fcc_schmid(orientation)
    # No noise clipping or identification of total recorded shear as plastic.
    nonzero = abs(q[:, 4]) > 0
    measured_schmid = np.median(q[nonzero, 6]/q[nonzero, 4]) if np.any(nonzero) else None
    first = q[q[:, 0] <= q[0, 0]+1.]
    last = q[q[:, 0] >= q[-1, 0]-1.]
    positive=(q[:,2]>0)&(q[:,4]>0)
    inferred=np.sqrt(4000*q[positive,2]/(np.pi*q[positive,4]))
    return dict(samples=len(q), duration_seconds=float(q[-1, 0]-q[0, 0]),
        median_sample_interval_seconds=float(np.median(np.diff(q[:, 0]))),
        duplicate_timestamps=int(np.sum(np.diff(q[:, 0]) == 0)),
        unnamed_extra_columns=q.shape[1]-7,
        unnamed_columns_max_absolute=None if q.shape[1]==7 else float(np.max(abs(q[:,7:]))),
        axial_stress_min_MPa=float(q[:, 4].min()), axial_stress_max_MPa=float(q[:, 4].max()),
        shear_stress_min_MPa=float(q[:, 6].min()), shear_stress_max_MPa=float(q[:, 6].max()),
        first_second_mean_shear_MPa=float(first[:, 6].mean()),
        last_second_mean_shear_MPa=float(last[:, 6].mean()),
        shear_relaxation_MPa=float(first[:, 6].mean()-last[:, 6].mean()),
        strain_recorded_initial=float(q[0, 3]), strain_recorded_final=float(q[-1, 3]),
        max_stress_from_load_discrepancy_MPa=float(np.max(abs(from_load-q[:, 4]))),
        diameter_used_for_unit_check_um=float(diameter_um),
        implied_diameter_median_um=float(np.median(inferred)) if len(inferred) else None,
        implied_diameter_range_um=[float(inferred.min()),float(inferred.max())] if len(inferred) else None,
        measured_shear_to_axial_ratio=None if measured_schmid is None else float(measured_schmid),
        geometric_max_schmid=factor['schmid_factor'],
        source_pin_spacing_m=None, proof_yield_MPa=None,
        collective_coordinate_mobility=None, model_physical_clock=None)


def audit(source_directory, output):
    source_directory, output = Path(source_directory), Path(output)
    if output.exists():
        raise FileExistsError('fresh report output required')
    metadata = json.loads((source_directory/'manifest.json').read_bytes())
    raw = (source_directory/'supplementary.zip').read_bytes()
    if hashlib.sha256(raw).hexdigest() != metadata['supplementary_sha256']:
        raise ValueError('supplementary source hash mismatch')
    outer = zipfile.ZipFile(io.BytesIO(raw)); inner = zipfile.ZipFile(io.BytesIO(outer.read('mmc2.zip')))
    specimens = {}
    for table in metadata['tables']:
        if table['label'] not in ('Table 1', 'Table 2'):
            continue
        for name, material, count, diameter, orientation, gauge in table['rows'][1:]:
            specimens[name.casefold()] = dict(sample=name, material=material,
                condition='as_cast' if table['label']=='Table 1' else 'annealed_500C_2h',
                declared_relaxations=int(count), diameter_um=float(diameter),
                orientation_cubic=[int(v) for v in re.findall(r'-?\d+', orientation)],
                gauge_length_um=float(gauge), source_table=table['label'])
    traces=[]; failures=[]; plot_data=[]
    for name in sorted(inner.namelist()):
        if not name.lower().endswith('.txt') or 'relaxcycle' not in name.lower():
            continue
        matches=[s for key,s in specimens.items() if key in [p.casefold() for p in name.split('/')]]
        try:
            if len(matches)!=1:
                raise ValueError('source file not uniquely matched to the published specimen table')
            specimen=matches[0]; values=parse_relaxation(inner.read(name))
            stats=summarize_trace(values,diameter_um=specimen['diameter_um'],orientation=specimen['orientation_cubic'])
            traces.append(dict(sample=specimen['sample'],condition=specimen['condition'],source_member=name,
                source_marked_not_considered=any(part.casefold()=='not considered' for part in name.split('/')),
                source_member_sha256=hashlib.sha256(inner.read(name)).hexdigest(),**stats))
            if specimen['sample']=='Al_20_r_10':
                # A few public data points for a portable validation plot. This
                # is explicit stride20 visualization, not a refitted curve.
                for row in values[::20]:
                    plot_data.append(dict(sample=specimen['sample'],source_member=name,
                        time_seconds=row[0]-values[0,0],axial_stress_MPa=row[4],shear_stress_MPa=row[6]))
        except ValueError as exc:
            failures.append(dict(source_member=name,reason=str(exc)))
    write_csv(output/'specimens.csv',list(specimens.values()))
    write_csv(output/'relaxation_summary.csv',traces)
    write_csv(output/'selected_measured_traces.csv',plot_data)
    counts=[dict(sample=s['sample'],declared=s['declared_relaxations'],
                 parsed=sum(r['sample']==s['sample'] for r in traces),
                 marked_not_considered=sum(r['sample']==s['sample'] and r['source_marked_not_considered'] for r in traces)) for s in specimens.values()]
    save_json(output/'source_audit.json',dict(completed=True,doi=metadata['doi'],
        source_sha256=metadata['supplementary_sha256'],source_bytes=len(raw),
        specimens=len(specimens),parsed_traces=len(traces),parse_failures=failures,counts=counts,
        source_marked_not_considered=sum(r['source_marked_not_considered'] for r in traces),
        source_exclusions_preserved=True,
        declared_count_mismatch_not_silently_repaired=True,
        specimen_data_scope='room-temperature monotonic tension with separate60s relaxation holds; not cyclic fatigue',
        units='measured mN,um,MPa,s; not a model-time conversion',
        yield_criterion_not_established=True,source_geometry_unknown=True,
        a_s_cell_mobility_available=False,production_physical_seconds_available=False))
    return traces,failures


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source',type=Path);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();rows,failures=audit(args.source,args.out)
    print(f'parsed physical measurements: {len(rows)}; explicit parse failures: {len(failures)}')
