"""Portable actual-source evidence; no potential fit or production clock."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from .line_drag_reference import load_line_drag_references, source_drag_pa_s
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json


def read_csv(path):
    with path.open(encoding='utf-8',newline='') as stream:return list(csv.DictReader(stream))


def save_figure(fig, path):
    # Stable text, not raster data or a regenerated numerical result.
    import io
    stream=io.StringIO();fig.savefig(stream,format='svg',metadata={'Date':None})
    path.write_text('\n'.join(line.rstrip() for line in stream.getvalue().splitlines())+'\n',encoding='utf-8')
    fig.savefig(path.with_suffix('.png'),dpi=125)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-cache',type=Path,required=True)
    parser.add_argument('--results',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();out=args.out
    if out.exists():raise FileExistsError('fresh report output required')
    out.mkdir(parents=True)
    kinetic=args.results/'crystalline_kinetics';wire=args.results/'wire_relaxation_audited'
    metadata=json.loads((kinetic/'source_metadata.json').read_bytes())
    assets=[]
    for name,expected in [('source_Al99.eam.alloy','eb0f0b204ea40787274efcf8e44d6de2'),
                          ('source_log.lammps','28a50d823517269ca0fef426c17a6874')]:
        raw=(args.source_cache/'zenodo_range_1024'/name).read_bytes()
        assert hashlib.md5(raw).hexdigest()==expected
        assets.append(dict(file=name,sha256=hashlib.sha256(raw).hexdigest(),published_md5_verified=True))
    # Derived coordinates are small; preserve them for offline numerical replay.
    raw=(args.source_cache/'zenodo_range_1024/plane_coordinates.npz').read_bytes()
    assert hashlib.sha256(raw).hexdigest()==metadata['projection_file_sha256']
    (out/'projected_coordinates.npz').write_bytes(raw)
    traces=read_csv(wire/'relaxation_summary.csv');specimens=read_csv(wire/'specimens.csv')
    review=[]
    for specimen in specimens:
        rows=[r for r in traces if r['sample']==specimen['sample']]
        diameter=float(specimen['diameter_um'])
        implied=np.array([float(r['implied_diameter_median_um']) for r in rows])
        review.append(dict(sample=specimen['sample'],source_table_diameter_um=diameter,
            implied_diameter_median_um=float(np.median(implied)),
            implied_diameter_min_um=float(implied.min()),implied_diameter_max_um=float(implied.max()),
            maximum_table_based_stress_discrepancy_MPa=max(float(r['max_stress_from_load_discrepancy_MPa']) for r in rows),
            source_declared_traces=int(specimen['declared_relaxations']),parsed_traces=len(rows),
            source_excluded_traces=sum(r['source_marked_not_considered']=='True' for r in rows)))
    write_csv(out/'experimental_source_consistency.csv',review)
    save_json(out/'evidence_scope.json',dict(completed=True,MD_source_assets=assets,
        projected_coordinate_sha256=hashlib.sha256(raw).hexdigest(),
        trajectory_full_checksum_verified=False,
        kinetic_source_temperature_K=300.,static_fit_reference_temperature_K=0.,
        raw_experimental_traces=len(traces),source_marked_not_considered=sum(r['source_marked_not_considered']=='True' for r in traces),
        duplicate_timestamps=sum(int(r['duplicate_timestamps']) for r in traces),
        files_with_unnamed_columns=sum(int(r['unnamed_extra_columns'])>0 for r in traces),
        source_notebooks_executed=False,source_discrepancies_repaired=False,
        physics_scope='actual source measurements/projections, not current LJ/Bessel model predictions',
        accepted_a_s_mobility=False,production_seconds_available=False,actual_yield_reproduced=False))

    line_rows=[]
    for key,source in load_line_drag_references().items():
        temps=(sorted(set([source['temperature_min_K'],300.,source['temperature_max_K']]))
               if source['quantity']=='B_line/T' else [source['temperature_K']])
        for T in temps:
            line_rows.append(dict(reference=key,temperature_K=T,B_line_Pa_s=source_drag_pa_s(key,T),
                method=source['method'],character=source['character'],doi=source['doi'],cell_mobility=None))
    write_csv(out/'physical_line_drag_references.csv',line_rows)

    modes=read_csv(kinetic/'correlation_modes.csv')
    fig,ax=plt.subplots(figsize=(7.5,4.5),layout='constrained')
    for name in ('first_half','full','eight_blocks','stride2'):
        rows=[r for r in modes if r['study']==name]
        t=np.array([float(r['lag_seconds'])*1e12 for r in rows])
        v=np.array([float(r['minimum_whitened_correlation_eigenvalue']) for r in rows])
        ax.plot(t,v,label=name)
        if name=='full':
            floor=np.array([float(r['empirical_error_floor']) for r in rows])
            ax.fill_between(t,v-floor,v+floor,alpha=.12,label='full empirical block+antisymmetry range')
    ax.axhline(0.,color='black',lw=.7)
    ax.set(title='Actual 300 K Al MD: not an accepted overdamped a/s clock',
           xlabel='Physical lag of source MD [ps]',ylabel='Minimum whitened correlation eigenvalue')
    ax.legend(fontsize=8);ax.grid(alpha=.2)
    save_figure(fig,out/'source_md_correlation.svg')

    values=read_csv(wire/'selected_measured_traces.csv')
    fig,ax=plt.subplots(figsize=(7.5,4.5),layout='constrained')
    for member in sorted(set(r['source_member'] for r in values)):
        rows=[r for r in values if r['source_member']==member]
        ax.plot([float(r['time_seconds']) for r in rows],
                [float(r['shear_stress_MPa']) for r in rows],label=Path(member).stem[-4:])
    ax.set(title='Measured Al_20_r_10: separate relaxation holds, NOT fatigue cycles',
           xlabel='Physical time within the experimental hold [s]',ylabel='Recorded resolved shear stress [MPa]')
    ax.legend(title='Source hold',ncol=4,fontsize=8);ax.grid(alpha=.2)
    save_figure(fig,out/'source_wire_relaxations.svg')
    print('actual-source evidence report complete',flush=True)


if __name__=='__main__':main()
