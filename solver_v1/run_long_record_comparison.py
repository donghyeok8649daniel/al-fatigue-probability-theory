"""Completed 1/2/5ns prefix and paired-dt mobility comparison; no certification."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json


def sampling_check(dense_source,coarse_source,output):
    """Same physical record, saved every5 versus25fs; not a new ensemble."""
    from .low_frequency_mobility import multitaper_spectrum,band_integral_proxy
    from .zero_frequency_kinetics import formal_zero_frequency_response
    dense_source,coarse_source,output=map(Path,(dense_source,coarse_source,output))
    if output.exists():raise FileExistsError('fresh output required')
    dm=json.loads((dense_source/'summary.json').read_bytes())
    cm=json.loads((coarse_source/'summary.json').read_bytes())
    if not dm['completed'] or not cm['completed'] or dm['restart_sha256']!=cm['restart_sha256']:
        raise ValueError('completed shared-initial-state records required')
    ratio=cm['frame_ps']/dm['frame_ps'];stride=int(round(ratio))
    if abs(ratio-stride)>1e-10 or stride<2:raise ValueError('integral sampling ratio required')
    with np.load(dense_source/'plane_coordinates.npz',allow_pickle=False) as z:
        excluded=int(round(25/dm['frame_ps']));q=z['coordinates_m'][excluded:-1]
        T=float(z['thermo'][excluded:-1,0].mean())
    if len(q)%stride:raise ValueError('matched FFT record durations required')
    coarse=q[::stride]
    with np.load(coarse_source/'plane_coordinates.npz',allow_pickle=False) as z:
        start=int(round(25/cm['frame_ps']));old=z['coordinates_m'][start:start+len(coarse)]
    records=[]
    for nw in (2.,3.):
        estimates=[multitaper_spectrum(x,frame*1e-12,nw=nw,tapers=2)
                   for x,frame in ((q,dm['frame_ps']),(coarse,cm['frame_ps']))]
        for lower,upper in ((.02,.04),(.04,.08),(.08,.16)):
            M=[]
            for e in estimates:
                b=band_integral_proxy(e,lower*1e12,upper*1e12)
                M.append(formal_zero_frequency_response(e['covariance_m2'],b['integral_proxy_m2_seconds'],1.380649e-23*T)['mobility_m2_per_J_second'])
            records.append(dict(nw=nw,lower_cycle_THz=lower,upper_cycle_THz=upper,
                dense_mobility=M[0],decimated_mobility=M[1],
                relative_diagonal_change=abs(np.diag(M[1])/np.diag(M[0])-1)))
    output.mkdir(parents=True)
    save_json(output/'summary.json',dict(completed=True,records=records,
        retained_ps=len(q)*dm['frame_ps'],temperature_K=T,
        matched_old_record_exact_equal=bool(np.array_equal(old,coarse)),
        matched_old_record_max_difference_m=float(np.max(abs(old-coarse))),
        dense_frame_ps=dm['frame_ps'],coarse_frame_ps=cm['frame_ps'],
        production_calibration_available=False))


def run(root,sources,output):
    root,sources,output=map(Path,(root,sources,output))
    if output.exists():raise FileExistsError('fresh output required')
    rows=[];invalid=[]
    for dt,prefix in ((.0025,'long_record_prefix_'),(.005,'long_record_dt5_prefix_')):
        for duration in (1000,2000,5000):
            report=json.loads((root/(prefix+str(duration))/'summary.json').read_bytes())
            if not report['completed']:raise ValueError('completed record study required')
            for r in report['records']:
                if not r['valid_finite_band']:
                    invalid.append(dict(dt_ps=dt,prefix_ps=duration,**r));continue
                M=np.asarray(r['response']['mobility_m2_per_J_second'])
                for axis,label in enumerate(('normal','direct110','transverse')):
                    rows.append(dict(dt_ps=dt,prefix_ps=duration,parts=r['parts'],block=r['block'],
                        nw=r['nw'],lower_cycle_THz=r['lower_cycle_THz'],upper_cycle_THz=r['upper_cycle_THz'],
                        axis=label,mobility_m2_per_J_s=M[axis,axis],
                        normalized_collective_energy_mobility_m2_per_J_s=144*M[axis,axis],
                        temperature_K=r['temperature_K']))
    comparisons=[]
    for band in (.001,.002,.005,.01,.02,.04,.08):
        for axis in ('normal','direct110','transverse'):
            groups=[[r for r in rows if r['dt_ps']==dt and r['prefix_ps']==5000
                     and r['lower_cycle_THz']==band and r['axis']==axis] for dt in (.0025,.005)]
            full=[[r['mobility_m2_per_J_s'] for r in g if r['parts']==1] for g in groups]
            if not all(full):continue
            comparison=dict(lower_cycle_THz=band,axis=axis,
                dt2p5_full_min=min(full[0]),dt2p5_full_max=max(full[0]),
                dt5_full_min=min(full[1]),dt5_full_max=max(full[1]),
                relative_dt_change_of_taper_mean=float(abs(np.mean(full[0])-np.mean(full[1]))/np.mean(full[0])))
            for tag,g in zip(('dt2p5','dt5'),groups):
                values=[r['mobility_m2_per_J_s'] for r in g]
                comparison[tag+'_all_blocks_min']=min(values)
                comparison[tag+'_all_blocks_max']=max(values)
            comparisons.append(comparison)
    old=sources/'nve_N6_serial_1000ps_dt2p5/plane_coordinates.npz'
    long=sources/'nve_N6_serial_5000ps_dt2p5/plane_coordinates.npz'
    thermal=sources/'nve_N6_thermal_1000ps_dt2p5/plane_coordinates.npz'
    with np.load(old,allow_pickle=False) as a,np.load(long,allow_pickle=False) as b,np.load(thermal,allow_pickle=False) as c:
        q=a['coordinates_m'];extended=b['coordinates_m'][:len(q)];observed=c['coordinates_m']
        reproducibility=dict(long_prefix_exact_equal=bool(np.array_equal(q,extended)),
            long_prefix_max_difference_m=float(np.max(abs(q-extended))),
            thermal_observation_exact_equal=bool(np.array_equal(q,observed)),
            thermal_observation_max_difference_m=float(np.max(abs(q-observed))))
    output.mkdir(parents=True)
    write_csv(output/'all_estimates.csv',rows)
    write_csv(output/'paired_dt_summary.csv',comparisons)
    save_json(output/'summary.json',dict(completed=True,comparisons=comparisons,
        invalid_bands=invalid,reproducibility=reproducibility,
        trajectory_sha256={p.parent.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (old,long,thermal)},
        nested_records_not_independent_replicas=True,
        ranges_are_sensitivity_not_confidence=True,production_clock_calibrated=False))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4.6),constrained_layout=True)
    for ax,axis in zip(axes,('normal','direct110')):
        for color,duration in enumerate((1000,2000,5000)):
            for dt,style in ((.0025,'-o'),(.005,'--s')):
                subset=sorted([r for r in rows if r['axis']==axis and r['parts']==1 and r['nw']==3
                    and r['dt_ps']==dt and r['prefix_ps']==duration],key=lambda r:r['lower_cycle_THz'])
                ax.loglog([np.sqrt(r['lower_cycle_THz']*r['upper_cycle_THz']) for r in subset],
                    [r['mobility_m2_per_J_s'] for r in subset],style,color='C'+str(color),
                    label=f'{duration/1000:g} ns, dt={dt*1000:g} fs',markersize=4)
        ax.set_title(axis+' plane-gap PMF');ax.set_xlabel('Band center [cycles/ps, source MD]')
        ax.set_ylabel('Finite-band mobility [m²/(J s)]');ax.grid(alpha=.25)
    axes[0].legend(fontsize=8)
    fig.suptitle('Record-length / timestep study: NW=3 shown; all tapers and blocks saved')
    fig.savefig(output/'record_length_comparison.png',dpi=160);plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root');p.add_argument('--sources',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();run(a.root,a.sources,a.out)
