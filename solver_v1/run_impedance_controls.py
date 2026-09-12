"""Raw-record reconstruction and figures for the v28 research calibration."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.signal.windows import dpss
from .impedance_mobility import scalar_low_band_drag
from .run_vector_registry_audit import save_json
from .run_low_stress_cyclic_diagnostic import write_csv


def audit_unavailable(results):
    """Retain every source band rejected as under-resolved; do not hide omissions."""
    results=Path(results); out=results/'unavailable_source_bands.csv'
    if out.exists(): raise FileExistsError('existing unavailable-band audit preserved')
    candidate=json.loads((results/'calibration_candidate.json').read_bytes())
    root=Path(__file__).resolve().parents[1]
    rows=[]; valid=0
    for source in candidate['source_inputs']:
        if 'long_record_' not in source['path']: continue
        raw=(root/source['path']).read_bytes()
        if hashlib.sha256(raw).hexdigest()!=source['sha256']: raise ValueError('changed source')
        data=json.loads(raw)
        for r in data['records']:
            if r['valid_finite_band']:
                valid+=1;continue
            rows.append(dict(source=source['path'],dt_ps=data['source_dt_ps'],
                prefix_ps=data['prefix_ps'],parts=r['parts'],block=r['block'],nw=r['nw'],
                lower_per_ps=r['lower_cycle_THz'],upper_per_ps=r['upper_cycle_THz'],
                reason=r['reason']))
    write_csv(out,rows)
    print(json.dumps(dict(valid_matrix_band_records=valid,unavailable_matrix_band_records=len(rows))))


def verify_raw(source, saved_summary, output):
    source,saved_summary,output=map(Path,(source,saved_summary,output))
    if output.exists(): raise FileExistsError('fresh verification file required')
    expected=json.loads(saved_summary.read_bytes())
    meta=json.loads((source/'summary.json').read_bytes())
    path=source/'plane_coordinates.npz'
    digest=hashlib.sha256(path.read_bytes()).hexdigest()
    if (not meta['completed'] or digest!=expected['source_projection_sha256']
            or meta['dt_ps']!=expected['source_dt_ps'] or meta['duration_ps']!=5000):
        raise ValueError('matched complete 5 ns raw record required')
    with np.load(path,allow_pickle=False) as z:
        t=z['time_seconds']; q=z['coordinates_m']; temperature=z['thermo'][:,0]
    dt=meta['frame_ps']*1e-12
    if not np.allclose(np.diff(t),dt,rtol=1e-10,atol=0): raise ValueError('nonuniform clock')
    skip=int(round(25e-12/dt));q=q[skip:];temperature=float(temperature[skip:].mean())
    q=q-q.mean(axis=0,keepdims=True)
    C=np.mean(q*q,axis=(0,1))
    windows=dpss(len(q),3.,Kmax=2,sym=False)
    # Independent scalar implementation; no cross-matrix inversion or one-sided doubling.
    spectra=np.zeros((len(q)//2+1,3))
    for window in windows:
        fft=np.fft.rfft(q*window[:,None,None],axis=0)
        spectra += dt*np.sum(fft.real**2+fft.imag**2,axis=1)/(np.dot(window,window)*q.shape[1]*2)
    freq=np.fft.rfftfreq(len(q),dt)
    checks=[]
    for r in expected['records']:
        if r['parts']!=1 or r['nw']!=3 or not r['valid_finite_band']: continue
        mask=(freq>=r['lower_cycle_THz']*1e12)&(freq<r['upper_cycle_THz']*1e12)
        K=spectra[mask].mean(axis=0)/2
        for axis in (0,1):
            reference_C=r['covariance_m2'][axis][axis]
            reference_K=r['band']['integral_proxy_m2_seconds'][axis][axis]
            errors=(abs(C[axis]/reference_C-1),abs(K[axis]/reference_K-1),
                    abs(temperature/r['temperature_K']-1))
            # Same saved discrete arithmetic, not a claim of physical accuracy.
            if max(errors)>2e-12: raise ValueError('raw/scalar spectral reconstruction failed')
            checks.append(dict(axis=axis,lower_per_ps=r['lower_cycle_THz'],
                upper_per_ps=r['upper_cycle_THz'],covariance_relative_difference=errors[0],
                spectrum_relative_difference=errors[1],temperature_relative_difference=errors[2],
                scalar_mobility_m2_J_s=1/scalar_low_band_drag(C[axis],K[axis],temperature)))
    save_json(output,dict(completed=True,dt_ps=meta['dt_ps'],source_sha256=digest,
        checks=checks,maximum_reconstruction_relative_difference=max(max(
            r['covariance_relative_difference'],r['spectrum_relative_difference'],
            r['temperature_relative_difference']) for r in checks),
        physical_convergence_certified=False,production_clock_calibrated=False))


def plot(results):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    results=Path(results);out=results/'mobility_controls.png'
    if out.exists(): raise FileExistsError('existing plot preserved')
    with (results/'low_band_residuals.csv').open() as stream: rows=list(csv.DictReader(stream))
    candidate=json.loads((results/'calibration_candidate.json').read_bytes())
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),constrained_layout=True)
    for axis,ax in enumerate(axes):
        rr=[r for r in rows if int(r['axis'])==axis]
        for lo in sorted(set(float(r['lower_per_ps']) for r in rr)):
            group=[r for r in rr if float(r['lower_per_ps'])==lo]
            hi=float(group[0]['upper_per_ps']); x=np.sqrt(lo*hi)
            values=[float(r['scalar_mobility_m2_J_s']) for r in group]
            ax.plot([x,x],[min(values),max(values)],color='0.6',lw=3)
            center=next(r for r in group if r['parts']=='1' and r['prefix_ps']=='5000'
                and r['dt_ps']=='0.0025' and r['nw']=='3.0')
            ax.scatter(x,float(center['scalar_mobility_m2_J_s']),color='C0',s=30,zorder=3)
        ax.axvspan(.02,.16,color='C0',alpha=.07,label='Three fit bands')
        ax.axhline(candidate['fits'][axis]['mobility_m2_per_J_s'],color='C1',ls='--',label='Constant-drag fit')
        ax.set(xscale='log',yscale='log',xlabel='Reference MD band frequency [cycles/ps]',
            ylabel='Scalar mobility proxy [m²/(J s)]',title=('Normal','Direct110 slip')[axis])
        ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Observed control spans, NOT confidence intervals or production calibration')
    fig.savefig(out,dpi=160);plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(__doc__);sub=p.add_subparsers(dest='action',required=True)
    v=sub.add_parser('verify');v.add_argument('--source',type=Path,required=True)
    v.add_argument('--saved-summary',type=Path,required=True);v.add_argument('--out',type=Path,required=True)
    v=sub.add_parser('plot');v.add_argument('--results',type=Path,required=True)
    v=sub.add_parser('unavailable');v.add_argument('--results',type=Path,required=True)
    args=p.parse_args()
    if args.action=='verify': verify_raw(args.source,args.saved_summary,args.out)
    elif args.action=='plot': plot(args.results)
    else: audit_unavailable(args.results)
