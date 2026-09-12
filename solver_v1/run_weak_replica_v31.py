"""Predeclared weak-drive independent-initialization reference MD campaign."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
from .run_low_frequency_forcing_v29 import sha
from .collective_forcing import harmonic_response
from .phase_resolution import null_lockin_windows
from .impedance_mobility import scalar_impedance,inverse_disk_bounds
from .run_vector_registry_audit import save_json
from .run_low_stress_cyclic_diagnostic import write_csv
from .low_frequency_mobility import multitaper_spectrum,band_integral_proxy


def internal_work_diagnostics(thermo,work=None):
    """One CSV schema; absent unforced work is unavailable, not zero loss."""
    keys=('max_work_residual_eV','rms_work_residual_eV','work_eV','internal_energy_change_eV')
    if work is None:
        return dict.fromkeys(keys)
    energy=np.asarray(thermo)[:,3]-np.asarray(thermo)[0,3]
    work=np.asarray(work)
    if work.shape!=energy.shape or not np.all(np.isfinite(energy+work)):
        raise ValueError('matching finite native-work and energy records required')
    residual=energy-work
    return dict(zip(keys,(float(np.max(abs(residual))),float(np.sqrt(np.mean(residual**2))),
                          float(work[-1]),float(energy[-1]))))


def definition(source):
    forces=source['force_scale_eV_A']
    return dict(seeds=[35461,49277],equilibration_ps=100.,dt_ps=.00125,frame_ps=.025,
        duration_ps=2000.,exclude_ps=100.,frequency_per_ps=.2,force_scale_eV_A=forces,
        potential_sha256=source['potential_sha256'],temperature_target_K=300.,
        initialization='distinct fixed Gaussian velocity seeds,100ps NVT, then unthermostatted NVE',
        statistical_independence_certified=False,force_fraction=1.,
        planned_sample_is_screen_not_precision_certificate=True,production_clock_calibrated=False)


def run(study,potential):
    study,potential=Path(study),Path(potential)
    if study.exists():raise FileExistsError('fresh campaign required; raw partial states preserved')
    source=json.loads(Path('results/low_frequency_forcing_v29/summary.json').read_bytes())['protocol']
    p=definition(source)
    if sha(potential)!=p['potential_sha256']:raise ValueError('source changed')
    study.mkdir(parents=True);save_json(study/'protocol.json',p)
    base=[sys.executable,'-X','utf8','-m','solver_v1.run_reference_thermostat_md',
          '--potential',str(potential),'--ensemble','nve','--repeats','6','--threads','1',
          '--dt-ps',str(p['dt_ps'])]
    def launch(name,extra):
        with (study/(name+'.log')).open('x') as log:
            subprocess.run(base+['--out',str(study/name)]+extra,stdout=log,stderr=subprocess.STDOUT,check=True)
        print('Completed '+name,flush=True)
    def initialize(seed):
        launch(f'seed{seed}_init',['--duration-ps','.025','--equilibration-seed',str(seed),
               '--equilibration-ps','100'])
    with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(initialize,p['seeds']))
    cases=[]
    for seed in p['seeds']:
        restart=study/f'seed{seed}_init'/'equilibrated.restart'
        cases.append(dict(name=f'seed{seed}_null',seed=seed,axis=None,sign=0,restart_sha256=sha(restart)))
        for axis in (0,1):
            for sign in (-1,1):
                cases.append(dict(name=f'seed{seed}_axis{axis}_sign{sign:+d}',seed=seed,axis=axis,
                    sign=sign,force_eV_A=sign*p['force_scale_eV_A'][axis],restart_sha256=sha(restart)))
    # Binding saved BEFORE measuring driven/unforced records.
    save_json(study/'cases.json',cases)
    def one(c):
        extra=['--restart',str(study/f'seed{c["seed"]}_init'/'equilibrated.restart'),
               '--duration-ps',str(p['duration_ps'])]
        if c['axis'] is not None:
            extra+=['--drive-axis',str(c['axis']),'--drive-force-eV-A',str(c['force_eV_A']),
                    '--drive-frequency-per-ps',str(p['frequency_per_ps']),'--timestep-work']
        launch(c['name'],extra)
    # Each record is a deterministic MD solve, not Monte Carlo fatigue counting.
    with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(one,cases))


def analyze(study,out):
    study,out=Path(study),Path(out)
    if out.exists():raise FileExistsError('fresh report required')
    p=json.loads((study/'protocol.json').read_bytes());cases=json.loads((study/'cases.json').read_bytes())
    fits=[];checks=[];null=[];pairs=[];fdt=[]
    for c in cases:
        d=study/c['name'];m=json.loads((d/'summary.json').read_bytes())
        if (not m['completed'] or m['restart_sha256']!=c['restart_sha256']
                or any(m[k]!=p[k] for k in ('dt_ps','frame_ps','duration_ps'))
                or m['ensemble']!='nve'):
            raise ValueError('incomplete or mismatched record')
        expected=None if c['axis'] is None else dict(axis=c['axis'],force_eV_A=c['force_eV_A'],frequency_per_ps=p['frequency_per_ps'])
        if m['conjugate_drive']!=expected:raise ValueError('drive mismatch')
        with np.load(d/'plane_coordinates.npz') as z:
            t=z['time_seconds']*1e12;q=z['coordinates_m']*1e10;thermo=z['thermo'];mask=t>=100
            row=dict(case=c['name'],seed=c['seed'],axis=c['axis'],trajectory_sha256=sha(d/'plane_coordinates.npz'),
                mean_temperature_K=float(thermo[mask,0].mean()),
                first100_K=float(thermo[t<100,0].mean()),last100_K=float(thermo[t>=1900,0].mean()))
            row.update(internal_work_diagnostics(thermo,
                z['internal_step_work_eV'] if c['axis'] is not None else None))
            checks.append(row)
            for parts in (1,2):
                width=1900/parts
                if c['axis'] is None:
                    for axis in (0,1):
                        for harmonic in (1,2,3):
                            null.extend(dict(seed=c['seed'],axis=axis,parts=parts,harmonic=harmonic,**r)
                                for r in null_lockin_windows(t[mask],q[mask,:,axis],harmonic*.2,width))
                    for block in range(parts):
                        keep=(t>=100+block*width)&(t<100+(block+1)*width)
                        temperature=float(thermo[keep,0].mean())
                        for nw in (2.,3.):
                            spec=multitaper_spectrum(q[keep]*1e-10,.025e-12,nw=nw)
                            for band in (.005,.01,.02):
                                proxy=band_integral_proxy(spec,(.2-band)*1e12,(.2+band)*1e12)
                                for axis in (0,1):
                                    loss=2*np.pi*.2e12*proxy['integral_proxy_m2_seconds'][axis,axis]/(1.380649e-23*temperature)*16.02176634
                                    fdt.append(dict(seed=c['seed'],axis=axis,parts=parts,block=block,nw=nw,
                                        half_band_cycles_ps=band,temperature_K=temperature,imag_A2_eV=float(-loss)))
                else:
                    for block in range(parts):
                        keep=(t>=100+block*width)&(t<100+(block+1)*width)
                        r=dict(**c,parts=parts,block=block,**harmonic_response(t[keep],q[keep,0,c['axis']],.2,c['force_eV_A']))
                        for harmonic in (2,3):
                            h=harmonic_response(t[keep],q[keep,0,c['axis']],harmonic*.2,1.)
                            r.update({f'h{harmonic}_real_A':h['real_A2_eV'],f'h{harmonic}_imag_A':h['imag_A2_eV']})
                        fits.append(r)
    for seed in p['seeds']:
        for axis in (0,1):
            for parts in (1,2):
                for block in range(parts):
                    rr=[r for r in fits if (r['seed'],r['axis'],r['parts'],r['block'])==(seed,axis,parts,block)]
                    if len(rr)!=2:raise ValueError('missing signed partner')
                    chi=sum(complex(r['real_A2_eV'],r['imag_A2_eV']) for r in rr)/2
                    nn=[r for r in null if (r['seed'],r['axis'],r['parts'],r['harmonic'])==(seed,axis,parts,1)]
                    radius=max(np.hypot(r['real_A'],r['imag_A']) for r in nn)/p['force_scale_eV_A'][axis]
                    imag_floor=max(abs(r['imag_A']) for r in nn)/p['force_scale_eV_A'][axis]
                    r=dict(seed=seed,axis=axis,parts=parts,block=block,real_A2_eV=chi.real,imag_A2_eV=chi.imag,
                        phase_deg=float(np.degrees(np.angle(chi))),complex_null=radius,imag_null=imag_floor,
                        loss_to_null=float(-chi.imag/imag_floor),**scalar_impedance(chi,.2))
                    for h in (2,3):
                        vals={x['sign']:complex(x[f'h{h}_real_A'],x[f'h{h}_imag_A']) for x in rr}
                        z=(vals[1]+vals[-1])/2 if h==2 else (vals[1]-vals[-1])/2
                        hn=[x for x in null if (x['seed'],x['axis'],x['parts'],x['harmonic'])==(seed,axis,parts,h)]
                        bound=max(np.hypot(x['real_A'],x['imag_A']) for x in hn)
                        r.update({f'parity_h{h}_A':abs(z),f'parity_h{h}_to_null':abs(z)/bound})
                    pairs.append(r)
    aggregate=[]
    for axis in (0,1):
        rr=[r for r in pairs if r['axis']==axis and r['parts']==1]
        zz=[complex(r['real_A2_eV'],r['imag_A2_eV']) for r in rr];chi=sum(zz)/len(zz)
        radius=max(r['complex_null'] for r in rr)+max(abs(z-chi) for z in zz)
        aggregate.append(dict(axis=axis,real_A2_eV=chi.real,imag_A2_eV=chi.imag,
            phase_deg=float(np.degrees(np.angle(chi))),empirical_radius=radius,
            **inverse_disk_bounds(chi,radius,.2),zero_frequency_certified=False))
    out.mkdir(parents=True)
    for name,rows in [('responses',fits),('checks',checks),('null',null),('pairs',pairs),('aggregate',aggregate),('fdt',fdt)]:
        write_csv(out/(name+'.csv'),rows)
    save_json(out/'summary.json',dict(completed=True,protocol=p,checks=checks,pairs=pairs,aggregate=aggregate,
        statistical_independence_certified=False,production_clock_calibrated=False))
    # save_json normalizes numpy scalar types, including inverse-disk booleans.
    print(json.dumps(json.loads((out/'summary.json').read_bytes())['aggregate'],indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(__doc__);p.add_argument('action',choices=('run','analyze'))
    p.add_argument('--study',type=Path,required=True);p.add_argument('--potential',type=Path);p.add_argument('--out',type=Path)
    a=p.parse_args()
    if a.action=='run':run(a.study,a.potential)
    else:analyze(a.study,a.out)
