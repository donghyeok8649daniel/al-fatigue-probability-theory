"""Predeclared low-frequency reference-MD forcing, not production calibration.

Measured coordinate/energy unchanged. A stronger experimental probe must pass
its amplitude, timestep, thermal and harmonic controls, not just detect lag.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
from scipy.integrate import trapezoid
from .collective_forcing import harmonic_response
from .phase_resolution import null_lockin_windows, loss_envelope
from .impedance_mobility import scalar_impedance, inverse_disk_bounds
from .low_frequency_mobility import multitaper_spectrum, band_integral_proxy
from .run_vector_registry_audit import save_json
from .run_low_stress_cyclic_diagnostic import write_csv


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def make_cases(force_scales):
    if len(force_scales)<2 or not np.all(np.isfinite(force_scales)) or min(force_scales)<=0:
        raise ValueError('positive thermal force scales required')
    cases=[]
    # Full-amplitude pairs first, then amplitude controls, then dt controls.
    for fraction,dt,duration in ((4.,.0025,1000.),(2.,.0025,500.),(4.,.00125,1000.)):
        for axis in (0,1):
            for sign in (-1,1):
                cases.append(dict(name=f'axis{axis}_fraction{fraction:g}_dt{dt:g}_sign{sign:+d}',
                    axis=axis,fraction=fraction,sign=sign,dt_ps=dt,duration_ps=duration,
                    frequency_per_ps=.2,force_eV_A=float(sign*fraction*force_scales[axis])))
    return cases


def prepare(source,restart,potential,study):
    source,restart,potential,study=map(Path,(source,restart,potential,study))
    if study.exists(): raise FileExistsError('fresh study required')
    meta=json.loads((source/'summary.json').read_bytes())
    if (not meta['completed'] or meta['duration_ps']!=5000 or meta['dt_ps']!=.0025
            or meta['repeats']!=6 or meta['ensemble']!='nve'
            or meta['restart_sha256']!=sha(restart)
            or meta['potential_md5']!=hashlib.md5(potential.read_bytes()).hexdigest()):
        raise ValueError('matched completed reference and restart required')
    with np.load(source/'plane_coordinates.npz',allow_pickle=False) as z:
        q=z['coordinates_m'][1000:];temperature=float(z['thermo'][1000:,0].mean())
    rms=np.sqrt(np.mean((q-q.mean(axis=0))**2,axis=(0,1)))
    forces=1.380649e-23*temperature/rms/1.602176634e-9
    study.mkdir(parents=True)
    save_json(study/'protocol.json',dict(cases=make_cases(forces),exclude_ps=100.,
        frame_ps=.025,frequency_per_ps=.2,force_scale_eV_A=forces,rms_m=rms,
        temperature_K=temperature,source_sha256=sha(source/'plane_coordinates.npz'),
        restart_sha256=sha(restart),potential_sha256=sha(potential),
        selection='v28 planning: lower frequency, 4RMS primary, 2RMS amplitude and dt/2 controls',
        protocol_precedes_driven_results=True,
        half_force_duration_is_storage_harmonic_control_not_equal_loss_precision=True,
        production_clock_calibrated=False))


def validate_case(case,protocol,directory):
    meta=json.loads((Path(directory)/'summary.json').read_bytes())
    if (not meta['completed'] or meta['restart_sha256']!=protocol['restart_sha256']
            or meta['duration_ps']!=case['duration_ps'] or meta['dt_ps']!=case['dt_ps']
            or meta['frame_ps']!=protocol['frame_ps'] or meta['ensemble']!='nve'
            or meta['conjugate_drive']!=dict(axis=case['axis'],force_eV_A=case['force_eV_A'],
                frequency_per_ps=case['frequency_per_ps'])):
        raise ValueError('incomplete/mismatched case '+case['name'])
    return meta


def run(study,restart,potential):
    study=Path(study);protocol=json.loads((study/'protocol.json').read_bytes())
    if sha(restart)!=protocol['restart_sha256'] or sha(potential)!=protocol['potential_sha256']:
        raise ValueError('changed restart/potential')
    def one(case):
        out=study/case['name']
        if out.exists():
            # Resume only independently validated COMPLETE cases; never overwrite partial data.
            validate_case(case,protocol,out)
            print('Preserved completed '+case['name'],flush=True);return
        cmd=[sys.executable,'-X','utf8','-m','solver_v1.run_reference_thermostat_md',
            '--potential',str(potential),'--out',str(out),'--ensemble','nve','--repeats','6',
            '--restart',str(restart),'--threads','1','--duration-ps',str(case['duration_ps']),
            '--dt-ps',str(case['dt_ps']),'--drive-axis',str(case['axis']),
            '--drive-force-eV-A',str(case['force_eV_A']),
            '--drive-frequency-per-ps',str(case['frequency_per_ps'])]
        with (study/(case['name']+'.log')).open('x',encoding='utf8') as log:
            subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
        validate_case(case,protocol,out)
        print('Completed '+case['name'],flush=True)
    with ThreadPoolExecutor(max_workers=2) as pool: list(pool.map(one,protocol['cases']))


def analyze(study,source,output):
    study,source,output=map(Path,(study,source,output))
    if output.exists(): raise FileExistsError('fresh analysis directory required')
    protocol=json.loads((study/'protocol.json').read_bytes())
    if sha(source/'plane_coordinates.npz')!=protocol['source_sha256']: raise ValueError('source changed')
    for case in protocol['cases']: validate_case(case,protocol,study/case['name'])
    with np.load(source/'plane_coordinates.npz',allow_pickle=False) as z:
        teq=z['time_seconds'][1000:]*1e12;qeq=z['coordinates_m'][1000:]
    f=protocol['frequency_per_ps'];null=[];fits=[];checks=[];fdt=[]
    durations=sorted(set((c['duration_ps']-protocol['exclude_ps'])/p for c in protocol['cases'] for p in (1,2)))
    for axis in (0,1):
        for duration in durations:
            for harmonic in (1,2,3):
                null.extend(dict(axis=axis,duration_ps=duration,harmonic=harmonic,**r) for r in
                    null_lockin_windows(teq,qeq[:,:,axis]*1e10,harmonic*f,duration))
    for parts in (1,2):
        width=len(qeq)//parts
        for block in range(parts):
            for nw in (2.,3.):
                spec=multitaper_spectrum(qeq[block*width:(block+1)*width],.025e-12,nw=nw)
                for halfband in (.005,.01,.02):
                    proxy=band_integral_proxy(spec,(f-halfband)*1e12,(f+halfband)*1e12)
                    for axis in (0,1):
                        loss=2*np.pi*f*1e12*proxy['integral_proxy_m2_seconds'][axis,axis]/(
                            1.380649e-23*protocol['temperature_K'])*16.02176634
                        fdt.append(dict(axis=axis,parts=parts,block=block,nw=nw,
                            halfband_per_ps=halfband,imag_A2_eV=float(-loss)))
    for case in protocol['cases']:
        directory=study/case['name'];clock=[]
        with (directory/'lammps.log').open(encoding='utf8') as log:
            for line in log:
                cols=line.split()
                if len(cols)!=7 or not cols[0].isdigit():continue
                try: clock.append(list(map(float,cols[:2])))
                except ValueError:continue
        clock=np.asarray(clock)
        if len(clock)==0 or clock[0,1]!=0 or abs(clock[-1,1]-case['duration_ps'])>1e-9:
            raise ValueError('engine clock endpoints mismatch')
        clock_error=float(np.max(abs(clock[:,1]-clock[:,0]*case['dt_ps'])))
        if clock_error>1e-9:raise ValueError('engine clock mismatch')
        with np.load(directory/'plane_coordinates.npz',allow_pickle=False) as z:
            t=z['time_seconds']*1e12;q=z['coordinates_m'][:,0,case['axis']]*1e10
            thermo=z['thermo'];power=z['external_power_eV_ps'];com=z['center_velocity_angstrom_ps']
            force_samples=z['conjugate_force_eV_A']
        if (t[0]!=0 or abs(t[-1]-case['duration_ps'])>1e-9 or not np.all(np.isfinite(q))
                or not np.allclose(np.diff(t),.025,atol=1e-10,rtol=0)):
            raise ValueError('invalid trajectory/sampling')
        F=case['force_eV_A'];omega=2*np.pi*f
        force_error=float(np.max(abs(force_samples-F*np.cos(omega*t))))
        if force_error>1e-10*abs(F):raise ValueError('saved conjugate forcing mismatch')
        work=F*(np.cos(omega*t[-1])*q[-1]-q[0])+trapezoid(F*omega*np.sin(omega*t)*q,t)
        checks.append(dict(case=case['name'],trajectory_sha256=sha(directory/'plane_coordinates.npz'),
            clock_error_ps=clock_error,force_sample_error_eV_A=force_error,
            mean_temperature_K=float(thermo[:,0].mean()),
            first100_temperature_K=float(thermo[t<100,0].mean()),
            last100_temperature_K=float(thermo[t>=t[-1]-100,0].mean()),
            internal_energy_change_eV=float(thermo[-1,3]-thermo[0,3]),
            work_parts_eV=float(work),work_power_eV=float(trapezoid(power,t)),
            parts_energy_residual_eV=float(thermo[-1,3]-thermo[0,3]-work),
            power_energy_residual_eV=float(thermo[-1,3]-thermo[0,3]-trapezoid(power,t)),
            max_com_change_A_ps=float(np.max(np.linalg.norm(com-com[0],axis=1)))))
        duration=case['duration_ps']-protocol['exclude_ps']
        for parts in (1,2):
            for block in range(parts):
                lo=100+duration*block/parts;hi=100+duration*(block+1)/parts
                mask=(t>=lo)&(t<hi)
                row=dict(**case,parts=parts,block=block,window_ps=duration/parts,
                    **harmonic_response(t[mask],q[mask],f,F))
                for k in (2,3):
                    hh=harmonic_response(t[mask],q[mask],k*f,1.)
                    row.update({f'harmonic{k}_real_A':hh['real_A2_eV'],f'harmonic{k}_imag_A':hh['imag_A2_eV']})
                fits.append(row)
        # Phase work and q*sin integral share data/algebra; not independent tests.
        start=int(round(protocol['exclude_ps']/protocol['frame_ps']))
        endpoint=F*(np.cos(omega*t[-1])*q[-1]-np.cos(omega*t[start])*q[start])
        steady_parts=endpoint+trapezoid(F*omega*np.sin(omega*t[start:])*q[start:],t[start:])
        steady_fit=next(r for r in fits if r['name']==case['name'] and r['parts']==1)
        periodic_work=-np.pi*F**2*steady_fit['imag_A2_eV']*f*duration
        identity_error=float(steady_parts-endpoint-periodic_work)
        if abs(identity_error)>1e-8:raise ValueError('integer-cycle phase/work identity mismatch')
        checks[-1].update(steady_endpoint_work_eV=float(endpoint),
            steady_parts_work_eV=float(steady_parts),fitted_periodic_work_eV=float(periodic_work),
            phase_work_identity_error_eV=identity_error,
            steady_internal_energy_change_eV=float(thermo[-1,3]-thermo[start,3]),
            steady_power_work_eV=float(trapezoid(power[start:],t[start:])))
    paired=[]
    for case in protocol['cases']:
        if case['sign']!=1:continue
        for parts in (1,2):
            for block in range(parts):
                rr=[r for r in fits if all(r[k]==case[k] for k in ('axis','dt_ps','fraction'))
                    and r['parts']==parts and r['block']==block]
                if len(rr)!=2:raise ValueError('incomplete signed pair')
                chi=sum(complex(r['real_A2_eV'],r['imag_A2_eV']) for r in rr)/2
                nn=[n for n in null if (n['axis'],n['duration_ps'],n['harmonic'])==
                    (case['axis'],rr[0]['window_ps'],1)]
                radius=max(np.hypot(n['real_A'],n['imag_A']) for n in nn)/case['force_eV_A']
                row=dict(axis=case['axis'],fraction=case['fraction'],dt_ps=case['dt_ps'],parts=parts,
                    block=block,window_ps=rr[0]['window_ps'],real_A2_eV=chi.real,imag_A2_eV=chi.imag,
                    phase_deg=float(np.degrees(np.angle(chi))),complex_null_radius_A2_eV=radius,
                    **loss_envelope(chi.imag,[n['imag_A'] for n in nn],case['force_eV_A']),
                    **{k:v for k,v in scalar_impedance(chi,f).items() if k!='production_clock_calibrated'},
                    **{k:v for k,v in inverse_disk_bounds(chi,radius,f).items()
                        if k not in ('production_clock_calibrated','confidence_interval')})
                for k in (2,3):
                    # Susceptibility fits include signed F; raw harmonic fits do not.
                    values={r['sign']:complex(r[f'harmonic{k}_real_A'],r[f'harmonic{k}_imag_A']) for r in rr}
                    z=(values[1]+values[-1])/2 if k==2 else (values[1]-values[-1])/2
                    hn=[n for n in null if (n['axis'],n['duration_ps'],n['harmonic'])==
                        (case['axis'],rr[0]['window_ps'],k)]
                    envelope=max(np.hypot(n['real_A'],n['imag_A']) for n in hn)
                    row.update({f'parity_harmonic{k}_A':abs(z),f'parity_harmonic{k}_null_A':envelope,
                                f'parity_harmonic{k}_to_null':abs(z)/envelope})
                paired.append(row)
    output.mkdir(parents=True)
    for name,rows in [('responses',fits),('paired_responses',paired),('work_checks',checks),
                      ('null_windows',null),('fdt_spectral',fdt)]: write_csv(output/(name+'.csv'),rows)
    save_json(output/'summary.json',dict(completed=True,protocol=protocol,paired=paired,checks=checks,
        independent_replicas=False,production_clock_calibrated=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(__doc__);sub=p.add_subparsers(dest='action',required=True)
    for action in ('prepare','run','analyze'):
        q=sub.add_parser(action);q.add_argument('--study',type=Path,required=True)
        if action in ('prepare','analyze'):q.add_argument('--source',type=Path,required=True)
        if action in ('prepare','run'):
            q.add_argument('--restart',type=Path,required=True);q.add_argument('--potential',type=Path,required=True)
        if action=='analyze':q.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if a.action=='prepare':prepare(a.source,a.restart,a.potential,a.study)
    elif a.action=='run':run(a.study,a.restart,a.potential)
    else:analyze(a.study,a.source,a.out)
