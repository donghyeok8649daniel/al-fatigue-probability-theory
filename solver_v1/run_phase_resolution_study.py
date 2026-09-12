"""v27: predeclared finite-frequency phase test, preserving unresolved v26.

prepare -> run -> analyze. Only reference MD; never changes production physics.
Frequencies are selected from equilibrium SNR planning BEFORE driven results.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
from scipy.integrate import trapezoid, simpson
from .collective_forcing import harmonic_response, susceptibility_from_covariance
from .phase_resolution import null_lockin_windows, loss_envelope, harmonic_content
from .low_frequency_mobility import multitaper_spectrum, band_integral_proxy
from .mode_kinetic_calibration import mode_correlations
from .run_vector_registry_audit import save_json
from .run_low_stress_cyclic_diagnostic import write_csv


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def prepare(source, restart, output):
    source, restart, output = map(Path, (source, restart, output))
    if output.exists(): raise FileExistsError('fresh study directory required')
    meta = json.loads((source/'summary.json').read_bytes())
    if (not meta['completed'] or meta['repeats'] != 6 or meta['ensemble'] != 'nve'
            or meta['dt_ps'] != .0025 or meta['frame_ps'] != .025
            or meta['duration_ps'] != 5000 or sha(restart) != meta['restart_sha256']):
        raise ValueError('matched v25 equilibrium/restart required')
    with np.load(source/'plane_coordinates.npz') as z:
        q = z['coordinates_m'][1000:]; temperature = float(z['thermo'][1000:, 0].mean())
    kT = 1.380649e-23*temperature
    rms = np.sqrt(((q-q.mean(axis=0))**2).mean(axis=(0, 1)))
    forces = kT/rms/1.602176634e-9
    spectrum = multitaper_spectrum(q, .025e-12)
    planning = []
    for frequency in (.05, .2, .5, 1., 2.):
        band = band_integral_proxy(spectrum, (frequency-.01)*1e12, (frequency+.01)*1e12)
        for axis in (0, 1):
            proxy = band['integral_proxy_m2_seconds'][axis, axis]
            omega = 2*np.pi*frequency*1e12
            loss = omega*proxy/kT
            required = 4*proxy*9/((forces[axis]*1.602176634e-9*loss)**2)
            planning.append(dict(axis=axis, frequency_per_ps=frequency,
                predicted_loss_A2_eV=float(loss*16.02176634),
                SNR3_duration_ps=float(required*1e12)))
    chosen = {axis: min(x['frequency_per_ps'] for x in planning
              if x['axis'] == axis and x['SNR3_duration_ps'] <= 360.) for axis in (0, 1)}
    cases = []
    for axis in (0, 1):
        for dt in (.0025, .00125):
            for fraction in ((.5, 1.) if dt == .0025 else (1.,)):
                for sign in (-1, 1):
                    cases.append(dict(name=f'axis{axis}_dt{dt:g}_fraction{fraction:g}_sign{sign:+d}',
                        axis=axis, dt_ps=dt, frequency_per_ps=chosen[axis],
                        fraction=fraction, sign=sign, force_eV_A=float(sign*fraction*forces[axis])))
    output.mkdir(parents=True)
    save_json(output/'protocol.json', dict(cases=cases, planning=planning,
        frequency_selection='lowest predeclared frequency with single-record planned SNR3 <=360ps',
        selection_precedes_driven_results=True, duration_ps=400., exclude_ps=40.,
        temperature_K=temperature, rms_m=rms, source_sha256=sha(source/'plane_coordinates.npz'),
        restart_sha256=sha(restart), force_scale_eV_A=forces,
        production_clock_calibrated=False))


def compare_responses(paired, fdt, cutoffs, null, protocol):
    comparisons = []
    for axis in (0,1):
        selected = [r for r in paired if r['axis']==axis and r['parts']==1]
        base = next(r for r in selected if r['dt_ps']==.0025 and r['fraction']==1.)
        fine = next(r for r in selected if r['dt_ps']==.00125)
        half = next(r for r in selected if r['fraction']==.5)
        f = base['frequency_per_ps']
        spectral = [r['imag_A2_eV'] for r in fdt if r['axis']==axis and r['frequency_per_ps']==f]
        storage = [r['real_A2_eV'] for r in cutoffs if r['axis']==axis and r['frequency_per_ps']==f and r['cutoff_ps']>=10]
        empirical = [r for r in null if r['axis']==axis and r['frequency_per_ps']==f and r['duration_ps']==360.]
        F = protocol['force_scale_eV_A'][axis]
        floor_complex = max(np.hypot(r['real_A'],r['imag_A']) for r in empirical)/F
        complex_of = lambda r: complex(r['real_A2_eV'],r['imag_A2_eV'])
        comparisons.append(dict(axis=axis,frequency_per_ps=f,
            static_covariance_chi_A2_eV=float(protocol['rms_m'][axis]**2/(1.380649e-23*protocol['temperature_K'])*16.02176634),
            storage_to_static_ratio=float(base['real_A2_eV']/(protocol['rms_m'][axis]**2/(1.380649e-23*protocol['temperature_K'])*16.02176634)),
            base_phase_deg=base['phase_deg'],fine_phase_deg=fine['phase_deg'],half_force_phase_deg=half['phase_deg'],
            base_loss_to_null_ratio=base['loss_to_null_ratio'],fine_loss_to_null_ratio=fine['loss_to_null_ratio'],
            base_loss_exceeds_observed_null=base['loss_exceeds_observed_null'],
            fine_loss_exceeds_observed_null=fine['loss_exceeds_observed_null'],
            fdt_imag_min_A2_eV=min(spectral),fdt_imag_max_A2_eV=max(spectral),
            fdt_real_min_A2_eV=min(storage),fdt_real_max_A2_eV=max(storage),
            dt_complex_change_A2_eV=abs(complex_of(base)-complex_of(fine)),
            amplitude_complex_change_A2_eV=abs(complex_of(base)-complex_of(half)),
            observed_complex_null_envelope_A2_eV=floor_complex,
            # Sum, not root-sum-square: no assumed independent noise cancellation.
            dt_observed_noise_comparison_A2_eV=2*floor_complex,
            amplitude_observed_noise_comparison_A2_eV=3*floor_complex,
            numerical_discretization_error_separated_from_thermal_noise=False,
            zero_frequency_mobility_validated=False))
    return comparisons


def run(study, restart, potential):
    study = Path(study); protocol = json.loads((study/'protocol.json').read_bytes())
    if sha(restart) != protocol['restart_sha256']: raise ValueError('restart mismatch')
    def one(case):
        out = study/case['name']
        if out.exists(): raise FileExistsError('existing case must not be overwritten')
        cmd = [sys.executable, '-X', 'utf8', '-m', 'solver_v1.run_reference_thermostat_md',
            '--potential', str(potential), '--out', str(out), '--ensemble', 'nve',
            '--repeats', '6', '--restart', str(restart), '--threads', '1',
            '--duration-ps', str(protocol['duration_ps']), '--dt-ps', str(case['dt_ps']),
            '--drive-axis', str(case['axis']), '--drive-force-eV-A', str(case['force_eV_A']),
            '--drive-frequency-per-ps', str(case['frequency_per_ps'])]
        with (study/(case['name']+'.log')).open('x', encoding='utf8') as log:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
        print('Completed '+case['name'], flush=True)
    with ThreadPoolExecutor(max_workers=2) as pool: list(pool.map(one, protocol['cases']))


def audit_harmonics(study, source, output):
    """Secondary audit; retains all individual and force-parity quadratures."""
    study,source,output=map(Path,(study,source,output))
    if output.exists(): raise FileExistsError('fresh harmonic audit directory required')
    protocol=json.loads((study/'protocol.json').read_bytes())
    if sha(source/'plane_coordinates.npz')!=protocol['source_sha256']:
        raise ValueError('equilibrium source mismatch')
    with np.load(source/'plane_coordinates.npz') as z:
        t=z['time_seconds'][1000:]*1e12;q=z['coordinates_m'][1000:]*1e10
    null=[];responses=[];paired=[]
    for axis in (0,1):
        f=next(c['frequency_per_ps'] for c in protocol['cases'] if c['axis']==axis)
        for k in (2,3):
            null.extend(dict(axis=axis,harmonic=k,**r) for r in
                null_lockin_windows(t,q[:,:,axis],k*f,360.))
    for c in protocol['cases']:
        meta=json.loads((study/c['name']/'summary.json').read_bytes())
        if (not meta['completed'] or meta['restart_sha256']!=protocol['restart_sha256']
                or meta['dt_ps']!=c['dt_ps'] or meta['duration_ps']!=400.
                or meta['conjugate_drive']!=dict(axis=c['axis'],force_eV_A=c['force_eV_A'],frequency_per_ps=c['frequency_per_ps'])):
            raise ValueError('mismatched completed harmonic case required')
        with np.load(study/c['name']/'plane_coordinates.npz') as z:
            t=z['time_seconds']*1e12;mask=(t>=40)&(t<400)
            q=z['coordinates_m'][:,0,c['axis']]*1e10
        for k in (2,3):
            r=harmonic_response(t[mask],q[mask],k*c['frequency_per_ps'],1.)
            responses.append(dict(**c,harmonic=k,real_A=r['real_A2_eV'],imag_A=r['imag_A2_eV']))
    for c in protocol['cases']:
        if c['sign']!=1:continue
        for k in (2,3):
            r=[x for x in responses if all(x[key]==c[key] for key in ('axis','dt_ps','fraction')) and x['harmonic']==k]
            plus=next(x for x in r if x['sign']==1);minus=next(x for x in r if x['sign']==-1)
            zp=complex(plus['real_A'],plus['imag_A']);zm=complex(minus['real_A'],minus['imag_A'])
            envelope=max(np.hypot(x['real_A'],x['imag_A']) for x in null if x['axis']==c['axis'] and x['harmonic']==k)
            expected=(zp+(-1)**k*zm)/2
            paired.append(dict(axis=c['axis'],dt_ps=c['dt_ps'],fraction=c['fraction'],harmonic=k,
                even_amplitude_A=abs((zp+zm)/2),odd_amplitude_A=abs((zp-zm)/2),
                individual_max_amplitude_A=max(abs(zp),abs(zm)),
                observed_null_envelope_A=float(envelope),expected_parity_amplitude_A=abs(expected),
                expected_parity_to_null_ratio=float(abs(expected)/envelope),confidence_interval=False))
    output.mkdir(parents=True)
    for name,rows in [('null_windows',null),('responses',responses),('paired',paired)]:
        write_csv(output/(name+'.csv'),rows)
    save_json(output/'summary.json',dict(completed=True,paired=paired,
        source_sha256=protocol['source_sha256'],production_clock_calibrated=False))


def analyze(study, source, output):
    study, source, output = map(Path, (study, source, output))
    if output.exists(): raise FileExistsError('fresh result directory required')
    protocol = json.loads((study/'protocol.json').read_bytes())
    if sha(source/'plane_coordinates.npz') != protocol['source_sha256']:
        raise ValueError('equilibrium source mismatch')
    with np.load(source/'plane_coordinates.npz') as z:
        tq = z['time_seconds'][1000:]*1e12; qeq = z['coordinates_m'][1000:]
    null = []; fdt = []; cutoffs = []; fits = []; checks = []
    frequencies = sorted({.05} | {x['frequency_per_ps'] for x in protocol['cases']})
    for f in frequencies:
        for axis in (0, 1):
            for duration in (180., 360.):
                null.extend(dict(axis=axis, frequency_per_ps=f, duration_ps=duration, **x)
                    for x in null_lockin_windows(tq, qeq[:, :, axis]*1e10, f, duration))
    for parts in (1, 2):
        width = len(qeq)//parts
        for block in range(parts):
            for nw in (2., 3.):
                spec = multitaper_spectrum(qeq[block*width:(block+1)*width], .025e-12, nw=nw)
                for f in frequencies:
                    for halfband in (.005, .01, .02):
                        proxy = band_integral_proxy(spec, (f-halfband)*1e12, (f+halfband)*1e12)
                        for axis in (0, 1):
                            loss = 2*np.pi*f*1e12*proxy['integral_proxy_m2_seconds'][axis, axis]/(1.380649e-23*protocol['temperature_K'])*16.02176634
                            fdt.append(dict(axis=axis, frequency_per_ps=f, parts=parts, block=block,
                                nw=nw, halfband_per_ps=halfband, imag_A2_eV=float(-loss)))
    covariance = mode_correlations(qeq, 1600).sum(axis=1)/qeq.shape[1]
    quadrature=[]
    for axis in (0,1):
        frequency=next(c['frequency_per_ps'] for c in protocol['cases'] if c['axis']==axis)
        omega=2*np.pi*frequency
        for stride in (1,2):
            tt=np.arange(0,1601,stride)*.025; cc=covariance[::stride,axis]
            for rule in (trapezoid,simpson):
                chi=(cc[0]-1j*omega*rule(cc*np.exp(-1j*omega*tt),x=tt))/(1.380649e-23*protocol['temperature_K'])*16.02176634
                quadrature.append(dict(axis=axis,frequency_per_ps=frequency,cutoff_ps=40.,
                    covariance_sampling_ps=.025*stride,quadrature=rule.__name__,
                    real_A2_eV=float(chi.real),imag_A2_eV=float(chi.imag),
                    coarsening_is_not_fine_sampling_certification=True))
    for f in frequencies:
        for axis in (0, 1):
            for cutoff in (2., 5., 10., 20., 40.):
                n = int(round(cutoff/.025))+1
                chi, truncated = susceptibility_from_covariance(np.arange(n)*.025e-12,
                    covariance[:n, axis], f*1e12, 1.380649e-23*protocol['temperature_K'])
                cutoffs.append(dict(axis=axis, frequency_per_ps=f, cutoff_ps=cutoff,
                    real_A2_eV=float(chi.real*16.02176634), imag_A2_eV=float(chi.imag*16.02176634),
                    endpoint_error_A2_eV=float(abs(chi-truncated)*16.02176634)))
    for case in protocol['cases']:
        directory = study/case['name']; meta = json.loads((directory/'summary.json').read_bytes())
        if (not meta['completed'] or meta['restart_sha256'] != protocol['restart_sha256']
                or meta['duration_ps'] != protocol['duration_ps'] or meta['dt_ps'] != case['dt_ps']
                or meta['frame_ps'] != .025 or meta['ensemble'] != 'nve'
                or meta['conjugate_drive'] != dict(axis=case['axis'], force_eV_A=case['force_eV_A'],
                    frequency_per_ps=case['frequency_per_ps'])):
            raise ValueError('incomplete or mismatched case '+case['name'])
        clock_rows=[]
        with (directory/'lammps.log').open(encoding='utf8') as log:
            for line in log:
                columns=line.split()
                if len(columns)!=7 or not columns[0].isdigit(): continue
                try: values=list(map(float,columns))
                except ValueError: continue
                clock_rows.append(values[:2])
        clock_rows=np.asarray(clock_rows)
        if (len(clock_rows)==0 or clock_rows[0,1]!=0
                or not np.isclose(clock_rows[-1,1],protocol['duration_ps'],atol=1e-9,rtol=0)):
            raise ValueError('logged engine clock endpoints do not match')
        clock_error=float(np.max(np.abs(clock_rows[:,1]-clock_rows[:,0]*case['dt_ps'])))
        if clock_error>1e-9: raise ValueError('logged engine and analysis clocks differ')
        with np.load(directory/'plane_coordinates.npz') as z:
            t = z['time_seconds']*1e12; q = z['coordinates_m'][:, 0, case['axis']]*1e10
            thermo = z['thermo']; power = z['external_power_eV_ps']; com = z['center_velocity_angstrom_ps']
        w = 2*np.pi*case['frequency_per_ps']; F = case['force_eV_A']
        work = float(F*np.cos(w*t[-1])*q[-1]-F*np.cos(w*t[0])*q[0]+trapezoid(F*w*np.sin(w*t)*q,t))
        checks.append(dict(case=case['name'], trajectory_sha256=sha(directory/'plane_coordinates.npz'),
            logged_clock_max_error_ps=clock_error,
            mean_temperature_K=float(thermo[:,0].mean()),
            first40_temperature_K=float(thermo[t<40,0].mean()), last40_temperature_K=float(thermo[t>=360,0].mean()),
            work_parts_eV=work, work_power_eV=float(trapezoid(power,t)),
            internal_energy_change_eV=float(thermo[-1,3]-thermo[0,3]),
            parts_energy_residual_eV=float(thermo[-1,3]-thermo[0,3]-work),
            power_energy_residual_eV=float(thermo[-1,3]-thermo[0,3]-trapezoid(power,t)),
            max_com_change_A_ps=float(np.max(np.linalg.norm(com-com[0],axis=1)))))
        for parts in (1,2):
            for block in range(parts):
                lo=40+360*block/parts;hi=40+360*(block+1)/parts;mask=(t>=lo)&(t<hi)
                harmonics=harmonic_content(t[mask],q[mask],case['frequency_per_ps'])
                fits.append(dict(**case, parts=parts, block=block,
                    **harmonic_response(t[mask],q[mask],case['frequency_per_ps'],F),
                    **{f'harmonic_{k}_amplitude_A':harmonics['amplitude_A'][k-1] for k in (1,2,3)},
                    joint_harmonic_residual_rms_A=harmonics['residual_rms_A']))
        start=int(round(protocol['exclude_ps']/.025))
        steady=next(r for r in fits if r['name']==case['name'] and r['parts']==1)
        steady_boundary=F*(np.cos(w*t[-1])*q[-1]-np.cos(w*t[start])*q[start])
        steady_parts=steady_boundary+trapezoid(F*w*np.sin(w*t[start:])*q[start:],t[start:])
        steady_internal=thermo[-1,3]-thermo[start,3]
        steady_power=trapezoid(power[start:],t[start:])
        checks[-1].update(
            fitted_periodic_work_eV=float(-np.pi*F**2*steady['imag_A2_eV']*case['frequency_per_ps']*360),
            steady_endpoint_work_eV=float(steady_boundary),
            steady_parts_work_eV=float(steady_parts), steady_power_work_eV=float(steady_power),
            steady_internal_energy_change_eV=float(steady_internal),
            steady_parts_energy_residual_eV=float(steady_internal-steady_parts),
            steady_power_energy_residual_eV=float(steady_internal-steady_power))
    paired = []
    for case in protocol['cases']:
        if case['sign'] != 1: continue
        for parts in (1,2):
            for block in range(parts):
                rr=[r for r in fits if all(r[k]==case[k] for k in ('axis','dt_ps','fraction')) and r['parts']==parts and r['block']==block]
                assert len(rr)==2
                chi=np.mean([r['real_A2_eV']+1j*r['imag_A2_eV'] for r in rr])
                nq=[r['imag_A'] for r in null if (r['axis'],r['frequency_per_ps'],r['duration_ps'])==(case['axis'],case['frequency_per_ps'],360/parts)]
                envelope=loss_envelope(chi.imag,nq,case['force_eV_A'])
                paired.append(dict(axis=case['axis'],frequency_per_ps=case['frequency_per_ps'],
                    dt_ps=case['dt_ps'],fraction=case['fraction'],parts=parts,block=block,
                    real_A2_eV=float(chi.real),imag_A2_eV=float(chi.imag),phase_deg=float(np.angle(chi)*180/np.pi),
                    sign_disagreement_A2_eV=float(abs(complex(rr[0]['real_A2_eV'],rr[0]['imag_A2_eV'])-complex(rr[1]['real_A2_eV'],rr[1]['imag_A2_eV']))/2),
                    **envelope))
    comparisons = compare_responses(paired, fdt, cutoffs, null, protocol)
    output.mkdir(parents=True)
    for name, rows in [('planning',protocol['planning']),('null_windows',null),('fdt_spectral',fdt),
                       ('fdt_cutoffs',cutoffs),('responses',fits),('paired_responses',paired),('work_checks',checks),
                       ('phase_comparison',comparisons),('covariance_quadrature_audit',quadrature)]:
        write_csv(output/(name+'.csv'),rows)
    save_json(output/'summary.json',dict(completed=True,protocol=protocol,paired=paired,checks=checks,comparisons=comparisons,
        no_independence_or_confidence_claim=True,production_clock_calibrated=False))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes=plt.subplots(2,2,figsize=(11,7),constrained_layout=True)
    for axis,label in enumerate(('normal','direct110')):
        f=next(c['frequency_per_ps'] for c in protocol['cases'] if c['axis']==axis)
        reference=[r for r in fdt if r['axis']==axis and r['frequency_per_ps']==f]
        for row,key in enumerate(('real_A2_eV','imag_A2_eV')):
            ax=axes[row,axis]
            values=([r[key] for r in cutoffs if r['axis']==axis and r['frequency_per_ps']==f and r['cutoff_ps']>=10]
                    if row==0 else [r[key] for r in reference])
            ax.axhspan(min(values),max(values),color='C1',alpha=.3,label='Equilibrium FDT sensitivity')
            for r in paired:
                if r['axis']!=axis or r['parts']!=1: continue
                x=r['fraction']+(0.04 if r['dt_ps']==.00125 else 0)
                component='real_A' if row==0 else 'imag_A'
                bound=max(abs(n[component]) for n in null if n['axis']==axis and n['frequency_per_ps']==f and n['duration_ps']==360.)
                bound/=protocol['force_scale_eV_A'][axis]*r['fraction']
                ax.errorbar(x,r[key],yerr=bound,fmt='s' if r['dt_ps']==.00125 else 'o',
                    color='C2' if r['dt_ps']==.00125 else 'C0',capsize=4)
            ax.set_title(f'{label}: {f:g} cycles/ps; '+('storage' if row==0 else 'loss quadrature'))
            ax.set_xlabel('Force / thermal-RMS scale (fine dt offset +0.04)')
            ax.set_ylabel('Susceptibility [Angstrom²/eV]')
            ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0));ax.grid(alpha=.2)
            if row==1: ax.axhline(0,color='black',lw=.8)
    axes[0,0].plot([],[],'o',color='C0',label='dt 2.5 fs, signed pair')
    axes[0,0].plot([],[],'s',color='C2',label='dt 1.25 fs, signed pair')
    axes[0,0].legend(fontsize=8)
    fig.suptitle('Reference MD phase validation — bars: observed null envelopes, NOT confidence intervals')
    fig.savefig(output/'phase_comparison.png',dpi=150);plt.close(fig)


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=('prepare','run','analyze','harmonics'));p.add_argument('--source')
    p.add_argument('--restart');p.add_argument('--potential');p.add_argument('--study');p.add_argument('--out')
    args=p.parse_args()
    if args.action=='prepare': prepare(args.source,args.restart,args.out)
    elif args.action=='run': run(args.study,args.restart,args.potential)
    elif args.action=='harmonics': audit_harmonics(args.study,args.source,args.out)
    else: analyze(args.study,args.source,args.out)
