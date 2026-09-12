"""Matched signed/amplitude reference-MD test, never a production clock.

Eight declared400ps runs at .05cycles/ps; amplitudes are .5 and1 thermal RMS
linear-response scales, determined BEFORE observing driven response.
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
from .collective_forcing import harmonic_response,susceptibility_from_covariance
from .mode_kinetic_calibration import mode_correlations
from .low_frequency_mobility import multitaper_spectrum,band_integral_proxy
from .run_vector_registry_audit import save_json
from .run_low_stress_cyclic_diagnostic import write_csv


def execute(source,restart,potential,output):
    source,restart,potential,output=map(Path,(source,restart,potential,output))
    if output.exists():raise FileExistsError('fresh study directory required')
    meta=json.loads((source/'summary.json').read_bytes())
    if (not meta['completed'] or meta['repeats']!=6 or meta['ensemble']!='nve'
            or meta['frame_ps']!=.025 or meta['dt_ps']!=.0025 or meta['duration_ps']!=5000.):
        raise ValueError('matched completed six-plane 5ns NVE reference required')
    if hashlib.sha256(restart.read_bytes()).hexdigest()!=meta['restart_sha256']:
        raise ValueError('matched initial state required')
    with np.load(source/'plane_coordinates.npz') as z:
        q=z['coordinates_m'][1000:];T=float(z['thermo'][1000:,0].mean())
    rms=np.sqrt(np.mean((q-q.mean(axis=0))**2,axis=(0,1)))
    # kBT/rms in newtons, converted to TOTAL generalized eV/Angstrom.
    scales=1.380649e-23*T/rms/(1.602176634e-19/1e-10)
    cases=[dict(name=f'axis{axis}_fraction{fraction:g}_sign{sign:+d}',axis=axis,
                fraction=fraction,sign=sign,force_eV_A=float(sign*fraction*scales[axis]))
           for axis in (0,1) for fraction in (.5,1.) for sign in (-1,1)]
    output.mkdir(parents=True)
    save_json(output/'protocol.json',dict(cases=cases,frequency_per_ps=.05,
        duration_ps=400.,dt_ps=.0025,temperature_K=T,rms_m=rms,
        source_projection_sha256=hashlib.sha256((source/'plane_coordinates.npz').read_bytes()).hexdigest(),
        restart_sha256=meta['restart_sha256'],amplitudes_chosen_before_response=True,
        production_clock_calibrated=False))
    def one(case):
        cmd=[sys.executable,'-X','utf8','-m','solver_v1.run_reference_thermostat_md',
            '--potential',str(potential),'--out',str(output/case['name']),
            '--ensemble','nve','--repeats','6','--restart',str(restart),'--threads','1',
            '--duration-ps','400','--dt-ps','.0025','--drive-axis',str(case['axis']),
            '--drive-force-eV-A',str(case['force_eV_A']),'--drive-frequency-per-ps','.05']
        with (output/(case['name']+'.log')).open('w',encoding='utf8') as log:
            subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
        print('Completed '+case['name'],flush=True)
    with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(one,cases))


def analyze(study,source,output):
    study,source,output=map(Path,(study,source,output))
    if output.exists():raise FileExistsError('fresh result directory required')
    protocol=json.loads((study/'protocol.json').read_bytes())
    if hashlib.sha256((source/'plane_coordinates.npz').read_bytes()).hexdigest()!=protocol['source_projection_sha256']:
        raise ValueError('equilibrium source mismatch')
    rows=[];checks=[];cross=[]
    frequency=protocol['frequency_per_ps']
    for case in protocol['cases']:
        directory=study/case['name'];meta=json.loads((directory/'summary.json').read_bytes())
        if (not meta['completed'] or meta['conjugate_drive']['force_eV_A']!=case['force_eV_A']
                or meta['conjugate_drive']['axis']!=case['axis']
                or meta['conjugate_drive']['frequency_per_ps']!=frequency
                or meta['duration_ps']!=protocol['duration_ps'] or meta['dt_ps']!=protocol['dt_ps']
                or meta['frame_ps']!=.025 or meta['restart_sha256']!=protocol['restart_sha256']):
            raise ValueError('completed matched drive required')
        with np.load(directory/'plane_coordinates.npz') as z:
            t=z['time_seconds']*1e12;all_q=z['coordinates_m'][:,0,:]*1e10
            q=all_q[:,case['axis']]
            thermo=z['thermo'];power=z['external_power_eV_ps'];com=z['center_velocity_angstrom_ps']
        force=case['force_eV_A']*np.cos(2*np.pi*frequency*t)
        force_dot=-case['force_eV_A']*2*np.pi*frequency*np.sin(2*np.pi*frequency*t)
        displacement_work=float(force[-1]*q[-1]-force[0]*q[0]-trapezoid(force_dot*q,t))
        checks.append(dict(case=case['name'],mean_temperature_K=float(thermo[:,0].mean()),
            work_eV=float(trapezoid(power,t)),internal_energy_change_eV=float(thermo[-1,3]-thermo[0,3]),
            work_energy_residual_eV=float(thermo[-1,3]-thermo[0,3]-trapezoid(power,t)),
            displacement_work_eV=displacement_work,
            displacement_work_energy_residual_eV=float(thermo[-1,3]-thermo[0,3]-displacement_work),
            max_center_velocity_change_A_ps=float(np.max(np.linalg.norm(com-com[0],axis=1))),
            trajectory_sha256=hashlib.sha256((directory/'plane_coordinates.npz').read_bytes()).hexdigest()))
        # Exclude first40ps =2 loading cycles; full, halves and quarters.
        for parts in (1,2,4):
            for block in range(parts):
                lo=40+360*block/parts;hi=40+360*(block+1)/parts
                mask=(t>=lo)&(t<hi)
                rows.append(dict(case=case['name'],axis=case['axis'],fraction=case['fraction'],
                    sign=case['sign'],parts=parts,block=block,
                    **harmonic_response(t[mask],q[mask],frequency,case['force_eV_A'])))
                if parts==1:
                    for response_axis in range(3):
                        cross.append(dict(case=case['name'],drive_axis=case['axis'],response_axis=response_axis,
                            **harmonic_response(t[mask],all_q[mask,response_axis],frequency,case['force_eV_A'])))
    paired=[]
    for axis in (0,1):
        for fraction in (.5,1.):
            for parts in (1,2,4):
                for block in range(parts):
                    r=[x for x in rows if (x['axis'],x['fraction'],x['parts'],x['block'])==(axis,fraction,parts,block)]
                    chi=np.array([x['real_A2_eV']+1j*x['imag_A2_eV'] for x in r])
                    paired.append(dict(axis=axis,fraction=fraction,parts=parts,block=block,
                        real_A2_eV=float(chi.mean().real),imag_A2_eV=float(chi.mean().imag),
                        signed_response_disagreement_A2_eV=float(abs(chi[0]-chi[1])/2)))
    with np.load(source/'plane_coordinates.npz') as z:q=z['coordinates_m'][1000:]
    loss=[]
    for parts in (1,2):
        width=len(q)//parts
        for block in range(parts):
            for nw in (2.,3.):
                spec=multitaper_spectrum(q[block*width:(block+1)*width],.025e-12,nw=nw,tapers=2)
                for lo,hi in ((.04,.08),(.04,.06),(.045,.055)):
                    band=band_integral_proxy(spec,lo*1e12,hi*1e12)
                    for axis in (0,1):
                        imaginary=-2*np.pi*frequency*1e12*band['integral_proxy_m2_seconds'][axis,axis]/(1.380649e-23*protocol['temperature_K'])
                        force_SI=next(c['force_eV_A'] for c in protocol['cases']
                            if c['axis']==axis and c['fraction']==1 and c['sign']==1)*1.602176634e-9
                        spectral=2*band['integral_proxy_m2_seconds'][axis,axis]
                        required=2*spectral/(force_SI*abs(imaginary)/3)**2
                        loss.append(dict(axis=axis,parts=parts,block=block,nw=nw,
                            lower_cycle_THz=lo,upper_cycle_THz=hi,
                            imag_A2_eV=float(imaginary*16.02176634),
                            planning_single_record_SNR3_seconds=float(required),
                            planning_is_not_confidence_interval=True))
    dt=.025;C=mode_correlations(q,1600).sum(axis=1)/q.shape[1]
    fdt=[]
    for axis in (0,1):
        for cutoff in (2.,5.,10.,20.,40.):
            n=int(round(cutoff/dt))+1
            chi,cut=susceptibility_from_covariance(np.arange(n)*dt*1e-12,C[:n,axis],
                frequency*1e12,1.380649e-23*protocol['temperature_K'])
            factor=1.602176634e-19/1e-20
            fdt.append(dict(axis=axis,cutoff_ps=cutoff,real_A2_eV=float(chi.real*factor),
                imag_A2_eV=float(chi.imag*factor),endpoint_difference_A2_eV=float(abs(chi-cut)*factor),
                static_A2_eV=float(C[0,axis]/(1.380649e-23*protocol['temperature_K'])*factor)))
    output.mkdir(parents=True)
    write_csv(output/'responses.csv',rows);write_csv(output/'paired_responses.csv',paired)
    write_csv(output/'fdt_cutoffs.csv',fdt);write_csv(output/'work_checks.csv',checks)
    write_csv(output/'cross_responses.csv',cross)
    write_csv(output/'fdt_spectral_loss.csv',loss)
    save_json(output/'summary.json',dict(completed=True,protocol=protocol,paired=paired,
        fdt=fdt,spectral_loss=loss,work_checks=checks,finite_microcanonical_canonical_difference_unresolved=True,
        time_blocks_not_independent=True,production_clock_calibrated=False))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(10,7),constrained_layout=True)
    for axis,label in enumerate(('normal','direct110')):
        for column,component in enumerate(('real_A2_eV','imag_A2_eV')):
            ax=axes[axis,column]
            predicted=([x[component] for x in fdt if x['axis']==axis] if column==0 else
                       [x[component] for x in loss if x['axis']==axis and x['parts']==1])
            ax.axhspan(min(predicted),max(predicted),alpha=.25,color='C1',label='Equilibrium FDT sensitivity')
            for fraction in (.5,1.):
                full=next(x for x in paired if x['axis']==axis and x['fraction']==fraction and x['parts']==1)
                values=[x[component] for x in paired if x['axis']==axis and x['fraction']==fraction]
                ax.plot([fraction,fraction],[min(values),max(values)],color='C0',lw=2)
                ax.plot(fraction,full[component],'o',color='C0')
            ax.set_title(label+(' in-phase' if column==0 else 'quadrature'))
            ax.set_xlabel('Force / thermal-RMS force scale');ax.set_xticks([.5,1.])
            ax.set_ylabel('Susceptibility [Angstrom²/eV]');ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
            ax.grid(alpha=.2)
    axes[0,0].plot([],[],'o-',color='C0',label='Signed pairs; full/half/quarter range')
    axes[0,0].legend(fontsize=8)
    fig.suptitle('Reference-MD forced response: ranges are sensitivity, not confidence intervals')
    fig.savefig(output/'response_comparison.png',dpi=150);plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=('run','analyze'));p.add_argument('--source',required=True)
    p.add_argument('--out',required=True);p.add_argument('--restart');p.add_argument('--potential');p.add_argument('--study')
    a=p.parse_args()
    if a.action=='run':execute(a.source,a.restart,a.potential,a.out)
    else:analyze(a.study,a.source,a.out)
