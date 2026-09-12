"""Internal-step reference-MD work audit; no changes to force or production clock."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
from scipy.integrate import trapezoid,cumulative_trapezoid
from .run_low_frequency_forcing_v29 import sha
from .run_vector_registry_audit import save_json
from .run_low_stress_cyclic_diagnostic import write_csv


def run(study,restart,potential):
    study=Path(study)
    if study.exists():raise FileExistsError('fresh audit required')
    source=json.loads(Path('results/low_frequency_forcing_v29/summary.json').read_bytes())['protocol']
    if sha(restart)!=source['restart_sha256'] or sha(potential)!=source['potential_sha256']:
        raise ValueError('matched restart/potential required')
    cases=[dict(name=f'axis{axis}_dt{dt:g}_{mode}',axis=axis,dt_ps=dt,
                frame_ps=dt if mode=='dense' else .025,
                duration_ps=2. if mode=='dense' else 25.,
                force_eV_A=4*source['force_scale_eV_A'][axis],mode=mode)
           for dt in (.0025,.00125,.000625) for axis in (0,1) for mode in ('coarse','dense')]
    study.mkdir(parents=True)
    save_json(study/'protocol.json',dict(cases=cases,restart_sha256=sha(restart),
        potential_sha256=sha(potential),frequency_per_ps=.2,production_clock_calibrated=False))
    def one(c):
        args=[sys.executable,'-X','utf8','-m','solver_v1.run_reference_thermostat_md',
              '--potential',str(potential),'--restart',str(restart),'--out',str(study/c['name']),
              '--ensemble','nve','--repeats','6','--threads','1','--drive-axis',str(c['axis']),
              '--drive-force-eV-A',str(c['force_eV_A']),'--drive-frequency-per-ps','.2',
              '--duration-ps',str(c['duration_ps']),'--dt-ps',str(c['dt_ps']),
              '--frame-ps',str(c['frame_ps']),'--timestep-work']
        with (study/(c['name']+'.log')).open('x') as stream:
            subprocess.run(args,stdout=stream,stderr=subprocess.STDOUT,check=True)
        print('Completed '+c['name'],flush=True)
    with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(one,cases))


def analyze(study,out):
    study,out=Path(study),Path(out)
    if out.exists():raise FileExistsError('existing report preserved')
    p=json.loads((study/'protocol.json').read_bytes());rows=[];bindings=[]
    for c in p['cases']:
        d=study/c['name'];m=json.loads((d/'summary.json').read_bytes())
        if (not m['completed'] or not m['internal_step_work'] or m['restart_sha256']!=p['restart_sha256']
            or any(m[k]!=c[k] for k in ('dt_ps','frame_ps','duration_ps'))):
            raise ValueError('incomplete/mismatched audit')
        bindings.append(dict(case=c['name'],trajectory_sha256=sha(d/'plane_coordinates.npz')))
    for c in p['cases']:
        if c['mode']!='coarse':continue
        d=study/c['name'];dense=study/c['name'].replace('coarse','dense')
        with np.load(d/'plane_coordinates.npz') as z,np.load(dense/'plane_coordinates.npz') as zd:
            t=z['time_seconds']*1e12;power=z['external_power_eV_ps'];work=z['internal_step_work_eV']
            q=z['coordinates_m'][:,0,c['axis']]*1e10;F=c['force_eV_A'];omega=2*np.pi*.2
            td=zd['time_seconds']*1e12;stride=round(.025/c['dt_ps']);n=len(td[::stride])
            dense_trap=trapezoid(zd['external_power_eV_ps'],td)
            energy=z['thermo'][:,3]-z['thermo'][0,3]
            residual=energy-work
            coarse_residual=energy-cumulative_trapezoid(power,t,initial=0.)
            rows.append(dict(axis=c['axis'],dt_ps=c['dt_ps'],duration_ps=c['duration_ps'],
                internal_energy_change_eV=float(z['thermo'][-1,3]-z['thermo'][0,3]),
                internal_step_work_eV=float(work[-1]),coarse_power_work_eV=float(trapezoid(power,t)),
                parts_work_eV=float(F*(np.cos(omega*t[-1])*q[-1]-q[0])+trapezoid(F*omega*np.sin(omega*t)*q,t)),
                work_energy_residual_eV=float(z['thermo'][-1,3]-z['thermo'][0,3]-work[-1]),
                max_work_energy_residual_eV=float(np.max(abs(residual))),
                rms_work_energy_residual_eV=float(np.sqrt(np.mean(residual**2))),
                max_coarse_power_energy_residual_eV=float(np.max(abs(coarse_residual))),
                dense_power_identity_error_eV=float(zd['internal_step_work_eV'][-1]-dense_trap),
                sparse_vs_dense_work_error_eV=float(work[n-1]-dense_trap),
                sparse_vs_dense_q_max_m=float(np.max(abs(z['coordinates_m'][:n]-zd['coordinates_m'][::stride])))))
    if max(abs(r['dense_power_identity_error_eV']) for r in rows)>1e-10:
        raise ValueError('native endpoint accumulation identity failed')
    if max(abs(r['sparse_vs_dense_work_error_eV']) for r in rows)>1e-9:
        raise ValueError('native block accumulation differs from dense record')
    out.mkdir(parents=True)
    write_csv(out/'work_convergence.csv',rows)
    save_json(out/'summary.json',dict(completed=True,protocol=p,bindings=bindings,rows=rows,
        production_clock_calibrated=False,work_is_not_loss_calibration=True))
    print(json.dumps(rows,indent=2))


def legacy_check(study,out,restart,potential):
    """Run old default (no internal work) and compare the same first2 ps."""
    study,out=Path(study),Path(out)
    if out.exists():raise FileExistsError('fresh legacy check required')
    p=json.loads((study/'protocol.json').read_bytes())
    if sha(restart)!=p['restart_sha256'] or sha(potential)!=p['potential_sha256']:
        raise ValueError('matched source required')
    out.mkdir(parents=True);rows=[]
    for c in p['cases']:
        if c['mode']!='coarse':continue
        directory=out/c['name']
        args=[sys.executable,'-X','utf8','-m','solver_v1.run_reference_thermostat_md',
              '--potential',str(potential),'--restart',str(restart),'--out',str(directory),
              '--ensemble','nve','--repeats','6','--threads','1','--drive-axis',str(c['axis']),
              '--drive-force-eV-A',str(c['force_eV_A']),'--drive-frequency-per-ps','.2',
              '--duration-ps','2','--dt-ps',str(c['dt_ps'])]
        subprocess.run(args,check=True)
        with np.load(directory/'plane_coordinates.npz') as z,np.load(study/c['name']/'plane_coordinates.npz') as ref:
            n=len(z['time_seconds'])
            qerr=float(np.max(abs(z['coordinates_m']-ref['coordinates_m'][:n])))
            eerr=float(np.max(abs(z['thermo']-ref['thermo'][:n])))
            if qerr>1e-20 or eerr>1e-8:raise ValueError('observational fix altered trajectory')
            rows.append(dict(axis=c['axis'],dt_ps=c['dt_ps'],coordinate_difference_m=qerr,
                thermo_difference_mixed_units=eerr,trajectory_sha256=sha(directory/'plane_coordinates.npz')))
    save_json(out/'comparison.json',dict(completed=True,rows=rows,production_clock_calibrated=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(__doc__);p.add_argument('action',choices=('run','analyze','legacy'))
    p.add_argument('--study',type=Path,required=True);p.add_argument('--restart',type=Path)
    p.add_argument('--potential',type=Path);p.add_argument('--out',type=Path)
    a=p.parse_args()
    if a.action=='run':run(a.study,a.restart,a.potential)
    elif a.action=='analyze':analyze(a.study,a.out)
    else:legacy_check(a.study,a.out,a.restart,a.potential)
