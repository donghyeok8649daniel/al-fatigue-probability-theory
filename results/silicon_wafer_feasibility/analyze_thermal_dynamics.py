"""Raw force correlations of actual constrained atom dynamics, without drag fitting."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid

from solver_v1.silicon_conditional_research import KB_EV_K
from solver_v1.silicon_thermal_dynamics import stationary_correlation
from .run_local_crack_audit import save_json


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--calculation',type=Path,required=True)
    args=parser.parse_args();out=args.calculation
    meta=json.loads((out/'summary.json').read_text(encoding='utf-8'))
    ensemble=Path(meta['ensemble_source'])
    ensemble_meta=json.loads((ensemble/'summary.json').read_text(encoding='utf-8'))
    rows=[];curves={}
    groups=sorted(set((r['reference'],r['temperature_K'],r['dt_ps']) for r in meta['records']))
    for name,temperature,dt in groups:
        records=sorted([r for r in meta['records'] if (r['reference'],r['temperature_K'],r['dt_ps'])==(name,temperature,dt)],
                       key=lambda r:r['seed'])
        matching=[r for r in ensemble_meta['records'] if r['reference']==name and r['temperature_K']==temperature]
        preparation_means=[json.loads((ensemble/(r['label']+'.json')).read_text(encoding='utf-8'))['statistics']['mean'][1]
                           for r in matching]
        center=float(np.mean(preparation_means))
        histories=[np.load(out/(r['label']+'.npz'),allow_pickle=False) for r in records]
        sample_step=float(histories[0]['time_ps'][1]-histories[0]['time_ps'][0])
        max_lag=min(round(10/sample_step),len(histories[0]['time_ps'])//4)
        times=np.arange(max_lag+1)*sample_step
        correlations=[];own_center=[];kinetic=[];early=[];late=[]
        for r,history in zip(records,histories):
            force=history['observations'][:,0]
            correlations.append(stationary_correlation(force,max_lag,center=center)/(KB_EV_K*temperature))
            own_center.append(stationary_correlation(force,max_lag)/(KB_EV_K*temperature))
            kinetic.append(float(2*history['kinetic_energy_eV'].mean()/(r['bath_dimension']*KB_EV_K)))
            half=len(force)//2
            early.append(stationary_correlation(force[:half],max_lag,center=center)/(KB_EV_K*temperature))
            late.append(stationary_correlation(force[-half:],max_lag,center=center)/(KB_EV_K*temperature))
        correlations=np.asarray(correlations);own_center=np.asarray(own_center)
        integrals=cumulative_trapezoid(correlations,times,axis=1,initial=0)
        own_integrals=cumulative_trapezoid(own_center,times,axis=1,initial=0)
        windows=[w for w in [.02,.05,.1,.2,.5,1.,2.,5.,10.] if w<=times[-1]]
        label=f'{name}_T{round(temperature)}_dt{dt:.5f}'
        curves.update({label+'_correlations':correlations,label+'_own_center_correlations':own_center,
            label+'_time_ps':times,label+'_integrals':integrals,label+'_early':np.asarray(early),label+'_late':np.asarray(late)})
        window_reports=[]
        for window in windows:
            index=round(window/sample_step)
            values=integrals[:,index]
            window_reports.append(dict(lag_ps=window,finite_integral_mean_eV_ps_A2=float(values.mean()),
                replica_sem=float(values.std(ddof=1)/np.sqrt(len(values))) if len(values)>1 else None,
                replica_values=values.tolist(),own_center_mean=float(own_integrals[:,index].mean())))
        rows.append(dict(label=label,reference=name,temperature_K=temperature,dt_ps=dt,
            independent_preparation_replicas=len(records),duration_per_replica_ps=records[0]['duration_ps'],
            mean_from_independent_preparation_eV_A=center,
            observed_force_means=[float(h['observations'][:,0].mean()) for h in histories],
            variance_over_kBT_replicas=correlations[:,0].tolist(),
            kinetic_temperatures_K=kinetic,finite_window_integrals=window_reports,
            early_late_kernel_max_difference=float(np.max(abs(np.mean(early,axis=0)-np.mean(late,axis=0)))),
            reflections=sum(r['reflection_count'] for r in records),
            contains_boundary_impacted_trajectory=any(r['reflection_count']>0 for r in records),
            trajectory_diagnostics=[dict(label=r['label'],reflections=r['reflection_count'],
                energy_residual_over_initial_kinetic=r['max_energy_error_over_initial_kinetic']) for r in records],
            max_energy_residual_over_initial_kinetic=max(r['max_energy_error_over_initial_kinetic'] for r in records),
            exact_nonlinear_Mori_kernel_certified=False,constant_drag_identified=False,physical_mobility=None))
        print(label,'Gamma0',correlations[:,0].mean(),'Tkin',np.mean(kinetic),
            'integrals',[(r['lag_ps'],r['finite_integral_mean_eV_ps_A2']) for r in window_reports],flush=True)
    np.savez_compressed(out/'correlations.npz',**curves)
    save_json(out/'analysis.json',dict(rows=rows,
        interpretation='nonlinear fixed-gap force correlations of actual mass-based atom trajectories; all records retained, including boundary-impacted records; not fitted overdamped mobility',
        overlapping_time_origins_are_independent=False,production_t0_seconds=None,
        input_summary_sha256=hashlib.sha256((out/'summary.json').read_bytes()).hexdigest(),
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))


if __name__=='__main__':main()
