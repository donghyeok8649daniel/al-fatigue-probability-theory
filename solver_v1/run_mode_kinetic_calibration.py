"""Fit actual Al MD modal decay with disjoint temporal validation.

No fit to fatigue life, yield, desired lag or crack probability. Every mode,
cutoff and deterministic start is retained, including unsuccessful fits.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from .mode_kinetic_calibration import (mode_correlations,fit_damped_correlation,
    damped_position_correlation,formal_modal_mobility)
from .run_vector_registry_audit import save_json


def run(source,output):
    source,output=Path(source),Path(output)
    if output.exists():raise FileExistsError('fresh modal calibration output required')
    meta=json.loads((source/'source_metadata.json').read_bytes())
    file=source/'plane_coordinates.npz'
    with np.load(file,allow_pickle=False) as raw:q=raw['coordinates_m'];t=raw['time_seconds']
    if (meta['doi']!='10.5281/zenodo.10014454' or meta['temperature_K']!=300
            or q.shape!=(8192,12,3) or t.shape!=(8192,)
            or not np.allclose(np.diff(t),meta['frame_seconds'],rtol=1e-10,atol=0)):
        raise ValueError('predeclared complete 8192-frame source study required')
    output.mkdir(parents=True)
    started=time.perf_counter();dt=meta['frame_seconds']*1e12
    save_json(output/'protocol.json',dict(source_metadata=meta,
        projection_sha256=hashlib.sha256(file.read_bytes()).hexdigest(),
        training_frames=[0,4096],validation_frames=[4096,8192],
        lag_counts=[126,254,510],mode_indices=[1,2,3,4,5,6],
        coordinates=['normal','direct110','transverse'],
        modes_are_independent_samples=False,
        fit_family='equilibrium underdamped position correlation',
        noise_model='equal normalized-correlation residuals; correlated lags',
        atomic_mass_used_to_set_solver_time=False))
    train=mode_correlations(q[:4096],510);valid=mode_correlations(q[4096:],510)
    blocks=np.array([mode_correlations(q[i*1024:(i+1)*1024],510) for i in range(8)])
    spectra=[]
    for part in (q[:4096],q[4096:]):
        z=np.fft.fft(part-part.mean(axis=0,keepdims=True),axis=1,norm='ortho')
        spectra.append(abs(np.fft.fft(z,axis=0))**2)
    omega=2*np.pi*np.fft.fftfreq(4096,dt)
    positive=np.flatnonzero(omega>0)
    summaries=[]
    for k in range(1,7):
        for axis,label in enumerate(('normal','direct110','transverse')):
            seed=float(omega[positive[np.argmax(spectra[0][positive,k,axis])]])
            seed_second=float(omega[positive[np.argmax(spectra[1][positive,k,axis])]])
            for lag in (126,254,510):
                times=np.arange(lag+1)*dt
                observed=train[:lag+1,k,axis]/train[0,k,axis]
                withheld=valid[:lag+1,k,axis]/valid[0,k,axis]
                fit=fit_damped_correlation(times,observed,seed)
                # Second-half refit is SENSITIVITY, not replacement for the
                # held-out prediction made with first-half parameters.
                refit=fit_damped_correlation(times,withheld,seed_second)
                block_curves=blocks[:,:lag+1,k,axis]/blocks[:,0,k,axis,None]
                empirical_rms=float(np.max(np.sqrt(np.mean((block_curves-withheld)**2,axis=1))))
                g=fit['damping_per_ps'];w=fit['damped_angular_frequency_rad_ps']
                mobility=formal_modal_mobility(train[0,k,axis],
                    fit['formal_integral_time_ps']*1e-12,meta['temperature_K'])
                row=dict(mode=k,coordinate=label,cutoff_ps=float(times[-1]),
                    damping_per_ps=g,envelope_time_ps=fit['envelope_time_ps'],
                    damped_angular_frequency_rad_ps=w,
                    formal_integral_time_ps=fit['formal_integral_time_ps'],
                    conditional_modal_mobility_m2_per_J_s=mobility,
                    mobility_is_phase_pooled_mode_not_cell=True,
                    second_half_damping_per_ps=refit['damping_per_ps'],
                    second_half_frequency_rad_ps=refit['damped_angular_frequency_rad_ps'],
                    training_rmse=fit['training_rmse'],
                    heldout_rmse=float(np.sqrt(np.mean((fit['prediction']-withheld)**2))),
                    heldout_overdamped_rmse=float(np.sqrt(np.mean((np.exp(-fit['overdamped_rate_per_ps']*times)-withheld)**2))),
                    block_shape_rms_sensitivity=empirical_rms,
                    validation_to_training_variance=float(valid[0,k,axis]/train[0,k,axis]),
                    singular_values=fit['singular_values_log_parameters'],
                    optimizer_success=fit['optimizer_success'],active_mask=fit['active_mask'],
                    production_clock_calibrated=False)
                summaries.append(row)
                save_json(output/f'mode_{k}_{label}_lag_{lag}.json',dict(
                    summary=row,times_ps=times,train_correlation=observed,
                    validation_correlation=withheld,training_fit=fit,
                    second_half_sensitivity_fit=refit,block_correlations=block_curves))
                print(json.dumps(row,default=lambda x:x.tolist()),flush=True)
    save_json(output/'summary.json',dict(completed=True,records=summaries,
        elapsed_seconds=time.perf_counter()-started,
        production_M_a_phys=None,production_M_s_phys=None,production_t0_seconds=None,
        fit_supplies_source_modal_parameters_not_a_production_calibration=True,
        source_thermostat_independence_tested=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source',type=Path);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();run(args.source,args.out)
