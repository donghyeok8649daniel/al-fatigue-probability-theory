"""Temporal holdout analysis of completed, independently generated Al99 MD.

No candidate selection by favorable damping, and no production clock output.
The first 25 ps are excluded in every ensemble by protocol, not by fit quality.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from .mode_kinetic_calibration import (mode_correlations, fit_damped_correlation,
                                     formal_modal_mobility)
from .run_reference_thermostat_md import SOURCE_MD5
from .run_vector_registry_audit import save_json


def analyze(source, output):
    source, output = Path(source), Path(output)
    meta = json.loads((source/'summary.json').read_bytes())
    if not meta['completed'] or meta['potential_md5'] != SOURCE_MD5:
        raise ValueError('completed hash-bound reference trajectory required')
    if output.exists():
        raise FileExistsError('fresh analysis output required')
    file = source/'plane_coordinates.npz'
    with np.load(file, allow_pickle=False) as raw:
        q=raw['coordinates_m'];times=raw['time_seconds'];thermo=raw['thermo']
    dt=meta['frame_ps']
    if not np.allclose(np.diff(times),dt*1e-12,rtol=1e-10,atol=0):
        raise ValueError('source sampling interval mismatch')
    excluded=int(round(25./dt))
    q=q[excluded:];thermo=thermo[excluded:];times=times[excluded:]
    if len(q)<8192:
        raise ValueError('at least 8192 retained frames required')
    output.mkdir(parents=True)
    started=time.perf_counter()
    midpoint=len(q)//2
    halves=(q[:midpoint],q[midpoint:])
    cutoffs=(256,512,1024)
    correlations=[mode_correlations(part,max(cutoffs)) for part in halves]
    # Cubic [111] transverse degeneracy supplies a symmetry-pooled control.
    # Retain both original polarizations; pooling is not an independence claim.
    correlations=[np.concatenate([c,c[:,:,1:3].mean(axis=2,keepdims=True)],axis=2)
                  for c in correlations]
    spectra=[];frequencies=[]
    for part in halves:
        modes=np.fft.fft(part-part.mean(axis=0,keepdims=True),axis=1,norm='ortho')
        power=abs(np.fft.fft(modes,axis=0))**2
        spectra.append(np.concatenate([power,power[:,:,1:3].mean(axis=2,keepdims=True)],axis=2))
        frequencies.append(2*np.pi*np.fft.fftfreq(len(part),dt))
    records=[];planes=q.shape[1]
    for mode in range(1,planes//2+1):
        for axis,label in enumerate(('normal','direct110','transverse','transverse_symmetry_pooled')):
            seeds=[]
            for freq,power in zip(frequencies,spectra):
                pos=np.flatnonzero(freq>0)
                seeds.append(float(freq[pos[np.argmax(power[pos,mode,axis])]]))
            for cutoff in cutoffs:
                lagtimes=np.arange(cutoff+1)*dt
                fits=[];observed=[]
                for corr,seed in zip(correlations,seeds):
                    obs=corr[:cutoff+1,mode,axis]/corr[0,mode,axis]
                    observed.append(obs)
                    fits.append(fit_damped_correlation(lagtimes,obs,seed))
                fit,other=fits
                variance=float(correlations[0][0,mode,axis])
                temperature=float(thermo[:midpoint,0].mean())
                M=formal_modal_mobility(variance,
                    fit['formal_integral_time_ps']*1e-12,temperature)
                # q_k=(exp(ik)-1)u_k: mobility is a contravariant tensor.
                weight=4*np.sin(np.pi*mode/planes)**2
                row=dict(mode=mode,coordinate=label,cutoff_ps=cutoff*dt,
                    damping_per_ps=fit['damping_per_ps'],
                    second_half_damping_per_ps=other['damping_per_ps'],
                    frequency_rad_ps=fit['damped_angular_frequency_rad_ps'],
                    second_half_frequency_rad_ps=other['damped_angular_frequency_rad_ps'],
                    envelope_time_ps=fit['envelope_time_ps'],
                    formal_integral_time_ps=fit['formal_integral_time_ps'],
                    training_variance_m2=variance,
                    variance_ratio_second_to_first=float(correlations[1][0,mode,axis]/variance),
                    temperature_used_K=temperature,
                    second_half_temperature_K=float(thermo[midpoint:,0].mean()),
                    conditional_difference_mode_mobility_m2_per_J_s=M,
                    conditional_plane_mode_mobility_m2_per_J_s=M/weight,
                    training_rmse=fit['training_rmse'],
                    heldout_rmse=float(np.sqrt(np.mean((fit['prediction']-observed[1])**2))),
                    heldout_overdamped_rmse=float(np.sqrt(np.mean((
                        np.exp(-fit['overdamped_rate_per_ps']*lagtimes)-observed[1])**2))),
                    fit_success=bool(fit['optimizer_success'] and other['optimizer_success']),
                    active_bounds=fit['active_mask'],
                    production_mobility=False)
                records.append(row)
                save_json(output/f'mode_{mode}_{label}_lag_{cutoff}.json',dict(
                    summary=row,lag_times_ps=lagtimes,observed_halves=observed,
                    training_fit=fit,validation_sensitivity_fit=other))
            print(json.dumps(dict(mode=mode,coordinate=label,completed=True)),flush=True)
    z=np.fft.fft(q-q.mean(axis=0,keepdims=True),axis=1,norm='ortho')
    C=np.einsum('tki,tkj->kij',z,z.conj())/len(q)
    variance=C.real.diagonal(axis1=1,axis2=2)
    corr=C[1:]/np.sqrt(variance[1:,:,None]*variance[1:,None,:])
    elapsed_ps=(times-times[0])*1e12
    energy=thermo[:,3]/meta['atom_count']
    slope=float(np.polyfit(elapsed_ps,energy,1)[0])
    sanitized={k:v for k,v in meta.items() if k!='initialization_commands'}
    save_json(output/'summary.json',dict(completed=True,source_metadata=sanitized,
        projection_sha256=hashlib.sha256(file.read_bytes()).hexdigest(),
        excluded_initial_ps=25.,retained_duration_ps=float(elapsed_ps[-1]),
        temporal_halves_are_not_independent_replicas=True,
        records=records,mean_temperature_K=float(thermo[:,0].mean()),
        energy_slope_eV_per_atom_ps=slope,energy_range_eV_per_atom=float(np.ptp(energy)),
        max_normalized_modal_offdiagonal=float(np.max(abs(corr-np.eye(3)))),
        elapsed_seconds=time.perf_counter()-started,
        conditional_modal_friction_is_not_production_cell_friction=True,
        physical_PDE_seconds_enabled=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('source',type=Path);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();analyze(a.source,a.out)
