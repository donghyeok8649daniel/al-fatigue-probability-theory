"""Undamped-mode window leakage control for measured low-frequency power.

Not a statistical noise floor or certificate of a fitted memory kernel.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from .low_frequency_mobility import multitaper_spectrum,band_integral_proxy
from .mode_kinetic_calibration import mode_correlations
from .run_vector_registry_audit import save_json


def run(source,analysis,spectrum,output):
    source,analysis,spectrum,output=map(Path,(source,analysis,spectrum,output))
    if output.exists():raise FileExistsError('fresh output required')
    meta=json.loads((source/'summary.json').read_bytes())
    fits=json.loads((analysis/'summary.json').read_bytes())
    actual=json.loads((spectrum/'summary.json').read_bytes())
    if not meta['completed'] or not fits['completed'] or not actual['completed']:
        raise ValueError('completed studies required')
    digest=hashlib.sha256((source/'plane_coordinates.npz').read_bytes()).hexdigest()
    if digest!=fits['projection_sha256'] or digest!=actual['source_projection_sha256']:
        raise ValueError('all analyses must bind the same completed trajectory')
    with np.load(source/'plane_coordinates.npz',allow_pickle=False) as z:
        q=z['coordinates_m'][int(round(25/meta['frame_ps'])):]
    C=mode_correlations(q,1)[0];N=q.shape[1];dt=meta['frame_ps']*1e-12
    time=np.arange(len(q))*dt
    records=[]
    for nw in (2.,3.):
        selected=[r for r in actual['records'] if r['parts']==1 and r['nw']==nw and r['valid_finite_band']]
        bound=np.zeros((len(selected),3))
        for axis,label in enumerate(('normal','direct110','transverse')):
            for mode in range(1,N//2+1):
                # Both temporal-half frequencies from the longest DHO window.
                row=max([r for r in fits['records'] if r['mode']==mode and r['coordinate']==label],key=lambda r:r['cutoff_ps'])
                mode_bounds=[]
                for key in ('frequency_rad_ps','second_half_frequency_rad_ps'):
                    w=row[key]*1e12
                    basis=np.stack([np.cos(w*time),np.sin(w*time)],axis=1)
                    e=multitaper_spectrum(basis[:,None,:],dt,nw=nw,tapers=2)
                    mode_bounds.append([2*C[mode,axis]*np.linalg.eigvalsh(
                        band_integral_proxy(e,r['lower_cycle_THz']*1e12,r['upper_cycle_THz']*1e12)['integral_proxy_m2_seconds'])[-1]
                        for r in selected])
                weight=1 if N%2==0 and mode==N//2 else 2
                bound[:,axis]+=weight/N*np.max(mode_bounds,axis=0)
        for index,r in enumerate(selected):
            measured=np.diag(r['band']['integral_proxy_m2_seconds'])
            records.append(dict(nw=nw,lower_cycle_THz=r['lower_cycle_THz'],
                upper_cycle_THz=r['upper_cycle_THz'],
                undamped_fitted_modes_worst_phase_integral_proxy=bound[index],
                measured_integral_proxy=measured,
                leakage_to_measured_ratio=bound[index]/measured))
    output.mkdir(parents=True)
    save_json(output/'summary.json',dict(completed=True,records=records,
        interpretation='Taper leakage of undamped controls at fitted modal frequencies; not complete spectral uncertainty.',
        production_calibration_available=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('source');p.add_argument('--analysis',required=True)
    p.add_argument('--spectrum',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();run(a.source,a.analysis,a.spectrum,a.out)
