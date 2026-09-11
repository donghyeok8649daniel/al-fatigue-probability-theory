"""Mode-dependent low-band friction: test locality before cell projection.

Conjugate standing-wave phases are symmetry-pooled, not independent samples.
This leaves phase-offdiagonal dynamics uncalibrated. No production clock.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from .low_frequency_mobility import multitaper_spectrum,band_integral_proxy
from .zero_frequency_kinetics import formal_zero_frequency_response
from .run_vector_registry_audit import save_json
from .run_reference_thermostat_md import SOURCE_MD5


def run(source,output):
    source,output=Path(source),Path(output)
    if output.exists():raise FileExistsError('fresh output required')
    meta=json.loads((source/'summary.json').read_bytes())
    if not meta['completed'] or meta['potential_md5']!=SOURCE_MD5:
        raise ValueError('completed source-bound MD required')
    dt=meta['frame_ps']*1e-12;excluded=int(round(25e-12/dt))
    with np.load(source/'plane_coordinates.npz',allow_pickle=False) as file:
        q=file['coordinates_m'][excluded:];temperature=file['thermo'][excluded:,0]
    modes=np.fft.fft(q-q.mean(axis=0,keepdims=True),axis=1,norm='ortho')
    N=q.shape[1];records=[]
    for parts in (1,2):
        width=len(q)//parts
        for block in range(parts):
            sl=slice(block*width,(block+1)*width)
            T=float(temperature[sl].mean())
            for k in range(1,N//2+1):
                z=modes[sl,k]
                standing=(z.real[:,None,:] if k==N/2 else
                          np.sqrt(2)*np.stack([z.real,z.imag],axis=1))
                e=multitaper_spectrum(standing,dt,nw=3.,tapers=2)
                for lower,upper in ((.02,.04),(.04,.08),(.08,.16)):
                    band=band_integral_proxy(e,lower*1e12,upper*1e12)
                    response=formal_zero_frequency_response(e['covariance_m2'],
                        band['integral_proxy_m2_seconds'],1.380649e-23*T)
                    records.append(dict(parts=parts,block=block,mode=k,
                        lower_cycle_THz=lower,upper_cycle_THz=upper,
                        temperature_K=T,covariance_m2=e['covariance_m2'],
                        phase_pooled_response=response))
    output.mkdir(parents=True)
    save_json(output/'summary.json',dict(completed=True,records=records,
        source_projection_sha256=hashlib.sha256((source/'plane_coordinates.npz').read_bytes()).hexdigest(),
        atoms_per_plane=meta['atoms_per_plane'],planes=N,
        source_ensemble=meta['ensemble'],source_dt_ps=meta['dt_ps'],
        phase_pooling_is_symmetry_assumption=True,
        zero_frequency_limit_certified=False,production_clock_calibrated=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('source');p.add_argument('--out',required=True)
    a=p.parse_args();run(a.source,a.out)
