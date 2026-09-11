"""Predeclared finite-band spectral study of completed reference MD.

Full record, halves and quarters; two taper widths and three low bands.
No favorable-band selection, automatic certification, or production clock.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from .low_frequency_mobility import multitaper_spectrum, band_integral_proxy
from .zero_frequency_kinetics import formal_zero_frequency_response
from .run_reference_thermostat_md import SOURCE_MD5
from .run_vector_registry_audit import save_json


def analyze(source, output, *, lower_limit_study=False, prefix_ps=None,
            long_record_study=False):
    source,output=Path(source),Path(output)
    meta=json.loads((source/'summary.json').read_bytes())
    if not meta['completed'] or meta['potential_md5']!=SOURCE_MD5:
        raise ValueError('completed source-bound MD required')
    if output.exists(): raise FileExistsError('fresh output required')
    file=source/'plane_coordinates.npz'
    with np.load(file,allow_pickle=False) as z:
        q=z['coordinates_m'];T=z['thermo'][:,0];times=z['time_seconds']
    dt=meta['frame_ps']*1e-12
    if not np.allclose(np.diff(times),dt,rtol=1e-10,atol=0):
        raise ValueError('uniform source sampling required')
    if prefix_ps is not None:
        steps=prefix_ps/meta['frame_ps']
        if (not np.isfinite(steps) or steps<1 or abs(steps-round(steps))>1e-8
                or round(steps)>=len(q)):
            raise ValueError('available positive frame-aligned prefix required')
        q=q[:int(round(steps))+1];T=T[:len(q)]
    excluded=int(round(25e-12/dt));q=q[excluded:];T=T[excluded:]
    records=[]
    for parts in (1,2,4):
        width=len(q)//parts
        for block in range(parts):
            lo=block*width;hi=lo+width
            temperature=float(T[lo:hi].mean())
            for nw in (2.,3.):
                estimate=multitaper_spectrum(q[lo:hi],dt,nw=nw,tapers=2)
                bands=((.02,.04),(.04,.08),(.08,.16))
                if lower_limit_study or long_record_study:
                    bands=((.005,.01),(.01,.02))+bands
                if long_record_study:
                    bands=((.001,.002),(.002,.005))+bands
                for lower,upper in bands:
                    base=dict(parts=parts,block=block,nw=nw,tapers=2,
                        duration_ps=width*dt*1e12,temperature_K=temperature,
                        lower_cycle_THz=lower,upper_cycle_THz=upper,
                        half_bandwidth_cycle_THz=estimate['half_bandwidth_hz']*1e-12)
                    try:
                        band=band_integral_proxy(estimate,lower*1e12,upper*1e12)
                        response=formal_zero_frequency_response(estimate['covariance_m2'],
                            band['integral_proxy_m2_seconds'],1.380649e-23*temperature)
                        records.append(dict(**base,valid_finite_band=True,band=band,
                            covariance_m2=estimate['covariance_m2'],response=response))
                    except ValueError as error:
                        records.append(dict(**base,valid_finite_band=False,reason=str(error)))
    output.mkdir(parents=True)
    save_json(output/'summary.json',dict(completed=True,records=records,
        source_projection_sha256=hashlib.sha256(file.read_bytes()).hexdigest(),
        followup_lower_frequency_study=lower_limit_study,
        followup_long_record_study=long_record_study,prefix_ps=prefix_ps,
        source_ensemble=meta['ensemble'],source_dt_ps=meta['dt_ps'],
        source_repeats=meta['repeats'],source_lattice_angstrom=meta['lattice_angstrom'],
        interpretation='Finite-band local plane-gap PMF response; positivity is not a zero-frequency plateau certificate.',
        zero_frequency_limit_certified=False,production_clock_calibrated=False))
    print(json.dumps(dict(completed=True,records=len(records),output=str(output))))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('source',type=Path);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--lower-limit-study',action='store_true')
    p.add_argument('--long-record-study',action='store_true')
    p.add_argument('--prefix-ps',type=float)
    a=p.parse_args();analyze(a.source,a.out,lower_limit_study=a.lower_limit_study,
                           prefix_ps=a.prefix_ps,long_record_study=a.long_record_study)
