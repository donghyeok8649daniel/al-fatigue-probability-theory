"""Test thermal-kinetic channels against low-frequency plane displacement.

Kinetic-energy correlation is not identification of a hydrodynamic heat mode;
potential energy/current and causality remain outside this diagnostic.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from .low_frequency_mobility import multitaper_spectrum,heldout_spectral_prediction
from .run_vector_registry_audit import save_json


def run(source,output, *, recorded_mass_amu=None,engine_mass_amu=None):
    source,output=Path(source),Path(output)
    if output.exists():raise FileExistsError('fresh output required')
    meta=json.loads((source/'summary.json').read_bytes())
    if not meta['completed'] or not meta['thermal_observables']:
        raise ValueError('completed thermal-observable MD required')
    file=source/'plane_coordinates.npz';dt=meta['frame_ps']*1e-12
    with np.load(file,allow_pickle=False) as z:
        excluded=int(round(25e-12/dt))
        q=z['coordinates_m'][excluded:]
        e=z['plane_internal_kinetic_eV_per_atom'][excluded:]
        raw=z['plane_kinetic_eV_per_atom'][excluded:]
        kinetic=z['thermo'][excluded:,2]
    original_discrepancy=np.max(abs(raw.mean(axis=1)*meta['atom_count']/kinetic-1))
    recorded_mass_amu=meta.get('plane_kinetic_postprocess_mass_amu',recorded_mass_amu)
    engine_mass_amu=meta.get('atomic_mass_amu',engine_mass_amu)
    if (recorded_mass_amu is None or engine_mass_amu is None
            or not np.isfinite(recorded_mass_amu+engine_mass_amu)
            or min(recorded_mass_amu,engine_mass_amu)<=0):
        raise ValueError('explicit actual and recorded mass required, including legacy data')
    correction=engine_mass_amu/recorded_mass_amu
    e=e*correction;raw=raw*correction
    discrepancy=np.max(abs(raw.mean(axis=1)*meta['atom_count']/kinetic-1))
    joined=np.concatenate([q,e[:,:,None],np.roll(e,-1,axis=1)[:,:,None]],axis=2)
    # One training scale for both halves: no target leakage by separate scaling.
    midpoint=len(joined)//2
    scale=joined[:midpoint].std(axis=(0,1))
    if np.min(scale)<=0:raise ValueError('nondegenerate channels required')
    joined=joined/scale
    rows=[]
    for negative_control in (False,True):
        spectra=[]
        for half in (joined[:midpoint],joined[midpoint:]):
            x=half.copy()
            if negative_control:x[:,:,3:]=np.roll(x[:,:,3:],len(x)//3,axis=0)
            spectra.append(multitaper_spectrum(x,dt,nw=3.,tapers=2))
        for lower,upper in ((.02,.04),(.04,.08),(.08,.16)):
            matrices=[]
            for spec in spectra:
                mask=(spec['frequency_hz']>=lower*1e12)&(spec['frequency_hz']<upper*1e12)
                if mask.sum()<3:raise ValueError('resolved bands required')
                matrices.append(spec['spectrum_m2_seconds'][mask].mean(axis=0))
            for axis,label in enumerate(('normal','direct110','transverse')):
                for train in (0,1):
                    rows.append(dict(axis=label,training_half=train,
                        lower_cycle_THz=lower,upper_cycle_THz=upper,
                        time_shift_negative_control=negative_control,
                        prediction=heldout_spectral_prediction(matrices[train],matrices[1-train],axis,[3,4])))
    output.mkdir(parents=True)
    save_json(output/'summary.json',dict(completed=True,records=rows,
        source_projection_sha256=hashlib.sha256(file.read_bytes()).hexdigest(),
        max_kinetic_energy_relative_normalization_error=float(discrepancy),
        original_kinetic_energy_relative_error=float(original_discrepancy),
        recorded_postprocess_mass_amu=recorded_mass_amu,engine_mass_amu=engine_mass_amu,
        kinetic_energy_postprocess_correction_factor=correction,
        energy_channels='internal kinetic energy per atom in the two adjacent planes',
        scaled_spectrum_units='seconds, dimensionless coordinates scaled by first-half standard deviation',
        thermal_causality_established=False,production_calibration_available=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('source');p.add_argument('--out',required=True)
    p.add_argument('--recorded-mass-amu',type=float);p.add_argument('--engine-mass-amu',type=float)
    a=p.parse_args();run(a.source,a.out,recorded_mass_amu=a.recorded_mass_amu,
                        engine_mass_amu=a.engine_mass_amu)
