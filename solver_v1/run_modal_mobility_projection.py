"""Calculate dimensional mobility from the completed source-mode fits.

Distinguishes a single periodic plane-gap PMF from a local intensive cell-law
hypothesis. The latter is not a certified production a/s mobility. No favorable
cutoff is selected: every predefined cutoff and both temporal halves are saved.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from .mode_kinetic_calibration import mode_correlations
from .run_vector_registry_audit import save_json


def scalar_mixture_mobility(variances, integral_times_s, *, multiplicities,
                            planes, temperature_K):
    C=np.asarray(variances,float);tau=np.asarray(integral_times_s,float)
    multiplicity=np.asarray(multiplicities,float)
    if (C.shape!=tau.shape or C.shape!=multiplicity.shape or C.ndim!=1
            or np.any(~np.isfinite(C+tau+multiplicity))
            or np.any(C<=0) or np.any(tau<=0) or np.any(multiplicity<=0)
            or planes<3 or not np.isfinite(temperature_K) or temperature_K<=0):
        raise ValueError('positive resolved modal covariance/time and explicit counts required')
    C0=float(multiplicity@C/planes)
    K=float(multiplicity@(C*tau)/planes)
    friction=1.380649e-23*temperature_K*K/C0**2
    return dict(variance_m2=C0,correlation_integral_m2_s=K,
                integral_time_s=K/C0,friction_J_s_per_m2=friction,
                mobility_m2_per_J_s=1/friction)


def run(source, fits, output):
    source,fits,output=Path(source),Path(fits),Path(output)
    if output.exists():raise FileExistsError('fresh output required')
    metadata=json.loads((source/'source_metadata.json').read_bytes())
    with np.load(source/'plane_coordinates.npz',allow_pickle=False) as data:
        q=data['coordinates_m']
    if q.shape!=(8192,12,3) or metadata['temperature_K']!=300:
        raise ValueError('matched completed source projection required')
    covariance=[mode_correlations(q[:4096],1),mode_correlations(q[4096:],1)]
    rows=[]
    for cutoff in (126,254,510):
        for half in (0,1):
            row=dict(cutoff_ps=cutoff*.025,temporal_half=half,
                source_temperature_K=300,atoms_per_plane=576,
                scalar_axis_approximation=True,production_calibrated=False)
            estimates={}
            for axis,label in enumerate(('normal','direct110','transverse')):
                values=[]
                for mode in range(1,7):
                    record=json.loads((fits/f'mode_{mode}_{label}_lag_{cutoff}.json').read_bytes())
                    fit=record['training_fit' if half==0 else 'second_half_sensitivity_fit']
                    values.append(fit['formal_integral_time_ps']*1e-12)
                estimate=scalar_mixture_mobility(covariance[half][0,1:7,axis],values,
                    multiplicities=[2,2,2,2,2,1],planes=12,temperature_K=300)
                # Exact rescaling if F_bar=F_plane/576 AND kBT_bar=kBT/576.
                # Interpreting F_bar/M_bar as a LOCAL production cell is the hypothesis.
                estimate['intensive_cell_mobility_hypothesis_m2_per_J_s']=576*estimate['mobility_m2_per_J_s']
                estimates[label]=estimate
            row['axes']=estimates
            row['slip_to_normal_mobility_ratio']=estimates['direct110']['mobility_m2_per_J_s']/estimates['normal']['mobility_m2_per_J_s']
            rows.append(row)
    output.mkdir(parents=True)
    save_json(output/'mobility.json',dict(completed=True,
        formula='Gamma=kBT integral(C)/C0^2; M=1/Gamma; diagonal scalar PMF approximation',
        units='M:m^2/(J s); Gamma:J s/m^2',records=rows,
        normalization='q_l is a difference of 576-atom periodic plane averages',
        conditional_on_fitted_DHO_zero_frequency_extrapolation=True,
        local_cell_transfer_unvalidated=True,
        energy_normalization_identity='F_bar=F_plane/Np and thermal_energy_bar=kBT/Np imply M_bar=Np*M_plane for the same stochastic q; local-cell transfer at unchanged kBT is unvalidated',
        source='10.5281/zenodo.10014454, actual 8192 frames',
        production_M_a_phys=None,production_M_s_phys=None,production_t0_seconds=None))
    for row in rows:
        print(json.dumps(row),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('source',type=Path);p.add_argument('--fits',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();run(a.source,a.fits,a.out)
