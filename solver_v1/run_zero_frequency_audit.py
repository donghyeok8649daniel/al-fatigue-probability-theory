"""Evaluate actual public Al data at every declared low-frequency cutoff."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from .zero_frequency_kinetics import correlation_integral_audit
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json
from .report_public_calibration_evidence import save_figure


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source',type=Path);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    if args.out.exists():raise FileExistsError('fresh result directory required')
    metadata=json.loads((args.source/'source_metadata.json').read_bytes())
    file=args.source/'plane_coordinates.npz'
    with np.load(file,allow_pickle=False) as raw:
        time,q=raw['time_seconds'],raw['coordinates_m']
    if len(q)!=metadata['frames'] or not np.allclose(np.diff(time),metadata['frame_seconds'],rtol=1e-10,atol=0):
        raise ValueError('source physical sampling mismatch')
    plans=[('first_half',q[:len(q)//2],62,4),('full_4_blocks',q,126,4),
           ('full_8_blocks',q,62,8)]
    rows=[];details={};summaries=[]
    for name,values,lag,blocks in plans:
        result=correlation_integral_audit(values,frame_seconds=metadata['frame_seconds'],
            max_lag=lag,blocks=blocks);details[name]=result
        for i,t in enumerate(result['cutoff_seconds']):
            rows.append(dict(study=name,cutoff_seconds=t,
                minimum_integral_eigenvalue_seconds=result['integral_eigenvalues_seconds'][i,0],
                maximum_integral_eigenvalue_seconds=result['integral_eigenvalues_seconds'][i,-1],
                empirical_floor_seconds=result['empirical_floor_seconds'][i],
                half_cutoff_change_seconds=result['half_cutoff_change_seconds'][i],
                positive_integral_above_floor=bool(result['positive_integral_above_floor'][i]),
                block_sensitivity_seconds=result['block_sensitivity_seconds'][i],
                quadrature_change_seconds=result['quadrature_change_seconds'][i],
                antisymmetric_seconds=result['antisymmetric_seconds'][i]))
        summaries.append(dict(study=name,frames=len(values),maximum_lag_seconds=result['cutoff_seconds'][-1],
            final_integral_eigenvalues_seconds=result['integral_eigenvalues_seconds'][-1],
            final_floor_seconds=result['empirical_floor_seconds'][-1],
            final_half_cutoff_change_seconds=result['half_cutoff_change_seconds'][-1],
            late_positive_fraction=float(np.mean(result['positive_integral_above_floor'][len(result['cutoff_seconds'])//2:]))))
    write_csv(args.out/'integral_cutoff_summary.csv',rows)
    save_json(args.out/'integral_details.json',details)
    save_json(args.out/'decision.json',dict(completed=True,source_doi=metadata['doi'],
        projection_sha256=hashlib.sha256(file.read_bytes()).hexdigest(),summaries=summaries,
        coordinate='adjacent periodic (111) plane means; NOT production a/s',
        thermostat=metadata['thermostat'],source_MD_time_is_seconds=True,
        no_cutoff_selected_to_manufacture_positive_friction=True,
        M_a_phys=None,M_s_phys=None,t0_seconds=None,production_clock_calibrated=False,
        next_requirement='longer matched-coordinate, thermostat-audited data; converged integral and appropriate low-frequency PMF'))
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(12,3.8),layout='constrained')
    for axis,(name,data) in zip(axes,details.items()):
        x=data['cutoff_seconds']*1e12
        for j in range(3):axis.plot(x,data['integral_eigenvalues_seconds'][:,j]*1e12,label=f'eigenvalue {j+1}')
        axis.fill_between(x,-data['empirical_floor_seconds']*1e12,data['empirical_floor_seconds']*1e12,
            alpha=.15,color='k',label='block + quadrature + antisymmetry')
        axis.axhline(0,color='k',lw=.5)
        axis.set(title=name,xlabel='MD lag cutoff [ps]',ylabel='Whitened integral [ps]');axis.grid(alpha=.2)
    axes[0].legend(fontsize=6)
    fig.suptitle('Actual Al MD: integral sensitivity, NOT a production physical clock')
    save_figure(fig,args.out/'low_frequency_integral.svg')
    print(json.dumps(summaries,default=lambda x:x.tolist(),indent=2),flush=True)


if __name__=='__main__':main()
