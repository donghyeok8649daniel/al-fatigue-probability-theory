"""Local classical Gaussian diagnostics of two same-grip stationary candidates.

These are labelled single-well quadratic approximations, not computed finite-T
Si free energies. Basin truncation, anharmonicity, quantum nuclei and basin
multiplicity are not evaluated. They provide no opening rate or physical clock.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.constants import Boltzmann,electron_volt


def main(args):
    if args.output.exists():raise ValueError('fresh harmonic diagnostic output required')
    with np.load(args.geometry) as d:free=d['free']
    records=[];arrays={}
    for name in ('loading8','return8'):
        source=args.results/('dense_'+name)
        summary=json.loads((source/'summary.json').read_text(encoding='utf-8'))
        if not summary['complete']:raise ValueError('source Hessian incomplete')
        with np.load(source/'raw_hessian.npz') as d:
            h=(d['hessian']+d['hessian'].T)/2;g=d['gradient'];basis=d['basis']
            positions=d['positions'];numbers=d['numbers'];energy=float(d['energy'])
        values,vectors=np.linalg.eigh(h)
        if values[0]<=0:raise ValueError('positive local Hessian required; no Gaussian well assigned')
        inverse=(vectors/values)@vectors.T
        minimum_shift=-inverse@g
        correction=-.5*float(g@inverse@g)
        covariance_diag=np.sum((basis@inverse)*basis,axis=1).reshape(-1,3)
        if np.min(covariance_diag)<-1e-10:raise ValueError('invalid quadratic covariance')
        arrays[name]=dict(positions=positions,numbers=numbers,basis=basis,covariance_diag=covariance_diag)
        records.append(dict(state=name,dimension=len(h),energy_eV=energy,log_determinant_in_eV_A2_units=float(np.log(values).sum()),
            minimum_curvature_eV_A2=float(values[0]),linearized_minimum_energy_correction_eV=correction,
            linearized_maximum_atom_displacement_A=float(np.linalg.norm((basis@minimum_shift).reshape(-1,3),axis=1).max()),
            maximum_gradient_eV_A=float(np.max(abs(g))),
            hessian_inverse_residual_norm=float(np.linalg.norm(h@inverse-np.eye(len(h)))),
            raw_hessian_sha256=hashlib.sha256((source/'raw_hessian.npz').read_bytes()).hexdigest()))
    first,last=records
    if (first['dimension']!=last['dimension'] or not np.array_equal(arrays['loading8']['basis'],arrays['return8']['basis'])
            or not np.array_equal(arrays['loading8']['numbers'],arrays['return8']['numbers'])):
        raise ValueError('quadratic integration coordinates differ')
    grip_difference=float(np.max(abs(arrays['loading8']['positions'][~free]-arrays['return8']['positions'][~free])))
    if grip_difference>1e-12:raise ValueError('the two states do not have the same fixed grips')
    delta_energy=last['energy_eV']-first['energy_eV']
    log_ratio=last['log_determinant_in_eV_A2_units']-first['log_determinant_in_eV_A2_units']
    correction=last['linearized_minimum_energy_correction_eV']-first['linearized_minimum_energy_correction_eV']
    temperatures=[];kb=Boltzmann/electron_volt
    for temperature in (10.,50.,150.,300.):
        kt=kb*temperature
        row=dict(temperature_K=temperature,kBT_eV=kt,delta_energy_eV=delta_energy,
            quadratic_entropy_term_eV=.5*kt*log_ratio,
            raw_centered_quadratic_delta_F_eV=delta_energy+.5*kt*log_ratio,
            linearly_center_corrected_quadratic_delta_F_eV=delta_energy+correction+.5*kt*log_ratio)
        for name,a in arrays.items():
            atom_variance=kt*a['covariance_diag'].sum(axis=1)
            row[name+'_Gaussian_free_atom_RMS_A']=float(np.sqrt(atom_variance[free].mean()))
            row[name+'_Gaussian_maximum_atom_RMS_A']=float(np.sqrt(atom_variance.max()))
        temperatures.append(row)
    result=dict(complete=True,states=records,temperature_diagnostics=temperatures,grip_coordinate_max_difference_A=grip_difference,
        delta_energy_eV=delta_energy,log_determinant_ratio=log_ratio,
        delta_linearized_center_energy_correction_eV=correction,
        formula='Delta F_quadratic = Delta E + kBT/2 * log(det H_return / det H_load); same dimension and coordinate measure',
        interpretation='single labelled local quadratic wells only; no finite basin integration or physical free-energy certification',
        omitted=['anharmonicity','basin boundaries and truncation','quantum nuclear effects','configurational basin multiplicity','transition saddle','kinetic mobility'],
        temperatures_are_not_MD_runs=True,new_potential_calls=0,new_DFT=0,new_MD=0,material_approved=False,physical_clock=None)
    args.output.mkdir(parents=True)
    (args.output/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(10.5,4.4),layout='constrained')
    t=[r['temperature_K'] for r in temperatures]
    axes[0].plot(t,[r['raw_centered_quadratic_delta_F_eV'] for r in temperatures],'o-',label='Quadratic Delta F')
    axes[0].axhline(delta_energy,color='#666',ls='--',label='Stored Delta E')
    axes[0].set(xlabel='Assumed temperature (K)',ylabel='Return8 minus loading8 (eV)');axes[0].legend(fontsize=8)
    for name in ('loading8','return8'):
        axes[1].plot(t,[r[name+'_Gaussian_maximum_atom_RMS_A'] for r in temperatures],'o-',label=name)
    axes[1].set(xlabel='Assumed temperature (K)',ylabel='Maximum Gaussian atom RMS (Angstrom)');axes[1].legend(fontsize=8)
    fig.suptitle('Same-grip local classical quadratic diagnostic\nBasin validity, anharmonicity and quantum effects not established',fontsize=11)
    fig.savefig(args.output/'local_harmonic_diagnostic.png',dpi=180);plt.close(fig)
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('results','geometry','output'):p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
