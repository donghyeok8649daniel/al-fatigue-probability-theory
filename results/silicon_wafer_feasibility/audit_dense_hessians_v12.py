"""Independent stored-matrix replay, units and stationary-residual diagnostics."""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh,solve

def main(args):
    if args.output.exists() and any(args.output.iterdir()):raise ValueError('fresh audit output required')
    args.output.mkdir(parents=True,exist_ok=True)
    rows=[];missing=[]
    for name in ('force100','loading8','return8','loading10'):
        folder=args.results/('dense_'+name)
        if not (folder/'summary.json').exists():missing.append(name);continue
        summary=json.loads((folder/'summary.json').read_text(encoding='utf-8'))
        if not summary['complete']:missing.append(name);continue
        protocol=json.loads((folder/'protocol.json').read_text(encoding='utf-8'))
        checks=json.loads((folder/'mode_checks.json').read_text(encoding='utf-8'))
        if (summary.get('independent_checked_directions')!=6 or len(checks)!=6
                or any(len(mode['force_differences'])!=2 or len(mode['energy_differences'])!=2 for mode in checks)):
            raise ValueError('full matrix exists but the six independent direction checks are incomplete: '+name)
        with np.load(folder/'raw_hessian.npz') as d:h=d['hessian'];basis=d['basis'];g=d['gradient']
        with np.load(folder/'spectrum.npz') as d:stored_values=d['eigenvalues'];vectors=d['eigenvectors']
        symmetric=(h+h.T)/2
        # scipy evr independently diagonalizes the saved matrix; original uses numpy eigh.
        values=eigh(symmetric,eigvals_only=True,driver='evr',check_finite=True)
        error=float(np.max(abs(values-stored_values)))
        residual=float(np.linalg.norm(symmetric@vectors-vectors*stored_values))
        raw_asymmetry=float(np.max(abs(h-h.T)))
        basis_error=float(np.max(abs(basis.T@basis-np.eye(len(h)))))
        if error>1e-9 or basis_error>1e-12:raise ValueError('independent eigensystem or basis replay failed')
        if abs(values[0]-summary['minimum_curvature_eV_A2'])>1e-9:raise ValueError('summary does not match raw matrix')
        force_error=max(check['error_norm_eV_A2'] for mode in checks for check in mode['force_differences'])
        force_half=max(mode['force_differences'][-1]['error_norm_eV_A2'] for mode in checks)
        energy_error=max(abs(check['curvature_eV_A2']-mode['rayleigh_eV_A2']) for mode in checks for check in mode['energy_differences'])
        row=dict(state=name,dimension=len(h),minimum_curvature_eV_A2=float(values[0]),
            maximum_curvature_eV_A2=float(values[-1]),negative_modes=int(np.sum(values<0)),
            basis_orthogonality_max_error=basis_error,independent_spectrum_max_error_eV_A2=error,
            eigen_residual_frobenius_eV_A2=residual,raw_asymmetry_max_eV_A2=raw_asymmetry,
            force_FD_max_error_norm_eV_A2=force_error,force_half_step_max_error_norm_eV_A2=force_half,
            energy_FD_max_error_eV_A2=energy_error,checked_directions=len(checks),
            reduced_gradient_norm_eV_A=float(np.linalg.norm(g)),reduced_gradient_max_eV_A=float(np.max(abs(g))),
            raw_hessian_sha256=hashlib.sha256((folder/'raw_hessian.npz').read_bytes()).hexdigest(),
            spectrum_sha256=hashlib.sha256((folder/'spectrum.npz').read_bytes()).hexdigest())
        if values[0]>0:
            newton=solve(symmetric,-g,assume_a='pos')
            cartesian=(basis@newton).reshape(-1,3)
            row['harmonic_stationarity_correction']=dict(coordinate_norm_A=float(np.linalg.norm(newton)),
                maximum_atom_displacement_A=float(np.linalg.norm(cartesian,axis=1).max()),
                predicted_energy_reduction_eV=float(-.5*g@newton),
                solve_residual_eV_A=float(np.linalg.norm(symmetric@newton+g)),
                status='linearized residual estimate only, not an actually relaxed state')
        if protocol['ensemble']=='force' and np.linalg.eigvalsh(symmetric[:-1,:-1])[0]>0:
            response=solve(symmetric[:-1,:-1],symmetric[:-1,-1],assume_a='pos')
            schur=float(symmetric[-1,-1]-symmetric[-1,:-1]@response)
            stiffness=schur*protocol['grip_extension_basis_norm']**2
            stored=summary['relaxed_grip_tangent']
            if abs(stiffness-stored['d_force_d_extension_eV_A2'])>1e-9:raise ValueError('grip Schur-complement conversion differs')
            row['relaxed_tangent']=stored
        rows.append(row)
    result=dict(complete=len(rows)==4,missing_or_incomplete=missing,states=rows,
        independent_method='scipy.linalg.eigh evr versus original numpy.linalg.eigh; raw matrix and basis replay',
        new_potential_calls=0,new_DFT=0,new_MD=0,material_approved=False,physical_clock=None,
        positive_computed_Hessian_does_not_certify_global_or_finite_temperature_stability=True)
    (args.output/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    if rows:
        fig,axes=plt.subplots(1,2,figsize=(11,4.6),layout='constrained')
        for row in rows:
            with np.load(args.results/('dense_'+row['state'])/'spectrum.npz') as d:values=d['eigenvalues']
            axes[0].plot(np.arange(1,min(25,len(values))+1),values[:25],'.-',label=row['state'])
        axes[0].axhline(0,color='#555',lw=.6);axes[0].set_xlabel('Ordered mode number');axes[0].set_ylabel('Curvature (eV / Angstrom^2)');axes[0].legend(fontsize=8)
        axes[1].bar(np.arange(len(rows)),[r['force_half_step_max_error_norm_eV_A2'] for r in rows],color='#087e8b')
        axes[1].set_xticks(np.arange(len(rows)),[r['state'] for r in rows],rotation=15);axes[1].set_yscale('log')
        axes[1].set_ylabel('Force FD check error norm (eV / Angstrom^2)')
        fig.suptitle('Full constrained Hessian and independent force differences\nStatic local diagnostic; no material, crack-initiation or kinetic certification')
        fig.savefig(args.output/'dense_hessian_checks.png',dpi=180);plt.close(fig)
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--results',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args())
