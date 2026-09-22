"""Matrix-free highest Gamma vibration, with full-Hessian and step controls.

The largest eigenpair is not a check of the lowest curvature or full stability.
Central force differences apply the mass-weighted Hessian to ARPACK vectors.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
from scipy.constants import electron_volt,atomic_mass,c
from scipy.sparse.linalg import LinearOperator,eigsh
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from results.silicon_wafer_feasibility.mace_force_only_v9 import force_only
from results.silicon_wafer_feasibility.run_mace_boron_v9 import MODEL_SHA256


def main(args):
    import torch
    from ase.io import read
    from mace.calculators import MACECalculator
    if hashlib.sha256(args.model.read_bytes()).hexdigest()!=MODEL_SHA256:raise ValueError('model mismatch')
    if len(args.states)!=len(args.labels):raise ValueError('label/state count mismatch')
    args.output.mkdir(parents=True,exist_ok=True)
    torch.set_num_threads(2);torch.set_num_interop_threads(1)
    calc=MACECalculator(model_paths=str(args.model),device='cpu',default_dtype='float64')
    rows=[];conversion=np.sqrt(electron_volt/1e-20/atomic_mass)/(2*np.pi*c*100)
    for index,(path,label) in enumerate(zip(args.states,args.labels)):
        output=args.output/(label+'_largest_mode.json')
        if output.exists():raise ValueError('refusing to overwrite completed mode')
        started=time.perf_counter();atoms=read(path);rootmass=np.sqrt(np.repeat(atoms.get_masses(),3))
        atoms.calc=calc;standard_force=atoms.get_forces();standard_energy=float(atoms.get_potential_energy())
        direct_energy,direct_force=force_only(calc,atoms)
        error=float(np.max(abs(standard_force-direct_force)))
        if error>1e-9 or abs(standard_energy-direct_energy)>1e-9:raise AssertionError('backend mismatch')
        evaluations=1;products=0
        def product(vector,step=1e-4):
            nonlocal evaluations,products
            v=np.asarray(vector).reshape(-1);norm=np.linalg.norm(v)
            if norm==0:return np.zeros_like(v)
            # The normalized mass-weighted direction fixes the displacement scale
            # independently of any ARPACK normalization or residual magnitude.
            displacement=(v/norm/rootmass).reshape(-1,3)
            plus,minus=atoms.copy(),atoms.copy()
            plus.positions+=step*displacement;minus.positions-=step*displacement
            value=-(force_only(calc,plus)[1]-force_only(calc,minus)[1]).ravel()/(2*step)/rootmass*norm
            evaluations+=2;products+=1
            if products%10==0:print('PRODUCT',label,products,flush=True)
            return value
        operator=LinearOperator((3*len(atoms),3*len(atoms)),matvec=product,dtype=float)
        seed=np.sin(np.arange(3*len(atoms))+.193);seed/=np.linalg.norm(seed)
        values,vectors=eigsh(operator,k=1,which='LA',v0=seed,ncv=20,tol=1e-7,maxiter=100)
        value=float(values[0]);v=vectors[:,0]
        if value<=0:raise AssertionError('highest mode not positive')
        checks=[]
        for step in [3e-4,1e-4,3e-5]:
            hv=product(v,step);rayleigh=float(v@hv)
            checks.append(dict(step_A_sqrt_u=step,Rayleigh_eV_A2_u=rayleigh,
                eigen_residual_L2_eV_A2_u=float(np.linalg.norm(hv-value*v)),
                Rayleigh_frequency_cm_minus1=float(conversion*np.sqrt(rayleigh))))
        reference=None
        if index==0 and args.reference_hessian:
            raw=np.load(args.reference_hessian)
            exact=float(raw['mass_weighted_projected_eigenvalues'][-1])
            reference=dict(full_Hessian_file=args.reference_hessian.name,
                full_Hessian_sha256=hashlib.sha256(args.reference_hessian.read_bytes()).hexdigest(),
                full_Hessian_eigenvalue=exact,absolute_eigenvalue_difference=abs(value-exact),
                full_Hessian_frequency_cm_minus1=float(conversion*np.sqrt(exact)))
            if abs(value-exact)>1e-6:raise AssertionError('full-Hessian control failed')
        record=dict(label=label,state_file=path.name,state_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            atoms=len(atoms),B_atoms=int(sum(atoms.numbers==5)),volume_A3=atoms.get_volume(),
            chemical_B_cm3=sum(atoms.numbers==5)/atoms.get_volume()*1e24,
            masses_u={str(int(z)):float(m) for z,m in zip(atoms.numbers,atoms.get_masses())},
            force_max_eV_A=float(np.max(np.linalg.norm(standard_force,axis=1))),
            highest_eigenvalue_eV_A2_u=value,highest_Gamma_frequency_cm_minus1=float(conversion*np.sqrt(value)),
            B_mass_weighted_mode_fraction=float(np.sum(v.reshape(-1,3)[atoms.numbers==5]**2)),
            independent_step_checks=checks,full_Hessian_control=reference,
            force_only_evaluations=evaluations,standard_ASE_evaluations=1,Hessian_vector_products=products,
            force_backend_max_difference_eV_A=error,elapsed_reported_s=time.perf_counter()-started,
            full_mechanical_stability_checked=False,material_calibrated=False,kinetic_calibrated=False,
            scope='highest Gamma eigenpair only; no full-q, lowest-curvature or finite-T certification')
        np.savez_compressed(args.output/(label+'_largest_mode.npz'),mass_weighted_eigenvector=v,
            eigenvalue=value,positions=atoms.positions,cell=atoms.cell.array,numbers=atoms.numbers,masses=atoms.get_masses())
        output.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
        rows.append(record);print('MODE',label,record['highest_Gamma_frequency_cm_minus1'],flush=True)
    (args.output/'summary.json').write_text(json.dumps(dict(cases=rows,complete=True,model_sha256=MODEL_SHA256,
        total_force_only_evaluations=sum(r['force_only_evaluations'] for r in rows),new_DFT_runs=0,new_MD_runs=0),indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--model',type=Path,required=True)
    p.add_argument('--states',type=Path,nargs='+',required=True);p.add_argument('--labels',nargs='+',required=True)
    p.add_argument('--reference-hessian',type=Path);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args())
