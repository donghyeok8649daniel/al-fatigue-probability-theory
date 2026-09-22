"""Full Cartesian force-difference Hessian of neutral B stationary candidates.

Only Gamma modes of the finite periodic supercell. Not a complete phonon
dispersion, finite-T free energy, carrier dynamics, or a probability clock.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
from scipy.linalg import null_space
from scipy.constants import electron_volt,atomic_mass,c
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from results.silicon_wafer_feasibility.run_mace_boron_v9 import MODEL_SHA256
from results.silicon_wafer_feasibility.mace_force_only_v9 import force_only


def main(args):
    import torch
    from ase.io import read
    from mace.calculators import MACECalculator
    if hashlib.sha256(args.model.read_bytes()).hexdigest()!=MODEL_SHA256:
        raise ValueError('model mismatch')
    torch.set_num_threads(2);torch.set_num_interop_threads(1)
    calc=MACECalculator(model_paths=str(args.model),device='cpu',default_dtype='float64')
    args.output.mkdir(parents=True,exist_ok=True)
    summary=[]
    names=args.cases or [p.name.removesuffix('_final.extxyz') for p in sorted(args.states.glob('*_final.extxyz'))]
    conversion=np.sqrt(electron_volt/1e-20/atomic_mass)/(2*np.pi*c*100)
    for name in names:
        record_path=args.output/(name+'_hessian.json')
        if record_path.exists():raise ValueError('refusing to overwrite Hessian evidence')
        start=time.perf_counter()
        atoms=read(args.states/(name+'_final.extxyz'));atoms.calc=calc
        force=atoms.get_forces();energy=float(atoms.get_potential_energy())
        backend_check=None
        if args.force_backend=='direct':
            direct_energy,direct_force=force_only(calc,atoms)
            backend_check=dict(energy_difference_eV=direct_energy-energy,
                               max_force_difference_eV_A=float(np.max(abs(direct_force-force))))
            if abs(direct_energy-energy)>1e-9 or np.max(abs(direct_force-force))>1e-9:
                raise AssertionError('force-only backend differs from ASE')
        def evaluate_force(a):
            if args.force_backend=='direct':return force_only(calc,a)[1]
            a.calc=calc
            return a.get_forces()
        step_hessian=1e-4
        hessian=np.empty((3*len(atoms),3*len(atoms)))
        for column in range(3*len(atoms)):
            plus,minus=atoms.copy(),atoms.copy()
            plus.positions.ravel()[column]+=step_hessian
            minus.positions.ravel()[column]-=step_hessian
            hessian[:,column]=-(evaluate_force(plus)-evaluate_force(minus)).ravel()/(2*step_hessian)
            if (column+1)%30==0:
                print('HESSIAN_COLUMN',name,column+1,len(hessian),flush=True)
        asym=float(np.max(abs(hessian-hessian.T)))
        # Remove precisely the translation subspace, not the lowest three
        # eigenvalues (which could include genuine unstable modes).
        trans=np.tile(np.eye(3),(len(atoms),1))/np.sqrt(len(atoms))
        basis=null_space(trans.T)
        values=np.linalg.eigvalsh(basis.T@(.5*(hessian+hessian.T))@basis)
        sqrt_mass=np.sqrt(np.repeat(atoms.get_masses(),3))
        dynamical=hessian/sqrt_mass[:,None]/sqrt_mass[None,:]
        mass_trans=sqrt_mass[:,None]*trans
        mass_basis=null_space(mass_trans.T)
        w2=np.linalg.eigvalsh(mass_basis.T@(.5*(dynamical+dynamical.T))@mass_basis)
        nB=int(sum(atoms.numbers==5))
        positive=w2[w2>0]
        top=(conversion*np.sqrt(positive[-3*nB:])).tolist()
        # Deterministic, nonrandom direction; full translation-projected vector.
        v=np.sin(np.arange(3*len(atoms))+0.31)
        v=basis@(basis.T@v);v/=np.linalg.norm(v)
        checks=[]
        for step in [1e-3,3e-4,1e-4]:
            plus,minus=atoms.copy(),atoms.copy()
            plus.positions+=step*v.reshape(-1,3);minus.positions-=step*v.reshape(-1,3)
            fd=-(evaluate_force(plus)-evaluate_force(minus)).ravel()/(2*step)
            checks.append(dict(step_A=step,max_abs_error_eV_A2=float(np.max(abs(fd-hessian@v))),
                               L2_error_eV_A2=float(np.linalg.norm(fd-hessian@v))))
        initial=read(args.states/(name+'_initial.extxyz'))
        df=(atoms.positions-initial.positions)@np.linalg.inv(atoms.cell.array)
        displacement=(df-np.rint(df))@atoms.cell.array
        b_indices=np.flatnonzero(atoms.numbers==5)
        b_distances=atoms.get_all_distances(mic=True)[np.ix_(b_indices,b_indices)]
        output=dict(structure=name,energy_eV=energy,atoms=len(atoms),B_atoms=nB,
            Hessian_method='central Cartesian force differences; one pair at a time',
            Hessian_step_A=step_hessian,force_evaluations=2*3*len(atoms)+7+(args.force_backend=='direct'),
            force_backend=args.force_backend,force_backend_comparison=backend_check,
            max_force_eV_A=float(np.max(np.linalg.norm(force,axis=1))),
            Hessian_max_asymmetry_eV_A2=asym,translation_residual_eV_A2=float(np.max(abs(hessian@trans))),
            Cartesian_projected_min_eigenvalue_eV_A2=float(values[0]),
            Cartesian_projected_first_10_eigenvalues_eV_A2=values[:10].tolist(),
            mass_weighted_min_eigenvalue=float(w2[0]),
            highest_3NB_Gamma_frequencies_cm_minus1=top,
            frequency_conversion_cm_minus1=conversion,
            independent_Hessian_vector_checks=checks,
            max_initial_to_final_minimum_image_displacement_A=float(np.max(np.linalg.norm(displacement,axis=1))),
            B_pair_distances_A=b_distances[np.triu_indices(nB,1)].tolist(),
            imposed_charge='none; neutral ML model',
            Hessian_status=('positive fixed-cell curvature and force tolerance passed'
                if values[0]>0 and np.max(np.linalg.norm(force,axis=1))<=1e-4
                else ('negative curvature; not a local minimum' if values[0]<=0
                      else 'positive curvature but force tolerance failed')),
            scope='Gamma supercell only; all cell strains excluded; finite-T and full-q stability unverified',
            elapsed_reported_s=time.perf_counter()-start)
        np.savez_compressed(args.output/(name+'_hessian.npz'),hessian=hessian,
            Cartesian_projected_eigenvalues=values,mass_weighted_projected_eigenvalues=w2)
        record_path.write_text(json.dumps(output,indent=2)+'\n',encoding='utf-8')
        summary.append(output)
        print('HESSIAN',name,values[0],'force',output['max_force_eV_A'],flush=True)
    (args.output/'hessian_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--model',type=Path,required=True)
    p.add_argument('--states',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--force-backend',choices=['ase','direct'],default='ase')
    p.add_argument('--cases',nargs='*');main(p.parse_args())
