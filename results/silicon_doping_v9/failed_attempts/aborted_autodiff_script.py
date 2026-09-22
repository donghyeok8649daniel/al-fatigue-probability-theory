"""Full Cartesian Hessian and independent force differences of neutral B minima.

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
        hessian=np.asarray(calc.get_hessian(atoms)).reshape(3*len(atoms),3*len(atoms))
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
            plus.calc=calc;minus.calc=calc
            fd=-(plus.get_forces()-minus.get_forces()).ravel()/(2*step)
            checks.append(dict(step_A=step,max_abs_error_eV_A2=float(np.max(abs(fd-hessian@v))),
                               L2_error_eV_A2=float(np.linalg.norm(fd-hessian@v))))
        initial=read(args.states/(name+'_initial.extxyz'))
        df=(atoms.positions-initial.positions)@np.linalg.inv(atoms.cell.array)
        displacement=(df-np.rint(df))@atoms.cell.array
        b_indices=np.flatnonzero(atoms.numbers==5)
        b_distances=atoms.get_all_distances(mic=True)[np.ix_(b_indices,b_indices)]
        output=dict(structure=name,energy_eV=energy,atoms=len(atoms),B_atoms=nB,
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
            Hessian_status='fixed-cell local minimum' if values[0]>0 else 'negative curvature; not a local minimum',
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
    p.add_argument('--cases',nargs='*');main(p.parse_args())
