"""Finite-cell neutral B controls at fixed lattice parameter, no new DFT.

Host-subtracted energies across sizes cancel the arbitrary B energy zero.
The chemical concentration of the periodic cell is not a measured Hall number.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from results.silicon_wafer_feasibility.mace_force_only_v9 import force_only
from results.silicon_wafer_feasibility.run_mace_boron_v9 import MODEL_SHA256


def main(args):
    import torch
    from ase import Atoms
    from ase.build import bulk
    from ase.calculators.calculator import Calculator,all_changes
    from ase.io import write
    from ase.optimize import LBFGS,FIRE
    from mace.calculators import MACECalculator
    if hashlib.sha256(args.model.read_bytes()).hexdigest()!=MODEL_SHA256:raise ValueError('model mismatch')
    if args.output.exists() and any(args.output.iterdir()):raise ValueError('fresh output required')
    args.output.mkdir(parents=True,exist_ok=True)
    torch.set_num_threads(2);torch.set_num_interop_threads(1)
    model=MACECalculator(model_paths=str(args.model),device='cpu',default_dtype='float64')
    class ForceCalculator(Calculator):
        implemented_properties=['energy','free_energy','forces']
        def __init__(self):super().__init__();self.evaluations=0
        def calculate(self,atoms=None,properties=None,system_changes=all_changes):
            super().calculate(atoms,properties,system_changes)
            e,f=force_only(model,atoms);self.evaluations+=1
            self.results=dict(energy=e,free_energy=e,forces=f)
    calc=ForceCalculator()
    source=json.loads(args.structures.read_text(encoding='utf-8'))
    base=json.loads(args.baseline.read_text(encoding='utf-8'))
    lattice=5.431
    perfect=bulk('Si','diamond',a=lattice,cubic=True)
    ebulk,_=force_only(model,perfect);ebulk/=len(perfect)
    rows=[]
    for size in args.repeats:
        for key in ['b_silicon_64_interstitial','2boron_int_cluster']:
            started=time.perf_counter();begin=calc.evaluations
            s=source[key];symbols=np.array(s['symbols']);cell=np.asarray(s['cell_A'])
            r=np.asarray(s['fractional_positions'])@cell
            # Confirm source Si host and the generated diamond are the SAME lattice.
            oldhost=bulk('Si','diamond',a=lattice,cubic=True).repeat((2,2,2))
            delta=r[symbols=='Si'][:,None,:]-oldhost.positions[None,:,:]
            delta-=np.rint(delta/(2*lattice))*(2*lattice)
            host_match=np.max(np.min(np.linalg.norm(delta,axis=-1),axis=1))
            if host_match>1e-9 or np.sum(symbols=='Si')!=64:raise AssertionError('source host lattice mismatch')
            b=r[symbols=='B'];length=size*lattice
            shift=np.rint((length/2-b.mean(axis=0))/lattice)*lattice
            boron=(b+shift)%length
            atoms=bulk('Si','diamond',a=lattice,cubic=True).repeat((size,size,size))
            host_atoms=len(atoms)
            atoms+=Atoms('B'*len(boron),positions=boron)
            atoms.calc=calc
            initial_energy=float(atoms.get_potential_energy());initial_force=atoms.get_forces().copy()
            atoms.calc=model
            standard_energy=float(atoms.get_potential_energy());standard_force=atoms.get_forces()
            backend_error=float(np.max(abs(standard_force-initial_force)))
            if abs(standard_energy-initial_energy)>1e-9 or backend_error>1e-9:raise AssertionError('force backend mismatch')
            atoms.calc=calc
            tag=f'{key}_{host_atoms}host'
            attempts=[]
            write(args.output/(tag+'_initial.extxyz'),atoms)
            opt=LBFGS(atoms,logfile=str(args.output/(tag+'_lbfgs.log')),maxstep=.1)
            ok=opt.run(fmax=1e-4,steps=250)
            attempts.append(dict(method='LBFGS',steps=opt.nsteps,converged=bool(ok),
                max_force_eV_A=float(np.max(np.linalg.norm(atoms.get_forces(),axis=1))),energy_eV=float(atoms.get_potential_energy())))
            write(args.output/(tag+'_lbfgs.extxyz'),atoms)
            if not ok:
                opt=FIRE(atoms,logfile=str(args.output/(tag+'_fire.log')),maxstep=.08,dt=.03)
                ok=opt.run(fmax=1e-4,steps=350)
                attempts.append(dict(method='FIRE',steps=opt.nsteps,converged=bool(ok),
                    max_force_eV_A=float(np.max(np.linalg.norm(atoms.get_forces(),axis=1))),energy_eV=float(atoms.get_potential_energy())))
            energy=float(atoms.get_potential_energy());force=atoms.get_forces().copy()
            base_case=next(r for r in base['cases'] if r['structure']==key)
            excess=energy-host_atoms*ebulk
            base_excess=base_case['final']['energy_eV']-64*ebulk
            row=dict(structure=key,host_atoms=host_atoms,B_atoms=len(boron),total_atoms=len(atoms),
                lattice_A=lattice,volume_A3=atoms.get_volume(),chemical_B_cm3=len(boron)/atoms.get_volume()*1e24,
                energy_eV=energy,pure_host_energy_per_atom_eV=ebulk,host_subtracted_energy_eV=excess,
                host_subtracted_change_from_64host_eV=excess-base_excess,
                B_chemical_potential_specified=False,absolute_formation_energy_eV=None,
                force_max_eV_A=float(np.max(np.linalg.norm(force,axis=1))),
                converged=bool(ok),attempts=attempts,source_host_match_A=float(host_match),
                source_B_shift_by_lattice_translation_A=shift.tolist(),backend_force_error_eV_A=backend_error,
                force_only_calculator_evaluations=calc.evaluations-begin,standard_ASE_control_evaluations=1,
                local_minimum_Hessian_checked=False,elapsed_reported_s=time.perf_counter()-started,
                status='neutral ML finite-periodic-cell control; no Si calibration')
            np.savez_compressed(args.output/(tag+'_state.npz'),positions=atoms.positions,cell=atoms.cell.array,
                                numbers=atoms.numbers,forces=force,energy=energy)
            write(args.output/(tag+'_final.extxyz'),atoms)
            (args.output/(tag+'_result.json')).write_text(json.dumps(row,indent=2)+'\n',encoding='utf-8')
            rows.append(row);print('DILUTE',json.dumps(row),flush=True)
    (args.output/'summary.json').write_text(json.dumps(dict(cases=rows,model_sha256=MODEL_SHA256,
        scope='fixed-volume neutral interstitial cell-size diagnostic; not Hall concentration or actual doping strength',
        pure_bulk_reference_eV_per_atom=ebulk,complete=True,new_DFT_runs=0,new_MD_runs=0),indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--model',type=Path,required=True)
    p.add_argument('--structures',type=Path,required=True);p.add_argument('--baseline',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--repeats',type=int,nargs='+',default=[3])
    main(p.parse_args())
