"""Follow BOTH signs of a detected negative fixed-cell Cartesian mode.

Preserve the stationary saddle and every failed minimization. This is a
downhill stability repair, not an activation-barrier or new DFT calculation.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
from scipy.linalg import eigh,null_space
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from results.silicon_wafer_feasibility.run_mace_boron_v9 import MODEL_SHA256
from results.silicon_wafer_feasibility.mace_force_only_v9 import force_only,force_calculator


def main(args):
    import torch
    from ase.io import read,write
    from ase.optimize import LBFGS,FIRE
    from mace.calculators import MACECalculator
    if hashlib.sha256(args.model.read_bytes()).hexdigest()!=MODEL_SHA256:raise ValueError('model mismatch')
    if args.output.exists() and any(args.output.iterdir()):raise ValueError('fresh output required')
    args.output.mkdir(parents=True,exist_ok=True)
    torch.set_num_threads(2);torch.set_num_interop_threads(1)
    model=MACECalculator(model_paths=str(args.model),device='cpu',default_dtype='float64')
    calc=force_calculator(model)
    atoms=read(args.state)
    h=np.load(args.hessian)['hessian']
    trans=np.tile(np.eye(3),(len(atoms),1))/np.sqrt(len(atoms));q=null_space(trans.T)
    values,vectors=eigh(q.T@(.5*(h+h.T))@q)
    if values[0]>=0:raise ValueError('no negative mode to follow')
    direction=q@vectors[:,0]
    initial_energy,initial_force=force_only(model,atoms)
    verification=[]
    for step in [1e-3,3e-4,1e-4,3e-5]:
        plus,minus=atoms.copy(),atoms.copy()
        plus.positions+=step*direction.reshape(-1,3);minus.positions-=step*direction.reshape(-1,3)
        ep,fp=force_only(model,plus);em,fm=force_only(model,minus)
        curvature=float(-direction@(fp-fm).ravel()/(2*step))
        verification.append(dict(step_A=step,force_curvature_eV_A2=curvature,
            energy_curvature_eV_A2=(ep-2*initial_energy+em)/step**2))
    if any(v['force_curvature_eV_A2']>=0 for v in verification):raise AssertionError('negative mode not reproducible')
    (args.output/'negative_mode_verification.json').write_text(json.dumps(dict(
        source=args.state.name,source_state_sha256=hashlib.sha256(args.state.read_bytes()).hexdigest(),
        source_hessian_sha256=hashlib.sha256(args.hessian.read_bytes()).hexdigest(),
        eigenvalue_eV_A2=float(values[0]),checks=verification),indent=2)+'\n',encoding='utf-8')
    rows=[]
    for sign,label in [(1,'plus'),(-1,'minus')]:
        started=time.perf_counter();before=calc.evaluations
        a=atoms.copy();a.positions+=sign*.05*direction.reshape(-1,3);a.calc=calc
        seed_energy=float(a.get_potential_energy())
        key=args.label+'_'+label
        write(args.output/(key+'_initial.extxyz'),a)
        attempts=[]
        opt=LBFGS(a,logfile=str(args.output/(key+'_lbfgs.log')),maxstep=.08)
        ok=opt.run(fmax=1e-5,steps=400)
        attempts.append(dict(method='LBFGS',steps=opt.nsteps,converged=bool(ok),
            force_max_eV_A=float(np.max(np.linalg.norm(a.get_forces(),axis=1))),energy_eV=float(a.get_potential_energy())))
        write(args.output/(key+'_lbfgs.extxyz'),a)
        if not ok:
            opt=FIRE(a,logfile=str(args.output/(key+'_fire.log')),maxstep=.05,dt=.02)
            ok=opt.run(fmax=1e-5,steps=500)
            attempts.append(dict(method='FIRE',steps=opt.nsteps,converged=bool(ok),
                force_max_eV_A=float(np.max(np.linalg.norm(a.get_forces(),axis=1))),energy_eV=float(a.get_potential_energy())))
        final_energy=float(a.get_potential_energy());final_force=a.get_forces().copy()
        write(args.output/(key+'_final.extxyz'),a)
        np.savez_compressed(args.output/(key+'_state.npz'),positions=a.positions,cell=a.cell.array,numbers=a.numbers,
                            energy=final_energy,forces=final_force)
        result=dict(structure=key,source_saddle=args.state.name,seed_mode_sign=sign,seed_displacement_norm_A=.05,
            saddle_energy_eV=initial_energy,seed_energy_eV=seed_energy,final_energy_eV=final_energy,
            energy_change_from_saddle_eV=final_energy-initial_energy,
            force_max_eV_A=float(np.max(np.linalg.norm(final_force,axis=1))),force_converged=bool(ok),
            attempts=attempts,calculator_evaluations=calc.evaluations-before,
            local_minimum_Hessian_checked=False,elapsed_reported_s=time.perf_counter()-started,
            status='negative-mode-following candidate; verify full Hessian before minimum claim')
        (args.output/(key+'_result.json')).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
        rows.append(result);print('MODE_FOLLOW',json.dumps(result),flush=True)
    (args.output/'summary.json').write_text(json.dumps(rows,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--model',type=Path,required=True)
    p.add_argument('--state',type=Path,required=True);p.add_argument('--hessian',type=Path,required=True)
    p.add_argument('--label',required=True);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args())
