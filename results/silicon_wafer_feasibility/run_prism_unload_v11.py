"""Unload a saved converged prism to a specified smaller grip displacement.

The old run remains immutable. A reference-affine predictor moves the free atoms,
then the new grips remain fixed during relaxation. This is a static branch check,
not a time history, fatigue cycle, crack detector, or proof of physical plasticity.
"""
from __future__ import annotations
import argparse, hashlib, json, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from solver_v1.silicon_initiation_research import prescribed_grips, grip_observables, connectivity_diagnostics
from results.silicon_wafer_feasibility.mace_force_only_v9 import force_calculator
from results.silicon_wafer_feasibility.run_mace_boron_v9 import MODEL_SHA256


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def main(args):
    import torch
    from ase import Atoms
    from ase.constraints import FixAtoms
    from ase.io import write
    from ase.optimize import LBFGSLineSearch, FIRE
    from mace.calculators import MACECalculator
    if hashlib.sha256(args.model.read_bytes()).hexdigest()!=MODEL_SHA256:
        raise ValueError('unexpected model weights')
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError('fresh output required')
    with np.load(args.geometry) as d:
        reference=d['reference'].copy(); lower=d['lower'].copy(); upper=d['upper'].copy(); free=d['free'].copy()
        area=float(d['area_A2']); gauge=float(d['gauge_length_A'])
    with np.load(args.state) as d:
        positions=d['positions'].copy(); numbers=d['numbers'].copy(); cell=d['cell'].copy()
        saved_energy=float(d['energy']); saved_force=d['forces'].copy()
    source_result=json.loads(args.result.read_text(encoding='utf-8'))
    if not source_result['converged'] or abs(saved_energy-source_result['energy_eV'])>1e-10:
        raise ValueError('matching converged source result required')
    if (not np.all(numbers==14) or positions.shape!=reference.shape
            or not np.array_equal(free,~(lower|upper))):
        raise ValueError('matching pure-Si geometry and masks required')
    source_extension=float(np.mean(positions[upper,2]-reference[upper,2])-np.mean(positions[lower,2]-reference[lower,2]))
    source_target=prescribed_grips(reference,lower=lower,upper=upper,extension_A=source_extension)
    extension=float(args.target_strain*gauge)
    if not np.isfinite(extension) or not 0<=extension<source_extension:
        raise ValueError('nonnegative smaller target strain required')
    target=prescribed_grips(reference,lower=lower,upper=upper,extension_A=extension)
    fixed=lower|upper
    if np.max(abs(positions[fixed]-source_target[fixed]))>1e-12:
        raise ValueError('saved state does not obey this fixed grip displacement')
    args.output.mkdir(parents=True,exist_ok=True)
    protocol=dict(model_sha256=MODEL_SHA256,source_raw_sha256=hashlib.sha256(args.state.read_bytes()).hexdigest(),
        source_result_sha256=hashlib.sha256(args.result.read_bytes()).hexdigest(),
        geometry_sha256=hashlib.sha256(args.geometry.read_bytes()).hexdigest(),
        source_extension_A=source_extension,source_strain=source_extension/gauge,
        extension_A=extension,strain=extension/gauge,fmax_eV_A=args.fmax,
        max_steps=args.max_steps,max_seconds=args.max_seconds,
        displacement_changed=True,new_initial_perturbation=False,initial_crack_seed_added=False,
        predictor='free atom z shift = clipped initial-reference axial fraction times change of grip extension',
        scope='static lower-displacement branch from a converged higher-displacement state, optimizer memory restarted',
        first_crack_event_certified=False,physical_clock=None)
    dump(args.output/'protocol.json',protocol)
    torch.set_num_threads(2);torch.set_num_interop_threads(1)
    model=MACECalculator(model_paths=str(args.model),device='cpu',default_dtype='float64')
    calc=force_calculator(model)
    atoms=Atoms(numbers=numbers,positions=positions,cell=cell,pbc=False)
    atoms.set_constraint(FixAtoms(mask=fixed));atoms.calc=calc
    started=time.perf_counter();checks=[];interruption=None;converged=False
    e0=float(atoms.get_potential_energy()); f0=atoms.get_forces(apply_constraint=False)
    replay=dict(energy_difference_eV=e0-saved_energy,force_difference_eV_A=float(abs(f0-saved_force).max()),
        fixed_coordinate_error_A=float(abs(atoms.positions[fixed]-positions[fixed]).max()))
    dump(args.output/'source_replay.json',replay)
    if abs(replay['energy_difference_eV'])>1e-8 or replay['force_difference_eV_A']>1e-8:
        raise ValueError('saved starting state does not reproduce with the pinned model')
    zlo,zhi=reference[lower,2].mean(),reference[upper,2].mean()
    fraction=np.clip((reference[:,2]-zlo)/(zhi-zlo),0.,1.)-.5
    predicted=positions.copy()
    predicted[free,2]+=fraction[free]*(extension-source_extension)
    predicted[fixed]=target[fixed]
    atoms.set_positions(predicted,apply_constraint=False)
    predictor_energy=float(atoms.get_potential_energy())
    predictor_force=atoms.get_forces(apply_constraint=False)
    np.savez_compressed(args.output/'predictor.npz',positions=atoms.positions,energy=predictor_energy,
        forces=predictor_force,numbers=numbers,cell=cell)
    write(args.output/'predictor.extxyz',atoms)
    try:
        for method in ('LBFGSLineSearch','FIRE'):
            opt=(LBFGSLineSearch(atoms,logfile=str(args.output/'linesearch.log'),maxstep=.08)
                 if method=='LBFGSLineSearch' else FIRE(atoms,logfile=str(args.output/'fire.log'),maxstep=.04,dt=.02))
            def checkpoint():
                if opt.nsteps%20==0:
                    np.savez_compressed(args.output/'checkpoint.npz',positions=atoms.positions,
                        numbers=numbers,cell=cell)
                    dump(args.output/'running.json',dict(method=method,step=opt.nsteps,
                        force_calls=calc.evaluations,elapsed_s=time.perf_counter()-started))
                if time.perf_counter()-started>=args.max_seconds:
                    raise TimeoutError('explicit unload wall-time budget')
            opt.attach(checkpoint,interval=1)
            converged=bool(opt.run(fmax=args.fmax,steps=args.max_steps))
            checks.append(dict(method=method,steps=opt.nsteps,converged=converged))
            if converged:break
    except TimeoutError as exc:
        interruption=str(exc);converged=False
    energy=float(atoms.get_potential_energy());force=atoms.get_forces(apply_constraint=False)
    observables=grip_observables(atoms.positions,force,lower=lower,upper=upper,area_A2=area)
    grip_error=float(abs(atoms.positions[fixed]-target[fixed]).max())
    if grip_error>1e-12:raise ValueError('fixed grips moved during unload relaxation')
    if converged and observables['free_force_max_eV_A']>args.fmax:
        raise ValueError('optimizer convergence disagrees with raw free forces')
    np.savez_compressed(args.output/'raw.npz',positions=atoms.positions,forces=force,energy=energy,numbers=numbers,cell=cell)
    write(args.output/'state.extxyz',atoms)
    row=dict(converged=converged,complete=converged,interruption=interruption,attempts=checks,
        strain=extension/gauge,extension_A=extension,energy_eV=energy,source_energy_eV=saved_energy,
        predictor_energy_eV=predictor_energy,
        energy_change_during_relaxation_eV=energy-predictor_energy,
        source_to_final_energy_change_eV=energy-saved_energy,
        energy_difference_is_physical_dissipation=False,grip_coordinate_error_A=grip_error,
        force_calls=calc.evaluations,elapsed_s=time.perf_counter()-started,**observables,
        connectivity=connectivity_diagnostics(atoms.positions,lower=lower,upper=upper),
        new_DFT=0,new_MD=0,first_crack_event_certified=False,initiation_probability=None,physical_clock=None)
    dump(args.output/'summary.json',row);print(json.dumps(row,indent=2))
    if not converged:raise RuntimeError('unload relaxation did not converge; raw state preserved')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('state','result','geometry','model','output'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--target-strain',type=float,required=True)
    p.add_argument('--fmax',type=float,default=.001);p.add_argument('--max-steps',type=int,default=600)
    p.add_argument('--max-seconds',type=float,default=3000);main(p.parse_args())
