"""Force-free and specified nominal-load states of the intact finite prism.

Minimize U(r,Delta)-F_external*Delta in the same explicit grip coordinate.
The zero imposed extension reference may carry surface-induced reaction force.
This calculation finds zero axial resultant separately. No fatigue rate or
automatic initiation classification follows from a converged local state.
"""
from __future__ import annotations
import argparse,hashlib,json,sys,time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from solver_v1.silicon_initiation_research import prescribed_grips,grip_observables,connectivity_diagnostics
from results.silicon_wafer_feasibility.mace_force_only_v9 import force_only
from results.silicon_wafer_feasibility.run_mace_boron_v9 import MODEL_SHA256


def dump(path,obj):path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def main(args):
    import torch
    from ase import Atoms
    from ase.io import write
    from mace.calculators import MACECalculator
    if hashlib.sha256(args.model.read_bytes()).hexdigest()!=MODEL_SHA256:raise ValueError('model changed')
    if args.output.exists() and any(args.output.iterdir()):raise ValueError('fresh output required')
    args.output.mkdir(parents=True,exist_ok=True)
    with np.load(args.geometry) as d:
        reference=d['reference'].copy();lower=d['lower'].copy();upper=d['upper'].copy();free=d['free'].copy()
        area=float(d['area_A2']);gauge=float(d['gauge_length_A'])
    fixed=lower|upper
    with np.load(args.state) as d:positions=d['positions'].copy();numbers=d['numbers'].copy();cell=d['cell'].copy()
    if positions.shape!=reference.shape or not np.all(numbers==14):raise ValueError('same pure-Si prism required')
    extension=float((positions[upper,2]-reference[upper,2]).mean()-(positions[lower,2]-reference[lower,2]).mean())
    exact=prescribed_grips(reference,lower=lower,upper=upper,extension_A=extension)
    if np.max(abs(positions[fixed]-exact[fixed]))>1e-8:raise ValueError('source has different rigid grip coordinates')
    atoms=Atoms(numbers=numbers,positions=positions,cell=cell,pbc=False)
    initial_graph=connectivity_diagnostics(positions,lower=lower,upper=upper)
    if not initial_graph[-1]['grip_connected']:raise ValueError('source reference lacks a connected grip-to-grip body')
    protocol=dict(model_sha256=MODEL_SHA256,input_state_sha256=hashlib.sha256(args.state.read_bytes()).hexdigest(),
        geometry_sha256=hashlib.sha256(args.geometry.read_bytes()).hexdigest(),area_A2=area,bulk_reference_gauge_A=gauge,
        initial_extension_A=extension,nominal_stresses_GPa=args.stresses,fmax_eV_A=args.fmax,
        energy='U - external_force * extension; grip centres move by -extension/2 and +extension/2',
        external_force_units='eV/Angstrom; nominal stress GPa = force/area_A2 * 160.2176634',
        extension_search_bounds_A=[-.25*gauge,.5*gauge],
        bounds_status='numerical geometry guard; active bound is failure, never a strength criterion',
        starting_graph_diagnostics=initial_graph,initial_crack_seed_added=False,
        scope='0 K local force-controlled research states; bare nanospecimen',physical_clock=None,
        crack_initiation_certified=False,production_enabled=False)
    dump(args.output/'protocol.json',protocol)
    torch.set_num_threads(2);torch.set_num_interop_threads(1)
    model=MACECalculator(model_paths=str(args.model),device='cpu',default_dtype='float64')
    started=time.perf_counter();calls=0;completed=[]
    x=np.r_[positions[free].ravel(),extension]
    for index,sigma in enumerate(args.stresses):
        folder=args.output/f'state_{index:03d}';folder.mkdir()
        target=float(sigma*area/160.2176634);case_calls=calls;cached={}
        def unpack(values):
            r=prescribed_grips(reference,lower=lower,upper=upper,extension_A=values[-1])
            r[free]=values[:-1].reshape(-1,3);return r
        def objective(values):
            nonlocal calls,cached
            atoms.positions[:]=unpack(values)
            energy,force=force_only(model,atoms);calls+=1
            conjugate=float((force[lower,2].sum()-force[upper,2].sum())/2)
            gradient=np.r_[-force[free].ravel(),conjugate-target]
            cached=dict(x=values.copy(),positions=atoms.positions.copy(),energy=energy,forces=force,
                enthalpy=energy-target*values[-1],gradient=gradient)
            if calls%10==0:
                np.savez_compressed(folder/'checkpoint.npz',**cached,numbers=numbers,cell=cell)
                dump(args.output/'running.json',dict(state=index,target_GPa=sigma,calls=calls,
                    extension_A=float(values[-1]),force_residual_eV_A=conjugate-target,
                    free_force_max_eV_A=float(np.linalg.norm(force[free],axis=1).max()),elapsed_s=time.perf_counter()-started))
            if time.perf_counter()-started>args.max_seconds:raise TimeoutError('explicit wall-time budget')
            return cached['enthalpy'],gradient
        interruption=None
        try:
            result=minimize(objective,x,method='L-BFGS-B',jac=True,
                bounds=[(None,None)]*(len(x)-1)+[(-.25*gauge,.5*gauge)],
                options=dict(maxiter=args.max_steps,maxfun=4*args.max_steps,maxls=40,ftol=0.,gtol=args.fmax/3,maxcor=25))
            x=result.x.copy();objective(x)
            optimizer=dict(success=bool(result.success),message=str(result.message),iterations=int(result.nit),evaluations=int(result.nfev))
        except TimeoutError as exc:
            interruption=str(exc);x=cached['x'].copy();optimizer=dict(success=False,message=interruption)
        atoms.positions[:]=unpack(x)
        observables=grip_observables(atoms.positions,cached['forces'],lower=lower,upper=upper,area_A2=area)
        residual=observables['conjugate_force_eV_A']-target
        active_bound=bool(x[-1]<=-.25*gauge+1e-8 or x[-1]>=.5*gauge-1e-8)
        converged=bool(observables['free_force_max_eV_A']<=args.fmax and abs(residual)<=args.fmax and not active_bound and interruption is None)
        row=dict(state=index,target_nominal_GPa=float(sigma),actual_nominal_GPa=observables['nominal_stress_GPa'],
            target_force_eV_A=target,axial_force_residual_eV_A=residual,extension_A=float(x[-1]),
            grip_centre_distance_A=gauge+float(x[-1]),energy_eV=cached['energy'],enthalpy_eV=float(cached['enthalpy']),
            converged=converged,active_guard_bound=active_bound,optimizer=optimizer,
            force_calls=calls-case_calls,elapsed_s=time.perf_counter()-started,**observables)
        np.savez_compressed(folder/'raw.npz',**cached,numbers=numbers,cell=cell)
        write(folder/'state.extxyz',atoms)
        dump(folder/'result.json',dict(**row,connectivity=connectivity_diagnostics(atoms.positions,lower=lower,upper=upper),
                                      crack_initiation_certified=False))
        completed.append(row);dump(args.output/'progress.json',dict(completed=completed,planned=len(args.stresses),calls=calls))
        print('FORCE_STATE',index,'target_GPa',sigma,'actual_GPa',observables['nominal_stress_GPa'],
            'extension',x[-1],'fmax',observables['free_force_max_eV_A'],'converged',converged,flush=True)
        if not converged:raise RuntimeError('force-controlled state not converged; prefix preserved')
    dump(args.output/'summary.json',dict(complete=True,states=completed,force_only_calls=calls,
        elapsed_s=time.perf_counter()-started,actual_new_DFT=0,actual_new_MD=0,physical_clock=None,initiation_probability=None))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--state',type=Path,required=True);p.add_argument('--geometry',type=Path,required=True)
    p.add_argument('--model',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--stresses',type=float,nargs='+',default=[0.,.1]);p.add_argument('--fmax',type=float,default=5e-4)
    p.add_argument('--max-steps',type=int,default=600);p.add_argument('--max-seconds',type=float,default=5400);main(p.parse_args())
