"""Low curvature modes for explicitly constrained intact-prism coordinates.

Fixed-displacement and axial-force control have different stability spaces.
No atomic mass or kinetic clock is inferred from static curvature. A truncated
Lanczos result is a numerical local diagnostic, not a global-minimum proof.
"""
from __future__ import annotations
import argparse,hashlib,json,sys,time
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import LinearOperator,eigsh,ArpackNoConvergence
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from results.silicon_wafer_feasibility.mace_force_only_v9 import force_only
from results.silicon_wafer_feasibility.run_mace_boron_v9 import MODEL_SHA256


def dump(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def main(args):
    import torch
    from ase import Atoms
    from mace.calculators import MACECalculator
    if hashlib.sha256(args.model.read_bytes()).hexdigest()!=MODEL_SHA256:raise ValueError('model changed')
    if args.output.exists() and any(args.output.iterdir()):raise ValueError('fresh output required')
    args.output.mkdir(parents=True,exist_ok=True)
    with np.load(args.state) as d:r0=d['positions'].copy();numbers=d['numbers'].copy();cell=d['cell'].copy()
    with np.load(args.geometry) as d:free=d['free'].copy();lower=d['lower'].copy();upper=d['upper'].copy()
    if len(r0)!=len(free) or np.any(lower&upper) or not np.array_equal(free,~(lower|upper)):
        raise ValueError('inconsistent constraint masks')
    atoms=Atoms(numbers=numbers,positions=r0,cell=cell,pbc=False)
    count=3*int(free.sum())+int(args.ensemble=='force');calls=0;started=time.perf_counter()
    torch.set_num_threads(2);torch.set_num_interop_threads(1)
    model=MACECalculator(model_paths=str(args.model),device='cpu',default_dtype='float64')
    protocol=dict(state_sha256=hashlib.sha256(args.state.read_bytes()).hexdigest(),
        geometry_sha256=hashlib.sha256(args.geometry.read_bytes()).hexdigest(),model_sha256=MODEL_SHA256,
        ensemble=args.ensemble,coordinates=count,free_atoms=int(free.sum()),step_A=args.step,
        requested_modes=args.modes,krylov_dimension=args.ncv,krylov_max_iterations=args.maxiter,
        force_control='free coordinates plus rigid grip extension; external dead load has zero second derivative',
        scope='local potential-energy curvature diagnostic; no finite-T or kinetic calibration',
        first_crack_event_certified=False,physical_clock=None)
    dump(args.output/'protocol.json',protocol)
    def displace(v):
        result=np.zeros_like(r0);result[free]=v[:3*int(free.sum())].reshape(-1,3)
        if args.ensemble=='force':result[lower,2]=-.5*v[-1];result[upper,2]=.5*v[-1]
        return result
    def reduced_gradient(force):
        v=-force[free].ravel()
        if args.ensemble=='force':v=np.r_[v,.5*(force[lower,2].sum()-force[upper,2].sum())]
        return v
    def evaluate(shift):
        nonlocal calls
        if time.perf_counter()-started>args.max_seconds or calls>=args.max_calls:
            raise TimeoutError('explicit curvature work budget')
        atoms.positions[:]=r0+shift;energy,force=force_only(model,atoms);calls+=1
        if calls%20==0:dump(args.output/'running.json',dict(calls=calls,elapsed_s=time.perf_counter()-started))
        return energy,reduced_gradient(force)
    e0,g0=evaluate(np.zeros_like(r0))
    def product(v,h=args.step):
        norm=np.linalg.norm(v)
        if norm==0:return np.zeros_like(v)
        shift=displace(v/norm)*h
        _,gp=evaluate(shift);_,gm=evaluate(-shift)
        return (gp-gm)*(norm/(2*h))
    rng=np.random.default_rng(4411);seed=rng.normal(size=count);seed/=np.linalg.norm(seed)
    operator=LinearOperator((count,count),matvec=product,dtype=float)
    converged=False;error=None;values=np.array([]);vectors=np.empty((count,0))
    try:
        values,vectors=eigsh(operator,k=args.modes,which='SA',v0=seed,ncv=args.ncv,
            tol=1e-5,maxiter=args.maxiter)
        converged=True
    except ArpackNoConvergence as exc:
        values=exc.eigenvalues;vectors=exc.eigenvectors;error='ARPACK no convergence; converged subset retained'
    except TimeoutError as exc:error=str(exc)
    records=[]
    # Preserve the solved modes before any independent refinement exceeds budget.
    np.savez_compressed(args.output/'modes.npz',values=values,vectors=vectors,gradient=g0,positions=r0,numbers=numbers,cell=cell)
    try:
        for i,value in enumerate(values):
            v=vectors[:,i];hv=product(v);half=product(v,args.step/2)
            row=dict(mode=i,ritz_curvature_eV_A2=float(value),
                residual_norm_eV_A2=float(np.linalg.norm(hv-value*v)),
                half_step_rayleigh_eV_A2=float(v@half),half_step_residual_eV_A2=float(np.linalg.norm(half-value*v)),
                finite_difference_change_norm_eV_A2=float(np.linalg.norm(half-hv)),
                projected_gradient_eV_A=float(g0@v),energy_checks=[])
            for amplitude in (1e-3,5e-4):
                ep,_=evaluate(displace(v)*amplitude);em,_=evaluate(-displace(v)*amplitude)
                row['energy_checks'].append(dict(amplitude_A=amplitude,even_curvature_eV_A2=float((ep+em-2*e0)/amplitude**2),
                    odd_slope_eV_A=float((ep-em)/(2*amplitude))))
            records.append(row);dump(args.output/'mode_checks.json',records)
    except TimeoutError as exc:error=str(exc)
    result=dict(lanczos_converged=converged,interruption=error,modes=records,
        energy_eV=float(e0),reduced_gradient_max_eV_A=float(abs(g0).max()),force_calls=calls,
        elapsed_s=time.perf_counter()-started,actual_new_MD=0,actual_new_DFT=0,
        positive_ritz_values_are_not_rigorous_global_stability_proof=True,
        global_minimum_certified=False,crack_initiation_certified=False,physical_clock=None)
    dump(args.output/'summary.json',result);print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--state',type=Path,required=True);p.add_argument('--geometry',type=Path,required=True)
    p.add_argument('--model',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--ensemble',choices=['displacement','force'],default='displacement')
    p.add_argument('--step',type=float,default=2e-4);p.add_argument('--modes',type=int,default=2)
    p.add_argument('--ncv',type=int,default=16);p.add_argument('--maxiter',type=int,default=24)
    p.add_argument('--max-calls',type=int,default=700);p.add_argument('--max-seconds',type=float,default=4500);main(p.parse_args())
