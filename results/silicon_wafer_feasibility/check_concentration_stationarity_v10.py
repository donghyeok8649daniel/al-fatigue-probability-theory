"""Matrix-free lowest Gamma curvature search for selected force-stationary cells.

Translations are explicitly penalized outside the zero-mean subspace. The
deterministic numerical starting vector is not a Monte Carlo sample.
No full-q, volume-relaxation or finite-T material stability certificate.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
from scipy.sparse.linalg import LinearOperator,eigsh,ArpackNoConvergence
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from results.silicon_wafer_feasibility.mace_force_only_v9 import force_only
from results.silicon_wafer_feasibility.run_mace_boron_v9 import MODEL_SHA256
from results.silicon_wafer_feasibility.run_concentration_cleavage_v10 import dump


def main(args):
    import torch
    from ase import Atoms
    from mace.calculators import MACECalculator
    if hashlib.sha256(args.model.read_bytes()).hexdigest()!=MODEL_SHA256:raise ValueError('model identity')
    if args.output.exists() and any(args.output.iterdir()):raise ValueError('fresh output required')
    args.output.mkdir(parents=True,exist_ok=True)
    torch.set_num_threads(2);torch.set_num_interop_threads(1)
    model=MACECalculator(model_paths=str(args.model),device='cpu',default_dtype='float64')
    raw=np.load(args.state);atoms=Atoms(numbers=raw['numbers'],positions=raw['positions'],cell=raw['cell'],pbc=True)
    e0,f0=force_only(model,atoms);calls=1;products=0;started=time.perf_counter()
    def product(v,step=1e-3,limited=True):
        nonlocal calls,products
        if limited and products>=args.max_products:raise RuntimeError('matrix-vector product budget exhausted')
        v=np.asarray(v).reshape(-1,3)
        mean=v.mean(axis=0,keepdims=True);p=v-mean;norm=np.linalg.norm(p)
        if norm==0:return np.broadcast_to(100.*mean,v.shape).ravel().copy()
        direction=p/norm
        plus,minus=atoms.copy(),atoms.copy()
        plus.positions+=step*direction;minus.positions-=step*direction
        fp=force_only(model,plus)[1];fm=force_only(model,minus)[1]
        calls+=2;products+=1
        h=-(fp-fm)/(2*step)*norm;h-=h.mean(axis=0,keepdims=True)
        if limited and products%20==0:print('LOWEST_PRODUCTS',products,flush=True)
        return (h+100.*mean).ravel()
    op=LinearOperator((3*len(atoms),)*2,matvec=product,dtype=float)
    seed=np.random.default_rng(74193).normal(size=(len(atoms),3));seed-=seed.mean(axis=0)
    seed=seed.ravel()/np.linalg.norm(seed)
    values=np.empty(0);vectors=np.empty((3*len(atoms),0));status='incomplete'
    try:
        values,vectors=eigsh(op,k=2,which='SA',v0=seed,ncv=args.krylov_dimension,tol=1e-5,maxiter=40)
        status='converged lowest Gamma Ritz pairs'
    except ArpackNoConvergence as error:
        values=error.eigenvalues;vectors=error.eigenvectors
        status='ARPACK iteration limit; partial pairs only'
    except RuntimeError as error:
        status=str(error)
    checks=[]
    for i,value in enumerate(values):
        vector=vectors[:,i]
        for step in [1e-3,3e-4,1e-4]:
            hv=product(vector,step,limited=False)
            checks.append(dict(mode=i,step_A=step,Rayleigh_eV_A2=float(vector@hv),
                eigen_residual_eV_A2=float(np.linalg.norm(hv-value*vector)),
                translation_component_L2=float(np.linalg.norm(vector.reshape(-1,3).mean(axis=0)))))
    np.savez_compressed(args.output/'modes.npz',eigenvalues=values,eigenvectors=vectors,
        positions=atoms.positions,cell=atoms.cell.array,numbers=atoms.numbers)
    result=dict(status=status,solver_converged=status.startswith('converged'),
        eigenvalues_eV_A2=values.tolist(),checks=checks,force_only_calculator_calls=calls,
        Hessian_vector_products=products,force_max_eV_A=float(np.linalg.norm(f0,axis=1).max()),
        state_sha256=hashlib.sha256(args.state.read_bytes()).hexdigest(),model_sha256=MODEL_SHA256,
        elapsed_s=time.perf_counter()-started,full_q_stability=False,material_calibrated=False,
        solver_parameters=dict(k=2,ncv=args.krylov_dimension,tolerance=1e-5,maxiter=40,
                               max_products=args.max_products,seed=74193),
        scope='selected fixed-cell Gamma curvature search, independent finite-difference step checks')
    dump(args.output/'summary.json',result);print(json.dumps(result),flush=True)
    if not result['solver_converged']:raise SystemExit(2)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--model',type=Path,required=True)
    p.add_argument('--state',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--max-products',type=int,default=240)
    p.add_argument('--krylov-dimension',type=int,default=28)
    main(p.parse_args())
