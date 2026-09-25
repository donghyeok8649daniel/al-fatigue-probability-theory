"""Full constrained static Hessian, using the same pinned MACE energy.

The grip coordinate is orthonormal in Cartesian displacement, not mass weighted.
An external dead force contributes a linear potential and zero Hessian. No
transition rate, global stability, first-crack event or physical clock follows.
"""
from __future__ import annotations
import argparse, hashlib, json, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from results.silicon_wafer_feasibility.mace_force_only_v9 import force_only
from results.silicon_wafer_feasibility.run_mace_boron_v9 import MODEL_SHA256


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def constraint_basis(free, lower, upper, ensemble):
    if np.any(lower & upper) or not np.array_equal(free, ~(lower | upper)):
        raise ValueError('constraint masks inconsistent')
    indices=np.flatnonzero(np.repeat(free, 3)); n=3*len(free)
    basis=np.zeros((n, len(indices)+(ensemble=='force')))
    basis[indices,np.arange(len(indices))]=1.
    grip_norm=float(np.sqrt((lower.sum()+upper.sum())/4))
    if ensemble=='force':
        if grip_norm==0: raise ValueError('nonempty grips required')
        grip=np.zeros((len(free),3));grip[lower,2]=-.5;grip[upper,2]=.5
        basis[:,-1]=grip.ravel()/grip_norm
    if not np.allclose(basis.T@basis,np.eye(basis.shape[1]),rtol=0,atol=1e-14):
        raise ValueError('basis is not orthonormal')
    return basis,grip_norm


def main(args):
    import torch
    from ase import Atoms
    from mace.calculators import MACECalculator
    if hashlib.sha256(args.model.read_bytes()).hexdigest()!=MODEL_SHA256:
        raise ValueError('model weights changed')
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError('fresh output folder required')
    args.output.mkdir(parents=True,exist_ok=True)
    with np.load(args.state) as d:
        r0=d['positions'].copy();numbers=d['numbers'].copy();cell=d['cell'].copy()
        saved_energy=float(d['energy']); saved_forces=d['forces'].copy()
    with np.load(args.geometry) as d:
        free=d['free'].copy();lower=d['lower'].copy();upper=d['upper'].copy()
        area=float(d['area_A2']);reference=d['reference'].copy();gauge=float(d['gauge_length_A'])
    basis,grip_norm=constraint_basis(free,lower,upper,args.ensemble)
    dimension=basis.shape[1];target=args.stress*area/160.2176634
    started=time.perf_counter();calls=0;vjp_rows=0
    torch.set_num_threads(2);torch.set_num_interop_threads(1)
    calc=MACECalculator(model_paths=str(args.model),device='cpu',default_dtype='float64')
    atoms=Atoms(numbers=numbers,positions=r0,cell=cell,pbc=False)
    protocol=dict(state_sha256=hashlib.sha256(args.state.read_bytes()).hexdigest(),
        geometry_sha256=hashlib.sha256(args.geometry.read_bytes()).hexdigest(),
        model_sha256=MODEL_SHA256,runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        ensemble=args.ensemble,dimension=dimension,
        basis='orthonormal Cartesian free displacements plus normalized rigid grip extension',
        grip_extension_basis_norm=grip_norm,external_nominal_stress_GPa=args.stress,
        energy='U - F*Delta; linear dead load has zero Hessian',method='autograd reduced force Jacobian, one row at a time',
        symmetrization='only after measuring and checking raw asymmetry',
        max_seconds=args.max_seconds,max_rows=args.max_rows,
        units='H: eV/Angstrom^2, coordinates: Angstrom in orthonormal Cartesian basis',
        no_atomic_mass_or_kinetic_time=True,new_DFT=0,new_MD=0,
        first_crack_certified=False,material_approved=False,physical_clock=None)
    dump(args.output/'protocol.json',protocol)
    batch=calc._clone_batch(calc._atoms_to_batch(atoms)).to_dict()
    out=calc.models[0](batch,training=True,compute_force=True,compute_virials=False,compute_stress=False)
    calls+=1
    energy=float(out['energy'].detach().cpu().reshape(-1)[0])
    force=out['forces']; positions=batch['positions']
    base_force=force.detach().cpu().numpy().copy()
    b=torch.as_tensor(basis,dtype=positions.dtype,device=positions.device)
    reduced=-b.T@force.reshape(-1)
    g0=reduced.detach().cpu().numpy().copy()
    if args.ensemble=='force': g0[-1]-=target/grip_norm
    replay=dict(saved_energy_difference_eV=energy-saved_energy,
        saved_force_max_difference_eV_A=float(np.max(abs(base_force-saved_forces))),
        reduced_enthalpy_gradient_max_eV_A=float(np.max(abs(g0))))
    dump(args.output/'source_replay.json',replay)
    if abs(energy-saved_energy)>1e-8 or replay['saved_force_max_difference_eV_A']>1e-8:
        raise ValueError('autograd graph energy/forces disagree with saved state')
    hessian=np.zeros((dimension,dimension));complete=False;interruption=None
    try:
        for row in range(dimension):
            if time.perf_counter()-started>args.max_seconds or row>=args.max_rows:
                raise TimeoutError('explicit direct-Hessian budget; partial rows retained')
            deriv=torch.autograd.grad(reduced[row],positions,retain_graph=True,create_graph=False)[0]
            hessian[row]=deriv.detach().cpu().numpy().ravel()@basis;vjp_rows=row+1
            if row%16==0 or row+1==dimension:
                np.savez_compressed(args.output/'checkpoint.npz',hessian_rows=hessian[:row+1],basis=basis,
                    positions=r0,numbers=numbers,gradient=g0)
                dump(args.output/'running.json',dict(completed_rows=row+1,planned_rows=dimension,
                    elapsed_s=time.perf_counter()-started))
                print('HESSIAN_ROW',row+1,dimension,'seconds',time.perf_counter()-started,flush=True)
        complete=True
    except TimeoutError as exc: interruption=str(exc)
    del reduced,force,out,batch,positions,b
    result=dict(complete=complete,interruption=interruption,completed_rows=vjp_rows,dimension=dimension,
        new_graph_evaluations=calls,autograd_vjp_rows=vjp_rows,elapsed_s=time.perf_counter()-started,
        new_DFT=0,new_MD=0,first_crack_certified=False,physical_clock=None)
    if not complete:
        np.savez_compressed(args.output/'partial.npz',hessian_rows=hessian[:vjp_rows],basis=basis,gradient=g0)
        dump(args.output/'summary.json',result);return
    raw_asym=float(np.max(abs(hessian-hessian.T)));norm=float(np.linalg.norm(hessian,ord='fro'))
    result.update(raw_max_asymmetry_eV_A2=raw_asym,raw_frobenius_norm_eV_A2=norm)
    np.savez_compressed(args.output/'raw_hessian.npz',hessian=hessian,basis=basis,positions=r0,
        numbers=numbers,cell=cell,forces=base_force,energy=energy,gradient=g0)
    if raw_asym>1e-9*max(norm,np.finfo(float).tiny):
        dump(args.output/'summary.json',dict(**result,validation_failed='raw Hessian asymmetry'));raise ValueError('Hessian asymmetry')
    symmetric=(hessian+hessian.T)/2
    eigenvalues,eigenvectors=np.linalg.eigh(symmetric)
    np.savez_compressed(args.output/'spectrum.npz',eigenvalues=eigenvalues,eigenvectors=eigenvectors)
    result.update(minimum_curvature_eV_A2=float(eigenvalues[0]),maximum_curvature_eV_A2=float(eigenvalues[-1]),
        negative_eigenvalues=int(np.sum(eigenvalues<0)),lowest_eight_eV_A2=eigenvalues[:8].tolist(),
        eigensystem_residual_frobenius=float(np.linalg.norm(symmetric@eigenvectors-eigenvectors*eigenvalues)),
        orthogonality_max_error=float(np.max(abs(eigenvectors.T@eigenvectors-np.eye(dimension)))))
    if args.ensemble=='force':
        # Coordinates are ordered so the leading block is the fixed-grip space.
        fixed_values=np.linalg.eigvalsh(symmetric[:-1,:-1])
        result['fixed_grip_at_same_state_lowest_eight_eV_A2']=fixed_values[:8].tolist()
        result['fixed_grip_at_same_state_negative_eigenvalues']=int(np.sum(fixed_values<0))
        if fixed_values[0]>0:
            free_response=np.linalg.solve(symmetric[:-1,:-1],symmetric[:-1,-1])
            schur=float(symmetric[-1,-1]-symmetric[-1,:-1]@free_response)
            extension=float((r0[upper,2]-reference[upper,2]).mean()-(r0[lower,2]-reference[lower,2]).mean())
            stiffness=schur*grip_norm**2
            result['relaxed_grip_tangent']=dict(schur_curvature_eV_A2=schur,
                d_force_d_extension_eV_A2=stiffness,
                d_nominal_stress_d_extension_GPa_A=stiffness/area*160.2176634,
                tangent_using_current_grip_distance_GPa=stiffness/area*160.2176634*(gauge+extension),
                convention='finite bare prism and rigid grips; current grip distance; not bulk Young modulus',
                free_response_solve_residual=float(np.linalg.norm(symmetric[:-1,:-1]@free_response-symmetric[:-1,-1])))
    checks=[]
    def evaluate(delta):
        nonlocal calls
        if time.perf_counter()-started>args.max_seconds: raise TimeoutError('independent-check budget')
        atoms.positions[:]=r0+delta
        e,f=force_only(calc,atoms);calls+=1
        return e,-basis.T@f.ravel()
    try:
        vectors=[('lowest_'+str(i),eigenvectors[:,i]) for i in range(min(args.check_modes,dimension))]
        rng=np.random.default_rng(12541)
        for i in range(2):
            v=rng.normal(size=dimension);v/=np.linalg.norm(v);vectors.append(('random_'+str(i),v))
        for name,v in vectors:
            displacement=(basis@v).reshape(r0.shape);hv=symmetric@v
            row=dict(name=name,rayleigh_eV_A2=float(v@hv),force_differences=[],energy_differences=[])
            for step in (2e-4,1e-4):
                ep,gp=evaluate(step*displacement);em,gm=evaluate(-step*displacement)
                numerical=(gp-gm)/(2*step)
                row['force_differences'].append(dict(step_A=step,
                    error_norm_eV_A2=float(np.linalg.norm(numerical-hv)),
                    rayleigh_eV_A2=float(v@numerical)))
            for step in (2e-3,1e-3):
                ep,_=evaluate(step*displacement);em,_=evaluate(-step*displacement)
                row['energy_differences'].append(dict(step_A=step,curvature_eV_A2=float((ep+em-2*energy)/step**2)))
            checks.append(row);dump(args.output/'mode_checks.json',checks)
    except TimeoutError as exc:result['independent_check_interruption']=str(exc)
    result.update(independent_checked_directions=len(checks),new_graph_evaluations=calls,
        elapsed_s=time.perf_counter()-started,local_hessian_is_not_global_or_material_validation=True)
    dump(args.output/'summary.json',result);print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('state','geometry','model','output'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--ensemble',choices=['force','displacement'],default='displacement')
    p.add_argument('--stress',type=float,default=0.)
    p.add_argument('--max-seconds',type=float,default=7200.)
    p.add_argument('--max-rows',type=int,default=100000)
    p.add_argument('--check-modes',type=int,default=4)
    main(p.parse_args())
