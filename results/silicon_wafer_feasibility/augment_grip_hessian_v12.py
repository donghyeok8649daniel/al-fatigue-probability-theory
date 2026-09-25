"""Add one rigid-grip coordinate to three completed fixed-grip Hessians.

Reuse the 648 x 648 matrices. Differentiate only the new grip-gradient row and
check the mixed column by independent force differences. This diagnoses local
force-ensemble stability at each state's own reaction, not a new relaxation.
"""
from __future__ import annotations
import argparse,hashlib,json,sys,time
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from results.silicon_wafer_feasibility.check_prism_dense_hessian_v12 import constraint_basis
from results.silicon_wafer_feasibility.mace_force_only_v9 import force_only
from results.silicon_wafer_feasibility.run_mace_boron_v9 import MODEL_SHA256


def dump(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def main(args):
    import torch
    from ase import Atoms
    from mace.calculators import MACECalculator
    if args.output.exists():raise ValueError('fresh output required')
    if hashlib.sha256(args.model.read_bytes()).hexdigest()!=MODEL_SHA256:raise ValueError('model changed')
    args.output.mkdir(parents=True)
    with np.load(args.geometry) as d:
        free=d['free'];lower=d['lower'];upper=d['upper'];reference=d['reference']
        area=float(d['area_A2']);gauge=float(d['gauge_length_A'])
    basis,norm=constraint_basis(free,lower,upper,'force')
    torch.set_num_threads(2);torch.set_num_interop_threads(1)
    calc=MACECalculator(model_paths=str(args.model),device='cpu',default_dtype='float64')
    started=time.perf_counter();calls=0;states=[]
    dump(args.output/'protocol.json',dict(model_sha256=MODEL_SHA256,runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        geometry_sha256=hashlib.sha256(args.geometry.read_bytes()).hexdigest(),
        method='one autograd rigid-grip row; reuse completed fixed-grip block; independent mixed-column force differences',
        force_ensemble='constant axial force matched to each stored state reaction; existing free-force residual retained',
        normalized_grip_coordinate_norm=norm,physical_clock=None,new_DFT=0,new_MD=0,first_crack_certified=False))
    for name in ('loading8','return8','loading10'):
        if time.perf_counter()-started>args.max_seconds:raise TimeoutError('augmentation time budget; prefix retained')
        source=args.results/('dense_'+name);outdir=args.output/name;outdir.mkdir()
        summary=json.loads((source/'summary.json').read_text(encoding='utf-8'))
        if not summary['complete']:raise ValueError('source fixed-grip Hessian incomplete')
        with np.load(source/'raw_hessian.npz') as d:
            fixed=d['hessian'];r0=d['positions'];numbers=d['numbers'];cell=d['cell']
            saved_energy=float(d['energy']);saved_forces=d['forces']
            if not np.array_equal(d['basis'],basis[:,:-1]):raise ValueError('fixed basis changed')
        atoms=Atoms(numbers=numbers,positions=r0,cell=cell,pbc=False)
        batch=calc._clone_batch(calc._atoms_to_batch(atoms)).to_dict()
        out=calc.models[0](batch,training=True,compute_force=True,compute_virials=False,compute_stress=False);calls+=1
        force=out['forces'];b=torch.as_tensor(basis,dtype=force.dtype,device=force.device)
        reduced=-b.T@force.reshape(-1);g=reduced.detach().cpu().numpy().copy()
        energy=float(out['energy'].detach().cpu().reshape(-1)[0]);base_force=force.detach().cpu().numpy().copy()
        if abs(energy-saved_energy)>1e-8 or np.max(abs(base_force-saved_forces))>1e-8:
            raise ValueError('stored source no longer replays')
        row=torch.autograd.grad(reduced[-1],batch['positions'],retain_graph=False,create_graph=False)[0]
        grip_row=row.detach().cpu().numpy().ravel()@basis
        del row,reduced,force,b,batch,out
        h=np.empty((len(basis.T),len(basis.T)))
        h[:-1,:-1]=(fixed+fixed.T)/2;h[-1,:]=grip_row;h[:-1,-1]=grip_row[:-1]
        checks=[]
        displacement=basis[:,-1].reshape(r0.shape)
        for step in (2e-4,1e-4):
            atoms.positions[:]=r0+step*displacement;ep,fp=force_only(calc,atoms)
            atoms.positions[:]=r0-step*displacement;em,fm=force_only(calc,atoms);calls+=2
            fd=-basis.T@(fp-fm).ravel()/(2*step)
            error=float(np.linalg.norm(fd-grip_row))
            checks.append(dict(step_A=step,mixed_column_error_norm_eV_A2=error,
                grip_energy_even_difference_eV_A2=float((ep+em-2*energy)/step**2)))
            if error>1e-5:raise ValueError('independent mixed-column force derivative failed')
        # Energy second differences use larger steps to control cancellation.
        energy_checks=[]
        for step in (2e-3,1e-3):
            atoms.positions[:]=r0+step*displacement;ep,_=force_only(calc,atoms)
            atoms.positions[:]=r0-step*displacement;em,_=force_only(calc,atoms);calls+=2
            energy_checks.append(dict(step_A=step,curvature_eV_A2=float((ep+em-2*energy)/step**2),
                error_eV_A2=float((ep+em-2*energy)/step**2-grip_row[-1])))
        values,vectors=np.linalg.eigh(h)
        fixed_values=np.linalg.eigvalsh(h[:-1,:-1])
        reaction=float(g[-1]*norm);g[-1]=0.
        extension=float((r0[upper,2]-reference[upper,2]).mean()-(r0[lower,2]-reference[lower,2]).mean())
        record=dict(state=name,complete=True,dimension=len(h),source_raw_hessian_sha256=hashlib.sha256((source/'raw_hessian.npz').read_bytes()).hexdigest(),
            source_energy_replay_error_eV=energy-saved_energy,source_force_replay_error_eV_A=float(np.max(abs(base_force-saved_forces))),
            matched_dead_force_eV_A=reaction,matched_nominal_stress_GPa=reaction/area*160.2176634,
            free_gradient_max_eV_A=float(np.max(abs(g[:-1]))),fixed_grip_minimum_curvature_eV_A2=float(fixed_values[0]),
            force_ensemble_minimum_curvature_eV_A2=float(values[0]),force_ensemble_negative_modes=int(np.sum(values<0)),
            independent_mixed_column_checks=checks,grip_energy_checks=energy_checks,
            eigen_residual_norm_eV_A2=float(np.linalg.norm(h@vectors-vectors*values)),
            new_relaxation_performed=False,material_approved=False,physical_clock=None)
        if fixed_values[0]>0:
            response=np.linalg.solve(h[:-1,:-1],h[:-1,-1])
            schur=float(h[-1,-1]-h[-1,:-1]@response);stiffness=schur*norm**2
            record.update(schur_curvature_eV_A2=schur,relaxed_axial_stiffness_eV_A2=stiffness,
                finite_prism_tangent_GPa=stiffness/area*160.2176634*(gauge+extension),
                tangent_convention='current grip distance and original area; not bulk Young modulus',
                free_response_residual_norm=float(np.linalg.norm(h[:-1,:-1]@response-h[:-1,-1])))
        np.savez_compressed(outdir/'raw_augmented_hessian.npz',hessian=h,basis=basis,gradient=g,positions=r0,
            eigenvalues=values,eigenvectors=vectors,grip_autograd_row=grip_row)
        dump(outdir/'summary.json',record);states.append(record)
        dump(args.output/'running.json',dict(completed=len(states),planned=3,elapsed_s=time.perf_counter()-started,new_graph_evaluations=calls))
    report=dict(complete=True,states=states,new_graph_evaluations=calls,new_grip_autograd_rows=3,
        reused_fixed_grip_rows=3*(basis.shape[1]-1),elapsed_s=time.perf_counter()-started,
        new_DFT=0,new_MD=0,physical_clock=None,first_crack_certified=False)
    dump(args.output/'summary.json',report);print(json.dumps(report,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('results','geometry','model','output'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--max-seconds',type=float,default=600.)
    main(p.parse_args())
