"""Static rigid-interface branch audit, never a first-crack initiation result.

The u=0 symmetry path is followed under normal traction. Full three-coordinate
curvature is checked; loss of this rigid local minimum is an ideal mechanical
diagnostic of an unapproved potential. No atomistic relaxation, thermal PMF,
mobility, probability, loading cycles or physical time is computed.
"""
import argparse,csv,hashlib,json,time
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from results.silicon_wafer_feasibility.run_static_probe import source_parameters
from .silicon_sw_anharmonic import AngularSW111BesselInterface
from .run_silicon_sw_material_fit import bulk_diagnostic,dump

EVA3_TO_GPA=1.602176634e-19/1e-30/1e9
EVA2_TO_JM2=1.602176634e-19/1e-20


def encode(jet):
    return dict(energy_eV=float(jet.value),gradient_eV_A=jet.gradient.tolist(),
                hessian_eV_A2=jet.hessian.tolist(),eigenvalues_eV_A2=np.linalg.eigvalsh(jet.hessian).tolist())


def main(args):
    if args.output.exists():raise ValueError('fresh output directory required')
    args.output.mkdir(parents=True);start=time.perf_counter()
    selected=json.loads(args.selected.read_text());reference=json.loads(args.reference.read_text())
    root=Path(__file__).resolve().parents[1];bindings=[]
    for relative in ['solver_v1/run_silicon_sw_interface_audit.py','solver_v1/silicon_sw_anharmonic.py',
        'solver_v1/silicon_bessel_reference.py','solver_v1/silicon_environment_research.py',
        'solver_v1/run_silicon_sw_material_fit.py']:
        raw=(root/relative).read_bytes();target=args.output/'source_snapshots'/relative
        target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(raw)
        bindings.append(dict(path=relative,sha256=hashlib.sha256(raw).hexdigest()))
    for path in [args.selected,args.reference]:
        raw=path.read_bytes();target=args.output/'input_snapshots'/path.name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(raw)
    dump(args.output/'source_manifest.json',bindings)
    p,_=source_parameters();beta=selected.get('angle_beta',selected.get('selected',{}).get('beta'))
    if beta is None:raise ValueError('selected angular beta required')
    label='selected_joint_candidate' if 'angle_beta' in selected else 'selected_angular_shape'
    configs=[(label,selected['parameters'],beta,reference['diamond_a0'])]
    if not args.candidate_only:
        original_bulk=bulk_diagnostic(p);configs.insert(0,('original_SW',p,0.,original_bulk['lattice_A']))
    roots=[];rows=[];calls=0;fd_checks=[];bessel=[]
    for label,parameters,beta,lattice in configs:
        for cut in ['shuffle','glide']:
            model=AngularSW111BesselInterface(parameters,lattice,angle_beta=beta,cut_kind=cut)
            area=float(model.geometry.atomic_cell_area);rc=parameters['sigma']*parameters['a']
            def calc(x):
                nonlocal calls
                calls+=1;return model.evaluate([x,0.,0.],method='direct')
            zero=calc(0.)
            if np.min(np.linalg.eigvalsh(zero.hessian))<=0:raise RuntimeError('initial rigid interface not stable')
            coarse=[];fine=[]
            for x in np.linspace(0.,rc,161):
                j=calc(float(x));eig=np.linalg.eigvalsh(j.hessian);fine.append((float(x),float(eig[0])))
                rows.append(dict(case=label,cut_kind=cut,opening_A=float(x),
                    energy_eV=float(j.value),work_J_m2=float(j.value/area*EVA2_TO_JM2),
                    normal_traction_GPa=float(j.gradient[0]/area*EVA3_TO_GPA),
                    lateral_force_maximum_eV_A=float(np.max(abs(j.gradient[1:]))),
                    minimum_curvature_eV_A2=float(eig[0]),normal_curvature_eV_A2=float(j.hessian[0,0])))
            for mesh,data in [('81',fine[::2]),('161',fine)]:
                brackets=[(left[0],right[0]) for left,right in zip(data,data[1:]) if left[1]>0 and right[1]<0]
                if not brackets:raise RuntimeError('no sampled first stability loss bracket')
                lo,hi=brackets[0];x=brentq(lambda a:np.linalg.eigvalsh(calc(float(a)).hessian)[0],lo,hi,xtol=2e-12)
                j=calc(x);eig,vec=np.linalg.eigh(j.hessian)
                roots.append(dict(case=label,cut_kind=cut,mesh_points=int(mesh),opening_A=float(x),
                    gap_A=float(model.gap),atomic_cell_area_A2=area,lattice_A=float(lattice),
                    normal_traction_GPa=float(j.gradient[0]/area*EVA3_TO_GPA),
                    critical_eigenvector=vec[:,0].tolist(),jet=encode(j),
                    interpretation='first full-curvature loss on sampled u=0 rigid symmetry branch; not crack initiation strength'))
            root=roots[-1]
            for x in [0.,.5*root['opening_A'],root['opening_A']]:
                q=np.array([x,0.,0.]);j=model.evaluate(q,method='direct');steps=[]
                for step in [2e-5,1e-5]:
                    h=np.column_stack([(model.evaluate(q+step*np.eye(3)[k],method='direct').gradient-
                        model.evaluate(q-step*np.eye(3)[k],method='direct').gradient)/(2*step) for k in range(3)])
                    calls+=6;steps.append(dict(step_A=step,maximum_hessian_error_eV_A2=float(np.max(abs(h-j.hessian)))))
                fd_checks.append(dict(case=label,cut_kind=cut,state_A=q.tolist(),checks=steps))
            separated=calc(rc+1.)
            root.update(zero=encode(zero),separated=encode(separated),
                rigid_separation_work_J_m2=float(separated.value/area*EVA2_TO_JM2))
            # One fresh Bessel jet at the branch loss, on each cut and model.
            bj=model.evaluate([root['opening_A'],0.,0.],method='bessel');dj=calc(root['opening_A'])
            bessel.append(dict(case=label,cut_kind=cut,bessel=encode(bj),direct=encode(dj),
                energy_error_eV=float(abs(bj.value-dj.value)),
                gradient_error_eV_A=float(np.max(abs(bj.gradient-dj.gradient))),
                hessian_error_eV_A2=float(np.max(abs(bj.hessian-dj.hessian)))))
            dump(args.output/'running.json',dict(completed_cases=len(bessel),planned_cases=len(configs)*2,elapsed_seconds=time.perf_counter()-start))
            print(json.dumps(dict(case=label,cut_kind=cut,traction_GPa=root['normal_traction_GPa'])),flush=True)
    with (args.output/'normal_branch.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    dump(args.output/'branch_roots.json',roots);dump(args.output/'finite_difference_checks.json',fd_checks)
    dump(args.output/'bessel_audit.json',bessel)
    if any(r['energy_error_eV']>2e-8 or r['gradient_error_eV_A']>2e-7 or r['hessian_error_eV_A2']>3e-6 for r in bessel):
        raise RuntimeError('fresh branch Bessel comparison failed')
    dump(args.output/'summary.json',dict(complete=True,cases=len(configs)*2,sampled_path_points=len(rows),
        direct_jet_calls=calls,bessel_jet_calls=len(bessel),elapsed_seconds=time.perf_counter()-start,
        definition='W=per-interface-cell excess static energy of two rigid half crystals, q=(opening,u1,u2)',
        unit_conversions=dict(eV_A3_to_GPa=EVA3_TO_GPA,eV_A2_to_J_m2=EVA2_TO_JM2),
        loading='only normal conjugate traction, W-A_atomic*T*opening; lateral force along u=0 checked',
        first_branch_loss='root of minimum eigenvalue of the full 3x3 rigid Hessian, using first sampled sign-change bracket',
        limitations=['no global minimum certification','no atomistic surface relaxation','no thermal free energy',
                     'no pre-existing crack is introduced, but a rigid plane is prescribed','no finite-size specimen mechanics'],
        new_DFT=0,new_MD=0,material_approved=False,first_initiation_validated=False,physical_clock_validated=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['selected','reference','output']:parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--candidate-only',action='store_true')
    main(parser.parse_args())
