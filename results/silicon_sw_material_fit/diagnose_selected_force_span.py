"""Post-selection fixed-shape force-span bounds, not a material refit.

No changes to selected parameters and no held-out-driven model selection.
Unconstrained and nonnegative force-only solutions do not retain elastic,
energy or stability gates. Bounds apply only to the locked three force bases.
"""
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
from scipy.optimize import nnls


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main(args):
    if args.output.exists():raise ValueError('fresh diagnostic output required')
    args.output.mkdir(parents=True)
    selection_path=args.cache/'final_selection.json';locked=digest(selection_path)
    selection=json.loads(selection_path.read_text());folder=args.cache/selection['selected_run']
    chosen=json.loads((folder/'selected_fit.json').read_text());alpha=np.array(chosen['amplitudes'])
    features=folder/'selected_features.npz';data=np.load(features);offset=data['offsets']
    # Frame metadata come from the original baseline, not inferred from forces.
    sys.path.insert(0,str(args.repo.resolve()))
    from results.silicon_wafer_feasibility.run_dft_material_audit import read_frames
    atoms,_=read_frames(args.archive)
    groups=[('dia','PW91'),('surface_001','PW91'),('surface_110','PW91'),
            ('surface_111','PW91'),('surface_111_3x3_das','PW91'),('surface_111_pandey','PBE')]
    rows=[]
    for kind,xc in groups:
        members=[i for i,a in enumerate(atoms) if str(a.info['config_type'])==kind and str(a.info.get('xc_functional','UNSPECIFIED'))==xc]
        if not members:continue
        indices=np.concatenate([np.arange(offset[i],offset[i+1]) for i in members])
        matrix=data['force_basis'][indices].reshape(-1,3);target=data['reference_force'][indices].ravel()
        free,_,rank,singular=np.linalg.lstsq(matrix,target,rcond=None);positive,_=nnls(matrix,target)
        def rmse(a):return float(np.linalg.norm(matrix@a-target)/np.sqrt(len(target)))
        rows.append(dict(config_type=kind,declared_xc=xc,frames=len(members),atoms=len(indices),components=len(target),
            Gram=(matrix.T@matrix).tolist(),linear=(matrix.T@target).tolist(),target_squared=float(target@target),
            locked_amplitudes=alpha.tolist(),locked_force_RMSE_eV_A=rmse(alpha),
            unconstrained_force_only_amplitudes=free.tolist(),unconstrained_force_only_RMSE_eV_A=rmse(free),
            nonnegative_force_only_amplitudes=positive.tolist(),nonnegative_force_only_RMSE_eV_A=rmse(positive),
            numerical_rank=int(rank),singular_values=singular.tolist(),
            force_only_solutions_preserve_material_gates=False,fixed_shape_only=True))
    if digest(selection_path)!=locked:raise RuntimeError('immutable selection changed during diagnosis')
    raw=Path(__file__).read_bytes();snap=args.output/'source_snapshots/diagnose_selected_force_span.py'
    snap.parent.mkdir(parents=True);snap.write_bytes(raw)
    save=lambda name,value:(args.output/name).write_bytes((json.dumps(value,indent=2,allow_nan=False)+'\n').encode())
    save('source_manifest.json',[dict(path='diagnose_selected_force_span.py',sha256=hashlib.sha256(raw).hexdigest())])
    save('group_force_bounds.json',rows)
    save('summary.json',dict(complete=True,groups=len(rows),selection_sha256=locked,feature_file_sha256=digest(features),
        archive_sha256=digest(args.archive),selected_fit_sha256=digest(folder/'selected_fit.json'),
        new_geometry_evaluations=0,new_DFT=0,new_MD=0,selected_parameters_changed=False,
        family_impossibility_proved=False,material_approved=False,
        interpretation='post-lock force-only lower bounds at a fixed basis; not admissible elastic/energy-constrained material fits'))
    print(json.dumps(rows))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['cache','archive','repo','output']:p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
