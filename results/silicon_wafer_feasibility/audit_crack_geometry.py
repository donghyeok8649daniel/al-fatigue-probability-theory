"""Spatial audit of same-grip Si structures; never assigns a crack label.

Pair-distance cutoffs are sensitivity controls, not chemical bond definitions.
Synthetic affine geometry controls are not new elastic material calculations.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def nonaffine(reference, changed, free):
    design = np.c_[reference[free], np.ones(int(free.sum()))]
    mapping, _, rank, _ = np.linalg.lstsq(design, changed[free], rcond=None)
    if rank != 4:
        raise ValueError('a full-rank three-dimensional reference required')
    residual = changed[free]-design@mapping
    return residual, mapping


def main(args):
    if args.output.exists():
        raise ValueError('fresh output required')
    root = args.root
    paths = [root/'results/silicon_initiation_v13'/f'dense_{name}'/'raw_hessian.npz'
             for name in ('loading8','return8')]
    states = []
    for path in paths:
        with np.load(path) as data:
            states.append({key:data[key].copy() for key in data.files})
    a, b = [state['positions'] for state in states]
    basis = states[0]['basis']
    if not np.array_equal(basis, states[1]['basis']):
        raise ValueError('different free coordinate basis')
    free = np.any(basis.reshape(len(a),3,-1) != 0, axis=(1,2))
    if not np.array_equal(a[~free], b[~free]):
        raise ValueError('grip locations must be identical')
    residual, mapping = nonaffine(a,b,free)
    affine_control = a@np.diag([1.,1.,1.1])+[.2,-.3,.1]
    control_residual, _ = nonaffine(a,affine_control,free)
    # Independent geometric invariances; no potential is called.
    rotation = np.array([[0.,-1.,0.],[1.,0.,0.],[0.,0.,1.]])
    rotated, _ = nonaffine(a@rotation+[1.,2.,3.],b@rotation+[1.,2.,3.],free)
    rigid_error = float(np.max(abs(np.linalg.norm(rotated,axis=1)-np.linalg.norm(residual,axis=1))))
    affine_error = float(np.max(abs(control_residual)))
    if rigid_error > 1e-11 or affine_error > 1e-11:
        raise ValueError('affine/rigid geometric control failed')
    i,j = np.triu_indices(len(a),1)
    da = np.linalg.norm(a[i]-a[j],axis=1)
    db = np.linalg.norm(b[i]-b[j],axis=1)
    dc = np.linalg.norm(affine_control[i]-affine_control[j],axis=1)
    rows, arrays = [], dict(positions_initial=a,positions_return=b,free=free,
                           nonaffine_displacement_A=residual,affine_mapping=mapping)
    for index, cutoff in enumerate([2.7,2.95,3.2,3.45,3.7]):
        lost = (da <= cutoff)&(db > cutoff)
        gained = (da > cutoff)&(db <= cutoff)
        control_lost = (da <= cutoff)&(dc > cutoff)
        arrays[f'lost_pairs_{index}'] = np.c_[i[lost],j[lost]]
        arrays[f'gained_pairs_{index}'] = np.c_[i[gained],j[gained]]
        rows.append(dict(distance_cutoff_A=cutoff,reference_pairs=int((da <= cutoff).sum()),
            lost_pairs=int(lost.sum()),gained_pairs=int(gained.sum()),
            changed_atoms=int(len(np.unique(np.r_[i[lost|gained],j[lost|gained]]))),
            affine_control_lost_pairs=int(control_lost.sum()),
            chemical_bond_or_crack_threshold=False))
    args.output.mkdir(parents=True)
    np.savez_compressed(args.output/'raw_spatial_diagnostics.npz',**arrays)
    summary = dict(source_sha256={str(path.relative_to(root)).replace('\\','/'):
            hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},
        free_atoms=int(free.sum()),fixed_atoms=int((~free).sum()),
        same_grip_max_difference_A=0.,
        nonaffine_rms_A=float(np.sqrt(np.mean(np.sum(residual*residual,axis=1)))),
        nonaffine_maximum_A=float(np.max(np.linalg.norm(residual,axis=1))),
        rigid_control_error_A=rigid_error,affine_control_error_A=affine_error,
        distance_sensitivity=rows,
        fresh_opposed_crack_surfaces_verified=False,cohesive_loss_verified=False,
        finite_time_nonrecovery_verified=False,first_crack_label=None,
        new_potential_calls=0,new_MD=0,new_DFT=0,
        interpretation='localized rearrangement diagnostics; no physical crack classification',
        runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    magnitude = np.linalg.norm(residual,axis=1)
    fig,axes = plt.subplots(1,2,figsize=(8.5,5),layout='constrained')
    for axis,positions,title in zip(axes,[a,b],['Initial 8%','Return 8%']):
        axis.scatter(positions[~free,0],positions[~free,2],s=12,c='grey',label='fixed grips')
        points = axis.scatter(positions[free,0],positions[free,2],s=16,c=magnitude,
                              vmin=0,vmax=magnitude.max(),cmap='viridis')
        axis.set(title=title,xlabel='x (Angstrom)',ylabel='z (Angstrom)',aspect='equal')
    fig.colorbar(points,ax=axes,label='Non-affine displacement between states (Angstrom)')
    fig.suptitle('Spatial projection of rearrangement; NOT crack probability',fontsize=11)
    fig.savefig(args.output/'spatial_rearrangement.png',dpi=150);plt.close(fig)
    print(json.dumps(summary,indent=2),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args())
