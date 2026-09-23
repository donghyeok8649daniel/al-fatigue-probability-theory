"""Geometry-only nonaffine diagnostics between converged intact-prism states.

A least-squares affine map of all free atoms removes uniform deformation and
translation. Its residual is a geometric descriptor, not plastic strain, a
crack detector, or a spatial decomposition of energy. No potential call is made.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from pathlib import Path
import numpy as np


def write_csv(path,records):
    with path.open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)


def fit_affine(before,after):
    design=np.column_stack((before,np.ones(len(before))))
    coefficients,_,rank,singular=np.linalg.lstsq(design,after,rcond=None)
    if rank!=4:raise ValueError('free-atom geometry does not span 3D affine coordinates')
    residual=after-design@coefficients
    return coefficients,residual,float(singular[0]/singular[-1])


def pairs_at(positions,cutoff):
    i,j=np.triu_indices(len(positions),1)
    distance=np.linalg.norm(positions[j]-positions[i],axis=1)
    return set(zip(i[distance<cutoff].tolist(),j[distance<cutoff].tolist()))


def main(args):
    if args.output.exists() and any(args.output.iterdir()):raise ValueError('fresh output required')
    args.output.mkdir(parents=True,exist_ok=True)
    with np.load(args.run/'geometry.npz') as data:
        reference=data['reference'];free=data['free'];fixed=~free
    ids=np.flatnonzero(free);initial=pairs_at(reference,2.8)
    touches_grip=np.zeros(len(reference),dtype=bool)
    for i,j in initial:
        if fixed[i] and free[j]:touches_grip[j]=True
        if fixed[j] and free[i]:touches_grip[i]=True
    states=[];manifest=[];excluded=[]
    for path in sorted(args.run.glob('state_*/result.json')):
        row=json.loads(path.read_text(encoding='utf-8'))
        if not row['converged']:
            excluded.append(dict(state=row['state'],reason='saved state is not force converged'));continue
        with np.load(path.parent/'raw.npz') as data:positions=data['positions'].copy()
        for graph in row['connectivity']:
            if len(pairs_at(positions,graph['cutoff_A']))!=graph['pair_count']:
                raise ValueError('direct all-pairs calculation differs from saved neighbor list')
        states.append((row,positions));manifest.append(dict(state=row['state'],file=path.parent.name+'/raw.npz',
            sha256=hashlib.sha256((path.parent/'raw.npz').read_bytes()).hexdigest()))
    if len(states)<2:raise ValueError('at least two converged states required')
    theta=.371
    rotation=np.array([[np.cos(theta),-np.sin(theta),0.],[np.sin(theta),np.cos(theta),0.],[0.,0.,1.]])
    shift=np.array([.47,-.21,.89])
    intervals=[];atom_rows=[];pair_rows=[];largest_change=None;max_invariance_error=0.
    for (a,ra),(b,rb) in zip(states,states[1:]):
        coefficients,residual,condition=fit_affine(ra[free],rb[free]);norm=np.linalg.norm(residual,axis=1)
        _,rotated,_=fit_affine(ra[free]@rotation+shift,rb[free]@rotation+shift)
        invariant_error=float(np.max(abs(np.linalg.norm(rotated,axis=1)-norm)))
        max_invariance_error=max(max_invariance_error,invariant_error)
        if invariant_error>1e-10:raise ValueError('geometric descriptor fails rigid-frame invariance')
        gram_error=float(np.max(abs(np.column_stack((ra[free],np.ones(len(ids)))).T@residual)))
        gram_relative=gram_error/max(float(np.linalg.norm(ra[free])*np.linalg.norm(rb[free])),1.)
        if gram_relative>1e-11:raise ValueError('least-squares orthogonality check failed')
        touches=touches_grip[free]
        row=dict(state_from=a['state'],state_to=b['state'],strain_from=a['strain'],strain_to=b['strain'],
            free_nonaffine_RMS_A=float(np.sqrt(np.mean(norm**2))),free_nonaffine_max_A=float(norm.max()),
            median_A=float(np.median(norm)),p90_A=float(np.quantile(norm,.9)),
            grip_neighbor_atoms=int(touches.sum()),other_free_atoms=int((~touches).sum()),
            grip_neighbor_RMS_A=float(np.sqrt(np.mean(norm[touches]**2))) if touches.any() else None,
            other_free_RMS_A=float(np.sqrt(np.mean(norm[~touches]**2))) if (~touches).any() else None,
            fraction_squared_residual_at_grip_neighbors=float(np.sum(norm[touches]**2)/np.sum(norm**2)),
            affine_design_condition=condition,rigid_frame_invariance_error_A=invariant_error,
            least_squares_orthogonality_relative_residual=gram_relative)
        intervals.append(row)
        for index,atom in enumerate(ids):
            atom_rows.append(dict(state_from=a['state'],state_to=b['state'],atom=int(atom),
                initial_grip_neighbor=bool(touches[index]),x_A=float(rb[atom,0]),y_A=float(rb[atom,1]),z_A=float(rb[atom,2]),
                residual_x_A=float(residual[index,0]),residual_y_A=float(residual[index,1]),residual_z_A=float(residual[index,2]),
                nonaffine_norm_A=float(norm[index])))
        for cutoff in (2.8,3.1,3.4,3.7):
            previous=pairs_at(ra,cutoff);current=pairs_at(rb,cutoff)
            for kind,pairs in [('formed',current-previous),('lost',previous-current)]:
                for i,j in sorted(pairs):
                    pair_rows.append(dict(state_from=a['state'],state_to=b['state'],cutoff_A=cutoff,change=kind,
                        atom_i=i,atom_j=j,contains_grip_atom=bool(fixed[i] or fixed[j]),
                        previous_length_A=float(np.linalg.norm(ra[i]-ra[j])),new_length_A=float(np.linalg.norm(rb[i]-rb[j]))))
        if largest_change is None or row['free_nonaffine_RMS_A']>largest_change[0]['free_nonaffine_RMS_A']:
            largest_change=(row,rb.copy(),norm.copy(),residual.copy())
    write_csv(args.output/'interval_nonaffine.csv',intervals);write_csv(args.output/'atom_nonaffine.csv',atom_rows)
    if pair_rows:write_csv(args.output/'neighbor_pair_changes.csv',pair_rows)
    summary=dict(converged_states=len(states),intervals=intervals,excluded_saved_states=excluded,
        maximum_rigid_frame_invariance_error_A=max_invariance_error,
        grip_neighbor_definition='free atom connected to a fixed grip atom by an initial lattice pair shorter than 2.8 A',
        affine_fit='unweighted least squares over all 216 free atoms, full 3D affine map plus translation',
        geometry_sha256=hashlib.sha256((args.run/'geometry.npz').read_bytes()).hexdigest(),source_raw_files=manifest,
        scope='geometric rearrangement only; not energy localization, plastic strain, crack classification or a transition barrier',
        new_potential_calls=0,new_DFT=0,new_MD=0,initiation_probability=None)
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    row,r,norm,residual=largest_change
    fig,axes=plt.subplots(1,3,figsize=(12,5.2),layout='constrained')
    for ax,axis in zip(axes[:2],(0,1)):
        scatter=ax.scatter(r[free,axis],r[free,2],c=norm,s=22,cmap='viridis',vmin=0,vmax=float(norm.max()))
        ax.scatter(r[fixed,axis],r[fixed,2],s=12,c='#b2b8c2')
        ax.set_aspect('equal');ax.set_xlabel(('x' if axis==0 else 'y')+' (Angstrom)');ax.set_ylabel('z (Angstrom)')
        ax.set_title('Free atom residual after best affine map',fontsize=9)
    fig.colorbar(scatter,ax=list(axes[:2]),label='Nonaffine displacement norm (Angstrom)',shrink=.8)
    labels=[f'{100*x["strain_from"]:.0f} -> {100*x["strain_to"]:.0f}%' for x in intervals]
    axes[2].plot(labels,[x['free_nonaffine_RMS_A'] for x in intervals],'o-',label='All free atoms')
    axes[2].plot(labels,[x['grip_neighbor_RMS_A'] for x in intervals],'s--',label='Initial grip neighbors')
    axes[2].plot(labels,[x['other_free_RMS_A'] for x in intervals],'^--',label='Other free atoms')
    axes[2].set_ylabel('Nonaffine RMS (Angstrom)');axes[2].tick_params(axis='x',rotation=20);axes[2].legend(fontsize=8)
    axes[2].set_title('Each interval refits its own affine map',fontsize=9)
    fig.suptitle(f'Largest saved rearrangement: {100*row["strain_from"]:.0f}% to {100*row["strain_to"]:.0f}% imposed strain\nGeometry diagnostic only; no first-crack or energy-localization claim',fontsize=12)
    fig.savefig(args.output/'prism_nonaffine.png',dpi=180);plt.close(fig)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args())
