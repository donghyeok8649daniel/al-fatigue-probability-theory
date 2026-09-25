"""Spatial participation of static constrained curvature modes.

The first-k projector, unlike individual vectors in a degenerate eigenspace,
is invariant under rotation within that chosen subspace. A gap at its boundary
is reported, not assumed. Cartesian participation is not a mass-weighted
phonon, crack probability, committor, reaction path, or physical transition rate.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from solver_v1.silicon_structure_diagnostics_v12 import pair_set,neighbor_lists


def participation(cartesian_modes,masks):
    modes=np.asarray(cartesian_modes,float)
    if modes.ndim!=3 or modes.shape[1]!=3 or modes.shape[2]<1:
        raise ValueError('atom x Cartesian x mode array required')
    count=modes.shape[2]
    weight=np.sum(modes*modes,axis=(1,2))/count
    if abs(weight.sum()-1)>1e-10:raise ValueError('normalized orthonormal modes required')
    effective=1/float(weight@weight)
    return dict(subspace_dimension=count,effective_atoms=effective,
        participation_ratio=effective/len(weight),maximum_atom_weight=float(weight.max()),
        **{name+'_weight':float(weight[mask].sum()) for name,mask in masks.items()}),weight


def main(args):
    if args.output.exists():raise ValueError('fresh localization output required')
    with np.load(args.geometry) as d:
        reference=d['reference'];free=d['free'];lower=d['lower'];upper=d['upper']
    neighbors=neighbor_lists(len(reference),pair_set(reference,2.8))
    coordination=np.array([len(n) for n in neighbors])
    grip_neighbors=free&np.array([any(not free[j] for j in n) for n in neighbors])
    masks=dict(free_initial_surface=free&(coordination<4),free_initial_core=free&(coordination==4),
        grips=~free,free_grip_neighbors=grip_neighbors,other_free=free&~grip_neighbors)
    # Only the first three masks form a partition. The last two are another
    # decomposition of the free atoms, not additional disjoint material regions.
    if np.any(sum(masks[k].astype(int) for k in ('free_initial_surface','free_initial_core','grips'))!=1):
        raise ValueError('reference labels do not partition atoms')
    rows=[];states=[];sources=[]
    for name in ('force100','loading8','return8','loading10'):
        folder=args.results/('dense_'+name)
        summary=json.loads((folder/'summary.json').read_text(encoding='utf-8'))
        if not summary['complete']:raise ValueError('all full matrices must be complete')
        with np.load(folder/'raw_hessian.npz') as d:basis=d['basis'];positions=d['positions']
        with np.load(folder/'spectrum.npz') as d:values=d['eigenvalues'];vectors=d['eigenvectors']
        cartesian=(basis@vectors[:,:8]).reshape(len(reference),3,8)
        state_weights={}
        for count in (1,4,8):
            metrics,weight=participation(cartesian[:,:,:count],masks)
            # A deterministic orthogonal change inside the chosen subspace
            # must leave its atom projector diagonal unchanged.
            q,_=np.linalg.qr(np.random.default_rng(431+count).normal(size=(count,count)))
            _,rotated=participation(cartesian[:,:,:count]@q,masks)
            rotation_error=float(np.max(abs(weight-rotated)))
            if rotation_error>1e-12:raise ValueError('subspace participation changed under basis rotation')
            row=dict(state=name,**metrics,minimum_curvature_eV_A2=float(values[0]),
                last_included_curvature_eV_A2=float(values[count-1]),
                boundary_eigen_gap_eV_A2=float(values[count]-values[count-1]),
                subspace_basis_rotation_error=rotation_error)
            rows.append(row);state_weights[count]=weight
        states.append((name,positions,state_weights))
        sources.append(dict(state=name,raw_hessian_sha256=hashlib.sha256((folder/'raw_hessian.npz').read_bytes()).hexdigest(),
            spectrum_sha256=hashlib.sha256((folder/'spectrum.npz').read_bytes()).hexdigest()))
    args.output.mkdir(parents=True)
    with (args.output/'participation.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    np.savez_compressed(args.output/'atom_weights.npz',**{name+'_first'+str(count):weight for name,positions,weights in states for count,weight in weights.items()})
    result=dict(complete=True,geometry_sha256=hashlib.sha256(args.geometry.read_bytes()).hexdigest(),sources=sources,
        region_atoms={name:int(mask.sum()) for name,mask in masks.items()},states=rows,
        region_partition='free_initial_surface + free_initial_core + grips; grip-neighbor masks separately overlap free regions',
        normalization='w_i=sum over Cartesian and first k mode components squared / k; sum_i w_i=1; effective_atoms=1/sum_i w_i^2',
        subspace='first1,4,8 by algebraic curvature; finite boundary eigen-gap recorded; no mass weighting',
        new_potential_calls=0,new_DFT=0,new_MD=0,crack_initiation_label=None,physical_clock=None)
    (args.output/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,4,figsize=(11.5,5.3),layout='constrained')
    maximum=max(weights[8].max() for name,positions,weights in states)
    for ax,(name,positions,weights) in zip(axes,states):
        dots=ax.scatter(positions[:,0],positions[:,2],c=weights[8],vmin=0,vmax=maximum,s=18,cmap='viridis')
        ax.set_title(name);ax.set_xlabel('x (Angstrom)');ax.set_aspect('equal')
    axes[0].set_ylabel('z (Angstrom)')
    fig.colorbar(dots,ax=axes,shrink=.7,label='First 8 mode Cartesian weight per atom')
    fig.suptitle('Static curvature subspace participation (x-z projection)\nNo mass weighting, fracture labels or rate calibration',fontsize=11)
    fig.savefig(args.output/'curvature_localization.png',dpi=180);plt.close(fig)
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('results','geometry','output'):p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
