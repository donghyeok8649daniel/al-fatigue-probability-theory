"""Independent replay of a same-displacement continuation, retaining its parent.

Relaxation energy changes at fixed grips are optimizer results, not a physical
time history, dissipation measurement, transition barrier, or first-crack label.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def all_pairs(r,cutoff):
    i,j=np.triu_indices(len(r),1)
    lengths=np.linalg.norm(r[i]-r[j],axis=1)
    return set(zip(i[lengths<cutoff].tolist(),j[lengths<cutoff].tolist()))


def main(args):
    if args.output.exists() and any(args.output.iterdir()):raise ValueError('fresh output required')
    protocol=json.loads((args.continuation/'protocol.json').read_text(encoding='utf-8'))
    row=json.loads((args.continuation/'summary.json').read_text(encoding='utf-8'))
    if digest(args.parent/'raw.npz')!=protocol['source_raw_sha256']:
        raise ValueError('continuation refers to a different parent')
    if digest(args.geometry)!=protocol['geometry_sha256']:
        raise ValueError('continuation geometry changed')
    with np.load(args.geometry) as d:
        reference=d['reference'];lower=d['lower'];upper=d['upper'];free=d['free'];area=float(d['area_A2'])
    with np.load(args.parent/'raw.npz') as d:parent_positions=d['positions'];parent_energy=float(d['energy'])
    with np.load(args.continuation/'raw.npz') as d:
        r=d['positions'];force=d['forces'];energy=float(d['energy'])
    previous=json.loads((args.previous/'result.json').read_text(encoding='utf-8'))
    if not previous['converged']:raise ValueError('previous lower-strain reference must be converged')
    with np.load(args.previous/'raw.npz') as d:r_previous=d['positions']
    fixed=lower|upper
    grip_difference=float(abs(r[fixed]-parent_positions[fixed]).max())
    if grip_difference>1e-12:raise ValueError('continuation changed the fixed boundary')
    axial=float((force[lower,2].sum()-force[upper,2].sum())/2)
    support=-force[fixed].sum(axis=0)
    support_torque=-np.cross(r[fixed],force[fixed]).sum(axis=0)
    observed=dict(energy_eV=energy,conjugate_force_eV_A=axial,
        nominal_stress_GPa=axial/area*160.2176634,
        free_force_max_eV_A=float(np.linalg.norm(force[free],axis=1).max()),
        reaction_force_residual_eV_A=float(np.linalg.norm(support)),
        reaction_torque_eV=float(np.linalg.norm(support_torque)),
        total_internal_force_eV_A=float(np.linalg.norm(force.sum(axis=0))),
        total_internal_torque_eV=float(np.linalg.norm(np.cross(r,force).sum(axis=0))))
    error=max(abs(observed[k]-row[k]) for k in observed)
    if error>1e-11:raise ValueError('stored observables differ from independent raw replay')
    graphs=[];changes=[]
    for stored in row['connectivity']:
        radius=stored['cutoff_A'];pairs=all_pairs(r,radius);old=all_pairs(reference,radius)
        before=all_pairs(r_previous,radius);indices=np.array(sorted(pairs),dtype=int).reshape(-1,2)
        graph=coo_matrix((np.ones(len(indices)),(indices[:,0],indices[:,1])),shape=(len(r),len(r)))
        count,labels=connected_components(graph,directed=False)
        bridge=bool(set(labels[lower])&set(labels[upper]))
        if len(pairs)!=stored['pair_count'] or count!=stored['components'] or bridge!=stored['grip_connected']:
            raise ValueError('independent all-pairs graph differs from stored neighbor diagnostics')
        graphs.append(dict(cutoff_A=radius,components=int(count),grip_connected=bridge,pairs=len(pairs),
            lost_from_initial_lattice=len(old-pairs),formed_from_initial_lattice=len(pairs-old),
            lost_from_previous_converged=len(before-pairs),formed_from_previous_converged=len(pairs-before)))
        for kind,values in [('lost',before-pairs),('formed',pairs-before)]:
            for i,j in sorted(values):
                changes.append(dict(cutoff_A=radius,change=kind,atom_i=i,atom_j=j,
                    contains_grip_atom=bool(fixed[i] or fixed[j]),
                    previous_distance_A=float(np.linalg.norm(r_previous[i]-r_previous[j])),
                    current_distance_A=float(np.linalg.norm(r[i]-r[j]))))
    # A geometric descriptor is meaningful even for an interrupted geometry, but
    # is explicitly separated from any claim of a converged equilibrium branch.
    design=np.column_stack((r_previous[free],np.ones(int(free.sum()))))
    fitted,_,rank,_=np.linalg.lstsq(design,r[free],rcond=None)
    if rank!=4:raise ValueError('affine reference does not span 3D')
    residual=r[free]-design@fitted;norm=np.linalg.norm(residual,axis=1)
    summary=dict(converged=bool(row['converged']),interruption=row['interruption'],
        strain=row['strain'],previous_converged_strain=previous['strain'],
        fixed_grip_difference_A=grip_difference,maximum_observable_replay_difference=error,
        force_identity_error_eV_A=float(np.linalg.norm(support-force[free].sum(axis=0))),
        torque_identity_error_eV=float(np.linalg.norm(support_torque-np.cross(r[free],force[free]).sum(axis=0))),
        **observed,energy_change_at_fixed_displacement_eV=energy-parent_energy,
        fixed_grip_external_work_eV=0.,
        previous_to_current_nonaffine_RMS_A=float(np.sqrt(np.mean(norm**2))),
        previous_to_current_nonaffine_max_A=float(norm.max()),graphs=graphs,
        relaxation_not_physical_dissipation=True,geometric_diagnostic_not_crack_detector=True,
        new_potential_calls=0,new_DFT=0,new_MD=0,initiation_probability=None,physical_clock=None,
        sources=[dict(role=label,sha256=digest(path)) for label,path in
            [('geometry',args.geometry),('parent_raw',args.parent/'raw.npz'),
             ('continuation_raw',args.continuation/'raw.npz'),('continuation_protocol',args.continuation/'protocol.json'),
             ('continuation_summary',args.continuation/'summary.json'),('previous_raw',args.previous/'raw.npz'),
             ('previous_result',args.previous/'result.json')]])
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    (args.output/'pair_changes.json').write_text(json.dumps(changes,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(11,5.4),layout='constrained')
    for ax,positions,title in zip(axes[:2],[r_previous,r],['Previous converged state','Continued state']):
        ax.scatter(positions[free,0],positions[free,2],s=18,c='#4f8da8')
        ax.scatter(positions[fixed,0],positions[fixed,2],s=18,c='#8b729f')
        ax.set_aspect('equal');ax.set_title(title);ax.set_xlabel('x (Angstrom)');ax.set_ylabel('z (Angstrom)')
    points=axes[2].scatter(r[free,0],r[free,2],s=22,c=norm,cmap='viridis')
    axes[2].set_aspect('equal');axes[2].set_title('Residual after affine map');axes[2].set_xlabel('x (Angstrom)')
    fig.colorbar(points,ax=axes[2],label='Geometric residual (Angstrom)')
    fig.suptitle(f'Same-boundary continuation: force convergence {row["converged"]}\nGeometric comparison only; no first-crack classification',fontsize=12)
    fig.savefig(args.output/'continuation_geometry.png',dpi=180);plt.close(fig)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('parent','continuation','previous','geometry','output'):
        p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
