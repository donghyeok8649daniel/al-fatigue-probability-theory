"""Independent force/constraint replay and geometric diagnostics of saved states.

No additional potential evaluation. Connectivity and reference-pair changes are
reported at several radii and are not labelled as physical crack initiation.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components


def pairs_at(r,cutoff):
    return {tuple(map(int,p)) for p in cKDTree(r).query_pairs(cutoff,output_type='ndarray')}


def main(args):
    if args.output.exists() and any(args.output.iterdir()):raise ValueError('fresh snapshot-audit output required')
    args.output.mkdir(parents=True,exist_ok=True)
    with np.load(args.run/'geometry.npz') as data:
        reference=data['reference'].copy();lower=data['lower'].copy();upper=data['upper'].copy();free=data['free'].copy()
        area=float(data['area_A2']);gauge=float(data['gauge_length_A'])
    fixed=lower|upper;records=[];graphs=[];manifest=[];latest=None
    cutoffs=(2.8,3.1,3.4,3.7)
    initial_pairs={c:pairs_at(reference,c) for c in cutoffs}
    for path in sorted(args.run.glob('state_*/result.json')):
        row=json.loads(path.read_text(encoding='utf-8'))
        with np.load(path.parent/'raw.npz') as d:
            r=d['positions'].copy();force=d['forces'].copy();energy=float(d['energy'])
        if not np.isfinite(energy) or not np.all(np.isfinite(r)) or not np.all(np.isfinite(force)):
            raise ValueError('nonfinite raw state')
        exact=reference.copy();exact[lower,2]-=row['extension_A']/2;exact[upper,2]+=row['extension_A']/2
        grip_error=float(np.max(abs(r[fixed]-exact[fixed])))
        if grip_error>1e-12:raise ValueError('rigid grips differ from prescribed boundary')
        conjugate=float(.5*(force[lower,2].sum()-force[upper,2].sum()))
        force_vector=force.sum(axis=0);torque_vector=np.cross(r,force).sum(axis=0)
        reaction=-force[fixed].sum(axis=0);free_resultant=force[free].sum(axis=0)
        reaction_torque=-np.cross(r[fixed],force[fixed]).sum(axis=0)
        free_torque=np.cross(r[free],force[free]).sum(axis=0)
        force_max=float(np.linalg.norm(force[free],axis=1).max())
        expected=dict(energy_eV=energy,conjugate_force_eV_A=conjugate,nominal_stress_GPa=conjugate/area*160.2176634,
            reaction_force_residual_eV_A=float(np.linalg.norm(reaction)),reaction_torque_eV=float(np.linalg.norm(reaction_torque)),
            total_internal_force_eV_A=float(np.linalg.norm(force_vector)),total_internal_torque_eV=float(np.linalg.norm(torque_vector)),
            free_force_max_eV_A=force_max)
        difference=max(abs(expected[k]-row[k]) for k in expected)
        if difference>1e-11:raise ValueError('saved observable differs from independent raw replay')
        records.append(dict(state=row['state'],strain=row['strain'],converged=row['converged'],
            extension_A=row['extension_A'],grip_coordinate_error_A=grip_error,
            maximum_observable_replay_difference=difference,
            reaction_free_force_identity_error_eV_A=float(np.linalg.norm(reaction-free_resultant)),
            reaction_free_torque_identity_error_eV=float(np.linalg.norm(reaction_torque-free_torque)),
            support_force_residual_relative_to_axial=float(np.linalg.norm(reaction)/abs(conjugate)) if conjugate else None,
            **expected))
        for cutoff in cutoffs:
            pairs=pairs_at(r,cutoff);old=initial_pairs[cutoff];new=pairs-old;lost=old-pairs
            indices=np.array(sorted(pairs),dtype=int).reshape(-1,2)
            graph=coo_matrix((np.ones(len(indices)),(indices[:,0],indices[:,1])),shape=(len(r),len(r)))
            count,labels=connected_components(graph,directed=False)
            graphs.append(dict(state=row['state'],cutoff_A=cutoff,pairs=len(pairs),formed_relative_to_initial_lattice=len(new),
                lost_relative_to_initial_lattice=len(lost),components=int(count),
                grip_connected=bool(set(labels[lower])&set(labels[upper])),
                largest_component_atoms=int(np.bincount(labels).max())))
        manifest.append(dict(state=row['state'],file=path.parent.name+'/raw.npz',
            sha256=hashlib.sha256((path.parent/'raw.npz').read_bytes()).hexdigest()))
        if row['converged']:latest=(r,row)
    if not records:raise ValueError('no saved result states')
    for filename,rows in [('force_and_constraint_replay.csv',records),('connectivity_radius_sensitivity.csv',graphs)]:
        with (args.output/filename).open('w',newline='',encoding='utf-8') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    # A trapezoidal force integral over the stored branch is a coarse diagnostic.
    # Its mismatch includes quadrature, residual relaxation and possible branch jumps.
    work=[]
    for a,b in zip(records,records[1:]):
        trapezoid=.5*(a['conjugate_force_eV_A']+b['conjugate_force_eV_A'])*(b['extension_A']-a['extension_A'])
        change=b['energy_eV']-a['energy_eV']
        work.append(dict(state_from=a['state'],state_to=b['state'],energy_change_eV=change,
            coarse_trapezoid_work_eV=trapezoid,residual_eV=change-trapezoid,
            meaning='not dissipation: includes finite load spacing, relaxation residual and any branch change'))
    summary=dict(states_verified=len(records),converged_states=sum(r['converged'] for r in records),
        maximum_observable_replay_difference=max(r['maximum_observable_replay_difference'] for r in records),
        maximum_grip_coordinate_error_A=max(r['grip_coordinate_error_A'] for r in records),
        geometry_sha256=hashlib.sha256((args.run/'geometry.npz').read_bytes()).hexdigest(),
        source_raw_files=manifest,coarse_work_diagnostics=work,new_potential_calls=0,new_DFT=0,new_MD=0,
        topology_is_first_crack_detector=False,initiation_probability=None,physical_clock=None)
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    if latest is not None:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.collections import LineCollection
        r,row=latest;cutoff=2.8;original=initial_pairs[cutoff];current=pairs_at(r,cutoff)
        fig,axes=plt.subplots(1,4,figsize=(12,6),layout='constrained')
        for ax,(positions,axis,title,bonds) in zip(axes,[(reference,0,'Initial lattice, x-z',original),
            (reference,1,'Initial lattice, y-z',original),(r,0,'Relaxed, x-z',current),(r,1,'Relaxed, y-z',current)]):
            segments=np.array([positions[list(pair)][:,[axis,2]] for pair in sorted(bonds)])
            ax.add_collection(LineCollection(segments,colors='#708090',linewidths=.55,alpha=.55,zorder=1))
            if positions is r:
                formed=current-original
                if formed:ax.add_collection(LineCollection([positions[list(pair)][:,[axis,2]] for pair in sorted(formed)],colors='#c26936',linewidths=1.4,zorder=2))
            ax.scatter(positions[free,axis],positions[free,2],s=15,c='#4f8da8',edgecolors='none',zorder=3)
            ax.scatter(positions[fixed,axis],positions[fixed,2],s=17,c='#8b729f',edgecolors='none',zorder=3)
            ax.set_aspect('equal');ax.autoscale_view();ax.set_title(title,fontsize=10)
            ax.set_xlabel(('x' if axis==0 else 'y')+' (Angstrom)');ax.set_ylabel('z (Angstrom)')
        fig.suptitle(f'Intact bare prism: {len(r)} atoms; last converged imposed strain {row["strain"]:.3f}\n'
            'Purple: rigid grips. Orange: new distance-neighbor pairs (2.8 A cutoff), not certified cracks.',fontsize=11)
        fig.savefig(args.output/'prism_geometry.png',dpi=180);plt.close(fig)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--run',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args())
