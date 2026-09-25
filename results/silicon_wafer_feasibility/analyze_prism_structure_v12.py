"""Spatial, cutoff-sensitive structural evidence, with no binary crack label."""
from __future__ import annotations
import argparse,csv,hashlib,json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from solver_v1.silicon_structure_diagnostics_v12 import pair_set,neighbor_lists,local_affine,edge_bottleneck,reference_sections


def dump(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def table(path,rows):
    if not rows:return
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def main(args):
    if args.output.exists() and any(args.output.iterdir()):raise ValueError('fresh output required')
    args.output.mkdir(parents=True,exist_ok=True)
    run=args.source/'intact_prism_360_linesearch_v2'
    with np.load(run/'geometry.npz') as d:
        reference=d['reference'].copy();free=d['free'].copy();lower=d['lower'].copy();upper=d['upper'].copy()
    cases=[(f'loading_{i*2}pct',run/f'state_{i:03d}'/'raw.npz',run/f'state_{i:03d}'/'result.json') for i in range(5)]
    cases+=[('loading_10pct',args.source/'prism_10pct_continuation'/'raw.npz',args.source/'prism_10pct_continuation'/'summary.json'),
        ('return_8pct',args.source/'prism_unload_10_to_8'/'raw.npz',args.source/'prism_unload_10_to_8'/'summary.json')]
    initial=pair_set(reference,2.8);neighbors=neighbor_lists(len(reference),initial)
    coord=np.array([len(x) for x in neighbors]);surface=coord<4;core=coord==4
    grip_neighbor=np.array([any(not free[j] for j in n) for n in neighbors])&free
    states=[];manifest=[];stats=[];atom_rows=[];section_rows=[];pair_rows=[];graph_rows=[]
    with np.load(cases[0][1]) as d:zero=d['positions'].copy()
    for name,path,summary_path in cases:
        summary=json.loads(summary_path.read_text(encoding='utf-8'))
        if not summary['converged']:raise ValueError('all descriptor states must be force converged')
        with np.load(path) as d:r=d['positions'].copy()
        states.append((name,r));manifest.append(dict(name=name,raw_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            summary_sha256=hashlib.sha256(summary_path.read_bytes()).hexdigest()))
        for radius in (2.8,4.2):
            local_neighbors=neighbor_lists(len(reference),pair_set(reference,radius))
            affine=local_affine(zero,r,local_neighbors)
            for row in affine:
                atom=row['atom'];eigen=row.pop('green_principal_strains')
                atom_rows.append(dict(state=name,neighborhood_cutoff_A=radius,**row,free=bool(free[atom]),initial_coordination=int(coord[atom]),
                    initial_surface=bool(surface[atom]),initial_grip_neighbor=bool(grip_neighbor[atom]),
                    x_A=float(r[atom,0]),y_A=float(r[atom,1]),z_A=float(r[atom,2]),
                    green_min=None if eigen is None else eigen[0],green_middle=None if eigen is None else eigen[1],
                    green_max=None if eigen is None else eigen[2]))
        for cutoff in (2.7,2.8,3.1,3.4,3.7):
            current=pair_set(r,cutoff);old=pair_set(zero,cutoff)
            lost=old-current;new=current-old
            graph=edge_bottleneck(len(r),current,lower,upper)
            graph_rows.append(dict(state=name,cutoff_A=cutoff,pairs=len(current),edge_disjoint_paths=graph['edge_disjoint_paths'],
                lost_from_zero=len(lost),formed_from_zero=len(new),
                lost_core_core=int(sum(core[i] and core[j] for i,j in lost)),
                lost_touching_initial_surface=int(sum(surface[i] or surface[j] for i,j in lost)),
                lost_touching_grip=int(sum(not free[i] or not free[j] for i,j in lost))))
            if cutoff==2.8:
                for row in reference_sections(reference,r,old,current,free):section_rows.append(dict(state=name,**row))
                for kind,pairs in [('lost',lost),('formed',new),('topological_min_cut',graph['cut_pairs'])]:
                    for i,j in sorted(pairs):
                        pair_rows.append(dict(state=name,change=kind,atom_i=i,atom_j=j,
                            initial_surface_i=bool(surface[i]),initial_surface_j=bool(surface[j]),
                            touches_grip=bool(not free[i] or not free[j]),
                            old_distance_A=float(np.linalg.norm(zero[i]-zero[j])),distance_A=float(np.linalg.norm(r[i]-r[j])),
                            midpoint_z_A=float((r[i,2]+r[j,2])/2)))
        for radius in (2.8,4.2):
            rows=[row for row in atom_rows if row['state']==name and row['neighborhood_cutoff_A']==radius]
            for label,mask in [('free_surface',free&surface),('free_core',free&core),('grip_neighbors',grip_neighbor),('other_free',free&~grip_neighbor)]:
                values=[row['D2_A2'] for row in rows if mask[row['atom']] and row['rank']==3 and row['residual_degrees_of_freedom']>0]
                stats.append(dict(state=name,neighborhood_cutoff_A=radius,region=label,atoms=int(mask.sum()),resolved_residual_atoms=len(values),
                    mean_D2_A2=float(np.mean(values)) if values else None,max_D2_A2=float(np.max(values)) if values else None))
    table(args.output/'atom_local_affine.csv',atom_rows);table(args.output/'region_affine.csv',stats)
    table(args.output/'graph_cutoff.csv',graph_rows);table(args.output/'reference_sections.csv',section_rows)
    table(args.output/'pair_changes.csv',pair_rows)
    dump(args.output/'summary.json',dict(states=len(states),source_states=manifest,
        geometry_sha256=hashlib.sha256((run/'geometry.npz').read_bytes()).hexdigest(),
        initial_surface_atoms=int(surface.sum()),initial_core_atoms=int(core.sum()),
        free_initial_surface_atoms=int((free&surface).sum()),free_initial_core_atoms=int((free&core).sum()),
        graph_cutoff_results=graph_rows,region_local_affine=stats,
        reference_neighbors='initial ideal 2.8 and 4.2 A neighborhoods; local map starts at relaxed zero-extension state',
        residual_resolution='rank3 and strictly positive residual degrees of freedom required for regional D2; three independent neighbors exactly interpolate and are excluded',
        surface_definition='ideal initial coordination < 4; includes geometrical edges, no oxidation/passivation',
        graph_cuts='unit-capacity topological bottlenecks; not a physical fracture path or strength',
        D2='least-squares relative-neighbor residual per neighbor; not plastic strain',
        stress_from_geometry_not_inferred=True,new_potential_calls=0,new_DFT=0,new_MD=0,
        crack_initiation_label=None,physical_clock=None))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    selected=[states[0],states[4],states[5],states[6]]
    fig,axes=plt.subplots(1,4,figsize=(13,6),layout='constrained')
    for ax,(name,r) in zip(axes,selected):
        old=pair_set(zero,2.8);new=pair_set(r,2.8)
        ax.add_collection(LineCollection([r[[i,j]][:,[0,2]] for i,j in new],colors='#bcc5ce',linewidths=.45))
        ax.add_collection(LineCollection([r[[i,j]][:,[0,2]] for i,j in old-new],colors='#d1495b',linewidths=1.3))
        ax.scatter(r[free&surface,0],r[free&surface,2],s=15,c='#e6ac41',label='Initial surface')
        ax.scatter(r[free&core,0],r[free&core,2],s=15,c='#2166a5',label='Initial core')
        ax.scatter(r[~free,0],r[~free,2],s=12,c='#555555',label='Fixed grips')
        ax.autoscale();ax.set_aspect('equal');ax.set_title(name.replace('_',' '));ax.set_xlabel('x (Angstrom)')
    axes[0].set_ylabel('z (Angstrom)');axes[-1].legend(fontsize=7,loc='upper right')
    fig.suptitle('Red: initial relaxed pairs now longer than 2.8 Angstrom\nGeometric descriptors only; no first-crack classification',fontsize=12)
    fig.savefig(args.output/'spatial_pair_changes.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4.6),layout='constrained')
    for name,r in selected:
        rows=[x for x in section_rows if x['state']==name]
        axes[0].plot([x['reference_midpoint_A'] for x in rows],[x['current_crossing_pairs'] for x in rows],'.-',label=name)
        axes[1].plot([x['reference_midpoint_A'] for x in rows],[x['lost_initial_pairs'] for x in rows],'.-',label=name)
    for ax in axes:ax.set_xlabel('Reference section z (Angstrom)');ax.legend(fontsize=8)
    axes[0].set_ylabel('Current crossing-pair count');axes[1].set_ylabel('Lost relaxed-zero crossing pairs')
    fig.suptitle('Reference-labelled sections, 2.8 Angstrom convention; not stress or crack area')
    fig.savefig(args.output/'sectional_pair_counts.png',dpi=180);plt.close(fig)
    print(json.dumps(dict(states=len(states),free_surface=int((free&surface).sum()),free_core=int((free&core).sum())),indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args())
