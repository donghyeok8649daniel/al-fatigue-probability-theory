"""All-pairs closest approach along a labelled same-grip straight chord.

The interpolation fraction is not time. This audits an initial path geometry,
not a minimum-energy path, crack classification, barrier, or probability.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar


def closest(before,after):
    i,j=np.triu_indices(len(before),1)
    initial=before[i]-before[j]
    velocity=(after[i]-after[j])-initial
    norm2=np.einsum('ij,ij->i',velocity,velocity)
    # Projection onto the closed segment is the exact geometric minimizer.
    t=np.clip(np.divide(-np.einsum('ij,ij->i',initial,velocity),norm2,
        out=np.zeros_like(norm2),where=norm2>0),0.,1.)
    distance=np.linalg.norm(initial+t[:,None]*velocity,axis=1)
    return i,j,t,distance,initial,velocity


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main(args):
    if args.output.exists():raise ValueError('fresh output required')
    with np.load(args.geometry) as d:
        free=d['free'].copy();reference=d['reference'].copy()
    snapshots=[];endpoints=[]
    for label,path in [('loading8',args.before),('return8',args.after)]:
        with np.load(path) as d:
            r=d['positions'].copy();z=d['numbers'].copy();f=d['forces'].copy()
            endpoints.append(dict(state=label,energy_eV=float(d['energy']),
                free_force_max_eV_A=float(np.linalg.norm(f[free],axis=1).max()),raw_sha256=sha(path)))
        if not np.all(z==14) or not np.all(np.isfinite(r)) or r.shape!=reference.shape:
            raise ValueError('finite same-specimen pure-Si geometry required')
        snapshots.append(r)
    before,after=snapshots;grip_error=float(np.max(abs(before[~free]-after[~free])))
    if grip_error>1e-12:raise ValueError('endpoint grips differ')
    i,j,t,distance,d0,dv=closest(before,after)
    reverse=closest(after,before)
    close_reverse=float(np.max(abs(reverse[3]-distance)))
    angle=.371;rot=np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1.]])
    transformed=closest(before@rot+np.array([.47,-.21,.89]),after@rot+np.array([.47,-.21,.89]))
    invariant_error=float(np.max(abs(transformed[3]-distance)))
    if max(close_reverse,invariant_error)>1e-10:raise ValueError('geometric invariance check failed')
    active=free[i]|free[j]
    # Independent bounded scalar minimization on the 20 shortest moving pairs.
    candidates=np.flatnonzero(active & (t>0) & (t<1))
    selected=candidates[np.argsort(distance[candidates])[:20]]
    independent=[]
    for k in selected:
        result=minimize_scalar(lambda x:float(np.linalg.norm((1-x)*(before[i[k]]-before[j[k]])+
            x*(after[i[k]]-after[j[k]]))**2),bounds=(0.,1.),method='bounded',options={'xatol':1e-12})
        error=abs(np.sqrt(result.fun)-distance[k])
        if not result.success or error>1e-9:raise ValueError('independent segment-distance minimization failed')
        independent.append(dict(atom_i=int(i[k]),atom_j=int(j[k]),distance_error_A=float(error)))
    rows=[]
    for k in np.argsort(distance)[:30]:
        rows.append(dict(atom_i=int(i[k]),atom_j=int(j[k]),closest_distance_A=float(distance[k]),
            interpolation_fraction=float(t[k]),before_distance_A=float(np.linalg.norm(d0[k])),
            after_distance_A=float(np.linalg.norm(d0[k]+dv[k])),both_free=bool(free[i[k]] and free[j[k]]),
            touches_grip=bool(not free[i[k]] or not free[j[k]])))
    fraction=np.linspace(0.,1.,257)
    minimum=np.array([np.linalg.norm(d0+x*dv,axis=1).min() for x in fraction])
    k=int(np.argmin(distance));movement=after[free]-before[free]
    result=dict(complete=True,atoms=len(before),free_atoms=int(free.sum()),all_pairs=len(i),
        endpoints=endpoints,geometry_sha256=sha(args.geometry),runner_sha256=sha(Path(__file__)),
        same_fixed_grip_max_difference_A=grip_error,
        labelled_free_atom_RMS_displacement_A=float(np.sqrt(np.mean(np.sum(movement**2,axis=1)))),
        labelled_free_atom_max_displacement_A=float(np.linalg.norm(movement,axis=1).max()),
        exact_continuous_chord_minimum_distance_A=float(distance[k]),
        closest_atom_pair=[int(i[k]),int(j[k])],closest_interpolation_fraction=float(t[k]),
        minimum_distance_at_start_A=float(minimum[0]),minimum_distance_at_end_A=float(minimum[-1]),
        sampled_257point_minimum_distance_A=float(minimum.min()),
        closest_pair_endpoint_distances_A=[float(np.linalg.norm(d0[k])),float(np.linalg.norm(d0[k]+dv[k]))],
        path_reversal_distance_error_A=close_reverse,rigid_frame_distance_error_A=invariant_error,
        independent_scalar_minimizations=independent,
        interpretation='unrelaxed labelled straight chord at identical fixed grips; compressed pairs motivate nonlinear path initialization, not a certified collision or transition barrier',
        interpolation_fraction_is_time=False,atom_permutation_performed=False,
        new_model_calls=0,new_DFT=0,new_MD=0,minimum_energy_path=False,
        first_crack_certified=False,physical_clock=None)
    args.output.mkdir(parents=True)
    with (args.output/'closest_pairs.csv').open('w',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    np.savez_compressed(args.output/'raw_geometry.npz',before=before,after=after,free=free,
        pair_i=i,pair_j=j,closest_fraction=t,closest_distance_A=distance,
        sampled_fraction=fraction,sampled_minimum_distance_A=minimum)
    (args.output/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(7,4.4),layout='constrained')
    ax.plot(fraction,minimum,color='#225c83',label='Shortest atom-pair distance')
    ax.scatter([t[k]],[distance[k]],color='#b53d41',zorder=3,label='Exact continuous minimum')
    ax.set(xlabel='Straight interpolation fraction (not time)',ylabel='Minimum pair distance (Angstrom)',
        title='Same 8% grips, different atomic arrangements\nGeometry only: no energy, barrier or crack label')
    ax.legend(fontsize=8);fig.savefig(args.output/'chord_geometry.png',dpi=170);plt.close(fig)
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    for name in ('before','after','geometry','output'):parser.add_argument('--'+name,type=Path,required=True)
    main(parser.parse_args())
