"""Compare two frozen Cartesian tangent spaces at exactly the same grips.

Geometric modal overlaps and finite-chord quadratic extrapolations only.
No new potential, relaxation, path, barrier, phonon, probability or clock.
"""
from __future__ import annotations
import argparse, csv, hashlib, json
from pathlib import Path
import numpy as np


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def main(args):
    if args.output.exists(): raise ValueError('Fresh output directory required')
    with np.load(args.geometry) as file:
        free = file['free']; reference = file['reference']
    states = []; pins = []
    for name in ('loading8', 'return8'):
        folder = args.results / ('dense_' + name)
        summary = json.loads((folder/'summary.json').read_text(encoding='utf-8'))
        protocol = json.loads((folder/'protocol.json').read_text(encoding='utf-8'))
        if summary.get('complete') is not True or summary.get('independent_checked_directions') != 6:
            raise ValueError('Completed matrix and six direction checks required')
        if protocol['ensemble'] != 'displacement' or protocol['geometry_sha256'] != sha(args.geometry):
            raise ValueError('Same fixed-grip geometry required')
        with np.load(folder/'raw_hessian.npz') as file:
            raw = {key: file[key].copy() for key in file.files}
        with np.load(folder/'spectrum.npz') as file:
            values, vectors = file['eigenvalues'], file['eigenvectors']
        for key in ('hessian','basis','positions','gradient','energy'):
            if not np.all(np.isfinite(raw[key])): raise ValueError('Nonfinite raw input: '+key)
        if not np.all(np.isfinite(values)) or not np.all(np.isfinite(vectors)):
            raise ValueError('Nonfinite spectrum')
        dimension = int(3*free.sum())
        if raw['hessian'].shape != (dimension,dimension) or raw['basis'].shape != (3*len(free),dimension):
            raise ValueError('Full free Cartesian basis required')
        if values.shape != (dimension,) or vectors.shape != (dimension,dimension):
            raise ValueError('Full eigensystem required')
        if values[0] <= 0 or np.any(np.diff(values) < 0): raise ValueError('Positive ordered spectrum required')
        h = (raw['hessian']+raw['hessian'].T)/2
        if np.max(abs(raw['basis'][np.repeat(~free,3)])) != 0:
            raise ValueError('Grips must have no tangent motion')
        b_error = float(np.max(abs(raw['basis'].T@raw['basis']-np.eye(dimension))))
        v_error = float(np.max(abs(vectors.T@vectors-np.eye(dimension))))
        eigen_error = float(np.linalg.norm(h@vectors-vectors*values))
        if max(b_error,v_error) > 1e-11 or eigen_error > 1e-8:
            raise ValueError('Source basis/eigensystem does not replay')
        states.append(dict(name=name,raw=raw,h=h,values=values,vectors=vectors,
                           checks=dict(basis_max_error=b_error,eigenvector_max_error=v_error,eigen_residual_norm=eigen_error)))
        pins.append(dict(state=name,**{filename+'_sha256':sha(folder/filename) for filename in
            ('raw_hessian.npz','spectrum.npz','summary.json','protocol.json')}))
    first,last = states
    if not np.array_equal(first['raw']['basis'],last['raw']['basis']) or not np.array_equal(first['raw']['numbers'],last['raw']['numbers']):
        raise ValueError('Coordinates or atom labels differ')
    if not np.all(first['raw']['numbers']==14) or first['raw']['positions'].shape != reference.shape:
        raise ValueError('Same pure-Si specimen required')
    grip_error = float(np.max(abs(first['raw']['positions'][~free]-last['raw']['positions'][~free])))
    if grip_error > 1e-12: raise ValueError('Different imposed grip positions')
    basis=first['raw']['basis']; v=first['vectors']; w=last['vectors']
    overlap = v.T@w
    completeness=max(float(np.max(abs(np.sum(overlap**2,axis=0)-1))),float(np.max(abs(np.sum(overlap**2,axis=1)-1))))
    if completeness > 1e-11: raise ValueError('Modal overlap completeness failed')
    chord=(last['raw']['positions']-first['raw']['positions']).ravel()
    d=basis.T@chord; length=float(np.linalg.norm(d))
    projection_error=float(np.max(abs(basis@d-chord)))
    if length <= 0 or projection_error > 1e-12: raise ValueError('Nonzero chord in the same free coordinates required')
    delta_energy=float(last['raw']['energy']-first['raw']['energy'])
    subspaces=[]; identity_max=0.; rotation_max=0.
    for count in (1,4,8,16):
        a,b=v[:,:count],w[:,:count]; block=a.T@b
        square=float(np.sum(block**2)); singular=np.linalg.svd(block,compute_uv=False)
        # Principal angles: permit rounding at 1 only after checking the excess.
        if np.max(singular)>1+1e-12: raise ValueError('Invalid subspace cosine')
        singular=np.minimum(singular,1.)
        angles=np.degrees(np.arccos(singular))
        projectors=float(np.linalg.norm(a@a.T-b@b.T,'fro')**2)
        identity_error=abs(projectors-(2*count-2*square)); identity_max=max(identity_max,identity_error)
        rng=np.random.default_rng(1709+count)
        qa,_=np.linalg.qr(rng.normal(size=(count,count)));qb,_=np.linalg.qr(rng.normal(size=(count,count)))
        rotated=float(np.sum(((a@qa).T@(b@qb))**2))
        rotation_error=abs(rotated-square);rotation_max=max(rotation_max,rotation_error)
        subspaces.append(dict(first_modes=count,mean_squared_subspace_overlap=square/count,
            principal_angles_degrees=angles.tolist(),projector_identity_error=identity_error,
            internal_basis_rotation_error=rotation_error,
            loading8_boundary_gap_eV_A2=float(first['values'][count]-first['values'][count-1]),
            return8_boundary_gap_eV_A2=float(last['values'][count]-last['values'][count-1]),
            loading8_chord_squared_projection_fraction=float(np.sum((a.T@d)**2)/length**2),
            return8_chord_squared_projection_fraction=float(np.sum((b.T@d)**2)/length**2)))
    if identity_max>1e-10 or rotation_max>1e-10: raise ValueError('Independent subspace identity failed')
    directions=[]
    for origin,source,target,matrix in [('loading8',first,last,overlap),('return8',last,first,overlap.T)]:
        for index in range(8):
            direction=source['vectors'][:,index]; target_rayleigh=float(direction@target['h']@direction)
            other_index=int(np.argmax(matrix[index]**2))
            directions.append(dict(direction_from=origin,mode=index+1,
                own_curvature_eV_A2=float(source['values'][index]),other_state_same_direction_curvature_eV_A2=target_rayleigh,
                best_overlap_other_mode=other_index+1,best_squared_overlap=float(matrix[index,other_index]**2),
                overlap_with_same_ordered_mode_squared=float(matrix[index,index]**2)))
    extrapolations=[]
    for state,vector,actual in [(first,d,delta_energy),(last,-d,-delta_energy)]:
        linear=float(state['raw']['gradient']@vector);quadratic=.5*float(vector@state['h']@vector)
        extrapolations.append(dict(from_state=state['name'],actual_endpoint_energy_difference_eV=actual,
            linear_term_eV=linear,quadratic_term_eV=quadratic,quadratic_predicted_difference_eV=linear+quadratic,
            prediction_minus_actual_eV=linear+quadratic-actual,
            unit_chord_local_curvature_eV_A2=2*quadratic/length**2))
    result=dict(complete=True,geometry_sha256=sha(args.geometry),sources=pins,
        runner_sha256=sha(Path(__file__)),dimension=len(d),same_fixed_grip_difference_A=grip_error,
        raw_eigensystem_checks=[dict(state=s['name'],**s['checks']) for s in states],
        overlap_completeness_error=completeness,free_coordinate_chord_norm_A=length,
        free_atom_chord_RMS_A=length/np.sqrt(free.sum()),chord_projection_error_A=projection_error,
        lowest_direction_squared_overlap=float(overlap[0,0]**2),subspaces=subspaces,
        first_eight_direction_comparisons=directions,finite_chord_quadratic_extrapolations=extrapolations,
        coordinate_convention='identical labelled free Cartesian coordinates; fixed grips; no mass weighting',
        limitations=['endpoint chord is not a transition path','quadratic extrapolation over a large chord is not a barrier',
            'approximate stationary states retain their gradient','ordered modes may mix; boundary gaps are reported',
            'mode weights are geometric and not probabilities'],
        new_potential_calls=0,new_DFT=0,new_MD=0,first_crack_certified=False,material_approved=False,physical_clock=None)
    args.output.mkdir(parents=True)
    np.savez_compressed(args.output/'modal_overlap.npz',squared_overlap=overlap**2,chord=d,
        loading8_values=first['values'],return8_values=last['values'])
    with (args.output/'directions.csv').open('w',newline='',encoding='utf-8') as file:
        writer=csv.DictWriter(file,fieldnames=list(directions[0]));writer.writeheader();writer.writerows(directions)
    (args.output/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    image=axes[0].imshow(overlap[:8,:8]**2,vmin=0,vmax=1,cmap='viridis',origin='lower')
    axes[0].set(xticks=range(8),xticklabels=range(1,9),yticks=range(8),yticklabels=range(1,9),
        xlabel='Returned 8% ordered direction',ylabel='Initial 8% ordered direction',title='Squared Cartesian overlap')
    fig.colorbar(image,ax=axes[0],shrink=.8,label='Geometric overlap (not probability)')
    for key,label in [('loading8','Initial 8%'),('return8','Returned 8%')]:
        axes[1].plot([r['first_modes'] for r in subspaces],
            [100*r[key+'_chord_squared_projection_fraction'] for r in subspaces],'o-',label=label)
    axes[1].set(xticks=[1,4,8,16],xlabel='Number of lowest-curvature directions',
        ylabel='Squared endpoint-chord projection (%)',title='Endpoint chord, not a reaction path')
    axes[1].legend(fontsize=8);axes[1].set_ylim(bottom=0)
    fig.suptitle('Same grips, different atomic states: static tangent comparison\nNo new model calls, transition barrier or initiation label',fontsize=11)
    fig.savefig(args.output/'same_grip_modes.png',dpi=180);plt.close(fig)
    print(json.dumps({key:result[key] for key in ('complete','lowest_direction_squared_overlap','subspaces',
         'finite_chord_quadratic_extrapolations','new_potential_calls')},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    for name in ('results','geometry','output'):parser.add_argument('--'+name,type=Path,required=True)
    main(parser.parse_args())
