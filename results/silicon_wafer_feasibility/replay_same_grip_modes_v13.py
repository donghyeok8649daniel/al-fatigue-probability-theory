"""Independent spectral/projector replay of the frozen same-grip comparison."""
from __future__ import annotations
import argparse, csv, hashlib, json
from pathlib import Path
import numpy as np
from scipy.linalg import subspace_angles


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def main(args):
    if args.output.exists(): raise ValueError('Preserve prior audit output')
    report=json.loads((args.comparison/'summary.json').read_text(encoding='utf-8'))
    if report.get('complete') is not True: raise ValueError('Completed comparison required')
    states={}; raw_errors={}
    def compare(name,a,b,tolerance=1e-9):
        a,b=np.asarray(a,dtype=float),np.asarray(b,dtype=float)
        if a.shape!=b.shape or not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
            raise ValueError('Shape/nonfinite comparison: '+name)
        error=float(np.max(abs(a-b)))
        if error>tolerance: raise ValueError(name+' mismatch: '+str(error))
        raw_errors[name]=error
    for pin in report['sources']:
        name=pin['state'];folder=args.results/('dense_'+name)
        for filename in ('raw_hessian.npz','spectrum.npz','summary.json','protocol.json'):
            if sha(folder/filename)!=pin[filename+'_sha256']: raise ValueError('Pinned input changed')
        with np.load(folder/'raw_hessian.npz') as file: raw={k:file[k].copy() for k in file.files}
        with np.load(folder/'spectrum.npz') as file: values=file['eigenvalues'];vectors=file['eigenvectors']
        states[name]=dict(raw=raw,values=values,vectors=vectors,cartesian=raw['basis']@vectors)
    if set(states)!={'loading8','return8'}: raise ValueError('Exact two-state schedule required')
    first,last=states['loading8'],states['return8']
    chord=(last['raw']['positions']-first['raw']['positions']).ravel()
    length2=float(chord@chord)
    with np.load(args.comparison/'modal_overlap.npz') as file:
        expected_overlap=(first['cartesian'].T@last['cartesian'])**2
        compare('full_cartesian_overlap',file['squared_overlap'],expected_overlap,1e-12)
        compare('chord_projection',first['raw']['basis']@file['chord'],chord,1e-12)
    compare('lowest_overlap',expected_overlap[0,0],report['lowest_direction_squared_overlap'],1e-12)
    if [x['first_modes'] for x in report['subspaces']]!=[1,4,8,16]: raise ValueError('Subspace schedule changed')
    for row in report['subspaces']:
        count=row['first_modes'];a=first['cartesian'][:,:count];b=last['cartesian'][:,:count]
        # Use SciPy QR/SVD principal-angle implementation in full Cartesian space.
        angles=np.sort(np.degrees(subspace_angles(a,b)))
        compare('angles_'+str(count),angles,row['principal_angles_degrees'],1e-6)
        p,q=a@a.T,b@b.T
        compare('projector_trace_'+str(count),np.einsum('ij,ji->',p,q)/count,row['mean_squared_subspace_overlap'],1e-12)
        compare('load_chord_'+str(count),chord@p@chord/length2,row['loading8_chord_squared_projection_fraction'],1e-12)
        compare('return_chord_'+str(count),chord@q@chord/length2,row['return8_chord_squared_projection_fraction'],1e-12)
    directions=report['first_eight_direction_comparisons']
    csv_rows=list(csv.DictReader((args.comparison/'directions.csv').open(encoding='utf-8')))
    if len(directions)!=16 or len(csv_rows)!=16: raise ValueError('Direction count changed')
    for row,csv_row in zip(directions,csv_rows):
        name=row['direction_from'];source=states[name];target=states['return8' if name=='loading8' else 'loading8']
        index=row['mode']-1;components=target['cartesian'].T@source['cartesian'][:,index]
        cross=float(target['values']@(components**2))
        key=name+'_'+str(index+1)
        compare('spectral_cross_curvature_'+key,cross,row['other_state_same_direction_curvature_eV_A2'])
        compare('csv_cross_curvature_'+key,cross,float(csv_row['other_state_same_direction_curvature_eV_A2']))
        compare('best_overlap_'+key,np.max(components**2),row['best_squared_overlap'],1e-12)
        if int(np.argmax(components**2))+1 != row['best_overlap_other_mode']: raise ValueError('Best direction differs')
    for row in report['finite_chord_quadratic_extrapolations']:
        name=row['from_state'];source=states[name];sign=1 if name=='loading8' else -1
        cartesian=sign*chord;normal=source['cartesian'].T@cartesian
        quadratic=.5*float(source['values']@(normal**2))
        linear=float(source['raw']['gradient']@(source['raw']['basis'].T@cartesian))
        actual=sign*float(last['raw']['energy']-first['raw']['energy'])
        compare('endpoint_energy_'+name,actual,row['actual_endpoint_energy_difference_eV'])
        compare('spectral_quadratic_'+name,quadratic,row['quadratic_term_eV'])
        compare('finite_chord_prediction_'+name,linear+quadratic,row['quadratic_predicted_difference_eV'])
    result=dict(complete=True,summary_sha256=sha(args.comparison/'summary.json'),runner_sha256=sha(Path(__file__)),
        independent_methods=['full Cartesian projectors','SciPy subspace_angles','spectral energy and Rayleigh sums'],
        checked_scalar_or_array_comparisons=len(raw_errors),maximum_absolute_difference=max(raw_errors.values()),
        per_comparison_errors=raw_errors,new_potential_calls=0,new_DFT=0,new_MD=0,
        scope='arithmetic/source integrity only; no new material, transition-path or first-crack certification')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='per_comparison_errors'},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    for name in ('results','comparison','output'):parser.add_argument('--'+name,type=Path,required=True)
    main(parser.parse_args())
