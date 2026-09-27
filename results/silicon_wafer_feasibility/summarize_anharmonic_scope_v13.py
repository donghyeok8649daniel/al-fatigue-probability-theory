"""Resolve nonlinear force residuals along and transverse to each probed mode.

Saved force arrays only. A small energy residual along a line does not bound
unprobed mode coupling or certify a finite-temperature basin integral.
"""
from pathlib import Path
import argparse,csv,hashlib,json
import numpy as np

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main(args):
    if args.output.exists():raise ValueError('fresh output required')
    probes=args.results/'anharmonic_probes'
    verified=json.loads((args.results/'anharmonic_replay.json').read_text())
    if not verified['complete'] or not verified['raw_replay_passed'] or verified['verified_points']!=72:
        raise ValueError('complete independently replayed72-point data required')
    if sha(probes/'summary.json')!=verified['reported_summary_sha256']:
        raise ValueError('summary changed after replay')
    for entry in verified['source_manifest']:
        if sha(args.results/entry['path'])!=entry['sha256']:raise ValueError('replayed raw input changed')
    rows=json.loads((probes/'points.json').read_text());computed=[];grouped={};errors=[]
    for state in ('loading8','return8'):
        with np.load(args.results/('dense_'+state)/'raw_hessian.npz') as f:
            basis=f['basis'];h=(f['hessian']+f['hessian'].T)/2;g0=f['gradient']
        with np.load(probes/(state+'_raw.npz')) as f:
            directions=f['directions'];schedule=f['evaluated'];forces=f['forces']
        selected=[r for r in rows if r['state']==state]
        if len(selected)!=len(schedule):raise ValueError('point count mismatch')
        for row,(mode,q),force in zip(selected,schedule,forces):
            if int(mode)!=row['mode'] or abs(q-row['coordinate_A'])>1e-12:raise ValueError('point ordering differs')
            v=basis.T@directions[int(mode)].ravel()
            if abs(v@v-1)>1e-12:raise ValueError('mode normalization differs')
            actual=-basis.T@force.ravel();linear=q*(h@v)
            residual=actual-g0-linear;longitudinal=float(v@residual)
            transverse=residual-longitudinal*v
            full=float(np.linalg.norm(residual));trans=float(np.linalg.norm(transverse))
            for a,b in ((full,row['nonlinear_full_gradient_norm_eV_A']),(longitudinal,row['nonlinear_directional_gradient_eV_A'])):
                if not np.isfinite(a) or not np.isclose(a,b,rtol=1e-9,atol=1e-10):raise ValueError('saved gradient diagnostic differs')
                errors.append(abs(a-b))
            error=abs(full*full-longitudinal*longitudinal-trans*trans)
            if error>1e-12:raise ValueError('orthogonal decomposition failed')
            point=dict(state=state,mode=int(mode),label=row['label'],sign=row['sign'],
                maximum_atom_displacement_A=row['maximum_atom_displacement_A'],
                nonlinear_energy_residual_eV=row['nonlinear_energy_residual_eV'],
                nonlinear_energy_residual_over_kBT300=row['nonlinear_energy_residual_over_kBT300'],
                full_gradient_residual_norm_eV_A=full,longitudinal_gradient_residual_eV_A=longitudinal,
                transverse_gradient_residual_norm_eV_A=trans,
                transverse_fraction_of_residual_norm=trans/full if full else None,
                linear_gradient_increment_norm_eV_A=float(np.linalg.norm(linear)),
                projection_identity_error_eV2_A2=error)
            computed.append(point);grouped.setdefault((state,row['label']),[]).append(point)
    groups=[]
    for (state,label),points in sorted(grouped.items()):
        worst=max(points,key=lambda p:p['full_gradient_residual_norm_eV_A'])
        groups.append(dict(state=state,label=label,points=len(points),
            maximum_absolute_energy_residual_eV=max(abs(p['nonlinear_energy_residual_eV']) for p in points),
            maximum_absolute_energy_residual_over_kBT300=max(abs(p['nonlinear_energy_residual_over_kBT300']) for p in points),
            largest_full_gradient_residual_point=worst))
    args.output.mkdir(parents=True)
    with (args.output/'points.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(computed[0]));writer.writeheader();writer.writerows(computed)
    result=dict(complete=True,points=len(computed),groups=groups,new_model_calls=0,new_DFT=0,new_MD=0,
        largest_saved_gradient_replay_difference_eV_A=max(errors),
        maximum_projection_identity_error_eV2_A2=max(p['projection_identity_error_eV2_A2'] for p in computed),
        sources=dict(points_sha256=sha(probes/'points.json'),raw_replay_sha256=sha(args.results/'anharmonic_replay.json'),
            script_sha256=sha(Path(__file__))),
        scope='three straight Cartesian directions per state; transverse force is force perpendicular to that coordinate, not crack flux or transition rate',
        finite_temperature_free_energy_certified=False,first_crack_certified=False)
    (args.output/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='groups'},indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--results',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args())
