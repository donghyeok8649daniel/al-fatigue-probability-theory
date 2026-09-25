"""Recompute line-probe summaries from stored absolute energies and full forces."""
import argparse, hashlib, json
from pathlib import Path
import numpy as np
from scipy.constants import Boltzmann,electron_volt


def main(args):
    output=args.result/'raw_replay.json'
    if output.exists():raise ValueError('new replay path required')
    summary=json.loads((args.result/'summary.json').read_text(encoding='utf-8'))
    if not summary['complete']:raise ValueError('completed probe set required')
    rows=json.loads((args.result/'points.json').read_text(encoding='utf-8'))
    protocol=json.loads((args.result/'protocol.json').read_text(encoding='utf-8'))
    errors=[];sources=[];pair_errors=[];count=0
    kt=Boltzmann/electron_volt*300
    for name in ('loading8','return8'):
        path=args.result/(name+'_raw.npz')
        with np.load(path) as raw,np.load(args.results/('dense_'+name)/'raw_hessian.npz') as hraw:
            if not np.array_equal(raw['positions'],hraw['positions']) or not np.array_equal(raw['basis'],hraw['basis']):
                raise ValueError('probe reference coordinates or basis changed')
            stored_source=next(s for s in protocol['sources'] if s['state']==name)
            if hashlib.sha256((args.results/('dense_'+name)/'raw_hessian.npz').read_bytes()).hexdigest()!=stored_source['raw_sha256']:
                raise ValueError('source Hessian hash changed')
            h=(hraw['hessian']+hraw['hessian'].T)/2;g=hraw['gradient'];b=raw['basis']
            for (mode,x),energy,force in zip(raw['evaluated'],raw['energies'],raw['forces']):
                mode=int(mode);direction=raw['directions'][mode];v=b.T@direction.ravel()
                if abs(v@v-1)>1e-12:raise ValueError('probe direction not normalized')
                delta=float(energy-raw['baseline_energy']);curvature=float(v@h@v)
                quadratic=float(x*g@v+.5*x*x*curvature)
                exact_gradient=-b.T@force.ravel();nonlinear=exact_gradient-g-x*(h@v)
                matches=[r for r in rows if r['state']==name and r['mode']==mode and abs(r['coordinate_A']-x)<1e-14]
                if len(matches)!=1:raise ValueError('raw point identity ambiguous')
                row=matches[0]
                scalars={
                    'actual_delta_energy_eV':delta,'linear_quadratic_delta_energy_eV':quadratic,
                    'nonlinear_energy_residual_eV':delta-quadratic,
                    'nonlinear_energy_residual_over_kBT300':(delta-quadratic)/kt,
                    'actual_directional_gradient_eV_A':float(exact_gradient@v),
                    'nonlinear_full_gradient_norm_eV_A':float(np.linalg.norm(nonlinear)),
                    'nonlinear_directional_gradient_eV_A':float(nonlinear@v),
                    'maximum_atom_displacement_A':float(abs(x)*np.linalg.norm(direction,axis=1).max())}
                errors.extend(abs(row[key]-value) for key,value in scalars.items());count+=1
            selected=[r for r in rows if r['state']==name]
            for row in selected:
                if row['sign']<0:continue
                partner=next(r for r in selected if r['mode']==row['mode'] and r['label']==row['label'] and r['sign']==-1)
                q=row['coordinate_A'];v=b.T@raw['directions'][row['mode']].ravel()
                symmetric=.5*(row['actual_delta_energy_eV']+partner['actual_delta_energy_eV'])-.5*q*q*(v@h@v)
                antisymmetric=.5*(row['actual_delta_energy_eV']-partner['actual_delta_energy_eV'])-q*(g@v)
                pair_errors.extend([abs(symmetric-row['symmetric_nonlinear_energy_eV']),abs(antisymmetric-row['antisymmetric_nonlinear_energy_eV'])])
        sources.append(dict(file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    maximum=max(errors+pair_errors)
    if maximum>1e-9 or count!=summary['actual_points']:raise ValueError('raw nonlinear replay failed')
    result=dict(complete=True,points=count,source_files=sources,maximum_scalar_replay_error=maximum,
        method='stored full forces and absolute energies; quadratic prediction recomputed with raw Hessian, not stored eigenvalue',
        new_potential_calls=0,new_DFT=0,new_MD=0,physical_clock=None)
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8');print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('results','result'):p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
