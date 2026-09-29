"""Independent raw arithmetic and SHA replay, with zero model evaluations."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.constants import Boltzmann, electron_volt


def main(root, result, output, cloud_override=None):
    if output.exists(): raise ValueError('fresh receipt required')
    errors={};checks=0
    def check(name,value,tolerance):
        nonlocal checks
        checks+=1;errors[name]=max(errors.get(name,0.),float(value))
        if not np.isfinite(value) or value>tolerance: raise ValueError((name,value,tolerance))
    input_record=json.loads((result/'input_manifest.json').read_text())
    manifest=input_record['sha256']
    for path,digest in manifest.items():
        allowed={digest,*input_record.get('accepted_line_ending_sha256',{}).get(path,[])}
        check('input_hash_mismatch',int(hashlib.sha256((root/path).read_bytes()).hexdigest() not in allowed),0)
    cloud=cloud_override or result/'full_dimensional_probes'
    for directory in (result,cloud):
        for name,digest in json.loads((directory/'artifact_manifest.json').read_text())['sha256'].items():
            check('artifact_hash_mismatch',int(hashlib.sha256((directory/name).read_bytes()).hexdigest()!=digest),0)
    summary=json.loads((cloud/'summary.json').read_text())
    rows=json.loads((cloud/'points.json').read_text())
    protocol=json.loads((cloud/'protocol.json').read_text())
    geometry=np.load(root/'results/silicon_initiation_v11/intact_prism_360_linesearch_v2/geometry.npz')
    free=geometry['free'];cache={};pairs={}
    for name in ('loading8','return8'):
        path=root/('results/silicon_initiation_v13/dense_'+name+'/raw_hessian.npz')
        check('cloud_source_hash_mismatch',int(hashlib.sha256(path.read_bytes()).hexdigest()!=protocol['source_sha256'][name]),0)
        d=np.load(path);h=(d['hessian']+d['hessian'].T)/2
        center=-np.linalg.solve(h,d['gradient'])
        cache[name]=(d,h,center)
    for i,row in enumerate(rows):
        raw=np.load(cloud/f'point_{i:03d}.npz');source,h,center=cache[row['state']]
        d=raw['reduced_displacement'];z=raw['normal_vector'];kt=Boltzmann/electron_volt*row['temperature_K']
        predicted=source['gradient']@d+.5*d@h@d
        expected=-.5*source['gradient']@np.linalg.solve(h,source['gradient'])+.5*kt*(z@z)
        check('Gaussian_energy_identity_eV',abs(predicted-expected),1e-10)
        check('position_lift_A',np.max(abs(raw['positions']-(source['positions']+(source['basis']@d).reshape(-1,3)))),1e-12)
        check('fixed_grip_displacement_A',np.max(abs(raw['positions'][~free]-source['positions'][~free])),0)
        e=float(raw['energy'])-float(source['energy'])
        residual=e-predicted
        check('quadratic_delta_replay_eV',abs(predicted-row['quadratic_delta_energy_eV']),1e-10)
        check('actual_energy_replay_eV',abs(e-row['actual_delta_energy_eV']),1e-8)
        check('energy_residual_replay_eV',abs(residual-row['energy_residual_eV']),1e-8)
        check('weight_log_replay',abs(-residual/kt-row['log_importance_weight']),1e-6)
        g=-source['basis'].T@raw['forces'].ravel()
        check('full_gradient_norm_replay_eV_A',abs(np.linalg.norm(g-source['gradient']-h@d)-row['full_gradient_residual_norm_eV_A']),1e-10)
        key=(row['state'],row['seed'],row['sample'],row['temperature_K'])
        pairs.setdefault(key,{})[row['sign']]=d
        generated=np.random.default_rng(row['seed']).standard_normal((4,len(d)))[row['sample']]
        check('normal_vector_seed_replay',np.max(abs(z-generated)),0)
    for key,pair in pairs.items():
        if len(pair)==2:
            check('antithetic_center_identity_A',np.max(abs(pair[-1]+pair[1]-2*cache[key[0]][2])),1e-10)
    for group in summary['groups']:
        logs=np.array([v['log_importance_weight'] for v in rows if v['state']==group['state'] and v['temperature_K']==group['temperature_K']])
        weights=np.exp(logs-logs.max());weights/=weights.sum()
        check('weight_ESS_replay',abs(1/(weights@weights)-group['weight_ESS_diagnostic']),1e-10)
    check('saved_point_count',abs(len(rows)-summary['points']),0)
    check('model_call_accounting',abs(summary['new_model_calls']+summary.get('reused_point_evaluations',0)-len(rows)-len(summary['source_replays'])),0)
    receipt=dict(complete=True,checks=checks,maximum_errors=errors,
                 cloud_complete=summary['complete'],points_replayed=len(rows),
                 new_model_evaluations=0,new_MD=0,
                 replay_does_not_validate_material_or_thermal_equilibrium=True)
    output.write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(receipt,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path('.'))
    p.add_argument('--result',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--cloud',type=Path)
    a=p.parse_args();main(a.root,a.result,a.output,a.cloud)
