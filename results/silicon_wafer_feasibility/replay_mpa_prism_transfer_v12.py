"""Independent algebraic replay of same-geometry cross-model prism observations."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np


def main(args):
    output=args.result/'raw_replay.json'
    if output.exists():raise ValueError('fresh transfer replay required')
    summary=json.loads((args.result/'summary.json').read_text(encoding='utf-8'))
    if not summary['complete']:raise ValueError('complete transfer audit required')
    maximum=0.;sources=[];energies={}
    for row in summary['states']:
        path=args.result/(row['state']+'_raw.npz')
        with np.load(path) as raw:
            r=raw['positions'];free=raw['free'];lower=raw['lower'];upper=raw['upper'];area=float(raw['area_A2'])
            expected={}
            for label in ('baseline','candidate'):
                f=raw[label+'_forces'];e=float(raw[label+'_energy'])
                expected[label+'_energy_eV']=e
                expected[label+'_free_force_max_eV_A']=float(np.linalg.norm(f[free],axis=1).max())
                expected[label+'_nominal_stress_GPa']=float((f[lower,2].sum()-f[upper,2].sum())*.5/area*160.2176634)
                energies[(row['state'],label)]=e
            f=raw['candidate_forces']
            expected['candidate_internal_force_norm_eV_A']=float(np.linalg.norm(f.sum(axis=0)))
            expected['candidate_internal_torque_norm_eV']=float(np.linalg.norm(np.cross(r,f).sum(axis=0)))
            expected['free_force_disagreement_component_RMS_eV_A']=float(np.sqrt(np.mean((f[free]-raw['baseline_forces'][free])**2)))
            maximum=max(maximum,max(abs(row[key]-value) for key,value in expected.items()))
        sources.append(dict(file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    for label in ('baseline','candidate'):
        delta=energies[('return8',label)]-energies[('loading8',label)]
        maximum=max(maximum,abs(delta-summary['same_grip_energy_comparison'][label+'_return8_minus_loading8_eV']))
    if maximum>1e-10:raise ValueError('same-geometry transfer raw replay failed')
    result=dict(complete=True,states=len(sources),maximum_scalar_replay_error=maximum,source_files=sources,
        new_model_calls=0,new_DFT=0,new_MD=0,physical_clock=None)
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8');print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--result',type=Path,required=True);main(p.parse_args())
