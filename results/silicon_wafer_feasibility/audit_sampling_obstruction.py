"""Cached-proposal audit; no new energy calls, equilibrium or physical clock."""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import numpy as np
from solver_v1.silicon_device_ensemble import FixedGripTarget
from solver_v1.silicon_sampling_score import nonlinear_correction_gradient


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def arrays(path):
    with np.load(path) as data: return {key:data[key].copy() for key in data.files}


def main(a):
    if a.output.exists(): raise ValueError('fresh output required')
    a.output.mkdir(parents=True)
    pilot=a.root/'results/silicon_crack_thermal_validation/nonlinear_pilot'
    protocol=json.loads((pilot/'protocol.json').read_text())
    if not json.loads((pilot/'summary.json').read_text())['complete']:
        raise ValueError('completed source pilot required')
    sources=[]
    for relative,expected in protocol['source_sha256'].items():
        if sha(a.root/relative)!=expected: raise ValueError('source changed')
        sources.append(arrays(a.root/relative))
    common=(sources[0]['positions']+sources[1]['positions'])/2
    results=[]; raw=[]
    for chain,source in enumerate(sources):
        t=FixedGripTarget(source['positions'],source['basis'],source['hessian'],source['gradient'],300.,common)
        files=sorted((pilot/f'chain_{chain}').glob('evaluation_*.npz'))
        if len(files)!=65: raise ValueError('source proposal record incomplete')
        data=[arrays(path) for path in files]
        projected=-t.basis.T@data[0]['forces'].ravel()
        np.testing.assert_allclose(projected,source['gradient'],rtol=0,atol=1e-8)
        scores=[nonlinear_correction_gradient(t,item['u'],item['forces']) for item in data]
        rng=np.random.default_rng(protocol['seeds'][chain]); current=0; rows=[]
        beta=protocol['proposal_scale']; persistence=np.sqrt(1-beta*beta)
        for i in range(1,len(data)):
            before,proposal=data[current],data[i]
            expected_u=persistence*before['u']+beta*rng.standard_normal(len(before['u']))
            np.testing.assert_allclose(proposal['u'],expected_u,rtol=0,atol=1e-12)
            dphi=float(proposal['correction']-before['correction'])
            accept=bool(np.log(rng.random())<min(0.,-dphi))
            du=proposal['u']-before['u']
            dr=proposal['positions'][t.free]-before['positions'][t.free]
            energy_delta=float(proposal['energy']-before['energy'])
            harmonic_delta=.5*t.reference.thermal_energy*(proposal['u']@proposal['u']-before['u']@before['u'])
            np.testing.assert_allclose((energy_delta-harmonic_delta)/t.reference.thermal_energy,dphi,rtol=0,atol=2e-11)
            rows.append([i,int(accept),dphi,energy_delta,harmonic_delta,
                np.linalg.norm(du),np.linalg.norm(dr),np.sqrt(np.mean(np.sum(dr*dr,axis=1))),
                np.linalg.norm(scores[current]),np.linalg.norm(scores[i]),dphi-scores[current]@du])
            if accept: current=i
        rows=np.asarray(rows);raw.append(rows)
        def segment(part):
            return dict(proposals=len(part),acceptance_fraction=float(part[:,1].mean()),
                mean_delta_correction=float(part[:,2].mean()),median_delta_correction=float(np.median(part[:,2])),
                mean_proposal_cartesian_rms_A=float(part[:,7].mean()),
                mean_current_score_norm=float(part[:,8].mean()),
                mean_absolute_linear_remainder=float(np.abs(part[:,10]).mean()))
        saved=arrays(pilot/f'chain_{chain}'/'chain.npz')
        np.testing.assert_array_equal(rows[16:,1].astype(bool),saved['accepted_at_saved_step'])
        results.append(dict(chain=chain,proposal_dimension=t.basis.shape[1],
            outside_proposals=sum(not bool(item['inside']) for item in data),
            warmup=segment(rows[:16]),retained_first_half=segment(rows[16:40]),
            retained_second_half=segment(rows[40:]),
            source_force_projection_error_eV_A=float(np.max(abs(projected-source['gradient'])))))
    np.savez_compressed(a.output/'proposal_diagnostics.npz',rows=np.asarray(raw))
    summary=dict(completed_utc=datetime.now(timezone.utc).isoformat(),chains=results,
        columns=['proposal','accepted','delta_phi','delta_energy_eV','delta_harmonic_eV',
                 'step_u_norm','step_cartesian_norm_A','step_atom_rms_A','current_score_norm',
                 'proposed_score_norm','nonlinear_score_linear_remainder'],
        new_potential_calls=0,new_MD=0,new_DFT=0,physical_clock=None,
        equilibrium_certified=False,material_approved=False,crack_probability_estimated=False,
        proposal_score_is_not_a_material_force_error=True,
        input_sha256={str(p.relative_to(a.root)).replace('\\','/'):sha(p)
            for p in sorted(pilot.rglob('*')) if p.is_file()},
        source_sha256={relative:sha(a.root/relative) for relative in (
            'solver_v1/silicon_sampling_score.py','solver_v1/silicon_device_ensemble.py',
            'results/silicon_wafer_feasibility/audit_sampling_obstruction.py')})
    (a.output/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(chains=results,new_potential_calls=0),indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path('.'));p.add_argument('--output',type=Path,required=True)
    main(p.parse_args())
