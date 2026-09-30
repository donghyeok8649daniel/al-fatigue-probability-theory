"""Independent arithmetic/hash replay of v15; never calls an atomistic model."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.constants import Boltzmann, electron_volt
from scipy.linalg import eigh

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def main(root, cad=None):
    r=root/'results/silicon_device_ensemble_v15'; count=0
    for path,digest in json.loads((r/'device_audit/input_manifest.json').read_text())['sha256'].items():
        assert sha(root/path)==digest, path
        count+=1
    audit=json.loads((r/'device_audit/device_sensitivity.json').read_text())
    with np.load(root/'results/silicon_initiation_v13/dense_return_zero/raw_hessian.npz') as z:
        h=z['hessian']; b=z['basis']
    v=np.zeros(len(h)); v[-1]=1/6
    values,vectors=eigh(h)
    compliance=float(np.sum((vectors.T@v)**2/values)); count+=1
    assert abs(audit['source_schur_stiffness_eV_A2']-1/compliance)<1e-10
    # Recompute physical extension directly from the rigid-grip entries.
    grip=b[:,-1].reshape(-1,3)
    assert abs((grip[:,2].max()-grip[:,2].min())-v[-1])<1e-14;count+=1
    kb=Boltzmann/electron_volt
    for row in audit['rows']:
        ratio=row['device_to_specimen_stiffness_ratio']
        var=0. if row['infinite_stiffness'] else kb*row['temperature_K']*compliance/(1+ratio)
        gain=1. if row['infinite_stiffness'] else ratio/(1+ratio)
        assert abs(var-row['extension_variance_A2'])<1e-12
        assert abs(gain-row['d_extension_d_actuator'])<1e-12
        count+=2
    facts=json.loads((r/'device_audit/cad_facts.json').read_text(encoding='utf-8'))
    g=facts['specimen_geometry'];area=np.pi*g['gauge_diameter']**2/4
    assert abs(area-facts['area_mm2'])<1e-12
    assert abs(100*area-facts['required_force_at_100MPa_N'])<1e-9
    assert facts['maximum_force_N'] is None and facts['measured_machine_stiffness_N_mm'] is None
    count+=3
    if cad:
        assert sha(cad)==facts['sha256'];count+=1
    protocol=json.loads((r/'nonlinear_pilot/protocol.json').read_text())
    for path,digest in protocol['source_sha256'].items():
        assert sha(root/path)==digest;count+=1
    for key,path in [('runner_sha256','results/silicon_wafer_feasibility/sample_nonlinear_ensemble_v15.py'),
                     ('target_module_sha256','solver_v1/silicon_device_ensemble.py'),
                     ('sampler_module_sha256','solver_v1/silicon_thermal_research.py')]:
        assert sha(root/path)==protocol[key];count+=1
    allocation=json.loads((r/'nonlinear_pilot/allocation_check.json').read_text())
    sample=json.loads((r/'nonlinear_pilot/summary.json').read_text())
    assert allocation['available_bytes']<allocation['required_bytes'] and not allocation['passed']
    assert sample['phase']=='preallocation' and sample['new_model_calls']==0
    assert sample['evaluations_recorded']==0 and not sample['complete']
    assert not list((r/'nonlinear_pilot').glob('chain_*'))
    count+=4
    print(json.dumps(dict(complete=True,checks=count,new_potential_evaluations=0,
        allocation_block_preserved=True,scope='stored results, source integrity and blocked attempt; not canonical sampling validation'),indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path('.'));p.add_argument('--cad',type=Path)
    a=p.parse_args();main(a.root,a.cad)
