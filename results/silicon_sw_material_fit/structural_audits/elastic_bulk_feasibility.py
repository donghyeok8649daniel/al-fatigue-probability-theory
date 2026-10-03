"""Exact two-variable quadratic minimum on the linear elastic bands.

The feasible domain relaxes away the C44 constraint. A certified bulk error
above its guard therefore establishes conflict for the more restricted model.
This says nothing about other SW shapes, longer range, or general Bessel families.
"""
from pathlib import Path
from itertools import combinations
import sys,json
import numpy as np
cache=Path(__file__).resolve().parent;root=cache.parents[2]/'aft-silicon-wafer';sys.path.insert(0,str(root))
from solver_v1.run_silicon_sw_guarded_fit import design
from results.silicon_wafer_feasibility.run_dft_material_audit import read_frames
frames,_=read_frames(root/'.cache/si-atomistic-v5/Si_PRX_GAP.zip')
fit=json.loads((cache/'results/fit.json').read_text());old=np.load(cache/'results/features.npz')
_,_,b,t=design(old['energy_basis'],old['force_basis'],fit['training_frames'],frames,
              old['offsets'],old['reference_energy'],old['reference_force'],fit['anchor_frame'])
b=np.c_[b[:,0]+b[:,1],b[:,2]];bound=float(np.sum((b@np.ones(2)-t)**2))
base=json.loads((cache/'results/bulk_diagnostic.json').read_text())['original']
c11,c12=base['C11_GPa'],base['C12_GPa'];pair=(c11+2*c12)/3;angle=(c11-c12)/3
c=np.array([[pair,2*angle],[pair,-angle]]);target=np.array([153.3,56.3])
a=np.vstack((c,-c,np.eye(2)));lower=np.r_[.9*target,-1.1*target,[0.,0.]]
gram=b.T@b;linear=b.T@t;records=[]
for number in range(3):
    for active in combinations(range(6),number):
        aa=a[list(active)]
        if number:
            kkt=np.block([[gram,-aa.T],[aa,np.zeros((number,number))]])
            try:solution=np.linalg.solve(kkt,np.r_[linear,lower[list(active)]])
            except np.linalg.LinAlgError:continue
            x,multiplier=solution[:2],solution[2:]
        else:x=np.linalg.solve(gram,linear);multiplier=np.array([])
        slack=a@x-lower
        if min(slack)<-1e-8 or (number and min(multiplier)<-1e-8):continue
        gradient=gram@x-linear-aa.T@multiplier
        value=float(np.sum((b@x-t)**2))
        # Half-objective Lagrange dual is quadratic; compute its global infimum.
        q=linear+aa.T@multiplier
        dual=float(t@t-q@np.linalg.solve(gram,q)+2*multiplier@lower[list(active)])
        records.append(dict(amplitudes=x.tolist(),active=list(active),multipliers=multiplier.tolist(),
             slack=slack.tolist(),bulk_squared=value,dual_lower_bound=dual,duality_gap=value-dual,
             kkt_max_abs=float(np.max(abs(gradient))),C11_C12_GPa=(c@x).tolist()))
assert records
best=min(records,key=lambda r:r['bulk_squared'])
result=dict(complete=True,scope='fixed original shapes, equilibrium, common pair amplitude; linear C11/C12 bands only',
             relative_elastic_band=.1,maximum_bulk_squared=bound,
             optimum=best,certified_incompatible=best['dual_lower_bound']>bound,
             C44_constraint_omitted=True,material_family_impossibility_claim=False)
(cache/'elastic_guard_results/linear_band_bulk_feasibility.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
