"""Frozen goal-guided numerical mesh; no potential, event or CDF change."""
from pathlib import Path
from functools import lru_cache
import json,hashlib,time,numpy as np
from scipy.integrate import cumulative_trapezoid
from mesh_plan import NORMAL_BREAKS,SLIP_BREAKS,SLIP_COUNTS
HERE=Path(__file__).resolve().parent
BASE=[4,42,40,42,60]
COUNTS=[[int(np.floor(n*1.25**level+.5)) for n in BASE] for level in range(5)]
FLOOR=.25

def core(level,source):
    x=np.array(source['normal_centers']); f=np.array(source['amax_event_flux'])
    dx=np.r_[np.diff(x)[0],(x[2:]-x[:-2])/2,np.diff(x)[-1]]; density=f/dx; parts=[]
    for k,(left,right,n) in enumerate(zip(NORMAL_BREAKS[:-1],NORMAL_BREAKS[1:],COUNTS[level])):
        if k<3: e=np.linspace(left,right,n+1)
        else:
            fine=np.linspace(left,right,4001); monitor=np.interp(fine,x,density,left=0.,right=0.)**(1/3)
            cum=cumulative_trapezoid(monitor,fine,initial=0.); assert cum[-1]>0
            equi=FLOOR*(fine-left)/(right-left)+(1-FLOOR)*cum/cum[-1]
            e=np.interp(np.linspace(0,1,n+1),equi,fine); e[0]=left; e[-1]=right
        assert np.diff(e).min()>0; parts.append(e[:-1])
    return np.r_[np.concatenate(parts),NORMAL_BREAKS[-1]]

@lru_cache(None)
def contract(): return json.loads((HERE/'MESH_V5_CONTRACT.json').read_text())

def mesh(level,axis,outer=1):
    if axis=='a':
        e=np.array(contract()['normal_core_edges'][str(level)])
        return np.r_[e,np.concatenate([np.linspace(l,r,outer+1)[1:] for l,r in zip([1.45,1.8,2.5,3.5,5.],[1.8,2.5,3.5,5.,6.])])]
    return np.r_[np.concatenate([np.linspace(l,r,n+1)[:-1] for l,r,n in zip(SLIP_BREAKS[:-1],SLIP_BREAKS[1:],SLIP_COUNTS[level])]),SLIP_BREAKS[-1]]

def validate():
    c=contract(); assert c['normal_counts']==COUNTS
    assert c['raw_CDF_tolerance']==.01 and c['outer_absolute_tolerance']==1e-12 and c['unresolved_budget']==1e-13
    assert c['planning_source_SHA256']==hashlib.sha256((HERE/'PLANNING_SENSITIVITY.json').read_bytes()).hexdigest()
    assert c['event_flux_source_SHA256']==hashlib.sha256((HERE/'EVENT_MESH_DIAGNOSTIC.json').read_bytes()).hexdigest()
    assert not c['potential_changed'] and not c['event_changed'] and not c['probability_output_correction']
    for k in range(5):
        e=mesh(k,'a'); assert e[0]==.7 and e[-1]==6 and np.diff(e).min()>0 and len(e)-1==sum(COUNTS[k])+5
    return c

if __name__=='__main__':
    path=HERE/'MESH_V5_CONTRACT.json'; assert not path.exists()
    assert not any((p/'progress.json').exists() for p in HERE.glob('r*_aug'))
    planning=json.loads((HERE/'PLANNING_SENSITIVITY.json').read_text())
    assert planning['planning_only'] and not planning['accepted_as_scientific_test']
    source=next(r for r in json.loads((HERE/'EVENT_MESH_DIAGNOSTIC.json').read_text())['rows'] if r['case']=='c1_0_o1_d1e-16_q8_pad1_aug')
    report=dict(mesh_family='goal_allocation_v5',normal_counts=COUNTS,normal_breaks=NORMAL_BREAKS,normal_core_edges={str(k):core(k,source).tolist() for k in range(5)},slip_counts=SLIP_COUNTS,slip_breaks=SLIP_BREAKS,refinement_ratio=1.25,uniform_monitor_fraction=FLOOR,planning_source_SHA256=hashlib.sha256((HERE/'PLANNING_SENSITIVITY.json').read_bytes()).hexdigest(),event_flux_source_SHA256=hashlib.sha256((HERE/'EVENT_MESH_DIAGNOSTIC.json').read_bytes()).hexdigest(),planning_source_not_an_error_certificate=True,counts_basis='Approximate finite-chain integrated goal sensitivity including absorption edges; 1/3 allocation with a fixed positive floor. Only a node-allocation heuristic.',frozen_UTC=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),before_any_v5_CDF=True,potential_changed=False,event_changed=False,probability_output_correction=False,raw_CDF_tolerance=.01,outer_absolute_tolerance=1e-12,unresolved_budget=1e-13)
    path.write_text(json.dumps(report,indent=2)); validate()
    from scipy import ndimage
    with np.load(HERE/'support/envelope_augmented.npz') as z:
        rho=z['density']; ee=[z[f'edges{k}'] for k in range(4)]
    rows=[]
    for a,s in [(0,0),(1,0),(0,1),(1,1),(2,0),(0,2)]:
        edges=[mesh(a,'a'),mesh(s,'s')]*2; counts=[]
        for k,e in enumerate(edges):
            x=(e[:-1]+e[1:])/2; index=np.clip(np.searchsorted(ee[k],x,side='right')-1,0,len(ee[k])-2)
            counts.append(np.bincount(index,minlength=rho.shape[k]))
        sizes={}
        for threshold in [1e-8,1e-9,1e-10]:
            active=ndimage.maximum_filter(rho>=threshold,size=3,mode='constant',cval=0).astype(float)
            for k in range(3,-1,-1): active=np.tensordot(active,counts[k],axes=(k,0))
            sizes[str(threshold)]=int(active)
        rows.append(dict(a_level=a,s_level=s,shape=[len(e)-1 for e in edges],retained_counts=sizes))
    (HERE/'GOAL_MESH_V5_CAPACITY.json').write_text(json.dumps(dict(rows=rows,capacity_only=True,unresolved_budget_not_verified=True),indent=2))
    print(json.dumps(rows,indent=2))
