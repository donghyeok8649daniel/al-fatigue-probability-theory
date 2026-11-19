"""Reactive-corridor mesh v2; actual G/event and all physical bounds fixed.

The first graded pilot was too coarse in the post-barrier normal corridor.
Its failed CDF and U diagnostics are preserved. This family allocates fixed
counts to fixed intervals and is frozen before any v2 probability result.
"""
import numpy as np
NORMAL_BREAKS=[.70,.74,.90,1.05,1.20,1.45]
SLIP_BREAKS=[-.08,-.04,.14,.25]
NORMAL_COUNTS=[[3,16,12,16,50],[4,20,15,20,63],[5,25,19,25,79],[6,31,24,31,99]]
SLIP_COUNTS=[[4,36,14],[5,45,18],[6,56,23],[8,70,29]]

def mesh(level,axis,outer=1):
    breaks=NORMAL_BREAKS if axis=='a' else SLIP_BREAKS
    counts=(NORMAL_COUNTS if axis=='a' else SLIP_COUNTS)[level]
    e=np.r_[np.concatenate([np.linspace(l,r,n+1)[:-1] for l,r,n in zip(breaks[:-1],breaks[1:],counts)]),breaks[-1]]
    if axis=='a':e=np.r_[e,np.concatenate([np.linspace(l,r,outer+1)[1:] for l,r in zip([1.45,1.8,2.5,3.5,5.],[1.8,2.5,3.5,5.,6.])])]
    return e

if __name__=='__main__':
    from pathlib import Path
    import json
    here=Path(__file__).parent;z=np.load(here/'support/envelope.npz');rho=z['density'];edges=[z[f'edges{k}'] for k in range(4)]
    rows=[]
    for ia,iss in [(0,0),(1,0),(0,1),(1,1),(2,1),(1,2)]:
        ee=[mesh(ia,'a'),mesh(iss,'s')]*2
        counts=[]
        for k,e in enumerate(ee):
            x=(e[:-1]+e[1:])/2;index=np.clip(np.searchsorted(edges[k],x,side='right')-1,0,len(edges[k])-2)
            counts.append(np.bincount(index,minlength=rho.shape[k]))
        estimates={}
        for threshold in [1e-12,1e-14,1e-16]:
            a=(rho>=threshold).astype(float)
            for k in range(3,-1,-1):a=np.tensordot(a,counts[k],axes=(k,0))
            estimates[str(threshold)]=int(a)
        rows.append(dict(a_level=ia,s_level=iss,shape=[len(e)-1 for e in ee],counts_estimate=estimates))
    report=dict(normal_breaks=NORMAL_BREAKS,normal_counts=NORMAL_COUNTS,slip_breaks=SLIP_BREAKS,slip_counts=SLIP_COUNTS,rows=rows,mask_density_threshold=1e-16,unresolved_acceptance_budget=1e-13,source_observable='support/envelope.npz; diagnostic only, U must be propagated')
    (here/'CORRIDOR_MESH_CONTRACT.json').write_text(json.dumps(report,indent=2));print(report)
