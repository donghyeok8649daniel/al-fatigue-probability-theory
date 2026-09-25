"""Nine conservative angular residual kernels, for a bounded material test.

For each O-Si-O, Si-O-Si and Si-Si-O central environment, use P0/P1/P2 of
the neighbor angle multiplied by two fixed short-range radial envelopes.
This extends the failed radial-only diagnostic. It is not a production model.
"""
from __future__ import annotations
import numpy as np
from .silicon_oxide_delta_v12 import correction_features,_rise


def angular_features(numbers,positions,cell,pbc):
    from ase import Atoms
    from ase.neighborlist import neighbor_list
    z=np.asarray(numbers,int);r=np.asarray(positions,float)
    atoms=Atoms(numbers=z,positions=r,cell=cell,pbc=pbc)
    i,j,d=neighbor_list('ijD',atoms,2.8,self_interaction=False)
    energy=np.zeros(9);force=np.zeros((len(z),3,9))
    for center in range(len(z)):
        if z[center] not in (8,14):continue
        keep=np.flatnonzero(i==center)
        js=j[keep];vectors=d[keep]
        for a in range(len(js)):
            for b in range(a+1,len(js)):
                za,zb=z[js[a]],z[js[b]]
                if z[center]==14 and za==8 and zb==8:offset=0
                elif z[center]==8 and za==14 and zb==14:offset=3
                elif z[center]==14 and {int(za),int(zb)}=={14,8}:offset=6
                else:continue
                u,v=vectors[a],vectors[b];ru,rv=np.linalg.norm(u),np.linalg.norm(v)
                if ru<=0 or rv<=0:raise ValueError('zero-distance angular neighbor')
                uhat,vhat=u/ru,v/rv;cos=float(uhat@vhat)
                radii=np.array([ru,rv]);transition,transition_prime=_rise(radii,2.2,2.8)
                gaussian=np.exp(-.5*((radii-1.8)/.45)**2)
                weight=gaussian*(1-transition)
                deriv=gaussian*(-transition_prime-(radii-1.8)/.45**2*(1-transition))
                p=np.array([1.,cos,.5*(3*cos*cos-1)])
                dp=np.array([0.,1.,3*cos])
                eu=deriv[0]*weight[1]*uhat[:,None]*p+weight[0]*weight[1]/ru*(vhat-cos*uhat)[:,None]*dp
                ev=weight[0]*deriv[1]*vhat[:,None]*p+weight[0]*weight[1]/rv*(uhat-cos*vhat)[:,None]*dp
                energy[offset:offset+3]+=weight[0]*weight[1]*p
                force[center,:,offset:offset+3]+=eu+ev
                force[js[a],:,offset:offset+3]-=eu
                force[js[b],:,offset:offset+3]-=ev
    return energy,force


def combined_features(numbers,positions,cell,pbc):
    pair_e,pair_f=correction_features(numbers,positions,cell,pbc)
    angle_e,angle_f=angular_features(numbers,positions,cell,pbc)
    return np.r_[pair_e,angle_e],np.concatenate((pair_f,angle_f),axis=2)
