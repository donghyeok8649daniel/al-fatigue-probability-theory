"""Minimal conservative Si-O radial residual candidate; never production.

Eight fixed smooth kernels test whether a purely radial residual can explain
the measured material error. No MACE weights, chemistry, charge or kinetics are
changed. A good force fit alone is not material or initiation validation.
"""
from __future__ import annotations
import numpy as np

CENTERS_A=np.linspace(1.4,3.0,8)
WIDTH_A=.24
SUPPORT_A=3.5


def _rise(r,lo,hi):
    t=np.clip((r-lo)/(hi-lo),0.,1.)
    value=t**3*(10-15*t+6*t*t)
    deriv=30*t*t*(1-t)**2/(hi-lo)
    deriv=np.where((r>lo)&(r<hi),deriv,0.)
    return value,deriv


def radial_kernels(distances):
    r=np.asarray(distances,float)
    if np.any(~np.isfinite(r)) or np.any(r<=0):raise ValueError('positive finite pair distances required')
    left,leftp=_rise(r,.6,1.0);right,rightp=_rise(r,2.8,SUPPORT_A)
    window=left*(1-right);windowp=leftp*(1-right)-left*rightp
    d=(r[:,None]-CENTERS_A[None,:])/WIDTH_A
    gaussian=np.exp(-.5*d*d)
    return gaussian*window[:,None],gaussian*(windowp[:,None]-d/WIDTH_A*window[:,None])


def correction_features(numbers,positions,cell,pbc):
    """Energy coefficients and physical FORCE coefficients with periodic images.

    The directed neighbor list is restricted to Si->O, counting each heteropair
    per periodic cell once. No factor 1/2 applies after selecting this direction.
    """
    from ase import Atoms
    from ase.neighborlist import neighbor_list
    numbers=np.asarray(numbers,int);r=np.asarray(positions,float)
    if r.shape!=(len(numbers),3):raise ValueError('matching numbers/positions required')
    atoms=Atoms(numbers=numbers,positions=r,cell=cell,pbc=pbc)
    i,j,delta=neighbor_list('ijD',atoms,SUPPORT_A,self_interaction=False)
    keep=(numbers[i]==14)&(numbers[j]==8)
    i=i[keep];j=j[keep];delta=delta[keep]
    energy=np.zeros(8);force=np.zeros((len(r),3,8))
    if len(i):
        distance=np.linalg.norm(delta,axis=1)
        value,derivative=radial_kernels(distance)
        energy=value.sum(axis=0)
        vector=delta[:,:,None]/distance[:,None,None]*derivative[:,None,:]
        np.add.at(force,i,vector);np.add.at(force,j,-vector)
    return energy,force
