"""Minimal scalar-environment modulation of SW angles, with Bessel representation.

E_i = E_pair,i + (rho_i/rho_ref)**nu * E_angle,i.
rho_i=sum_j exp(gamma*sigma/(r_ij-r_cut)) is a model environment, not measured
electron density. rho_ref is fixed at a stated cubic bulk lattice and nu>=0.
No new radial dictionary or free density scale. Ordinary SW is exactly nu=0.
This separate research family has no material, initiation, or kinetic approval.
"""
from __future__ import annotations
import numpy as np
from .silicon_environment_research import (Jet,DiamondCell,_site_energy,_state,_finish,
                                           _value,_dot,_exp,_parameters)
from .silicon_bessel_reference import SW111BesselInterface,site_energy


def _screen_contract(power,reference):
    if (not np.isfinite(power) or power<0 or not np.isfinite(reference) or reference<=0):
        raise ValueError('nonnegative exponent and fixed positive density reference required')


def cubic_density_reference(parameters,lattice):
    """Exact original-support per-site rho in the chosen perfect cubic bulk."""
    model=DiamondCell(parameters,lattice);p=model.p;rc=p['sigma']*p['a'];values=[]
    for sites in model.environments:
        radii=np.array([np.linalg.norm(v) for v,_ in sites]);r=radii[radii<rc]
        values.append(float(np.sum(np.exp(p['gamma']*p['sigma']/(r-rc)))))
    if abs(values[0]-values[1])>1e-12 or values[0]<=0:
        raise ValueError('cubic density reference inconsistent between diamond basis sites')
    return values[0]


def screened_site(vectors,parameters,power,reference):
    """Explicit centered triples and differentiated scalar density product."""
    _screen_contract(power,reference);p=_parameters(parameters)
    original=_site_energy(vectors,p,'direct')
    if not power:return original
    pair=_site_energy(vectors,dict(p,**{'lambda':0.}),'direct')
    rho=0.
    for v in vectors:
        radius=_dot(v,v)**.5
        if _value(radius)<p['sigma']*p['a']:
            rho+=_exp(p['gamma']*p['sigma']/(radius-p['sigma']*p['a']))
    if _value(rho)==0:return pair
    if _value(rho)<0:raise ValueError('negative environment density is unresolved')
    return pair+(rho/reference)**power*(original-pair)


class ScreenedDiamondCell(DiamondCell):
    def __init__(self,parameters,lattice,*,density_power,density_reference,image_shell=2):
        super().__init__(parameters,lattice,image_shell=image_shell)
        _screen_contract(density_power,density_reference)
        self.power,self.reference=float(density_power),float(density_reference)

    def evaluate(self,state,representation='direct',derivatives=True):
        if representation!='direct':raise ValueError('screened diamond explicit representation required')
        q=_state(state,9,derivatives);numeric=np.asarray(state,float)
        f_numeric=np.eye(3)+np.array([[numeric[0],numeric[5]/2,numeric[4]/2],
            [numeric[5]/2,numeric[1],numeric[3]/2],[numeric[4]/2,numeric[3]/2,numeric[2]]])
        if np.linalg.det(f_numeric)<=0:raise ValueError('orientation-preserving cell required')
        lower=(np.linalg.svd(self.cell@f_numeric.T,compute_uv=False)[-1]*(self.shell+1)
               -np.linalg.norm(self.lattice/4*np.ones(3)@f_numeric.T)-np.linalg.norm(numeric[6:]))
        if lower<=self.p['a']*self.p['sigma']:raise ValueError('screened diamond image shell insufficient')
        f=[[1+q[0],q[5]/2,q[4]/2],[q[5]/2,1+q[1],q[3]/2],[q[4]/2,q[3]/2,1+q[2]]]
        result=0.;cutoff=self.p['a']*self.p['sigma']
        for sites in self.environments:
            vectors=[]
            for base,sign in sites:
                if np.linalg.norm(base@f_numeric.T+sign*numeric[6:])>=cutoff:continue
                vectors.append([sum(f[a][b]*base[b] for b in range(3))+sign*q[6+a] for a in range(3)])
            result+=screened_site(vectors,self.p,self.power,self.reference)
        return _finish(result,9,derivatives)


class ScreenedSW111BesselInterface(SW111BesselInterface):
    """Same plane Hankel integrals; apply screening only after all-site moments sum."""
    def __init__(self,parameters,lattice,*,density_power,density_reference,**kwargs):
        super().__init__(parameters,lattice,**kwargs)
        _screen_contract(density_power,density_reference)
        self.power,self.reference=float(density_power),float(density_reference)

    def _screened_moments(self,moments):
        original=site_energy(moments,self.p)
        if not self.power:return original
        pair,rho=moments.jets()[:2];pair=.5*pair
        if rho.value<=0:raise ValueError('nonpositive Bessel environment; resolution unresolved')
        return pair+(rho/self.reference)**self.power*(original-pair)

    def _direct_screened_site(self,center,planes,q,active):
        coordinates=[Jet.variable(v,3,i) for i,v in enumerate(q)];vectors=[]
        for plane in planes:
            sign=int(plane[2]>self.cut)-int(center[2]>self.cut) if active else 0
            dz,delta=plane[2]-center[2],plane[3]-center[3]
            current_d,current_delta=dz+sign*q[0],delta+sign*q[1:]
            for v in self.plane.direct_vectors(current_d,current_delta,exclude_origin=plane[:2]==center[:2]):
                base=v-sign*q[[1,2,0]]
                vectors.append([base[0]+sign*coordinates[1],base[1]+sign*coordinates[2],base[2]+sign*coordinates[0]])
        return screened_site(vectors,self.p,self.power,self.reference)

    def evaluate(self,state,*,method='bessel'):
        q=np.asarray(state,float)
        if q.shape!=(3,) or not np.all(np.isfinite(q)) or q[0]+self.gap<=0:
            raise ValueError('finite three-coordinate state must keep planes ordered')
        if method not in ['bessel','direct']:raise ValueError('method must be bessel or direct')
        planes,cutoff=self._planes(q[0]),self.p['a']*self.p['sigma'];result=Jet.constant(0.,3)
        for center in planes:
            opposite=[p for p in planes if (p[2]>self.cut)!=(center[2]>self.cut)]
            sign=1 if center[2]<self.cut else -1
            if all(abs(p[2]-center[2])>=cutoff and abs(p[2]-center[2]+sign*q[0])>=cutoff for p in opposite):continue
            if method=='bessel':
                now=self._screened_moments(self._center_moments(center,planes,q,True))
                bulk=self._screened_moments(self._center_moments(center,planes,np.zeros(3),False))
            else:
                now=self._direct_screened_site(center,planes,q,True)
                bulk=self._direct_screened_site(center,planes,np.zeros(3),False)
            result+=now-bulk
        return result
