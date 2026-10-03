"""One angular shape parameter, conservative SW and rank-0..4 Bessel moments.

theta=cos(angle)-c0; f(theta)=theta**2*(1-beta*theta)**2.
At tetrahedral theta=0, f, f' and f'' agree with original SW. Bulk curvature
is therefore preserved on the four-neighbor branch, not on arbitrary states.
The beta bounds keep a single angular zero and monotone cost on either side.
This is a separate unapproved research family; no probability/clock gate.
"""
from __future__ import annotations
from functools import lru_cache
from math import factorial
import numpy as np
from scipy.special import jv
from .silicon_environment_research import (Jet,DiamondCell,_site_energy,_state,
                                          _finish,_value,_dot,_exp,_parameters)
from .silicon_bessel_reference import SWPlaneBessel,SW111BesselInterface,PlaneMoments


MONOMIALS=tuple((a,b,m-a-b) for m in range(5) for a in range(m+1) for b in range(m-a+1))
MULTIPLICITIES=np.array([factorial(sum(t))/np.prod([factorial(a) for a in t]) for t in MONOMIALS])
RADIAL_POWERS=tuple(sorted(set((a+b,c) for a,b,c in MONOMIALS)))
RADIAL_INDEX=[RADIAL_POWERS.index((a+b,c)) for a,b,c in MONOMIALS]
MODES=np.arange(-4,5)
PHI=2*np.pi*np.arange(16)/16
ANGULAR_COEFFICIENTS=np.array([np.exp(-1j*np.outer(MODES,PHI))@
                             (np.cos(PHI)**a*np.sin(PHI)**b)/16 for a,b,_ in MONOMIALS])
# These exact finite Fourier coefficients have parity/degree zeros. Set only
# known analytic zeros, not small physical energy/force values or fit weights.
for k,(a,b,_) in enumerate(MONOMIALS):
    forbidden=(abs(MODES)>a+b)|((MODES-(a+b))%2!=0)
    ANGULAR_COEFFICIENTS[k,forbidden]=0.


def beta_contract(beta,c):
    if not np.isfinite(beta) or not np.isfinite(c) or abs(c)>1:
        raise ValueError('finite angular shape and preferred cosine required')
    if max(beta*(-1-c),beta*(1-c))>.5+2e-16:
        raise ValueError('angular shape must retain monotone single-zero cost')


def polynomial_coefficients(beta,c):
    beta_contract(beta,c)
    theta=np.polynomial.Polynomial([-c,1.])
    coefficients=(theta**2*(1-beta*theta)**2).coef
    return np.pad(coefficients,(0,5-len(coefficients)))


def angular_site(vectors,parameters,beta,representation='direct'):
    """Independent centered triples or full Cartesian tensor contractions."""
    p=_parameters(parameters);beta_contract(beta,p['costheta0'])
    if representation not in ['direct','moments']:raise ValueError('direct or moments required')
    if beta==0 and representation=='direct':return _site_energy(vectors,p,'direct')
    pair=0.;neighbors=[];cutoff=p['sigma']*p['a']
    for v in vectors:
        square=_dot(v,v)
        if _value(square)<=1e-20:raise ValueError('coincident atoms in angular model')
        r=square**.5
        if _value(r)>=cutoff:continue
        z=p['sigma']/(r-cutoff)
        pair+=.5*p['A']*p['epsilon']*(p['B']*(p['sigma']/r)**p['p']-(p['sigma']/r)**p['q'])*_exp(z)
        neighbors.append((_exp(p['gamma']*z),[a/r for a in v]))
    if representation=='direct':
        angular=0.
        for j,(w,n) in enumerate(neighbors):
            for weight,m in neighbors[j+1:]:
                theta=_dot(n,m)-p['costheta0']
                angular+=w*weight*theta**2*(1-beta*theta)**2
    else:
        coefficients=polynomial_coefficients(beta,p['costheta0'])
        contraction=0.
        for axes,multiplicity in zip(MONOMIALS,MULTIPLICITIES):
            t=sum(w*np.prod([n[k]**axes[k] for k in range(3)]) for w,n in neighbors)
            contraction+=coefficients[sum(axes)]*multiplicity*t*t
        diagonal=sum(w*w for w,_ in neighbors)
        angular=.5*(contraction-diagonal*np.sum(coefficients))
    return pair+p['lambda']*p['epsilon']*angular


class AngularDiamondCell(DiamondCell):
    def __init__(self,parameters,lattice,*,angle_beta,image_shell=2):
        super().__init__(parameters,lattice,image_shell=image_shell)
        beta_contract(angle_beta,self.p['costheta0']);self.beta=float(angle_beta)

    def evaluate(self,state,representation='direct',derivatives=True):
        q=_state(state,9,derivatives);numeric=np.asarray(state,float)
        fn=np.eye(3)+np.array([[numeric[0],numeric[5]/2,numeric[4]/2],
            [numeric[5]/2,numeric[1],numeric[3]/2],[numeric[4]/2,numeric[3]/2,numeric[2]]])
        if np.linalg.det(fn)<=0:raise ValueError('orientation-preserving cell required')
        lower=(np.linalg.svd(self.cell@fn.T,compute_uv=False)[-1]*(self.shell+1)
               -np.linalg.norm(self.lattice/4*np.ones(3)@fn.T)-np.linalg.norm(numeric[6:]))
        if lower<=self.p['a']*self.p['sigma']:raise ValueError('angular diamond image shell insufficient')
        f=[[1+q[0],q[5]/2,q[4]/2],[q[5]/2,1+q[1],q[3]/2],[q[4]/2,q[3]/2,1+q[2]]]
        result=0.;cutoff=self.p['sigma']*self.p['a']
        for sites in self.environments:
            vectors=[]
            for base,sign in sites:
                if np.linalg.norm(base@fn.T+sign*numeric[6:])>=cutoff:continue
                vectors.append([sum(f[a][b]*base[b] for b in range(3))+sign*q[6+a] for a in range(3)])
            result+=angular_site(vectors,self.p,self.beta,representation)
        return _finish(result,9,derivatives)


def zero_moments():
    return PlaneMoments(np.zeros(37),np.zeros((37,3)),np.zeros((37,3,3)))


def combine_moments(parts):
    result=zero_moments()
    for m in parts:
        result=PlaneMoments(result.value+m.value,result.gradient+m.gradient,result.hessian+m.hessian)
    return result


def moment_energy(moments,parameters,beta):
    p=_parameters(parameters);coeff=polynomial_coefficients(beta,p['costheta0'])
    pair,diagonal,*t=moments.jets()
    contraction=sum(coeff[sum(axes)]*mult*v*v for axes,mult,v in zip(MONOMIALS,MULTIPLICITIES,t))
    return .5*pair+.5*p['lambda']*p['epsilon']*(contraction-diagonal*np.sum(coeff))


class AngularSWPlaneBessel(SWPlaneBessel):
    """J0..J4 Hankel integrals for fixed symmetric Cartesian moments.

37 statistics are pair, diagonal, and 35 rank-0..4 components. They are not
independent material parameters. Tensor multiplicities restore full norms.
The complete G/-G disk is retained; imaginary residuals fail explicitly.
"""
    @lru_cache(maxsize=64)
    def _coefficients(self,distance):
        d=Jet.variable(float(distance),1,0);p=self.p;rc=p['a']*p['sigma']
        maximum=np.sqrt(rc**2-distance**2);rho=maximum*(self.abscissa+1)/2
        weights=np.pi*maximum*self.weights*rho
        channels=np.empty((17,3,self.nodes))
        for k,rxy in enumerate(rho):
            r=(d*d+rxy*rxy)**.5;z=p['sigma']/(r-rc);w=(p['gamma']*z).exp()
            pair=p['A']*p['epsilon']*(p['B']*(p['sigma']/r)**p['p']-(p['sigma']/r)**p['q'])*z.exp()
            t,u=rxy/r,d/r
            values=[pair,w*w]+[w*t**a*u**b for a,b in RADIAL_POWERS]
            for j,value in enumerate(values):channels[j,:,k]=value.value,value.gradient[0],value.hessian[0,0]
        output=np.empty((len(self.waves),3,17,5))
        for order in range(5):
            output[:,:,:,order]=np.einsum('gn,cdn,n->gdc',jv(order,self.waves[:,None]*rho),channels,weights)
        return output

    def evaluate(self,distance,delta):
        d,delta=float(distance),np.asarray(delta,float)
        if not np.isfinite(d) or d==0 or delta.shape!=(2,) or not np.all(np.isfinite(delta)):
            raise ValueError('nonzero finite d and two-dimensional delta required')
        if abs(d)>=self.p['a']*self.p['sigma']:return zero_moments()
        return self._evaluate_cached(d,float(delta[0]),float(delta[1]))

    @lru_cache(maxsize=128)
    def _evaluate_cached(self,d,dx,dy):
        coefficients=self._coefficients(d);groups=[(np.zeros((1,2)),coefficients[0])]
        groups.extend((s.vectors,c) for s,c in zip(self.shells,coefficients[1:]))
        total=zero_moments()
        for vectors,c in groups:
            psi=np.arctan2(vectors[:,1],vectors[:,0]);f=np.zeros((3,len(vectors),37),complex)
            f[:,:,:2]=c[:,None,:2,0]
            for k,(axes,radial) in enumerate(zip(MONOMIALS,RADIAL_INDEX)):
                a,b,_=axes
                for m,angular in zip(MODES,ANGULAR_COEFFICIENTS[k]):
                    if abs(m)>a+b or (m-a-b)%2:continue
                    signed=(-1)**abs(int(m)) if m<0 else 1
                    f[:,:,k+2]+=c[:,None,radial+2,abs(m)]*angular*((-1j)**int(m))*signed*np.exp(1j*m*psi)
            phase=np.exp(1j*(vectors@np.array([dx,dy])))/self.geometry.atomic_cell_area
            value=np.einsum('vk,v->k',f[0],phase);gradient=np.zeros((37,3),complex);hessian=np.zeros((37,3,3),complex)
            gradient[:,0]=np.einsum('vk,v->k',f[1],phase)
            hessian[:,0,0]=np.einsum('vk,v->k',f[2],phase)
            for a in range(2):
                gradient[:,a+1]=np.einsum('vk,v->k',f[0],phase*1j*vectors[:,a])
                hessian[:,0,a+1]=np.einsum('vk,v->k',f[1],phase*1j*vectors[:,a])
                hessian[:,a+1,0]=hessian[:,0,a+1]
                for b in range(2):hessian[:,a+1,b+1]=np.einsum('vk,v->k',f[0],-phase*vectors[:,a]*vectors[:,b])
            scale=max(1.,np.max(abs(value)),np.max(abs(gradient)),np.max(abs(hessian)))
            if max(np.max(abs(value.imag)),np.max(abs(gradient.imag)),np.max(abs(hessian.imag)))>1e-11*scale:
                raise FloatingPointError('rank-four reciprocal symmetry left nonreal moment')
            total=combine_moments([total,PlaneMoments(value.real,gradient.real,hessian.real)])
        for array in (total.value,total.gradient,total.hessian):array.setflags(write=False)
        return total

    def intraplane(self):
        values=np.zeros(37);p=self.p;rc=p['a']*p['sigma']
        for v in self.direct_vectors(0.,[0.,0.],exclude_origin=True):
            r=np.linalg.norm(v);n=v/r;z=p['sigma']/(r-rc);w=np.exp(p['gamma']*z)
            values[0]+=p['A']*p['epsilon']*(p['B']*(p['sigma']/r)**p['p']-(p['sigma']/r)**p['q'])*np.exp(z)
            values[1]+=w*w
            values[2:]+=[w*np.prod(n**np.array(axes)) for axes in MONOMIALS]
        return PlaneMoments(values,np.zeros((37,3)),np.zeros((37,3,3)))


class AngularSW111BesselInterface(SW111BesselInterface):
    def __init__(self,parameters,lattice,*,angle_beta,cut_kind='shuffle',shell_index=96,nodes=512):
        super().__init__(parameters,lattice,cut_kind=cut_kind,shell_index=shell_index,nodes=nodes)
        beta_contract(angle_beta,self.p['costheta0']);self.beta=float(angle_beta)
        self.plane=AngularSWPlaneBessel(parameters,lattice,shell_index=shell_index,nodes=nodes)

    def _center_moments(self,center,planes,q,active):
        parts=[];side=center[2]>self.cut
        for plane in planes:
            sign=int(plane[2]>self.cut)-int(side) if active else 0
            if plane[:2]==center[:2]:parts.append(self.plane.intraplane())
            else:parts.append(self.plane.evaluate(plane[2]-center[2]+sign*q[0],plane[3]-center[3]+sign*q[1:]).mapped(sign))
        return combine_moments(parts)

    def _direct_site(self,center,planes,q,active):
        coordinates=[Jet.variable(v,3,i) for i,v in enumerate(q)];vectors=[]
        for plane in planes:
            sign=int(plane[2]>self.cut)-int(center[2]>self.cut) if active else 0
            dz,delta=plane[2]-center[2],plane[3]-center[3]
            for v in self.plane.direct_vectors(dz+sign*q[0],delta+sign*q[1:],exclude_origin=plane[:2]==center[:2]):
                base=v-sign*q[[1,2,0]]
                vectors.append([base[0]+sign*coordinates[1],base[1]+sign*coordinates[2],base[2]+sign*coordinates[0]])
        return angular_site(vectors,self.p,self.beta,'direct')

    def evaluate(self,state,*,method='bessel'):
        q=np.asarray(state,float)
        if q.shape!=(3,) or not np.all(np.isfinite(q)) or q[0]+self.gap<=0:
            raise ValueError('finite three-coordinate state must keep planes ordered')
        if method not in ['bessel','direct']:raise ValueError('method must be bessel or direct')
        planes,rc=self._planes(q[0]),self.p['a']*self.p['sigma'];result=Jet.constant(0.,3)
        for center in planes:
            sign=1 if center[2]<self.cut else -1
            opposite=[p for p in planes if (p[2]>self.cut)!=(center[2]>self.cut)]
            if all(abs(p[2]-center[2])>=rc and abs(p[2]-center[2]+sign*q[0])>=rc for p in opposite):continue
            if method=='bessel':
                now=moment_energy(self._center_moments(center,planes,q,True),self.p,self.beta)
                bulk=moment_energy(self._center_moments(center,planes,np.zeros(3),False),self.p,self.beta)
            else:
                now=self._direct_site(center,planes,q,True);bulk=self._direct_site(center,planes,np.zeros(3),False)
            result+=now-bulk
        return result
