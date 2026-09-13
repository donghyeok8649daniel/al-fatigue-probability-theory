"""STF rank-four infinite environment, separate research hypothesis.

Q4 is NONZERO in cubic FCC. Never drop its bulk background or apply the
odd-moment zero-background formula. No existing model or pair term is altered.
"""
from functools import lru_cache
import itertools
import math
import numpy as np
from scipy.special import gamma, gammaincc
from .angular_environment_reference import _radial
from .fcc111_geometry import fcc111_geometry_from_b
from .fcc111_lattice_sum import triangular_reciprocal_shells
from .tail_constrained_material import normalized_amplitude
from .vector_interface_reference import _square_jet

INDICES=tuple(itertools.product(range(3),repeat=4))


def stf_fourth(raw):
    raw=np.asarray(raw)
    if raw.shape[-4:]!=(3,3,3,3):
        raise ValueError('four Cartesian tensor indices required')
    prefix=tuple(range(raw.ndim-4));offset=len(prefix)
    value=sum(raw.transpose(prefix+tuple(offset+i for i in p))
              for p in itertools.permutations(range(4)))/24
    trace=np.einsum('...iikl->...kl',value);twice=np.trace(trace,axis1=-2,axis2=-1)
    out=value.copy()
    for i,j,k,l in INDICES:
        out[...,i,j,k,l]-=((i==j)*trace[...,k,l]+(i==k)*trace[...,j,l]
            +(i==l)*trace[...,j,k]+(j==k)*trace[...,i,l]
            +(j==l)*trace[...,i,k]+(k==l)*trace[...,i,j])/7
        out[...,i,j,k,l]+=((i==j)*(k==l)+(i==k)*(j==l)+(i==l)*(j==k))*twice/35
    return out


@lru_cache(maxsize=1)
def fourth_basis():
    P=stf_fourth(np.eye(81).reshape(81,3,3,3,3)).reshape(81,81)
    values,vectors=np.linalg.eigh(P)
    if np.count_nonzero(values>.5)!=9 or np.max(abs(P@P-P))>1e-13:
        raise ArithmeticError('STF rank-four projector failed')
    return vectors[:,values>.5]


@lru_cache(maxsize=32)
def xy_derivatives(indices):
    # D_G(P(G) F^(t)(|G|²)) = P_G F^(t) + 2G P F^(t+1).
    terms={(0,0,0):1.}
    for direction in indices:
        new={}
        for (px,py,t),c in terms.items():
            p=[px,py]
            if p[direction]:
                v=p.copy();v[direction]-=1;key=(*v,t)
                new[key]=new.get(key,0.)+c*p[direction]
            p[direction]+=1;key=(*p,t+1)
            new[key]=new.get(key,0.)+2*c
        terms=new
    return tuple((*key,c) for key,c in terms.items())


def fourth_fourier(d,vectors,k,normal_order=0):
    vectors=np.asarray(vectors);Q=np.sqrt(k*k+np.sum(vectors*vectors,axis=1))
    radial={(t,z):_radial(d,Q,k,t,normal_order,z) for z in range(5) for t in range(5-z)}
    values=[]
    for indices in INDICES:
        xy=tuple(i for i in indices if i!=2);z=4-len(xy)
        values.append((1j)**len(xy)*sum(c*vectors[:,0]**px*vectors[:,1]**py*radial[t,z]
            for px,py,t,c in xy_derivatives(xy)))
    return np.asarray(values).T@fourth_basis()


class RankFourEnvironment:
    def __init__(self,decay,*,tolerance=2e-12):
        if not np.isfinite(decay+tolerance) or min(decay,tolerance)<=0:
            raise ValueError('positive finite decay/tolerance required')
        self.k=float(decay);self.C=normalized_amplitude(decay)
        self.geometry=fcc111_geometry_from_b(1.);self.h=self.geometry.h111
        self.tolerance=tolerance

    @lru_cache(maxsize=256)
    def plane(self,d,delta,k=None):
        k=self.k if k is None else k
        if d<=0 or k<=0 or not np.all(np.isfinite([d,k,*delta])):
            raise ValueError('positive separation/decay and finite registry required')
        result=np.zeros((10,9));small=0
        shells=[np.zeros((1,2))]+[s.vectors for s in triangular_reciprocal_shells(self.geometry,32)]
        for count,vectors in enumerate(shells):
            gx,gy=vectors.T
            c0,c1,c2=[fourth_fourier(d,vectors,k,n) for n in (0,1,2)]
            cs=np.stack([c0,c1,1j*gx[:,None]*c0,1j*gy[:,None]*c0,c2,
                1j*gx[:,None]*c1,1j*gy[:,None]*c1,-gx[:,None]**2*c0,
                -(gx*gy)[:,None]*c0,-gy[:,None]**2*c0],axis=1)
            cs*=self.C/self.geometry.atomic_cell_area
            result+=np.real(np.einsum('g,gic->ic',np.exp(1j*(vectors@delta)),cs))
            last=float(np.max(np.sum(abs(cs),axis=0)))
            small=small+1 if count and last<self.tolerance*.1 else 0
            if small>=3:
                return result,dict(shells=count,last_shell_envelope=last)
        raise ArithmeticError('rank-four reciprocal convergence not reached')

    def same_plane(self,k):
        g=self.geometry;cover=g.b/np.sqrt(3)
        for radius in range(4,100,2):
            low=radius-2*cover
            if k*low<=4:
                continue
            tail=2*np.pi*self.C/g.atomic_cell_area*(cover*gamma(5)*gammaincc(5,k*low)/k**5
                                                       +gamma(6)*gammaincc(6,k*low)/k**6)
            if tail<self.tolerance*.1:
                extent=int(np.ceil(2*radius/g.b))+1
                m,n=np.meshgrid(np.arange(-extent,extent+1),np.arange(-extent,extent+1),indexing='ij')
                r=g.b*np.sqrt(m*m+m*n+n*n);keep=(r>0)&(r<=radius)
                S4=float(np.sum(self.C*r[keep]**4*np.exp(-k*r[keep])))
                # Sixfold in-plane symmetry fixes degree-four angular averages.
                raw=np.zeros(81)
                for j,indices in enumerate(INDICES):
                    if 2 not in indices:
                        nx=indices.count(0)
                        raw[j]=S4*({0:3/8,2:1/8,4:3/8}.get(nx,0.))
                return raw@fourth_basis(),tail
        raise ArithmeticError('rank-four same-plane tail not resolved')

    @lru_cache(maxsize=24)
    def bulk(self,stretch=1.):
        if not np.isfinite(stretch) or stretch<=0:
            raise ValueError('positive stretch required')
        k=self.k*stretch;same,tail=self.same_plane(k)
        # value, normal-strain first/second, engineering-shear first/second.
        out=np.zeros((5,9));out[0]=same;small=0
        for layer in range(1,81):
            d=layer*self.h
            jet,_=self.plane(d,tuple(self.geometry.abc_shift(layer)),k)
            term=2*np.array([jet[0],d*jet[1],d*d*jet[4],d*jet[2],d*d*jet[7]])
            out+=term;last=float(np.max(abs(term)))
            small=small+1 if last<self.tolerance*.1 else 0
            if small>=3:
                return out*stretch**4,dict(layers=layer,last_layer=last,same_plane_tail=tail)
        raise ArithmeticError('rank-four bulk layer convergence not reached')

    def bulk_column(self,step=4e-4):
        Q,_=self.bulk();v,n,nn,s,ss=Q
        samples={i:np.sum(self.bulk(1+i*step)[0][0]**2) for i in (-2,-1,1,2)}
        hydro=(-samples[2]+16*samples[1]-30*(v@v)+16*samples[-1]-samples[-2])/(12*step**2)
        return np.array([2*v@n,-v@v,hydro,2*(n@n+v@nn),2*(s@s+v@ss)])

    @lru_cache(maxsize=192)
    def interface_jet(self,state):
        a,x,y=map(float,state)
        if a<=0 or not np.all(np.isfinite(state)):
            raise ValueError('finite state and positive opening required')
        changes=[];small=0;bulk=self.bulk()[0][0]
        for layer in range(1,81):
            delta=self.geometry.abc_shift(layer)
            old,_=self.plane(layer*self.h,tuple(delta))
            active,_=self.plane(a+(layer-1)*self.h,tuple(delta+[x,y]))
            active=active.copy();active[0]-=old[0]
            changes.append(active);last=float(np.max(abs(active)))
            small=small+1 if last<self.tolerance*.1 else 0
            if small>=3:
                depths=np.cumsum(np.asarray(changes)[::-1],axis=0)[::-1]
                full=depths.copy();full[:,0]+=bulk
                jet=2*_square_jet(full)
                jet[0]=2*np.sum(2*depths[:,0]@bulk+np.sum(depths[:,0]**2,axis=1))
                return jet,dict(layers=layer,last_layer=last,bulk_moment_norm=np.linalg.norm(bulk))
        raise ArithmeticError('rank-four interface convergence not reached')


@lru_cache(maxsize=1)
def polynomial_weights():
    powers=sorted({tuple(indices.count(i) for i in range(3)) for indices in INDICES})
    weights=np.zeros((len(powers),9));basis=fourth_basis()
    for j,indices in enumerate(INDICES):
        weights[powers.index(tuple(indices.count(i) for i in range(3)))]+=basis[j]
    return tuple(powers),weights


def neighbor_fourth_jets(R,k,C):
    """Analytic value/Jacobian/Hessian of STF(R^4) C exp(-k|R|)."""
    R=np.asarray(R,float)
    if R.ndim!=2 or R.shape[1]!=3 or np.any(~np.isfinite(R)) or not np.isfinite(k+C) or min(k,C)<=0:
        raise ValueError('nonzero finite three-dimensional neighbor vectors required')
    r=np.linalg.norm(R,axis=1)
    if np.any(r<=0):raise ValueError('self excluded from neighbor jets')
    powers,weights=polynomial_weights()
    P=np.zeros((len(R),9));dP=np.zeros((len(R),3,9));HP=np.zeros((len(R),3,3,9))
    def monomial(p):return np.prod(R**np.asarray(p)[None,:],axis=1)
    for alpha,weight in zip(powers,weights):
        P+=monomial(alpha)[:,None]*weight
        for i in range(3):
            if not alpha[i]:continue
            p=list(alpha);p[i]-=1
            dP[:,i]+=alpha[i]*monomial(p)[:,None]*weight
            for j in range(3):
                if not p[j]:continue
                t=p.copy();t[j]-=1
                HP[:,i,j]+=alpha[i]*p[j]*monomial(t)[:,None]*weight
    f=C*np.exp(-k*r);e=R/r[:,None];df=-k*f[:,None]*e
    HF=f[:,None,None]*((k*k+k/r)[:,None,None]*e[:,:,None]*e[:,None,:]
                       -(k/r)[:,None,None]*np.eye(3))
    value=f[:,None]*P
    first=f[:,None,None]*dP+df[:,:,None]*P[:,None,:]
    second=(f[:,None,None,None]*HP+HF[:,:,:,None]*P[:,None,None,:]
            +df[:,:,None,None]*dP[:,None,:,:]+df[:,None,:,None]*dP[:,:,None,:])
    return value,first,second


@lru_cache(maxsize=8)
def _bulk_neighbor_data(model,bulk):
    return neighbor_fourth_jets(bulk.R,model.k,model.C)


def rank_four_bulk_column(model,bulk,q_cubic):
    """Controlled direct Bloch reference, including NONZERO bulk Q4.

The local/interface energy remains Poisson-summed. This direct radius is
accompanied by an absolute infinite-lattice tail, not silently canonicalized.
"""
    if bulk.stretch!=1:
        raise ValueError('pristine cubic reference required')
    value,first,second=_bulk_neighbor_data(model,bulk)
    q=bulk.wavevector(q_cubic);qr=bulk.R@q;qn=np.linalg.norm(q)
    Q0=value.sum(axis=0)
    L=np.einsum('n,nic->ic',np.sin(qr),first)
    K=np.einsum('n,nijc->ijc',2*np.sin(qr/2)**2,second)
    H=2*L@L.T+4*np.einsum('c,ijc->ij',Q0,K)
    tail=bulk.tail;k=model.k;C=model.C
    def Jtail(p):return C*(4*np.sqrt(3)*tail.exponential(k,3+p)+k*tail.exponential(k,4+p))
    def Htail(p):return C*(36*tail.exponential(k,2+p)+(8*np.sqrt(3)+2)*k*tail.exponential(k,3+p)
                              +k*k*tail.exponential(k,4+p))
    dQ=C*tail.exponential(k,4)
    dL=min(Jtail(0),qn*Jtail(1));dK=min(2*Htail(0),.5*qn*qn*Htail(2))
    bound=2*(2*np.linalg.norm(L)*dL+dL*dL)+4*(np.linalg.norm(Q0)*dK+np.linalg.norm(K)*dQ+dQ*dK)
    return H,float(bound)
