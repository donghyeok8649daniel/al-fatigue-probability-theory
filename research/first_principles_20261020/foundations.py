"""Independent, dimensionless verification tools, not a production energy model.

The half-stack is a declared 2D zigzag lattice, not calibrated aluminium.
All nonzero row reciprocal modes are retained up to an exposed truncation.
The algebraic layer zero mode is summed with Hurwitz zeta, not dropped.
"""
from __future__ import annotations

import math
import numpy as np
from scipy.special import gamma, kv, zeta

EPSILON = 1.0
SIGMA = 0.82
B = 1.0


def row_zero(a, p, b=B):
    a = np.asarray(a, dtype=float)
    value = math.sqrt(math.pi)*gamma(p-.5)/gamma(p)/b*a**(1-2*p)
    zeros = np.zeros_like(value)
    return np.stack((value, (1-2*p)*value/a, zeros,
                     (1-2*p)*(-2*p)*value/a**2, zeros, zeros), axis=-1)


def row_reciprocal(a, s, p, *, modes=64, b=B):
    """S, S_a, S_s, S_aa, S_as, S_ss for the manuscript staggered row."""
    a, s = np.broadcast_arrays(np.asarray(a, float), np.asarray(s, float))
    if np.any(a <= 0) or p <= .5 or b <= 0:
        raise ValueError('Need a,b>0 and p>1/2')
    m = np.arange(1, modes+1, dtype=float)
    wave = 2*np.pi*m/b
    aa = a[..., None]
    phase = wave*(b/2-s[..., None])
    def coefficients(exponent):
        nu = exponent-.5
        return 4*math.sqrt(math.pi)/(b*gamma(exponent))*(wave/(2*aa))**nu*kv(nu, aa*wave)
    f = coefficients(p)
    fp1 = coefficients(p+1)
    fp2 = coefficients(p+2)
    da = -2*p*aa*fp1
    daa = -2*p*fp1+4*p*(p+1)*aa**2*fp2
    c, sn = np.cos(phase), np.sin(phase)
    nz = np.stack((np.sum(f*c, axis=-1), np.sum(da*c, axis=-1),
                   np.sum(f*wave*sn, axis=-1), np.sum(daa*c, axis=-1),
                   np.sum(da*wave*sn, axis=-1), np.sum(-f*wave**2*c, axis=-1)), axis=-1)
    return row_zero(a, p, b)+nz


def row_direct(a, s, p, *, radius=8192, b=B):
    """Elementary real-space reference; no Bessel/Fourier routines."""
    dtype = np.longdouble
    x = (np.arange(-radius, radius+1, dtype=dtype)+dtype('.5'))*dtype(str(b))-dtype(str(s))
    a = dtype(str(a));r2 = a*a+x*x
    t, t1, t2 = r2**(-p), r2**(-p-1), r2**(-p-2)
    return np.asarray((np.sum(t), -2*p*a*np.sum(t1), 2*p*np.sum(x*t1),
        -2*p*np.sum(t1)+4*p*(p+1)*a*a*np.sum(t2),
        -4*p*(p+1)*a*np.sum(x*t2),
        -2*p*np.sum(t1)+4*p*(p+1)*np.sum(x*x*t2)), dtype=float)


def row_tail_value_bound(p, radius, b=B):
    """Uniform s in [-b/2,b/2]; deliberately conservative elementary bound."""
    return 2*b**(-2*p)/(2*p-1)*(radius-1)**(1-2*p)


def lj_row(a, s, *, modes=64, b=B, sigma=SIGMA, epsilon=EPSILON):
    return 4*epsilon*(sigma**12*row_reciprocal(a,s,6,modes=modes,b=b)
                      -sigma**6*row_reciprocal(a,s,3,modes=modes,b=b))


class ZigzagCut:
    """Pair-only cut excess per lateral repeat b, for two rigid half-stacks.

    Reference row l: (jb+l*b/2, l*h). Upper half l>=1 is displaced by
    (-s,delta). Thus layer k uses W_row(k*h+delta,s+(1-k)*b/2), weight k.
    This is a 2D example, not an area-normalized three-dimensional crystal.
    """
    def __init__(self, h, *, layers=16, modes=48, sigma=SIGMA, epsilon=EPSILON, b=B):
        self.h=float(h);self.layers=layers;self.modes=modes
        self.sigma=sigma;self.epsilon=epsilon;self.b=b
        self.k=np.arange(1,layers+1,dtype=float)
        self.base_a=self.k*self.h
        self.base_s=(1-self.k)*self.b/2
        self._base={p:row_reciprocal(self.base_a,self.base_s,p,modes=modes,b=b)
                      -row_zero(self.base_a,p,b) for p in (3,6)}

    def _power_zero_cut(self, delta, p):
        r=2*p-1;x=delta/self.h
        if x <= -1:raise ValueError('First active gap must be positive')
        c=math.sqrt(math.pi)*gamma(p-.5)/gamma(p)/self.b
        value=c*self.h**(-r)*(zeta(r-1,1+x)-x*zeta(r,1+x)-zeta(r-1,1.0))
        da=-r*c*self.h**(-r-1)*(zeta(r,1+x)-x*zeta(r+1,1+x))
        daa=r*(r+1)*c*self.h**(-r-2)*(zeta(r+1,1+x)-x*zeta(r+2,1+x))
        return np.array([value,da,0.,daa,0.,0.])

    def features(self, delta, s):
        values={}
        for p in (3,6):
            d=self.base_a+delta
            nz=row_reciprocal(d,self.base_s+s,p,modes=self.modes,b=self.b)-row_zero(d,p,self.b)
            # Only the fixed reference VALUE is subtracted. Its q derivatives
            # are zero; subtracting reference gradients would be an error.
            nz[:,0]-=self._base[p][:,0]
            values[p]=self._power_zero_cut(delta,p)+np.sum(self.k[:,None]*nz,axis=0)
        return 4*self.epsilon*(self.sigma**12*values[6]-self.sigma**6*values[3])

    def separation_limit(self):
        def power(p):
            r=2*p-1;c=math.sqrt(math.pi)*gamma(p-.5)/gamma(p)/self.b
            return c*self.h**(-r)*zeta(r-1,1.0)+np.sum(self.k*self._base[p][:,0])
        return -4*self.epsilon*(self.sigma**12*power(6)-self.sigma**6*power(3))


def triangular_basis(b=1.0):
    return np.array([[b,b/2],[0,math.sqrt(3)*b/2]])


def plane_reciprocal(d, delta, p, *, shell=12, b=1.0):
    """Independent 2D Poisson power sum and full (d,delta_x,delta_y) derivatives."""
    basis=triangular_basis(b);area=float(np.linalg.det(basis))
    reciprocal=2*np.pi*np.linalg.inv(basis).T
    m,n=np.meshgrid(np.arange(-shell,shell+1),np.arange(-shell,shell+1),indexing='ij')
    indices=np.stack((m.ravel(),n.ravel()),axis=1)
    indices=indices[np.any(indices!=0,axis=1)]
    g=indices@reciprocal.T;length=np.linalg.norm(g,axis=1)
    phase=g@np.asarray(delta,float)
    def hat(exponent):
        return 2*np.pi/gamma(exponent)*(length/(2*d))**(exponent-1)*kv(exponent-1,d*length)
    f=hat(p);f1=hat(p+1);f2=hat(p+2)
    fd=-2*p*d*f1;fdd=-2*p*f1+4*p*(p+1)*d*d*f2
    zero=np.pi/(p-1)*d**(2-2*p)
    grad=np.r_[((2-2*p)*zero/d+np.sum(fd*np.cos(phase))), -np.sum((f*np.sin(phase))[:,None]*g,axis=0)]/area
    hessian=np.zeros((3,3));hessian[0,0]=((2-2*p)*(1-2*p)*zero/d**2+np.sum(fdd*np.cos(phase)))/area
    hessian[0,1:]=hessian[1:,0]=-np.sum((fd*np.sin(phase))[:,None]*g,axis=0)/area
    hessian[1:,1:]=-np.einsum('i,ij,ik->jk',f*np.cos(phase),g,g)/area
    return (zero+np.sum(f*np.cos(phase)))/area,grad,hessian


def plane_direct(d, delta, p, *, shell=256, b=1.0):
    basis=triangular_basis(b)
    m,n=np.meshgrid(np.arange(-shell,shell+1,dtype=np.longdouble),np.arange(-shell,shell+1,dtype=np.longdouble),indexing='ij')
    xy=np.stack((m.ravel(),n.ravel()),axis=1)@basis.astype(np.longdouble).T+np.asarray(delta,np.longdouble)
    rr=np.column_stack((np.full(len(xy),d,dtype=np.longdouble),xy))
    r2=np.sum(rr*rr,axis=1);value=np.sum(r2**(-p))
    grad=-2*p*np.sum(rr*r2[:,None]**(-p-1),axis=0)
    hessian=-2*p*np.sum(r2**(-p-1))*np.eye(3,dtype=np.longdouble)+4*p*(p+1)*np.einsum('i,ij,ik->jk',r2**(-p-2),rr,rr)
    return float(value),np.asarray(grad,float),np.asarray(hessian,float)


def plane_tail_bound(p, shell, delta, b=1.0):
    c=math.sqrt(3)*b/2-np.linalg.norm(delta)/(shell+1)
    if c<=0 or p<=1:raise ValueError('Tail bound domain')
    return 8*c**(-2*p)*shell**(2-2*p)/(2*p-2)


def finite_crossing_counts(lower_rows,upper_rows):
    counts={}
    for low in range(1-lower_rows,1):
        for high in range(1,upper_rows+1):
            counts[high-low]=counts.get(high-low,0)+1
    return counts


def finite_crossing_formula(k,lower_rows,upper_rows):
    return max(0,min(0,upper_rows-k)-max(1-lower_rows,1-k)+1)
