"""N=4 fixed-end TwoRowLJ; an independent algebraic Poisson-Bessel evaluator.

Coordinates are (...,4)=(a1,s1,a2,s2). The six upper pairs are counted once.
Half-integer K is evaluated algebraically, with 16 harmonics on a>=0.70.
No pair cutoff, energy fit, event threshold, clock calibration or closure.
"""
from math import factorial, gamma, pi, sqrt
import numpy as np

A0 = 0.7713438268704838
THETA = 0.02
MOBILITY = np.array([1., .05, 1., .05])
FORCES = (1., 4.5, -4.5, 1.)
DURATIONS = (.1, .3, .3, .1)
CYCLES = 6
SIGMA = .82


def power_jet(p, a, s, modes=16):
    a, s = np.broadcast_arrays(np.asarray(a, float), np.asarray(s, float))
    if np.any(a < .7-1e-14):
        raise ValueError('Algebraic kernel is restricted to a>=0.70')
    exponent = 1-2*p
    value = sqrt(pi)*gamma(p-.5)/gamma(p)*a**exponent
    da = exponent*value/a
    daa = exponent*(exponent-1)*value/a**2
    ds, das, dss = [np.zeros_like(a) for _ in range(3)]
    for m in range(1, modes+1):
        c = 2*pi*m
        phase = c*(.5-s)
        cc, ss = np.cos(phase), np.sin(phase)
        coeff, ca, caa = [np.zeros_like(a) for _ in range(3)]
        common = 2*pi**p*m**(p-1)/gamma(p)*np.exp(-c*a)
        for k in range(p):
            power = p+k
            term = common*factorial(p-1+k)/(factorial(k)*factorial(p-1-k)*(2*c)**k)*a**(-power)
            coeff += term
            ca -= term*(c+power/a)
            caa += term*((c+power/a)**2+power/a**2)
        value += coeff*cc
        da += ca*cc
        ds += coeff*c*ss
        daa += caa*cc
        das += ca*c*ss
        dss -= coeff*c*c*cc
    return value, da, ds, daa, das, dss


def local_jet(a, s):
    p6, p3 = power_jet(6, a, s), power_jet(3, a, s)
    return tuple(4*(SIGMA**12*x-SIGMA**6*y) for x, y in zip(p6, p3))


def local(a, s):
    # Energy-only form of the same half-integer Bessel identity.
    a, s = np.broadcast_arrays(np.asarray(a, float), np.asarray(s, float))
    if np.any(a < .7-1e-14):
        raise ValueError('a outside declared domain')
    z3, z6 = np.zeros_like(a), np.zeros_like(a)
    for m in range(1, 17):
        v = 1/(2*pi*m*a)
        term = np.exp(-2*pi*m*a)*np.cos(2*pi*m*(.5-s))
        z3 += m**2*term*(1+v*(3+3*v))
        z6 += m**5*term*(1+v*(15+v*(105+v*(420+v*(945+945*v)))))
    s3 = 3*pi/8*a**-5+pi**3*a**-3*z3
    s6 = 63*pi/256*a**-11+pi**6/60*a**-6*z6
    return 4*(SIGMA**12*s6-SIGMA**6*s3)


def phi(r):
    z = SIGMA/np.asarray(r)
    return 4*(z**12-z**6)


def energy(q, force):
    q = np.asarray(q, float)
    a, s, aa, ss = np.moveaxis(q, -1, 0)
    e = 2*local(A0, 0.)+local(a, s)+local(aa, ss)+phi(3.)
    e += phi(np.hypot(a-A0, 1+s))+phi(np.hypot(a-A0, 2-s))
    e += phi(np.hypot(aa-A0, 2+ss))+phi(np.hypot(aa-A0, 1-ss))
    e += phi(np.hypot(aa-a, 1+ss-s))
    return e-force*((a-A0)+.2*s+(aa-A0)+.2*ss)


def jet(q, force):
    """Analytic full gradient/Hessian; fixed-coordinate rows are eliminated."""
    q = np.asarray(q, float)
    if q.shape != (4,):
        raise ValueError('jet requires a single four-coordinate point')
    a = np.array([A0, q[0], q[2], A0])
    s = np.array([0., q[1], q[3], 0.])
    g, h = np.zeros(4), np.zeros((4, 4))
    for atom in (1, 2):
        _, da, ds, daa, das, dss = local_jet(a[atom], s[atom])
        ix = slice(2*(atom-1), 2*atom)
        g[ix] = [da-force, ds-.2*force]
        h[ix, ix] = [[daa, das], [das, dss]]
    for i in range(4):
        for j in range(i+1, 4):
            d = np.array([a[j]-a[i], j-i+s[j]-s[i]])
            r = np.linalg.norm(d)
            dp = 4*(-12*SIGMA**12/r**13+6*SIGMA**6/r**7)
            ddp = 4*(156*SIGMA**12/r**14-42*SIGMA**6/r**8)
            block = (ddp-dp/r)*np.outer(d, d)/r**2+dp/r*np.eye(2)
            for atom, sign in ((i, -1), (j, 1)):
                if atom in (1, 2):
                    ix = slice(2*(atom-1), 2*atom)
                    g[ix] += sign*dp*d/r
                    for other, osign in ((i, -1), (j, 1)):
                        if other in (1, 2):
                            jx = slice(2*(other-1), 2*other)
                            h[ix, jx] += sign*osign*block
    return float(energy(q, force)), g, h


def energy_tail_bound(amin=.7):
    m = 17
    c = 2*pi*amin
    v = 1/(c*m)
    t3 = pi**3*amin**-3*m**2*np.exp(-c*m)*(1+v*(3+3*v))/(1-np.exp(-c)*(18/17)**2)
    t6 = pi**6/60*amin**-6*m**5*np.exp(-c*m)*(1+v*(15+v*(105+v*(420+v*(945+945*v)))))/(1-np.exp(-c)*(18/17)**5)
    return 4*(SIGMA**12*t6+SIGMA**6*t3)
