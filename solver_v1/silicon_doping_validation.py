"""Independent mathematical references for the v8 audit, NOT silicon models.

The half-space calculation starts from displacement equilibrium and a Schur
subspace, not the crack solver's compliance quartic or Airy coefficients.
The two-orbital model is exactly soluble synthetic electronic data.
"""
from __future__ import annotations

import itertools
import numpy as np
from scipy.linalg import schur
from scipy.special import expit, logsumexp


def fourier_surface_compliance(c11, c12, c44):
    """(111)[11-2] plane strain, positive Fourier wave number, upper half-space.

    Surface traction = -sigma[:, y]. With u=A exp(ik(x+p*y)),
    (Q+p*(R+R.T)+p*p*T)A=0. Impedance Z=-i*B/A; Y=inv(Z).
    Mode-I energy coefficient H=Re(Y_yy)/2, units GPa^-1.
    This validation helper accepts stable cubic constants only.
    """
    if not np.all(np.isfinite([c11, c12, c44])) or min(c11-c12, c11+2*c12, c44) <= 0:
        raise ValueError('stable finite cubic constants required')
    crystal = np.zeros((3, 3, 3, 3))
    for i, j, k, l in itertools.product(range(3), repeat=4):
        if i == j == k == l:
            crystal[i, j, k, l] = c11
        elif i == j and k == l:
            crystal[i, j, k, l] = c12
        elif (i == k and j == l) or (i == l and j == k):
            crystal[i, j, k, l] = c44
    frame = np.array([[1, 1, -2], [1, 1, 1], [1, -1, 0]], float)
    frame /= np.linalg.norm(frame, axis=1)[:, None]
    rotated = np.einsum('ai,bj,ck,dl,ijkl->abcd', frame, frame, frame, frame, crystal)
    q = rotated[:2, 0, :2, 0]
    r = rotated[:2, 0, :2, 1]
    t = rotated[:2, 1, :2, 1]
    companion = np.block([[np.zeros((2, 2)), np.eye(2)],
                          [-np.linalg.solve(t, q), -np.linalg.solve(t, r+r.T)]])
    triangular, basis, count = schur(companion, output='complex', sort=lambda p: p.imag > 0)
    if count != 2:
        raise ValueError('two decaying modes required')
    a, pa = basis[:2, :2], basis[2:, :2]
    b = r.T@a+t@pa
    y = 1j*np.linalg.solve(b.T, a.T).T
    residual = np.linalg.norm(companion@basis[:, :2]-basis[:, :2]@triangular[:2, :2])
    if (not np.allclose(y, y.conj().T, atol=1e-11, rtol=1e-10)
            or np.linalg.eigvalsh(y).min() <= 0):
        raise ValueError('surface compliance must be Hermitian positive definite')
    return dict(compliance=y, H_GPa_inv=float(y[1, 1].real/2), subspace_residual=float(residual))


def two_orbital_data(q):
    """Two spin-resolved, independently occupied levels; eV and arbitrary q."""
    q = np.asarray(q, float)
    if q.shape != (2,):
        raise ValueError('two coordinates required')
    h = np.array([[[.8, .1], [.1, .4]], [[.3, -.08], [-.08, .7]]])
    linear = np.array([[.4, -.15], [-.25, .3]])
    e = np.array([-.03, .04])+linear@q+.5*np.einsum('i,aij,j->a', q, h, q)
    g = linear+np.einsum('aij,j->ai', h, q)
    return e, g, h


def two_orbital_canonical(q, kt):
    """Canonical N=0,1,2: internal entropy included once in the N=1 sector."""
    e, g, h = two_orbital_data(q)
    log_z = logsumexp(-e/kt)
    w = np.exp(-e/kt-log_z)
    mean = w@g
    d = g-mean
    energies = np.array([0., -kt*log_z, sum(e)])
    gradients = np.array([np.zeros(2), mean, g.sum(axis=0)])
    hessians = np.array([np.zeros((2, 2)),
                         np.einsum('a,aij->ij', w, h)-np.einsum('a,ai,aj->ij', w, d, d)/kt,
                         h.sum(axis=0)])
    return energies, gradients, hessians


def two_orbital_product(q, mu, kt):
    """Independent Fermi product Z=product(1+exp((mu-e_i)/kT))."""
    e, g, h = two_orbital_data(q)
    f = expit((mu-e)/kt)
    susceptibility = f*(1-f)/kt
    hessian = np.einsum('a,aij->ij', f, h)-np.einsum('a,ai,aj->ij', susceptibility, g, g)
    return dict(grand_potential_eV=float(-kt*np.logaddexp(0., (mu-e)/kt).sum()),
                gradient=f@g, hessian=hessian, mean_excess_electrons=float(sum(f)),
                dmean_dmu=float(sum(susceptibility)),
                dmean_dq=-(susceptibility@g), level_occupancy=f)
