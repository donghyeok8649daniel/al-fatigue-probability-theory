"""Static, long-wavelength interface elasticity derived from bulk C_ijkl.

Research only. C must come from the selected infinite-lattice energy, not an
extra fit. This is the continuum limit, NOT an atomically exact finite-q
kernel, a probability generator, or a physical kinetic calibration.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.linalg import schur


def cubic_elastic_tensor(c11, c12, c44):
    """Cubic Cartesian tensor; inputs/outputs have the SAME stress units."""
    values = np.array([c11, c12, c44], dtype=float)
    if not np.all(np.isfinite(values)) or min(c11-c12, c11+2*c12, c44) <= 0:
        raise ValueError("positive cubic elastic energy is required")
    eye = np.eye(3)
    result = (c12*np.einsum("ij,kl->ijkl", eye, eye)
              + c44*(np.einsum("ik,jl->ijkl", eye, eye)
                     + np.einsum("il,jk->ijkl", eye, eye)))
    for i in range(3):
        result[i, i, i, i] += c11-c12-2*c44
    return result


def rotate_elastic_tensor(tensor, basis_rows):
    """Rows are the new unit vectors expressed in the old Cartesian basis."""
    basis = np.asarray(basis_rows, dtype=float)
    if basis.shape != (3, 3) or not np.allclose(basis@basis.T, np.eye(3), atol=1e-12):
        raise ValueError("an orthonormal basis is required")
    return np.einsum("ia,jb,kc,ld,abcd->ijkl", basis, basis, basis, basis, tensor)


@dataclass(frozen=True)
class SurfaceImpedance:
    upper: np.ndarray
    lower: np.ndarray
    jump_per_wave_number: np.ndarray
    hermitian_residual: float
    riccati_residual: float


def halfspace_impedance(tensor, direction):
    r"""Eliminate two identical elastic half spaces for exp(i |k| t.x).

    Coordinates 0,1 are in-plane and 2 is the interface normal. Returns Z
    in stress units, so the jump stiffness is |k| * jump_per_wave_number.
    The surface energy is 1/2 d_hat^H K d_hat (Fourier normalization separate).
    Ordered invariant subspaces avoid a degenerate-eigenvector inverse in
    the isotropic limit. No inertial mass or frequency appears.
    """
    c = np.asarray(tensor, dtype=float)
    t = np.asarray(direction, dtype=float)
    if c.shape != (3, 3, 3, 3) or not np.all(np.isfinite(c)):
        raise ValueError("finite 3x3x3x3 elastic tensor required")
    if t.shape != (2,) or not np.all(np.isfinite(t)) or np.linalg.norm(t) == 0:
        raise ValueError("nonzero finite in-plane direction required")
    t = np.r_[t/np.linalg.norm(t), 0.]
    scale = float(np.max(np.abs(c)))
    if scale == 0:
        raise ValueError("zero elastic tensor")
    c = c/scale
    # Exact minor and major symmetries, not an arbitrary symmetrization of C.
    if max(np.max(abs(c-c.swapaxes(0, 1))), np.max(abs(c-c.swapaxes(2, 3))),
           np.max(abs(c-c.transpose(2, 3, 0, 1)))) > 2e-12:
        raise ValueError("elastic tensor symmetries violated")
    kelvin_basis=[]
    for i,j in ((0,0),(1,1),(2,2),(1,2),(0,2),(0,1)):
        mode=np.zeros((3,3))
        if i==j:
            mode[i,j]=1.
        else:
            mode[i,j]=mode[j,i]=1/np.sqrt(2.)
        kelvin_basis.append(mode)
    kelvin=np.einsum("aij,ijkl,bkl->ab",kelvin_basis,c,kelvin_basis)
    if np.linalg.eigvalsh(kelvin)[0]<=0:
        raise ValueError("nonpositive bulk elastic energy")
    a = np.einsum("ijkl,j,l->ik", c, t, t)
    b = np.einsum("ikl,l->ik", c[:, 2, :, :], t)
    d = c[:, 2, :, 2]
    if np.linalg.eigvalsh(d)[0] <= 0:
        raise ValueError("nonpositive normal acoustic block")
    inv_d = np.linalg.inv(d)
    state = np.block([[-1j*inv_d@b, inv_d],
                      [a-b.T@inv_d@b, -1j*b.T@inv_d]])
    impedances = []
    errors = []
    riccati = []
    for upper in (True, False):
        _, vectors, count = schur(state, output="complex",
            sort=(lambda z: z.real < 0) if upper else (lambda z: z.real > 0))
        if count != 3:
            raise ValueError("three decaying elastic modes are required in each half space")
        u = vectors[:3, :3]
        traction = vectors[3:, :3]
        z = np.linalg.solve(u.T, traction.T).T * (-1 if upper else 1)
        errors.append(float(np.max(abs(z-z.conj().T))))
        if errors[-1] > 1e-9:
            raise FloatingPointError("non-Hermitian surface impedance")
        z = (z+z.conj().T)/2  # remove recorded roundoff only
        if np.linalg.eigvalsh(z)[0] <= 0:
            raise ValueError("unstable static elastic half space")
        # T=traction/displacement; invariant graph obeys a matrix Riccati identity.
        T = -z if upper else z
        residual = state[3:, :3]+state[3:, 3:]@T-T@state[:3, :3]-T@state[:3, 3:]@T
        riccati.append(float(np.max(abs(residual))))
        impedances.append(z*scale)
    upper, lower = impedances
    compliance = np.linalg.inv(upper)+np.linalg.inv(lower)
    jump = np.linalg.inv(compliance)
    return SurfaceImpedance(upper, lower, jump, max(errors), max(riccati))


def interface_kernel(tensor, wave_vector):
    """3-vector displacement-jump stiffness; q=0 is FREE rigid translation.

    With tensor [Pa] and wave_vector [1/m], returns [Pa/m]. Do not add the
    local interface Hessian here: it belongs to the misfit energy only once.
    """
    wave = np.asarray(wave_vector, dtype=float)
    if wave.shape != (2,) or not np.all(np.isfinite(wave)):
        raise ValueError("finite two-component wave vector required")
    magnitude = np.linalg.norm(wave)
    if magnitude == 0:
        return np.zeros((3, 3), dtype=complex)
    return magnitude*halfspace_impedance(tensor, wave).jump_per_wave_number


def screw_elastic_factor(tensor, slip_index=0, variation_index=1):
    """Independent anti-plane reduction: sqrt(A D - B^2), same units as C.

    Only valid when the designated slip polarization decouples by symmetry.
    For FCC(111), m=e1=[1,-1,0]/sqrt(2), variation=e2 has this mirror symmetry.
    """
    c = np.asarray(tensor, dtype=float)
    m, x, z = slip_index, variation_index, 2
    others = [i for i in range(3) if i != m]
    cross = [abs(c[m, j, k, l]) for j in (x, z) for k in others for l in (x, z)]
    if max(cross) > 1e-11*np.max(abs(c)):
        raise ValueError("anti-plane polarization does not decouple")
    A, B, D = c[m, x, m, x], c[m, z, m, x], c[m, z, m, z]
    if min(A, D, A*D-B*B) <= 0:
        raise ValueError("unstable anti-plane energy")
    return float(np.sqrt(A*D-B*B))
