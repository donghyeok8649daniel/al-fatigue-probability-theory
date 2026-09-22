"""Research bridge from anisotropic crack loading to pure-Si atomic boundaries.

No strength fit, probability generator, material approval, or physical clock.
Length: Angstrom in the atomic field; m in specimen formulas. Elasticity: GPa.
K: MPa sqrt(m). G: J/m^2 per projected crack area (two new surfaces).

The plane-strain field follows the anisotropic Airy-potential construction
of Sih, Paris & Irwin (1965), DOI 10.1007/BF00186854. See the accompanying
v6 note for the derivation and independent constitutive/J-integral checks.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .silicon_crack_research import diamond_crack_strip


SI_111_FRAME = np.array([[1., 1., -2.], [1., 1., 1.], [1., -1., 0.]])
SI_111_FRAME /= np.linalg.norm(SI_111_FRAME, axis=1)[:, None]
K_MPA_M_TO_GPA_A = 100.  # 1e-3 GPa/MPa * sqrt(1e10 A/m)
G_GPA_A_TO_J_M2 = .1


def _positive(value, name):
    value = np.asarray(value, float)
    if not np.all(np.isfinite(value)) or np.any(value <= 0):
        raise ValueError(name + ' must be finite and positive')
    return value


def specimen_K(stress_MPa, crack_length_m, *, geometry_factor):
    """K=Y sigma sqrt(pi a); Y must be explicitly supplied for the geometry.

    For an infinite plate with a central through crack, a is HALF crack length
    and Y=1. This function does not infer Y from a mesh or a notch radius.
    """
    stress = np.asarray(stress_MPa, float)
    if not np.all(np.isfinite(stress)) or np.any(stress < 0):
        raise ValueError('tensile stress must be finite and nonnegative')
    return stress*_positive(geometry_factor, 'geometry factor')*np.sqrt(
        np.pi*_positive(crack_length_m, 'crack length in m'))


def specimen_stress(K_MPa_sqrt_m, crack_length_m, *, geometry_factor):
    k = np.asarray(K_MPa_sqrt_m, float)
    if not np.all(np.isfinite(k)) or np.any(k < 0):
        raise ValueError('mode-I K must be finite and nonnegative')
    return k/(_positive(geometry_factor, 'geometry factor')*np.sqrt(
        np.pi*_positive(crack_length_m, 'crack length in m')))


def cubic_tensor(c11, c12, c44, frame=SI_111_FRAME):
    """Rotate a cubic stiffness with ROWS of frame expressed in cubic axes."""
    if not np.all(np.isfinite([c11, c12, c44])) or min(c44, c11-c12, c11+2*c12) <= 0:
        raise ValueError('stable finite cubic stiffness required')
    frame = np.asarray(frame, float)
    if (frame.shape != (3, 3) or not np.all(np.isfinite(frame))
            or not np.allclose(frame@frame.T, np.eye(3), atol=1e-12, rtol=0)
            or np.linalg.det(frame) < 0):
        raise ValueError('right-handed orthonormal row frame required')
    eye = np.eye(3)
    return (c12*np.einsum('ij,kl->ijkl', eye, eye)
            + c44*(np.einsum('ik,jl->ijkl', eye, eye)+np.einsum('il,jk->ijkl', eye, eye))
            + (c11-c12-2*c44)*np.einsum('ia,ja,ka,la->ijkl', frame, frame, frame, frame))


class AnisotropicModeI:
    """Plane strain, straight crack x<0, with decoupled anti-plane elasticity.

    Engineering strain ordering is (xx, yy, 2xy). Arbitrary orientations that
    couple the in-plane field to anti-plane displacement are rejected.
    Exactly isotropic material is evaluated with its repeated-root limit.
    """
    def __init__(self, c11, c12, c44, *, frame=SI_111_FRAME):
        c = cubic_tensor(c11, c12, c44, frame)
        pairs = [(0, 0), (1, 1), (0, 1)]
        if max(abs(c[i, j, k, 2]) for i, j in pairs for k in (0, 1)) > 1e-10*max(c11, c44):
            raise ValueError('this orientation requires a coupled three-component crack field')
        self.stiffness = np.array([[c[i, j, k, l] for k, l in pairs] for i, j in pairs])
        self.compliance = np.linalg.inv(self.stiffness)
        self.isotropic = abs(c11-c12-2*c44) <= 1e-11*max(abs(c11), abs(c44))
        if self.isotropic:
            self.mu = float(c44)
            self.nu = float(c12/(2*(c12+c44)))
            self.kappa = 3-4*self.nu
            self.H_GPa_inv = (1-self.nu)/(2*self.mu)
            self.roots = np.array([1j, 1j])
            self.root_residual = 0.
            return
        s = self.compliance
        polynomial = np.array([s[0, 0], -2*s[0, 2], 2*s[0, 1]+s[2, 2], -2*s[1, 2], s[1, 1]])
        roots = np.roots(polynomial/polynomial[0])
        p = roots[roots.imag > 1e-8]
        if len(p) != 2 or abs(p[0]-p[1]) < 1e-7:
            raise ValueError('anisotropic roots not sufficiently resolved')
        self.roots = p
        self.root_residual = float(np.max(abs(np.polyval(polynomial, p)))/max(abs(polynomial)))
        # Complex Airy coefficients enforce sum(d)=1, sum(p*d)=0. Hence
        # sigma_yy(x,0)=K/sqrt(2*pi*x), sigma_xy(x,0)=0 exactly.
        self.weights = np.linalg.solve(np.array([np.ones(2), p]), [1., 0.])
        displacement_vectors = np.array([s[0, 0]*p*p+s[0, 1]-s[0, 2]*p,
                                         s[0, 1]*p+s[1, 1]/p-s[1, 2]])
        self.amplitudes = displacement_vectors*self.weights
        self.H_GPa_inv = float(.5*np.real(1j*np.sum(self.amplitudes[1])))
        if not self.H_GPa_inv > 0 or self.root_residual > 1e-9:
            raise RuntimeError('crack compliance/root validation failed')

    def field(self, xy_A, K_MPa_sqrt_m):
        """Return displacement A, gradient du_i/dx_j, and tensile stress GPa.

        Crack-face side is selected by signed y. The tip itself is singular
        and rejected; no continuum tip stress is sampled as an atomic force.
        """
        xy = np.asarray(xy_A, float)
        k = float(K_MPa_sqrt_m)*K_MPA_M_TO_GPA_A
        if xy.shape[-1:] != (2,) or not np.all(np.isfinite(xy)) or not np.isfinite(k) or k < 0:
            raise ValueError('finite xy coordinates and nonnegative K required')
        r = np.linalg.norm(xy, axis=-1)
        if np.any(r <= 0):
            raise ValueError('continuum crack field is singular at its tip')
        if self.isotropic:
            theta = np.arctan2(xy[..., 1], xy[..., 0])
            unit = np.stack([np.cos(theta/2), np.sin(theta/2)], axis=-1)
            angular = (self.kappa-np.cos(theta))[..., None]*unit
            derivative = (np.sin(theta)[..., None]*unit
                          + .5*(self.kappa-np.cos(theta))[..., None]*np.stack([-unit[..., 1], unit[..., 0]], axis=-1))
            scale = k/(2*self.mu*np.sqrt(2*np.pi))
            u = scale*np.sqrt(r)[..., None]*angular
            radial = scale*.5/np.sqrt(r)[..., None]*angular
            tangential = scale/np.sqrt(r)[..., None]*derivative
            grad = np.stack([radial*np.cos(theta)[..., None]-tangential*np.sin(theta)[..., None],
                             radial*np.sin(theta)[..., None]+tangential*np.cos(theta)[..., None]], axis=-1)
        else:
            z = xy[..., 0, None]+xy[..., 1, None]*self.roots
            # Preserve negative-zero y on the lower crack face (complex
            # addition can otherwise erase its signed imaginary zero).
            z = z.real+1j*z.imag
            z.imag = xy[..., 1, None]*self.roots.imag
            root = np.sqrt(z)
            u = 2*k/np.sqrt(2*np.pi)*np.real(np.einsum('in,...n->...i', self.amplitudes, root))
            dx = k/np.sqrt(2*np.pi)*np.real(np.einsum('in,...n->...i', self.amplitudes, 1/root))
            dy = k/np.sqrt(2*np.pi)*np.real(np.einsum('in,...n->...i', self.amplitudes, self.roots/root))
            grad = np.stack([dx, dy], axis=-1)
        strain = np.stack([grad[..., 0, 0], grad[..., 1, 1], grad[..., 0, 1]+grad[..., 1, 0]], axis=-1)
        stress = np.einsum('ij,...j->...i', self.stiffness, strain)
        return u, grad, stress

    def energy_release_J_m2(self, K_MPa_sqrt_m):
        k = np.asarray(K_MPa_sqrt_m, float)
        if not np.all(np.isfinite(k)) or np.any(k < 0):
            raise ValueError('nonnegative finite mode-I K required')
        return self.H_GPa_inv*(K_MPA_M_TO_GPA_A*k)**2*G_GPA_A_TO_J_M2

    def griffith_K(self, separation_work_J_m2):
        """Energetic equality only; not lattice-trapping K+ or measured K_Ic."""
        return np.sqrt(_positive(separation_work_J_m2, 'separation work')/
                       (self.H_GPa_inv*K_MPA_M_TO_GPA_A**2*G_GPA_A_TO_J_M2))

    def contour_J(self, K_MPa_sqrt_m, *, radius_A, points=4096):
        """Independent circular J integral from stress, strain and du/dx."""
        radius = float(_positive(radius_A, 'contour radius'))
        if type(points) is not int or points < 64:
            raise ValueError('at least 64 angular integration points required')
        theta = -np.pi+(np.arange(points)+.5)*2*np.pi/points
        normal = np.stack([np.cos(theta), np.sin(theta)], axis=-1)
        _, grad, stress = self.field(radius*normal, K_MPa_sqrt_m)
        strain = np.stack([grad[:, 0, 0], grad[:, 1, 1], grad[:, 0, 1]+grad[:, 1, 0]], axis=-1)
        energy = .5*np.sum(stress*strain, axis=-1)
        traction = np.stack([stress[:, 0]*normal[:, 0]+stress[:, 2]*normal[:, 1],
                              stress[:, 2]*normal[:, 0]+stress[:, 1]*normal[:, 1]], axis=-1)
        integrand = energy*normal[:, 0]-np.sum(traction*grad[:, :, 0], axis=-1)
        return float(np.mean(integrand)*2*np.pi*radius*G_GPA_A_TO_J_M2)


@dataclass
class AtomicCrackBoundary:
    reference: np.ndarray
    fixed: np.ndarray
    bonds: np.ndarray
    period_A: float
    radius_A: float
    grip_width_A: float
    tip_offset_A: float

    def displaced(self, elasticity, K_MPa_sqrt_m):
        r = self.reference.copy()
        r[:, :2] += elasticity.field(r[:, :2], K_MPa_sqrt_m)[0]
        return r


def atomic_crack_boundary(lattice_A, *, radius_A, front_repeats=4, grip_width_A=8.):
    """Cylindrical (111)[11-2] boundary layer, periodic along [1-10].

    No interactions are removed. Only outer atoms have prescribed K-field
    positions. All three coordinates of interior atoms remain independent.
    The nominal tip lies halfway between neighboring cross-plane bond columns.
    """
    lattice = float(_positive(lattice_A, 'lattice'))
    radius = float(_positive(radius_A, 'radius'))
    grip = float(_positive(grip_width_A, 'grip width'))
    if grip >= radius/2:
        raise ValueError('grip must be thinner than half the domain radius')
    nx = int(np.ceil((2*radius+3*lattice)/(lattice*np.sqrt(1.5))))
    ny = int(np.ceil((2*radius+3*lattice)/(lattice*np.sqrt(3))))
    strip = diamond_crack_strip(lattice, nx=nx, ny=ny, nz=front_repeats, grip_width=grip)
    r = strip.positions.copy()
    columns = np.unique(np.round(r[strip.crossing_bonds].mean(axis=1)[:, 0], 8))
    tips = .5*(columns[:-1]+columns[1:])
    tip = float(tips[np.argmin(abs(tips))])
    r[:, 0] -= tip
    distance = np.linalg.norm(r[:, :2], axis=1)
    selected = distance < radius
    index = np.full(len(r), -1, int)
    index[selected] = np.arange(np.sum(selected))
    bonds = strip.crossing_bonds[np.all(selected[strip.crossing_bonds], axis=1)]
    fixed = distance[selected] >= radius-grip
    if not np.any(fixed) or np.all(fixed):
        raise ValueError('boundary must contain fixed and free atoms')
    return AtomicCrackBoundary(r[selected], fixed, index[bonds], strip.front_period,
                               radius, grip, tip)
