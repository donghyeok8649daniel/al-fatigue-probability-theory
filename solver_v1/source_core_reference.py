"""Mishin Al99 core TARGET evaluator; NEVER an LJ/Bessel production model.

The published source has a finite cutoff. Its finite row sum is exact for
that source and is not substituted for the research hybrid's infinite LJ.
Same geometry, per-atom counting, vector degrees of freedom and boundary
protocol provide an independent static comparator, not a fitted yield law.
"""
from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from .isolated_screw_core import IsolatedScrewCore
from .periodic_plane_covariance import SourceEAMBlochHessian


def source_elastic_tensor(reference, length_scale_m):
    """Source affine energy tensor, eV/L0^3, including actual prestress.

    C_iJkL=(.5 sum_R V_ik R_J R_L + F'' D_iJ D_kL)/Omega,
    V=Hess(phi+2F'_bulk f), D_iJ=sum_R f_i R_J.
    This is the same source's long-wave Hessian, not fitted Al moduli.
    """
    if not np.isfinite(length_scale_m) or length_scale_m <= 0:
        raise ValueError('positive fixed physical length scale required')
    operator = SourceEAMBlochHessian(reference)
    R = operator.R
    D = np.einsum('ni,nj->ij', operator.density_gradient, R)
    tensor = .5*np.einsum('nik,nj,nl->ijkl', operator.pair_hessian, R, R)
    tensor += operator.second*np.einsum('ij,kl->ijkl', D, D)
    volume = reference.geometry.lattice_constant**3/4
    return tensor/volume*(length_scale_m/1e-10)**3


class SourceRowKernel:
    """All actual source neighbors in each infinite-periodic atomic row."""
    def __init__(self, rows):
        self.rows = rows
        self.parity = np.ones(2)

    def evaluate(self, vectors, *, order=1):
        vectors = np.asarray(vectors, float)
        if vectors.ndim != 2 or vectors.shape[1] != 3 or np.any(~np.isfinite(vectors)) or order not in (1, 2):
            raise ValueError('finite vector rows and derivative order 1/2 required')
        source, length = self.rows.source, self.rows.length_angstrom
        shifted = vectors.copy()
        shifted[:, 0] -= self.rows.b*np.floor(shifted[:, 0]/self.rows.b+.5)
        images = int(np.ceil(source.r[-1]/(self.rows.b*length)))+1
        atomic = np.repeat(shifted[:, None, :], 2*images+1, axis=1)*length
        atomic[:, :, 0] += self.rows.b*length*np.arange(-images, images+1)
        radius = np.linalg.norm(atomic, axis=-1)
        if np.any(radius == 0):
            raise ValueError('source core contains coincident atomic sites')
        inside = radius < source.r[-1]
        r = radius[inside]; unit = atomic[inside]/r[:, None]
        rr = unit[:, :, None]*unit[:, None, :]
        z, zp, zpp = source._rphi(r), source._rphi(r, 1), source._rphi(r, 2)
        radial = [(z/r, zp/r-z/r**2, zpp/r-2*zp/r**2+2*z/r**3),
                  (source._rho(r), source._rho(r, 1), source._rho(r, 2))]
        values, gradients, hessians = [], [], []
        for value, first, second in radial:
            v = np.zeros(radius.shape); g = np.zeros(radius.shape+(3,))
            v[inside] = value; g[inside] = first[:, None]*unit*length
            values.append(v.sum(axis=1)); gradients.append(g.sum(axis=1))
            if order == 2:
                h = np.zeros(radius.shape+(3, 3))
                h[inside] = (second[:, None, None]*rr
                    +(first/r)[:, None, None]*(np.eye(3)-rr))*length**2
                hessians.append(h.sum(axis=1))
        return dict(value=np.stack(values, axis=1), gradient=np.stack(gradients, axis=1),
            hessian=np.stack(hessians, axis=1) if order == 2 else None,
            modes_used=0, maximum_last_mode_envelope=0.,
            source_cutoff_exact=True, reciprocal_Bessel_model=False,
            published_cutoff_angstrom=source.cutoff, effective_support_angstrom=source.r[-1])


class MishinScrewCoreReference(IsolatedScrewCore):
    """Static target-only core; source F after each atom's full density."""
    def __init__(self, source, far_field, *, length_scale_m, **kwargs):
        if not np.isfinite(length_scale_m) or length_scale_m <= 0:
            raise ValueError('positive fixed length scale required')
        self.source = source
        self.length_angstrom = length_scale_m/1e-10
        super().__init__(source, far_field, **kwargs)
        self.production_eligible = False
        self.reference_only = True

    def _make_rows(self, source):
        b, h = source.geometry.b/self.length_angstrom, source.h/self.length_angstrom
        return SimpleNamespace(b=b, h=h, d=np.sqrt(3)*b/2, rho_bulk=source._rho_bulk,
            source=source, length_angstrom=self.length_angstrom)

    def _make_kernel(self, tolerance):
        return SourceRowKernel(self.rows)

    def _site_response(self, channels, *, second):
        density = self.rows.rho_bulk+channels[:, 1]
        energy = .5*channels[:, 0]+self.source.F(density)-self.source.F(self.rows.rho_bulk)
        weights = np.column_stack([np.full(len(density), .5), self.source._F(density, 1)])
        apply = None
        if second:
            curvature = self.source._F(density, 2)
            def apply(delta):
                return np.column_stack([np.zeros(len(delta)), curvature*delta[:, 1]])
        return energy, weights, apply

    def _assemble(self, field, *, second=False):
        full = self.full_field(field)
        transform = np.array([[self.rows.d, self.rows.d/3], [0., self.rows.h]])
        alpha = np.linalg.svd(transform, compute_uv=False).min()
        displacement_bound = 2*float(np.max(np.linalg.norm(full[:, 1:], axis=1)))
        if alpha*(self.ring+1)-displacement_bound <= self.source.r[-1]/self.length_angstrom:
            raise ValueError('source neighborhood cannot enclose all actual cutoff neighbors; enlarge ring')
        return super()._assemble(field, second=second)
