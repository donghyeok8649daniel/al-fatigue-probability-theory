"""Current analytic research site law on infinite Bessel atomic rows.

This is a representation bridge, NOT a refit or production registration.
Keep every existing v19/v20 site term, including the extra radial quadrupole,
quartic rank3 invariant, even saturation, and density/rank1 coupling. Each
atom sums its whole environment BEFORE this local nonlinear function.
"""
from __future__ import annotations

import numpy as np

from .coordination_screening import CoordinationScreenedInterface, screening_factor
from .isolated_screw_core import IsolatedScrewCore
from .vector_fcc_rows import (
    VectorRowKernel, exponential_radial, symmetric_monomials, _radial_polynomial,
)


class CurrentMaterialRowKernel(VectorRowKernel):
    """22 infinite-row channels: pair,rho,Q1,Q2_even,Q3,Q2_odd.

    The inherited 17 channels retain their exact arithmetic. The added five
    components are the existing odd-range STF2 needed to form the Eg radial
    combination. They are NOT a new environmental contribution/parameter.
    """
    def __init__(self, rows, *, tolerance=2e-12, max_modes=40):
        super().__init__(rows, tolerance=tolerance, max_modes=max_modes)
        self.parity = np.r_[self.parity, np.ones(5)]

    def _mode(self, vectors, mode, order):
        old, oldg, oldh, envelope = super()._mode(vectors, mode, order)
        x, y, z = vectors.T
        g = 2*np.pi*mode/self.rows.b
        moment = self.rows.surface.angular  # unchanged rank3/odd radial kernel
        H, Ha, Haa = exponential_radial(np.hypot(y, z), g, moment.kappa)
        factor = (1 if mode == 0 else 2)*moment.amplitude/self.rows.b
        exponents, basis = symmetric_monomials(2)
        vals, grads, hess = [], [], []
        for nx, ny, nz in exponents:
            value, grad, tensor = _radial_polynomial(
                H[:, nx], Ha[:, nx], Haa[:, nx], y, z, ny, nz, order)
            coefficient = factor*1j**nx
            value = value*coefficient
            grad = grad*coefficient
            vals.append(value)
            grads.append(np.column_stack([1j*g*value, grad]))
            if order == 2:
                tensor = tensor*coefficient
                block = np.zeros((len(x), 3, 3), complex)
                block[:, 0, 0] = -g*g*value
                block[:, 0, 1:] = 1j*g*grad
                block[:, 1:, 0] = 1j*g*grad
                block[:, 1, 1] = tensor[:, 0]
                block[:, 1, 2] = block[:, 2, 1] = tensor[:, 1]
                block[:, 2, 2] = tensor[:, 2]
                hess.append(block)
        extra = np.stack(vals, axis=1)@basis
        eg = np.einsum('nci,ck->nki', np.stack(grads, axis=1), basis)
        bound = np.maximum(np.max(abs(extra), axis=1), np.max(abs(eg), axis=(1, 2)))
        phase = np.exp(1j*g*x)
        eh = None
        if order == 2:
            eh = np.einsum('ncij,ck->nkij', np.stack(hess, axis=1), basis)
            bound = np.maximum(bound, np.max(abs(eh), axis=(1, 2, 3)))
            eh = np.concatenate([oldh, np.real(eh*phase[:, None, None, None])], axis=1)
        return (np.concatenate([old, np.real(extra*phase[:, None])], axis=1),
                np.concatenate([oldg, np.real(eg*phase[:, None, None])], axis=1),
                eh, np.maximum(envelope, bound))


class CurrentMaterialSiteLaw:
    """Analytic energy, gradient and Hessian in the 22 summed site channels.

    Channels are changes from the fixed cubic infinite reference. All Q_bulk
    vanish there, but all actual noncubic site moments below remain present.
    Energy is eV per atom; channels use the SAME frozen density/radial gauges
    as CoordinationScreenedInterface. No per-plane nonlinear embedding.
    """
    def __init__(self, model):
        if not isinstance(model, CoordinationScreenedInterface):
            raise TypeError('explicit current CoordinationScreenedInterface required')
        self.model = model
        self.coefficients = model.coefficients.copy()
        self.embedding = model.face.bulk.embedding
        self.rho_bulk = float(model.face.rho_bulk)
        self.rho_ref = float(self.embedding.rho_ref)
        if not np.isclose(self.rho_bulk/self.rho_ref, 1., rtol=2e-12, atol=0):
            raise ValueError('this bridge requires the unchanged cubic density reference')
        self.maps = {}
        for name, start, size in [('Q1', 2, 3), ('Q2', 5, 5), ('Q3', 10, 7)]:
            self.maps[name] = np.eye(22)[start:start+size]
        self.maps['QE'] = (np.eye(22)[17:22]-model.moment.eta*np.eye(22)[5:10])/model.moment.normalization

    def evaluate(self, channels, *, second=True):
        z = np.asarray(channels, float)
        if z.ndim != 2 or z.shape[1] != 22 or np.any(~np.isfinite(z)):
            raise ValueError('finite (sites,22) summed environment changes required')
        rho = self.rho_bulk+z[:, 1]
        if np.any(rho <= 0):
            raise ValueError('nonpositive actual site density; no clipping')
        c = self.coefficients
        energy = .5*z[:, 0]+self.embedding.value(rho)-self.embedding.value(self.rho_bulk)
        grad = np.zeros_like(z)
        grad[:, 0], grad[:, 1] = .5, self.embedding.first_derivative(rho)
        hess = np.zeros((len(z), 22, 22)) if second else None
        if second:
            hess[:, 1, 1] = self.embedding.second_derivative(rho)

        def invariant(name):
            B = self.maps[name]
            Q = z@B.T
            return np.sum(Q*Q, axis=1), 2*Q@B, 2*(B.T@B)

        def add_norm(name, scale, alpha=0., normalization=1., quartic=False):
            nonlocal energy, grad, hess
            I, Ip, Ipp = invariant(name)
            if quartic:
                value, first, second_value = I*I, 2*I, np.full_like(I, 2.)
            else:
                k = alpha/normalization
                den = 1+k*I
                value, first, second_value = I/den, den**-2, -2*k*den**-3
            energy += scale*value
            grad += scale*first[:, None]*Ip
            if second:
                hess += scale*(first[:, None, None]*Ipp
                    + second_value[:, None, None]*Ip[:, :, None]*Ip[:, None, :])

        add_norm('Q3', c[5])
        add_norm('Q3', c[9], quartic=True)
        add_norm('Q2', c[7], self.model.alpha, self.model.even_reference_curvature)
        add_norm('QE', c[8], self.model.alpha)
        I, Ip, Ipp = invariant('Q1')
        x = rho/self.rho_ref
        g, gx, gxx = screening_factor(x, self.model.screening, law=self.model.law)
        energy += c[6]*g*I
        grad += c[6]*g[:, None]*Ip
        grad[:, 1] += c[6]*gx*I/self.rho_ref
        if second:
            hess += c[6]*g[:, None, None]*Ipp
            hess[:, 1, :] += c[6]*gx[:, None]*Ip/self.rho_ref
            hess[:, :, 1] += c[6]*gx[:, None]*Ip/self.rho_ref
            hess[:, 1, 1] += c[6]*gxx*I/self.rho_ref**2
        return energy, grad, hess


class CurrentMaterialScrewCore(IsolatedScrewCore):
    """Full current per-atom law on the existing isolated-row geometry.

    Fixed boundary and transverse neighborhood remain explicit approximations.
    Constructor accepts a research candidate, not material/kinetic adoption.
    """
    def __init__(self, model, far_field, **kwargs):
        self.current_material = model
        self.site_law = CurrentMaterialSiteLaw(model)
        super().__init__(model.base.surface, far_field, **kwargs)

    def _make_kernel(self, tolerance):
        return CurrentMaterialRowKernel(self.rows, tolerance=tolerance)

    def _site_response(self, channels, *, second):
        energy, weights, hessian = self.site_law.evaluate(channels, second=second)
        apply = None if hessian is None else lambda dz: np.einsum('nij,nj->ni', hessian, dz)
        return energy, weights, apply
