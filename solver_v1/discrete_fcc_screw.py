"""Infinite-row LJ/Bessel static FCC screw Hessian (research, not a PDE).

The line direction is e1. Atoms along each line are summed by the SAME
infinite one-dimensional Poisson identity as the reduced LJ reference; the
transverse row/layer indices remain discrete. No physical neighbor cutoff,
empirical core length, kinetic mass, or probability aggregation area is used.

This is an exact harmonic reduction for row-uniform anti-plane displacement,
with explicitly converged reciprocal/transverse series. It is not a nonlinear
dislocation core, nor a prescription for adding a second local GSF stiffness.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product
import math

import numpy as np
from scipy.special import kv

from .angular_environment_reference import traceless_second, traceless_third
from .lattice_bessel import _mode_coefficient_and_deda


def exponential_row_transform_derivatives(radius, wave_number, kappa):
    """H, H_G, H_GG, H_GGG for H=2 A kappa K1(A Q)/Q.

    Q=sqrt(kappa^2+G^2), t=G^2 and
    d_t^j H=2 A kappa (-A/2)^j Q^(-j-1) K_(j+1)(A Q).
    These are analytic derivatives, not finite differences.
    """
    radius, wave_number = np.broadcast_arrays(np.asarray(radius, float),
                                              np.asarray(wave_number, float))
    if np.any(radius <= 0) or not np.isfinite(kappa) or kappa <= 0:
        raise ValueError("positive row radius and exponential decay required")
    q = np.sqrt(kappa**2+wave_number**2)
    derivatives = [2*radius*kappa*(-radius/2)**j/q**(j+1)*kv(j+1, radius*q)
                   for j in range(4)]
    return (derivatives[0], 2*wave_number*derivatives[1],
            2*derivatives[1]+4*wave_number**2*derivatives[2],
            12*wave_number*derivatives[2]+8*wave_number**3*derivatives[3])


def row_exponential_slip_derivatives(y, z, delta, *, b, kappa, amplitude=1.,
                                     tolerance=1e-13, max_modes=96):
    """Sum f_xx and d_x[f STF(R^rank)] (rank 1,2,3) over an INFINITE row.

    R=(n*b+delta,y,z). Multiplication by x^p transforms to i^p H^(p),
    and differentiation to iG. STF is applied after summing raw tensors.
    G=0 contributes zero to every returned derivative. A same-line bond is
    omitted by the caller because both atoms undergo the same displacement.
    """
    y, z, delta = np.broadcast_arrays(np.asarray(y, float), np.asarray(z, float),
                                      np.asarray(delta, float))
    if b <= 0 or tolerance <= 0 or max_modes < 2 or not np.isfinite(amplitude):
        raise ValueError("invalid infinite-row controls")
    radius = np.hypot(y, z)
    density_xx = np.zeros(y.shape)
    # p is the number of x factors, up to rank 3. Other factors y,z are
    # multiplied before STF. This avoids a huge symbolic expression.
    polynomial_derivatives = np.zeros(y.shape+(4,))
    small = 0
    for m in range(1, max_modes+1):
        g = 2*np.pi*m/b
        hd = exponential_row_transform_derivatives(radius, g, kappa)
        phase = np.exp(1j*g*delta)
        prefactor = 2*amplitude/b
        xx = -prefactor*g*g*hd[0]
        density_xx += xx*phase.real
        envelope = float(np.max(abs(xx)))
        for p in range(4):
            coefficient = prefactor*1j*g*1j**p*hd[p]
            polynomial_derivatives[..., p] += (coefficient*phase).real
            envelope = max(envelope, float(np.max(abs(coefficient)*
                           np.maximum(1., radius)**(3-p))))
        if envelope < tolerance:
            small += 1
            if small == 2:
                break
        else:
            small = 0
    else:
        raise ArithmeticError("infinite exponential-row derivative modes did not converge")
    tensors = {}
    for rank in (1, 2, 3):
        raw = np.empty(y.shape+(3,)*rank)
        for indices in product(range(3), repeat=rank):
            raw[(Ellipsis,)+indices] = (polynomial_derivatives[..., indices.count(0)] *
                                        y**indices.count(1)*z**indices.count(2))
        tensors[rank] = (raw if rank == 1 else
                         traceless_second(raw) if rank == 2 else traceless_third(raw))
    return dict(density_xx=density_xx, moments=tensors, modes=m,
                last_mode_absolute_envelope=envelope)


def row_lj_slip_curvature(radius, delta, *, b, epsilon, sigma,
                          tolerance=1e-13, max_modes=96):
    """d_delta^2 sum_n phi_LJ(sqrt(A^2+(nb+delta)^2)).

    Reuses the verified LJ Bessel mode coefficient, not a new pair potential.
    The zero mode is registry independent and vanishes on differentiation.
    """
    radius, delta = np.broadcast_arrays(np.asarray(radius, float), np.asarray(delta, float))
    if np.any(radius <= 0) or b <= 0 or epsilon <= 0 or sigma <= 0 or tolerance <= 0:
        raise ValueError("invalid row/pair parameters")
    total = np.zeros(radius.shape)
    small = 0
    for m in range(1, max_modes+1):
        g = 2*np.pi*m/b
        c6 = _mode_coefficient_and_deda(radius, 6, m, b)[0]
        c3 = _mode_coefficient_and_deda(radius, 3, m, b)[0]
        term = -4*epsilon*g*g*(sigma**12*c6-sigma**6*c3)
        total += term*np.cos(g*delta)
        envelope = float(np.max(4*epsilon*g*g*(sigma**12*abs(c6)+sigma**6*abs(c3))))
        if envelope < tolerance:
            small += 1
            if small == 2:
                return total, m, envelope
        else:
            small = 0
    raise ArithmeticError("infinite LJ-row curvature modes did not converge")


@dataclass(frozen=True)
class RowSumDiagnostics:
    rings_used: int
    validation_rings: int
    rows_used: int
    reciprocal_modes: int
    last_ring_absolute_envelope: float
    omitted_validation_ring_envelope: float
    tail_status: str = "absolute-shell refinement estimate; not a rigorous infinite bound"


class FCCScrewRowHessian:
    """Exact anti-plane bulk symbol and explicitly defined layer-jump compliance.

    Rows (j,l): y=d*(j+l/3), z=h*l, delta_x=b*(j+l)/2 modulo b,
    d=sqrt(3)*b/2. The displacement is along e1 and constant along that row.
    A Bloch state has phase q_y*d*j+theta*l. This keeps ABC phases explicitly.
    q_z=(theta-q_y*d/3)/h. The physical symbol is K(q_y,theta), eV/L0^2.

    Scalar density has zero first variation for EVERY complete row at perfect
    registry, so F'' does not contribute in this subspace. F' DOES contribute
    through phi_eff=phi_LJ+2 F'(rho_bulk) f. Angular STF moments are per atom,
    summed over all rows before squaring. No plane-wise embedding is used.
    """
    def __init__(self, bulk, *, D1=0., D2=0., D3=0., angular_decay=None,
                 angular_amplitude=None, tolerance=1e-11, max_ring=16, fixed_ring=None):
        self.bulk = bulk
        self.b = float(bulk.geometry.b)
        self.h = float(bulk.a0)
        self.d = math.sqrt(3)*self.b/2
        self.amplitudes = {1: float(D1), 2: float(D2), 3: float(D3)}
        if tolerance <= 0 or max_ring < 4 or (fixed_ring is not None and not 1 <= fixed_ring < max_ring):
            raise ValueError("invalid transverse-series controls")
        if D2 and not np.isclose(self.h/self.b, math.sqrt(2/3), rtol=2e-8):
            raise ValueError("Q_bulk=0 requires the actual cubic reference")
        rho_bulk = bulk.full_environment_density(self.h, 0.).value
        first = float(bulk.embedding.first_derivative(rho_bulk))
        j, l = np.meshgrid(np.arange(-max_ring, max_ring+1),
                           np.arange(-max_ring, max_ring+1), indexing="ij")
        nonself = (j != 0) | (l != 0)
        j, l = j[nonself], l[nonself]
        y, z = self.d*(j+l/3), self.h*l
        delta = self.b/2*((j+l) % 2)
        radius = np.hypot(y, z)
        pair, pair_modes, _ = row_lj_slip_curvature(radius, delta, b=self.b,
            epsilon=bulk.p.epsilon, sigma=bulk.p.sigma_lj, tolerance=tolerance/100)
        density = row_exponential_slip_derivatives(y, z, delta, b=self.b,
            kappa=bulk.density_params.kappa, amplitude=bulk.density_params.C_rho,
            tolerance=tolerance/100)
        effective = pair+2*first*density["density_xx"]
        moments = {rank: np.zeros((len(j), 3**rank)) for rank in (1, 2, 3)}
        modes = max(pair_modes, density["modes"])
        if any(self.amplitudes.values()):
            if angular_decay is None or angular_amplitude is None:
                raise ValueError("angular decay AND normalization must be supplied from the same energy")
            angular = row_exponential_slip_derivatives(y, z, delta, b=self.b,
                kappa=angular_decay, amplitude=angular_amplitude, tolerance=tolerance/100)
            moments = {r: angular["moments"][r].reshape(len(j), 3**r) for r in (1, 2, 3)}
            modes = max(modes, angular["modes"])
        ring = np.maximum(abs(j), abs(l))
        # Absolute row envelopes, BEFORE Bloch cancellation. Successive square
        # rings have increasing minimum transverse distance. Nonzero row modes
        # decay exponentially ~exp(-2*pi*A/b). Report a refinement estimate,
        # not an unproved rigorous remainder bound.
        # Bound a SYMBOL difference in consistent stiffness units, not a sum
        # of unlike pair and moment-gradient coefficient magnitudes. With
        # ||L_r||<=2 S_r and omitted ||Delta L_r||<=2 T_r,
        # |Delta K_r| <=16 |D_r| S_total,r T_r (finite validation set).
        # The still-unseen infinite remainder is assessed by ring refinement.
        envelope = 2*abs(effective)
        for rank in (1, 2, 3):
            norm = np.linalg.norm(moments[rank], axis=1)
            envelope += 16*abs(self.amplitudes[rank])*sum(norm)*norm
        ring_sum = np.array([sum(envelope[ring == r]) for r in range(1, max_ring+1)])
        selected = fixed_ring
        if selected is None:
            for r in range(4, max_ring-1):
                if np.max(ring_sum[r-2:r]) < tolerance:
                    selected = r
                    break
        if selected is None:
            raise ArithmeticError("transverse infinite-row series did not converge; increase max_ring")
        keep = ring <= selected
        self.j, self.l = j[keep], l[keep]
        self.effective_pair = effective[keep]
        self.moments = {r: moments[r][keep] for r in (1, 2, 3)}
        self.diagnostics = RowSumDiagnostics(selected, max_ring, int(sum(keep)), modes,
                                             float(ring_sum[selected-1]), float(sum(ring_sum[selected:])))

    @classmethod
    def from_surface(cls, surface, **controls):
        return cls(surface.interface.bulk, D1=surface.vector_amplitude_ev,
                   D2=surface.quadrupole_amplitude_ev, D3=surface.amplitude_ev,
                   angular_decay=surface.angular.kappa,
                   angular_amplitude=surface.angular.amplitude, **controls)

    def symbol(self, q_y, theta, *, components=False):
        """K for e1 displacement; theta is the logical layer Bloch phase."""
        q_y, theta = np.broadcast_arrays(np.asarray(q_y, float), np.asarray(theta, float))
        if not np.all(np.isfinite(q_y)) or not np.all(np.isfinite(theta)):
            raise ValueError("wave numbers must be finite")
        phase = q_y[..., None]*self.d*self.j+theta[..., None]*self.l
        pair = (2*np.sin(phase/2)**2)@self.effective_pair
        angular = np.zeros_like(pair)
        for rank in (1, 2, 3):
            if self.amplitudes[rank]:
                collective = np.expm1(1j*phase)@self.moments[rank]
                angular += 2*self.amplitudes[rank]*np.sum(abs(collective)**2, axis=-1)
        total = pair+angular
        return dict(total=total, effective_pair=pair, angular=angular) if components else total

    def rigid_step_stiffness(self):
        """Uniform rigid upper/lower translation, exactly the GSF Hessian.

        A cross-plane separation l occurs l times across a cut. Angular
        derivative at depth r sums layers l>=r+1 BEFORE taking its norm.
        Factor 4 = two half crystals times d^2(D ||Q||^2)/ds^2.
        """
        positive = self.l > 0
        pair = float(sum(self.l[positive]*self.effective_pair[positive]))
        angular = 0.
        for rank in (1, 2, 3):
            for depth in range(int(max(self.l))):
                collective = self.moments[rank][self.l > depth].sum(axis=0)
                angular += 4*self.amplitudes[rank]*float(collective@collective)
        return dict(total=pair+angular, effective_pair=pair, angular=angular)

    def jump_stiffness(self, q_y, *, layer_phase_points=2048, constraint="logical_rows"):
        """Infinite-stack Schur stiffness for s_j=u_(j,1)-u_(j,0).

        C(q_y)=integral |exp(i theta)-1|^2/K(q_y,theta) dtheta/(2*pi).
        This constraint pairs LOGICAL rows, whose physical y positions differ
        by d/3. It is NOT silently an interpolated same-y continuum slip.
        The midpoint trapezoid avoids the removable acoustic zero at (0,0).
        Must be independently refined, especially near q_y=0.
        """
        if layer_phase_points < 32 or layer_phase_points % 2:
            raise ValueError("use an even layer-phase quadrature of at least 32 points")
        if constraint not in ("logical_rows", "spectral_same_y"):
            raise ValueError("explicit logical_rows or spectral_same_y constraint required")
        q_y = float(q_y)
        theta = -np.pi+2*np.pi*(np.arange(layer_phase_points)+.5)/layer_phase_points
        # Chunks bound memory and avoid changing the numerical quadrature.
        pieces = [self.symbol(q_y, block) for block in np.array_split(theta, max(1, len(theta)//256))]
        symbol = np.concatenate(pieces)
        if np.any(symbol <= 0) or not np.all(np.isfinite(symbol)):
            raise ValueError("unstable/undefined anti-plane bulk mode; cannot invert compliance")
        # spectral_same_y translates the upper layer by -d/3 in the Fourier
        # interpolation convention. It is a DIFFERENT finite-q observable,
        # not an algebraic correction to the actual adjacent-atom jump.
        offset = 0. if constraint == "logical_rows" else self.d/3
        compliance = float(np.mean(4*np.sin((theta-q_y*offset)/2)**2/symbol))
        return dict(stiffness=1/compliance, compliance=compliance,
                    layer_phase_points=layer_phase_points, minimum_sampled_bulk_symbol=float(min(symbol)),
                    constraint=constraint)

    def long_wave_jump_slope(self, tensor, stiffness_zero, *, constraint="logical_rows"):
        """Coefficient of |q_y| in the discrete K_jump, from the acoustic pole.

        tensor is C in eV/L0^3 and the ACTUAL stacked plane frame. This slope
        generally is NOT A_atomic_cell*mu_bar/2: an atomic gap constraint is
        not a mathematical displacement discontinuity between cut halfspaces.
        Its microscopic compliance and the ABC phase enter explicitly.
        """
        if constraint not in ("logical_rows", "spectral_same_y"):
            raise ValueError("unknown layer-jump constraint")
        c = np.asarray(tensor, float)
        a, b, d = c[0, 1, 0, 1], c[0, 2, 0, 1], c[0, 2, 0, 2]
        if d <= 0 or a*d <= b*b or stiffness_zero <= 0:
            raise ValueError("stable acoustic tensor and positive zero-wave stiffness required")
        volume = self.bulk.geometry.atomic_cell_area*self.h
        curvature = volume*d/self.h**2
        alpha = self.d/3-self.h*b/d
        beta = self.h*math.sqrt(a*d-b*b)/d
        offset = 0. if constraint == "logical_rows" else self.d/3
        slope = stiffness_zero**2*(beta**2-(alpha-offset)**2)/(2*curvature*beta)
        return dict(slope=slope, pole_shift=alpha, pole_width_per_q=beta,
                    layer_acoustic_curvature=curvature, constraint=constraint)

    def converged_jump_stiffness(self, q_y, *, constraint="logical_rows", initial_points=256,
                                 max_points=65536, absolute_tolerance=1e-10, relative_tolerance=1e-10):
        """Deterministic doubling; an unresolved acoustic pole is not certified.

        Two successive differences must be below the declared tolerance. This
        is an observed quadrature-convergence certificate, not a kinetic one.
        """
        if absolute_tolerance <= 0 or relative_tolerance < 0 or max_points < 4*initial_points:
            raise ValueError("invalid layer-phase refinement controls")
        previous = None
        small = 0
        history = []
        points = initial_points
        while points <= max_points:
            result = self.jump_stiffness(q_y,layer_phase_points=points,constraint=constraint)
            history.append(result)
            change = math.inf if previous is None else abs(result["stiffness"]-previous)
            if change <= absolute_tolerance+relative_tolerance*abs(result["stiffness"]):
                small += 1
                if small == 2:
                    return dict(**result, observed_quadrature_change=change,
                                quadrature_converged=True, refinement_history=history)
            else:
                small = 0
            previous = result["stiffness"]
            points *= 2
        raise ArithmeticError(f"layer-phase quadrature unresolved at q_y={q_y:g}; last change={change:g}")
