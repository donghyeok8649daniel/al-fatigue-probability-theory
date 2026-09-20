"""Classical harmonic conditional measure and inertial Si research diagnostics.

No production energy, probability generator or overdamped clock is supplied.
Atomic mass fixes the time of the explicitly integrated Newtonian reference,
NOT a conversion of the existing overdamped a/s model time.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.constants import atomic_mass, electron_volt, Boltzmann
from scipy.linalg import cho_factor, cho_solve, eigh
from scipy.sparse import csr_matrix, eye, issparse


KB_EV_K = Boltzmann / electron_volt
AMU_EV_PS2_A2 = atomic_mass * 1e4 / electron_volt
SI_MASS_AMU = 28.0855  # Explicit reference mass, not a dopant/carrier model.


@dataclass
class ConditionalHarmonic:
    bath_hessian: np.ndarray
    cross: np.ndarray
    bare_curvature: np.ndarray
    relaxed_curvature: np.ndarray
    bath_response: np.ndarray
    relaxed_lift: np.ndarray
    q_metric: np.ndarray
    logdet_bath: float


def conditional_harmonic(hessian, bath_basis, q_lift):
    """Integrate the Gaussian on a linear, orthogonal Cartesian constraint.

    R = B z + D q, B.T B=I, B.T D=0; D need not be normalized.
    The coarea measure for a physical relative gap differs by a CONSTANT.
    Its q dependence is zero and it cancels only for the same constraint family.
    An unstable/zero conditional mode is rejected, never regularized or abs'ed.
    """
    h = csr_matrix(hessian, dtype=float)
    b = csr_matrix(bath_basis, dtype=float)
    d = np.asarray(q_lift, dtype=float)
    if d.ndim == 1:
        d = d[:, None]
    n = h.shape[0]
    if (h.shape != (n, n) or d.ndim != 2 or d.shape[0] != n
            or b.shape != (n, n-d.shape[1]) or not d.shape[1] or not b.shape[1]
            or not np.isfinite(h.data).all() or not np.isfinite(b.data).all()
            or not np.isfinite(d).all()):
        raise ValueError('finite complete Cartesian constraint basis required')
    asymmetry = h-h.T
    if np.max(np.abs(asymmetry.data), initial=0) > 1e-10:
        raise ValueError('symmetric Hessian required')
    orthogonality = b.T@b-eye(b.shape[1], format='csr')
    if (np.max(np.abs(orthogonality.data), initial=0) > 1e-10
            or np.max(np.abs(b.T@d), initial=0) > 1e-10):
        raise ValueError('orthonormal bath orthogonal to q lift required')
    metric = d.T@d
    cho_factor(metric, lower=True)  # reject dependent retained coordinates
    bath = (b.T@h@b).toarray()
    factor = cho_factor(bath, lower=True, check_finite=True)
    cross = np.asarray(b.T@h@d)
    bare = np.asarray(d.T@h@d)
    response = cho_solve(factor, cross)
    return ConditionalHarmonic(bath, cross, bare, bare-cross.T@response,
        response, np.asarray(d-b@response), metric,
        float(2*np.log(np.diag(factor[0])).sum()))


def harmonic_free_energy_difference(energy_difference, logdet, reference_logdet, temperature):
    """Same-constraint classical Gaussian difference, not a validated PMF."""
    t = np.asarray(temperature, dtype=float)
    if np.any(~np.isfinite(t)) or np.any(t < 0):
        raise ValueError('finite nonnegative temperature required')
    if not np.all(np.isfinite([energy_difference, logdet, reference_logdet])):
        raise ValueError('finite energy and log determinants required')
    return energy_difference+.5*KB_EV_K*t*(logdet-reference_logdet)


def harmonic_memory(reduction, time_ps, *, mass_amu=SI_MASS_AMU):
    """Exact matrix memory for a finite, undamped, equal-mass harmonic bath.

    m_q q'' + K q + integral Gamma(t-u) q'(u) du = noise,
    m_q = m D.T D, K = Hqq-Hqz Hzz^-1 Hzq.
    Gamma(t) = sum_j c_j c_j.T/lambda_j cos(sqrt(lambda_j/m) t).
    Conditional canonical preparation gives <noise(t) noise(0).T>=kBT Gamma(t).
    This finite cosine sum has no certified constant-drag/Markov limit.
    """
    times = np.asarray(time_ps, float)
    if (times.ndim != 1 or not np.isfinite(times).all() or np.any(times < 0)
            or not np.isfinite(mass_amu) or mass_amu <= 0):
        raise ValueError('finite times and positive atomic mass required')
    values, vectors = eigh(reduction.bath_hessian, check_finite=True)
    if values[0] <= 0:
        raise ValueError('positive conditional spectrum required')
    couplings = vectors.T@reduction.cross
    weights = np.einsum('ji,jk->jik', couplings, couplings)/values[:, None, None]
    omega = np.sqrt(values/(mass_amu*AMU_EV_PS2_A2))
    kernel = np.einsum('tj,jik->tik', np.cos(np.outer(times, omega)), weights)
    return dict(time_ps=times, kernel_eV_A2=kernel, eigenvalues_eV_A2=values,
        eigenvectors=vectors, angular_frequency_ps=omega, weights_eV_A2=weights,
        q_mass_eV_ps2_A2=mass_amu*AMU_EV_PS2_A2*reduction.q_metric)


def harmonic_release(hessian, displacement, observation, time_ps, *, mass_amu=SI_MASS_AMU):
    """Exact full Cartesian stable harmonic release from zero velocity."""
    h = hessian.toarray() if issparse(hessian) else np.asarray(hessian, float)
    x, c, times = np.asarray(displacement, float), np.asarray(observation, float), np.asarray(time_ps, float)
    if (h.ndim != 2 or h.shape[0] != h.shape[1] or x.shape != (len(h),)
            or c.shape[-1] != len(h) or not np.isfinite(h).all()
            or not np.allclose(h, h.T, rtol=0, atol=1e-10)
            or not np.isfinite(x).all() or not np.isfinite(c).all()
            or times.ndim != 1 or not np.isfinite(times).all() or np.any(times < 0)
            or not np.isfinite(mass_amu) or mass_amu <= 0):
        raise ValueError('finite stable Cartesian system required')
    values, vectors = eigh(h)
    if values[0] <= 0:
        raise ValueError('stable full Hessian required for bounded release')
    omega = np.sqrt(values/(mass_amu*AMU_EV_PS2_A2))
    weights = (c@vectors)*(vectors.T@x)
    return np.cos(np.outer(times, omega))@weights.T


def velocity_verlet(model, positions, fixed, *, dt_ps, steps, observation,
                    velocities=None, mass_amu=SI_MASS_AMU, stride=1, progress=None):
    """Fresh conservative trajectory of free atoms with truly fixed grips.

    Total site energy differences are used to avoid extensive-energy cancellation.
    The optional progress callback receives completed steps, not a physical clock
    for the production probability solver. No thermostat or force rescaling.
    """
    r = np.asarray(positions, float).copy()
    fixed = np.asarray(fixed, bool)
    v = np.zeros_like(r) if velocities is None else np.asarray(velocities, float).copy()
    if (r.ndim != 2 or r.shape[1] != 3 or fixed.shape != (len(r),) or v.shape != r.shape
            or not np.isfinite(r).all() or not np.isfinite(v).all() or np.all(fixed)
            or np.any(v[fixed] != 0) or not np.isfinite(dt_ps) or dt_ps <= 0
            or int(steps) != steps or steps < 1 or int(stride) != stride or stride < 1
            or steps % stride or not np.isfinite(mass_amu) or mass_amu <= 0):
        raise ValueError('valid positions, fixed grips, mass and integral step counts required')
    mass = mass_amu*AMU_EV_PS2_A2
    free = ~fixed
    initial_grips = r[fixed].copy()
    current = model.evaluate(r)
    sites0 = current.site_energy.copy()
    kinetic0 = .5*mass*np.sum(v[free]**2)
    times, observables, potentials, kinetics, residuals = [], [], [], [], []
    max_residual = 0.
    for step in range(int(steps)+1):
        if step:
            v[free] -= .5*dt_ps*current.gradient[free]/mass
            r[free] += dt_ps*v[free]
            current = model.evaluate(r)
            v[free] -= .5*dt_ps*current.gradient[free]/mass
        potential = float(np.sum(current.site_energy-sites0))
        kinetic = float(.5*mass*np.sum(v[free]**2))
        residual = potential+kinetic-kinetic0
        if not np.isfinite(residual):
            raise FloatingPointError('nonfinite conservative trajectory')
        max_residual = max(max_residual, abs(residual))
        if step % stride == 0:
            times.append(step*dt_ps)
            observables.append(np.asarray(observation(r), float))
            potentials.append(potential)
            kinetics.append(kinetic)
            residuals.append(residual)
        if progress is not None and step and step % 1000 == 0:
            progress(step)
    return dict(time_ps=np.asarray(times), observations=np.asarray(observables),
        potential_difference_eV=np.asarray(potentials), kinetic_energy_eV=np.asarray(kinetics),
        energy_residual_eV=np.asarray(residuals), max_energy_residual_eV=max_residual,
        max_fixed_displacement_A=float(np.max(np.abs(r[fixed]-initial_grips), initial=0)),
        final_positions_A=r, final_velocities_A_ps=v, dt_ps=float(dt_ps), steps=int(steps),
        mass_amu=float(mass_amu), stride=int(stride), production_t0_seconds=None)
