"""Conditional classical ensembles and conservative reference dynamics.

Sampling is a numerical integration/preparation tool, never a fracture counter
or a clock for the production probability equation. Finite configuration
domains and sampling diagnostics must accompany every free-energy estimate.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import cholesky, solve_triangular
from scipy.special import logsumexp

from .silicon_conditional_research import KB_EV_K


@dataclass
class GaussianReference:
    """x = mean + sqrt(kBT) L^-T u, with u standard normal and H=L L.T."""

    mean: np.ndarray
    lower: np.ndarray
    temperature_K: float

    @classmethod
    def from_hessian(cls, hessian, temperature_K, mean=None):
        h = np.asarray(hessian, float)
        if (h.ndim != 2 or h.shape[0] != h.shape[1] or not len(h)
                or not np.isfinite(h).all() or not np.allclose(h, h.T, atol=1e-10, rtol=0)
                or not np.isfinite(temperature_K) or temperature_K <= 0):
            raise ValueError('positive temperature and symmetric finite Hessian required')
        center = np.zeros(len(h)) if mean is None else np.asarray(mean, float)
        if center.shape != (len(h),) or not np.isfinite(center).all():
            raise ValueError('one finite mean per coordinate required')
        return cls(center.copy(), cholesky(h, lower=True), float(temperature_K))

    @property
    def thermal_energy(self):
        return KB_EV_K*self.temperature_K

    @property
    def logdet_hessian(self):
        return float(2*np.log(np.diag(self.lower)).sum())

    def transform(self, standard_normal):
        u = np.asarray(standard_normal, float)
        if u.ndim not in (1, 2) or u.shape[-1] != len(self.mean) or not np.isfinite(u).all():
            raise ValueError('finite row vectors of reference dimension required')
        return self.mean+np.sqrt(self.thermal_energy)*solve_triangular(
            self.lower.T, u.T, lower=False, check_finite=False).T

    def whiten(self, coordinates):
        x = np.asarray(coordinates, float)
        if x.ndim not in (1, 2) or x.shape[-1] != len(self.mean) or not np.isfinite(x).all():
            raise ValueError('finite row vectors of reference dimension required')
        return (x-self.mean)@self.lower/np.sqrt(self.thermal_energy)


def exponential_reweight(residual_eV, temperature_K, *, inside=None, observables=None):
    """Finite-domain partition / UNTRUNCATED Gaussian reference partition.

    Outside-domain proposals retain their place in the denominator N with zero
    target weight. Discarding them and dividing by accepted N changes the free
    energy. ESS and maximum weight diagnose overlap; neither certifies it.
    """
    residual = np.asarray(residual_eV, float)
    if (residual.ndim != 1 or not len(residual) or not np.isfinite(residual).all()
            or not np.isfinite(temperature_K) or temperature_K <= 0):
        raise ValueError('finite residuals and positive temperature required')
    mask = np.ones(len(residual), bool) if inside is None else np.asarray(inside, bool)
    if mask.shape != residual.shape or not np.any(mask):
        raise ValueError('at least one proposal in the declared domain required')
    logw = np.where(mask, -residual/(KB_EV_K*temperature_K), -np.inf)
    normalizer = logsumexp(logw)
    weights = np.exp(logw-normalizer)
    result = dict(correction_eV=float(-KB_EV_K*temperature_K*(normalizer-np.log(len(logw)))),
        effective_sample_size=float(1/(weights@weights)), sample_count=len(residual),
        maximum_normalized_weight=float(weights.max()), domain_fraction=float(mask.mean()),
        normalized_weights=weights)
    if observables is not None:
        values = np.asarray(observables, float)
        if values.shape[0] != len(residual) or not np.isfinite(values).all():
            raise ValueError('finite observables at every proposal required')
        result['mean'] = np.einsum('i,i...->...', weights, values)
    return result


def pcn_sample(evaluate, initial, *, proposal_scale, draws, warmup, rng, stride=1,
               progress=None):
    """Metropolis-corrected Gaussian-preserving proposals in whitened space.

    evaluate(u) returns (Phi(u), observables); target is exp(-Phi) N(0,I).
    Phi=+infinity denotes the exterior of a finite, explicitly declared domain.
    Proposal count has NO physical-time interpretation. Repeated states after
    rejection are retained. No adaptation or per-sample retuning occurs.
    """
    current = np.asarray(initial, float).copy()
    if (current.ndim != 1 or not len(current) or not np.isfinite(current).all()
            or not 0 < proposal_scale <= 1 or int(draws) != draws or draws < 1
            or int(warmup) != warmup or warmup < 0 or int(stride) != stride or stride < 1):
        raise ValueError('finite initial vector and valid sampling settings required')
    potential, observation = evaluate(current)
    if not np.isfinite(potential):
        raise ValueError('initial sample must lie in the finite target domain')
    persistence = np.sqrt(1-proposal_scale**2)
    samples, values, potentials, accepted = [], [], [], []
    count, accepted_warmup = 0, 0
    for iteration in range(int(warmup+draws*stride)):
        proposed = persistence*current+proposal_scale*rng.standard_normal(len(current))
        proposed_potential, proposed_observation = evaluate(proposed)
        if np.isnan(proposed_potential) or proposed_potential == -np.inf:
            raise FloatingPointError('invalid target correction')
        accept = np.log(rng.random()) < min(0., potential-proposed_potential)
        if accept:
            current, potential, observation = proposed, proposed_potential, proposed_observation
        if iteration < warmup:
            accepted_warmup += int(accept)
        else:
            count += int(accept)
            if (iteration-warmup+1) % stride == 0:
                samples.append(current.copy())
                values.append(np.asarray(observation, float).copy())
                potentials.append(float(potential))
                accepted.append(bool(accept))
        if progress is not None and (iteration+1) % 500 == 0:
            progress(iteration+1)
    return dict(samples=np.asarray(samples), observations=np.asarray(values),
        potentials=np.asarray(potentials), accepted_at_saved_step=np.asarray(accepted),
        acceptance_fraction=count/(draws*stride),
        warmup_acceptance_fraction=accepted_warmup/warmup if warmup else None,
        proposal_scale=float(proposal_scale), warmup=int(warmup), stride=int(stride),
        final_rng_state=rng.bit_generator.state, physical_time_ps=None)


def block_statistics(values, *, blocks=16):
    """Batch-mean standard error diagnostic; no automatic convergence claim."""
    x = np.asarray(values, float)
    if (x.ndim < 1 or len(x) < 2*blocks or blocks < 2 or int(blocks) != blocks
            or not np.isfinite(x).all()):
        raise ValueError('finite history with at least two samples per block required')
    chunks = np.array_split(x, blocks)
    means = np.asarray([chunk.mean(axis=0) for chunk in chunks])
    return dict(mean=x.mean(axis=0), block_means=means,
        block_standard_error=means.std(axis=0, ddof=1)/np.sqrt(blocks),
        first_half_mean=x[:len(x)//2].mean(axis=0),
        second_half_mean=x[len(x)//2:].mean(axis=0), block_count=int(blocks))


def split_rhat(chains):
    """Classical split-Rhat diagnostic, not rank normalized or a proof of mixing."""
    x = np.asarray(chains, float)
    if x.ndim < 2 or x.shape[0] < 2 or x.shape[1] < 8 or not np.isfinite(x).all():
        raise ValueError('at least two finite chains with eight draws required')
    half = x.shape[1]//2
    split = np.concatenate([x[:, :half], x[:, -half:]], axis=0)
    within = split.var(axis=1, ddof=1).mean(axis=0)
    between = half*split.mean(axis=1).var(axis=0, ddof=1)
    if np.any(within <= 0):
        raise ValueError('nonzero within-chain variance required')
    return np.sqrt(((half-1)*within/half+between/half)/within)


def bounded_force_control(coordinates, bath_gradient, response, temperature_K, halfwidth):
    """Zero-mean integration-by-parts control for the SAME fixed bath box.

    g_i = response_i (1-x_i^2/R^2) vanishes on each respective face. Therefore
    E[g.grad U-kBT div(g)]=0 exactly for exp(-U/kBT) on (-R,R)^d, including a
    nonzero boundary density. Subtract from raw dU/dq to reduce variance.
    This identity needs the FULL declared ensemble, not an unannounced basin.
    """
    x, grad, v = (np.asarray(a, float) for a in (coordinates, bath_gradient, response))
    if (x.ndim != 1 or grad.shape != x.shape or v.shape != x.shape
            or not np.isfinite(x).all() or not np.isfinite(grad).all() or not np.isfinite(v).all()
            or not np.isfinite(halfwidth) or halfwidth <= 0 or np.any(abs(x) >= halfwidth)
            or not np.isfinite(temperature_K) or temperature_K <= 0):
        raise ValueError('finite vectors inside a positive fixed box required')
    return float((v*(1-(x/halfwidth)**2))@grad
                 +2*KB_EV_K*temperature_K*(v@x)/halfwidth**2)


def harmonic_split_proposal(evaluate, position, momentum, step, steps, initial_evaluation=None):
    """Reversible symplectic split of H=(u.u+p.p)/2+Phi(u).

    The harmonic rotation is exact. An attempted excursion outside the target
    domain rejects the entire proposal; accepted paths have valid reversed paths.
    The integration coordinate and momentum are sampling auxiliaries, not atoms.
    """
    u, p = np.asarray(position, float).copy(), np.asarray(momentum, float).copy()
    if u.ndim != 1 or p.shape != u.shape or not np.isfinite([step]).all() or step <= 0 or steps < 1:
        raise ValueError('matching vectors and positive integration controls required')
    potential, gradient, observation = evaluate(u) if initial_evaluation is None else initial_evaluation
    if not np.isfinite(potential):
        raise ValueError('start inside the declared sampling domain')
    cosine, sine = np.cos(step), np.sin(step)
    for _ in range(steps):
        p -= .5*step*gradient
        u, p = cosine*u+sine*p, -sine*u+cosine*p
        potential, gradient, observation = evaluate(u)
        if potential == np.inf:
            return None
        if not np.isfinite(potential) or not np.isfinite(gradient).all():
            raise FloatingPointError('nonfinite HMC potential or gradient')
        p -= .5*step*gradient
    return u,p,potential,gradient,observation


def harmonic_hmc_sample(evaluate, initial, *, draws, warmup, rng, step=.2, steps=(6,10)):
    """Metropolized harmonic-reference HMC with a full momentum refresh.

    The exact target is exp(-Phi(u))N(0,I) on the declared domain. A random,
    state-independent integral step count avoids choosing a resonance. All
    rejections are retained; no discretization bias is accepted as physics.
    """
    u = np.asarray(initial,float).copy()
    if (u.ndim != 1 or not np.isfinite(u).all() or int(draws)!=draws or draws<1
            or int(warmup)!=warmup or warmup<0 or len(steps)!=2
            or any(int(n)!=n or n<1 for n in steps) or steps[1]<steps[0]):
        raise ValueError('valid finite HMC state and proposal budget required')
    current = evaluate(u)
    if not np.isfinite(current[0]):
        raise ValueError('initial HMC point must be in the declared domain')
    samples, observations, potentials, accepted, energy_errors = [],[],[],[],[]
    counts = [0,0]
    outside = 0
    for iteration in range(warmup+draws):
        momentum = rng.standard_normal(len(u))
        length = int(rng.integers(steps[0],steps[1]+1))
        initial_energy = .5*(u@u+momentum@momentum)+current[0]
        proposal = harmonic_split_proposal(evaluate,u,momentum,step,length,current)
        if proposal is None:
            accept, error = False, np.inf
            outside += 1
        else:
            proposed,new_momentum,potential,gradient,observation = proposal
            error = float(.5*(proposed@proposed+new_momentum@new_momentum)+potential-initial_energy)
            accept = np.log(rng.random()) < min(0.,-error)
            if accept:
                u, current = proposed,(potential,gradient,observation)
        counts[int(iteration>=warmup)] += int(accept)
        if iteration>=warmup:
            samples.append(u.copy()); observations.append(np.asarray(current[2],float).copy())
            potentials.append(float(current[0])); accepted.append(bool(accept)); energy_errors.append(error)
    return dict(samples=np.asarray(samples),observations=np.asarray(observations),potentials=np.asarray(potentials),
        accepted_at_saved_step=np.asarray(accepted),hamiltonian_errors=np.asarray(energy_errors),
        acceptance_fraction=counts[1]/draws,warmup_acceptance_fraction=counts[0]/warmup if warmup else None,
        domain_rejected_trajectories=outside,step=step,steps=steps,
        final_rng_state=rng.bit_generator.state,physical_time_ps=None)
