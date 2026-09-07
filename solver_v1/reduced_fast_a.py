"""Experimental finite-temperature elimination of the fast normal coordinate.

This module does not replace :mod:`solver_v1.probability_pde_2d`.  It provides
fixed-normal-domain potential-of-mean-force calculations, moving-boundary
diagnostics, a discrete fast-process quasi-stationary distribution, and an
experimental conservative slow-``s`` equation for validation against the full
two-dimensional reference.

The effective configurational free energy is denoted ``F_eff`` in Python and
``mathcal F_eff`` in documentation.  It is not a geometric area.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Callable

import numpy as np
from scipy import sparse
from scipy.optimize import root_scalar
from scipy.sparse.linalg import spsolve

from .model import TwoRowLJ
from .probability_pde_2d import CyclicLoad2D, _bernoulli
from .dynamics_diagnostics import stable_normal_equilibrium


@lru_cache(maxsize=16)
def _gauss_legendre_rule(order: int) -> tuple[np.ndarray, np.ndarray]:
    if order < 16:
        raise ValueError("normal quadrature order must be at least 16")
    nodes, weights = np.polynomial.legendre.leggauss(int(order))
    nodes.setflags(write=False)
    weights.setflags(write=False)
    return nodes, weights


def _generalized_energy_derivatives(
    model: TwoRowLJ,
    a,
    s,
    force: float,
    *,
    energy_offset: float = 0.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    aa, ss = np.broadcast_arrays(np.asarray(a, dtype=float), np.asarray(s, dtype=float))
    energy, deda, deds = model.local_energy_gradient_array(aa, ss)
    generalized = (
        energy
        - float(force)
        * ((aa - model.a0) + model.p.chi_axial_projection * ss)
        + float(energy_offset)
    )
    return (
        np.asarray(generalized, dtype=float),
        np.asarray(deda, dtype=float) - float(force),
        np.asarray(deds, dtype=float)
        - float(force) * model.p.chi_axial_projection,
    )


@dataclass(frozen=True)
class FastAConditional:
    """Quadrature representation of equilibrium in ``a`` conditional on ``s``."""

    s: float
    force: float
    a_lower: float
    a_upper: float
    log_partition: float
    F_eff: float
    mean_a: float
    mean_dG_ds: float
    boundary_correction_s: float
    dF_eff_ds: float
    nodes: np.ndarray
    probability_weights: np.ndarray


def conditional_fast_a(
    model: TwoRowLJ,
    s: float,
    force: float,
    *,
    a_lower: float,
    a_upper: float,
    quadrature_order: int = 128,
    lower_s_derivative: float = 0.0,
    upper_s_derivative: float = 0.0,
    energy_offset: float = 0.0,
) -> FastAConditional:
    r"""Evaluate ``rho_eq(a|s,f)`` and ``mathcal F_eff=-kT log Z_a``.

    ``probability_weights`` are normalized quadrature masses and therefore sum
    to one.  When either integration limit depends on ``s``, the returned
    derivative includes the exact Leibniz boundary contribution.
    """

    if model.p.kT <= 0.0:
        raise ValueError("conditional equilibrium requires kT > 0")
    if not np.isfinite(a_lower) or not np.isfinite(a_upper) or a_upper <= a_lower:
        raise ValueError("normal integration bounds must be finite and ordered")
    base_nodes, base_weights = _gauss_legendre_rule(int(quadrature_order))
    half_width = 0.5 * (float(a_upper) - float(a_lower))
    midpoint = 0.5 * (float(a_upper) + float(a_lower))
    nodes = midpoint + half_width * base_nodes
    weights = half_width * base_weights
    generalized, _, dG_ds = _generalized_energy_derivatives(
        model, nodes, float(s), float(force), energy_offset=energy_offset
    )
    reference = float(np.min(generalized))
    scaled = weights * np.exp(
        np.clip(-(generalized - reference) / model.p.kT, -745.0, 0.0)
    )
    scaled_partition = float(np.sum(scaled))
    if not np.isfinite(scaled_partition) or scaled_partition <= 0.0:
        raise FloatingPointError("normal conditional partition is not positive")
    probability_weights = scaled / scaled_partition
    log_partition = -reference / model.p.kT + np.log(scaled_partition)
    F_eff = -model.p.kT * log_partition
    mean_a = float(np.sum(probability_weights * nodes))
    mean_dG_ds = float(np.sum(probability_weights * dG_ds))

    boundary_energy, _, _ = _generalized_energy_derivatives(
        model,
        np.array([a_lower, a_upper], dtype=float),
        float(s),
        float(force),
        energy_offset=energy_offset,
    )
    # Ratios exp(-G(boundary)/kT)/Z are evaluated in the same shifted scale.
    boundary_density = np.exp(
        np.clip(-(boundary_energy - reference) / model.p.kT, -745.0, 700.0)
    ) / scaled_partition
    boundary_correction = -model.p.kT * (
        float(boundary_density[1]) * float(upper_s_derivative)
        - float(boundary_density[0]) * float(lower_s_derivative)
    )
    return FastAConditional(
        s=float(s),
        force=float(force),
        a_lower=float(a_lower),
        a_upper=float(a_upper),
        log_partition=float(log_partition),
        F_eff=float(F_eff),
        mean_a=mean_a,
        mean_dG_ds=mean_dG_ds,
        boundary_correction_s=float(boundary_correction),
        dF_eff_ds=float(mean_dG_ds + boundary_correction),
        nodes=nodes,
        probability_weights=probability_weights,
    )


def _normal_stationary_roots(
    model: TwoRowLJ,
    s: float,
    force: float,
    *,
    scan_points: int = 1000,
) -> tuple[tuple[float, float], ...]:
    lower = float(model.p.a_min * (1.0 + 1.0e-10))
    upper = float(model.p.a_max)
    scan = np.linspace(lower, upper, int(scan_points))
    values = model.local_deda_array(scan, float(s)) - float(force)
    roots: list[tuple[float, float]] = []
    for index in range(scan.size - 1):
        left = float(values[index])
        right = float(values[index + 1])
        if not np.isfinite(left) or not np.isfinite(right) or left * right > 0.0:
            continue
        if left == 0.0:
            root_value = float(scan[index])
        else:
            solved = root_scalar(
                lambda a: model.local_deda(float(a), float(s)) - float(force),
                bracket=(float(scan[index]), float(scan[index + 1])),
                method="brentq",
                xtol=5.0e-15,
                rtol=4.0 * np.finfo(float).eps,
            )
            if not solved.converged:
                continue
            root_value = float(solved.root)
        stiffness = float(model.local_hessian(root_value, float(s))[0, 0])
        if not roots or abs(root_value - roots[-1][0]) > 1.0e-10:
            roots.append((root_value, stiffness))
    return tuple(roots)


def bound_basin_conditional(
    model: TwoRowLJ,
    s: float,
    force: float,
    *,
    a_lower: float | None = None,
    quadrature_order: int = 192,
) -> FastAConditional:
    """Return truncated Gibbs statistics below the exact opening saddle.

    The upper-bound derivative follows from differentiating
    ``G_a(a_dagger,s;f)=0``.  This object is not asserted to be the absorbing
    fast process's quasi-stationary distribution.
    """

    roots = _normal_stationary_roots(model, float(s), float(force))
    unstable = [(a, stiffness) for a, stiffness in roots if stiffness < 0.0]
    if not unstable:
        raise ValueError("no finite opening saddle exists at this state")
    a_dagger, saddle_stiffness = max(unstable, key=lambda item: item[0])
    hessian = model.local_hessian(a_dagger, float(s))
    upper_s_derivative = -float(hessian[0, 1]) / float(saddle_stiffness)
    lower = (
        float(model.p.a_min * (1.0 + 1.0e-8))
        if a_lower is None
        else float(a_lower)
    )
    return conditional_fast_a(
        model,
        float(s),
        float(force),
        a_lower=lower,
        a_upper=a_dagger,
        quadrature_order=quadrature_order,
        upper_s_derivative=upper_s_derivative,
    )


def laplace_F_eff(model: TwoRowLJ, s: float, force: float) -> dict[str, float]:
    r"""Low-temperature interior-minimum approximation to ``mathcal F_eff``.

    Repository coordinates are nondimensionalized by the lattice spacing, so
    ``da`` is a dimensionless integration measure.  Up to an additive constant
    independent of ``s`` and ``force``, the one-dimensional Laplace formula is

    ``F_eff = G(a*,s;f) + kT/2 log[G_aa(a*,s;f)/(2*pi*kT)]``.
    """

    equilibrium = stable_normal_equilibrium(model, float(s), float(force))
    if equilibrium is None:
        raise ValueError("no stable normal point exists for Laplace approximation")
    generalized, _, _ = _generalized_energy_derivatives(
        model, equilibrium.a, float(s), float(force)
    )
    correction = 0.5 * model.p.kT * np.log(
        equilibrium.normal_stiffness / (2.0 * np.pi * model.p.kT)
    )
    return {
        "a_star": equilibrium.a,
        "G_star": float(generalized),
        "normal_stiffness": equilibrium.normal_stiffness,
        "F_eff_laplace": float(generalized + correction),
        "correction": float(correction),
    }


def effective_profile(
    model: TwoRowLJ,
    s_values: np.ndarray,
    force: float,
    *,
    a_lower: float,
    a_upper: float,
    quadrature_order: int = 96,
) -> dict[str, np.ndarray]:
    """Vector profile of fixed-domain effective free energy and moments."""

    s_values = np.asarray(s_values, dtype=float).reshape(-1)
    rows = [
        conditional_fast_a(
            model,
            float(s),
            float(force),
            a_lower=float(a_lower),
            a_upper=float(a_upper),
            quadrature_order=quadrature_order,
        )
        for s in s_values
    ]
    return {
        "s": s_values,
        "F_eff": np.asarray([row.F_eff for row in rows], dtype=float),
        "mean_a": np.asarray([row.mean_a for row in rows], dtype=float),
        "dF_eff_ds": np.asarray([row.dF_eff_ds for row in rows], dtype=float),
        "mean_dG_ds": np.asarray([row.mean_dG_ds for row in rows], dtype=float),
        "boundary_correction_s": np.zeros(s_values.shape, dtype=float),
    }


def truncated_bound_profile(
    model: TwoRowLJ,
    s_values: np.ndarray,
    force: float,
    *,
    a_lower: float,
    a_upper: float,
    quadrature_order: int = 96,
) -> dict[str, np.ndarray]:
    """Profile from instantaneous bound-basin truncated Gibbs conditionals.

    This metastable object is distinct from the absorbing fast-process QSD and
    therefore carries no validated crack probability.  It is exposed only for
    comparison in regimes where its QSD discrepancy is acceptably small.
    """

    s_values = np.asarray(s_values, dtype=float).reshape(-1)
    if s_values.size < 3:
        raise ValueError("truncated bound profile requires at least three s points")
    if force > 0.0:
        _, saddles, bound = model.opening_saddle_batch(
            s_values[:, None], float(force)
        )
        saddles = saddles[:, 0]
        bound = bound[:, 0]
        if not np.all(bound):
            raise ValueError("opening stable basin is absent for part of the s grid")
        upper_limits = np.minimum(saddles, float(a_upper))
    else:
        upper_limits = np.full(s_values.shape, float(a_upper))
    if np.any(upper_limits <= float(a_lower)):
        raise ValueError("truncated normal basin contains an empty integration domain")
    rows = [
        conditional_fast_a(
            model,
            float(s),
            float(force),
            a_lower=float(a_lower),
            a_upper=float(upper),
            quadrature_order=quadrature_order,
        )
        for s, upper in zip(s_values, upper_limits)
    ]
    F_eff = np.asarray([row.F_eff for row in rows], dtype=float)
    mean_dG_ds = np.asarray([row.mean_dG_ds for row in rows], dtype=float)
    dF_eff_ds = np.gradient(F_eff, s_values, edge_order=2)
    return {
        "s": s_values,
        "F_eff": F_eff,
        "mean_a": np.asarray([row.mean_a for row in rows], dtype=float),
        # The SG drift uses differences of F_eff itself.  These arrays make the
        # moving-limit correction visible rather than reporting the fixed-limit
        # identity where it does not apply.
        "dF_eff_ds": dF_eff_ds,
        "mean_dG_ds": mean_dG_ds,
        "boundary_correction_s": dF_eff_ds - mean_dG_ds,
        "a_upper": upper_limits,
    }


@dataclass(frozen=True)
class FastAQSD:
    """Discrete QSD of the killed fast normal SG generator."""

    a: np.ndarray
    density: np.ndarray
    truncated_gibbs_density: np.ndarray
    escape_rate: float
    saddle: float
    l1_from_truncated_gibbs: float


def fast_a_qsd(
    model: TwoRowLJ,
    s: float,
    force: float,
    *,
    a_lower: float,
    a_upper: float,
    n_a: int = 121,
) -> FastAQSD:
    """Compute the dominant QSD of the discrete absorbing fast-``a`` process."""

    roots = _normal_stationary_roots(model, float(s), float(force))
    unstable = [a for a, stiffness in roots if stiffness < 0.0]
    if not unstable:
        raise ValueError("QSD diagnostic requires a finite opening saddle")
    saddle = float(max(unstable))
    edges = np.linspace(float(a_lower), float(a_upper), int(n_a) + 1)
    a = 0.5 * (edges[:-1] + edges[1:])
    da = float(edges[1] - edges[0])
    generalized, _, _ = _generalized_energy_derivatives(
        model, a, float(s), float(force)
    )
    beta = 1.0 / model.p.kT
    diffusivity = model.p.mobility_a * model.p.kT
    psi = beta * (generalized[1:] - generalized[:-1])
    rate_lr = (diffusivity / da**2) * _bernoulli(psi)
    rate_rl = (diffusivity / da**2) * _bernoulli(-psi)
    generator = np.zeros((a.size, a.size), dtype=float)
    for index in range(a.size - 1):
        generator[index + 1, index] += rate_lr[index]
        generator[index, index] -= rate_lr[index]
        generator[index, index + 1] += rate_rl[index]
        generator[index + 1, index + 1] -= rate_rl[index]
    intact = a < saddle
    killed = generator[np.ix_(intact, intact)]
    eigenvalues, eigenvectors = np.linalg.eig(killed)
    principal = int(np.argmax(eigenvalues.real))
    eigenvalue = complex(eigenvalues[principal])
    vector = np.asarray(eigenvectors[:, principal].real, dtype=float)
    if float(np.sum(vector)) < 0.0:
        vector = -vector
    vector = np.maximum(vector, 0.0)
    density = vector / (float(np.sum(vector)) * da)

    bound_energy = generalized[intact]
    reference = float(np.min(bound_energy))
    truncated = np.exp(
        np.clip(-(bound_energy - reference) / model.p.kT, -745.0, 0.0)
    )
    truncated /= float(np.sum(truncated) * da)
    l1_error = float(np.sum(np.abs(density - truncated)) * da)
    return FastAQSD(
        a=a[intact],
        density=density,
        truncated_gibbs_density=truncated,
        escape_rate=float(max(0.0, -eigenvalue.real)),
        saddle=saddle,
        l1_from_truncated_gibbs=l1_error,
    )


@dataclass(frozen=True)
class ReducedFastAParams:
    """Numerical controls for the experimental slow-configurational PDE."""

    n_s: int = 91
    s_wells: int = 3
    a_lower_factor: float = 0.82
    a_upper: float = 1.60
    quadrature_order: int = 96
    max_dt: float = 1.0e-3
    record_interval: float = 1.0e-3
    negative_tolerance: float = 2.0e-12
    normal_domain: str = "fixed"


def _reduced_s_grid(model: TwoRowLJ, params: ReducedFastAParams) -> tuple[np.ndarray, float]:
    if params.n_s < 5:
        raise ValueError("n_s must be at least 5")
    if params.s_wells < 1 or params.s_wells % 2 == 0:
        raise ValueError("s_wells must be a positive odd integer")
    half_span = 0.5 * params.s_wells * model.p.b
    edges = np.linspace(-half_span, half_span, params.n_s + 1)
    return 0.5 * (edges[:-1] + edges[1:]), float(edges[1] - edges[0])


def _profile_for_domain(
    model: TwoRowLJ,
    s: np.ndarray,
    force: float,
    a_lower: float,
    params: ReducedFastAParams,
) -> dict[str, np.ndarray]:
    if params.normal_domain == "fixed":
        return effective_profile(
            model,
            s,
            force,
            a_lower=a_lower,
            a_upper=params.a_upper,
            quadrature_order=params.quadrature_order,
        )
    if params.normal_domain == "truncated_bound":
        return truncated_bound_profile(
            model,
            s,
            force,
            a_lower=a_lower,
            a_upper=params.a_upper,
            quadrature_order=params.quadrature_order,
        )
    raise ValueError("normal_domain must be 'fixed' or 'truncated_bound'")


def _sg_generator_s(F_eff: np.ndarray, model: TwoRowLJ, ds: float) -> sparse.csr_matrix:
    beta = 1.0 / model.p.kT
    diffusivity = model.p.mobility_s * model.p.kT
    psi = beta * (F_eff[1:] - F_eff[:-1])
    rate_lr = (diffusivity / ds**2) * _bernoulli(psi)
    rate_rl = (diffusivity / ds**2) * _bernoulli(-psi)
    index = np.arange(F_eff.size, dtype=np.int64)
    left = index[:-1]
    right = index[1:]
    rows = np.concatenate((right, left, left, right))
    columns = np.concatenate((left, left, right, right))
    values = np.concatenate((rate_lr, -rate_lr, rate_rl, -rate_rl))
    return sparse.coo_matrix(
        (values, (rows, columns)), shape=(F_eff.size, F_eff.size)
    ).tocsr()


def _reduced_observables(
    density_s: np.ndarray,
    s: np.ndarray,
    ds: float,
    mean_a: np.ndarray,
    model: TwoRowLJ,
    force: float,
) -> dict[str, float]:
    mass = float(np.sum(density_s) * ds)
    if mass <= 0.0:
        raise FloatingPointError("reduced configurational density has zero mass")
    weights = density_s * ds / mass
    nwell = model.well_index(s)
    xi = s - model.p.b * nwell
    normal = float(np.sum(weights * (mean_a - model.a0) / model.a0))
    intrawell = float(
        np.sum(weights * model.p.chi_axial_projection * xi / model.a0)
    )
    plastic = float(
        np.sum(
            weights
            * model.p.chi_axial_projection
            * model.p.b
            * nwell
            / model.a0
        )
    )
    return {
        "force": float(force),
        "probability_mass": mass,
        "mean_s": float(np.sum(weights * s)),
        "normal_strain": normal,
        "intrawell_strain": intrawell,
        "plastic_strain": plastic,
        "strain": normal + intrawell + plastic,
        "mean_well_index": float(np.sum(weights * nwell)),
        "plastic_well_activity": float(np.sum(weights * np.abs(nwell))),
    }


def reduced_equilibrium(
    model: TwoRowLJ,
    force: float,
    *,
    params: ReducedFastAParams = ReducedFastAParams(),
    principal_well_only: bool = True,
) -> dict[str, np.ndarray | float]:
    """Return the stationary fixed-domain reduced Gibbs state at one force."""

    s, ds = _reduced_s_grid(model, params)
    a_lower = max(model.p.a_min * 1.001, params.a_lower_factor * model.a0)
    profile = _profile_for_domain(model, s, float(force), a_lower, params)
    mask = np.abs(s) < 0.5 * model.p.b if principal_well_only else np.ones_like(
        s, dtype=bool
    )
    reference = float(np.min(profile["F_eff"][mask]))
    density_s = np.zeros_like(s)
    density_s[mask] = np.exp(
        np.clip(-(profile["F_eff"][mask] - reference) / model.p.kT, -745.0, 0.0)
    )
    density_s /= float(np.sum(density_s) * ds)
    values = _reduced_observables(
        density_s, s, ds, profile["mean_a"], model, float(force)
    )
    return {**values, "s": s, "ds": ds, "density_s": density_s, **profile}


def run_reduced_fast_a(
    model: TwoRowLJ,
    *,
    params: ReducedFastAParams = ReducedFastAParams(),
    load: CyclicLoad2D = CyclicLoad2D(),
    preload_force: float = 0.0,
    record_callback: Callable[[dict[str, float]], None] | None = None,
) -> dict[str, np.ndarray | str]:
    """Evolve the experimental fixed-normal-domain slow-``s`` PDE.

    This solver deliberately has no reduced crack probability.  It is valid
    only for strain/configurational comparisons where full-model opening loss
    is negligible.
    """

    if model.p.n_cells != 1:
        raise ValueError("reduced fast-a solver requires n_cells=1")
    if params.max_dt <= 0.0 or params.record_interval <= 0.0:
        raise ValueError("reduced time controls must be positive")
    if params.normal_domain not in {"fixed", "truncated_bound"}:
        raise ValueError("normal_domain must be 'fixed' or 'truncated_bound'")
    s, ds = _reduced_s_grid(model, params)
    a_lower = max(model.p.a_min * 1.001, params.a_lower_factor * model.a0)
    if params.a_upper <= a_lower:
        raise ValueError("a_upper must exceed the reduced normal lower bound")
    initial_profile = _profile_for_domain(
        model, s, float(preload_force), a_lower, params
    )
    principal = np.abs(s) < 0.5 * model.p.b
    reference = float(np.min(initial_profile["F_eff"][principal]))
    density_s = np.zeros_like(s)
    density_s[principal] = np.exp(
        np.clip(
            -(initial_profile["F_eff"][principal] - reference) / model.p.kT,
            -745.0,
            0.0,
        )
    )
    density_s /= float(np.sum(density_s) * ds)

    record_names = (
        "time",
        "force",
        "probability_mass",
        "mass_balance_residual",
        "minimum_density",
        "mean_s",
        "normal_strain",
        "intrawell_strain",
        "plastic_strain",
        "strain",
        "mean_well_index",
        "plastic_well_activity",
    )
    records = {name: [] for name in record_names}
    time = 0.0
    next_record = 0.0
    minimum_density = float(np.min(density_s))

    def append_record(now: float, force: float, profile: dict[str, np.ndarray]) -> None:
        values = _reduced_observables(
            density_s, s, ds, profile["mean_a"], model, force
        )
        snapshot = {
            "time": float(now),
            **values,
            "mass_balance_residual": float(values["probability_mass"] - 1.0),
            "minimum_density": float(minimum_density),
        }
        for name in records:
            records[name].append(float(snapshot[name]))
        if record_callback is not None:
            record_callback(dict(snapshot))

    profile = initial_profile
    while True:
        force = load.value(time)
        if time + 1.0e-14 >= next_record:
            profile = _profile_for_domain(model, s, force, a_lower, params)
            append_record(time, force, profile)
            next_record += params.record_interval
        if time >= load.duration - 1.0e-14:
            break
        dt = min(params.max_dt, load.duration - time)
        next_time = time + dt
        next_force = load.value(next_time)
        next_profile = _profile_for_domain(
            model, s, next_force, a_lower, params
        )
        generator = _sg_generator_s(next_profile["F_eff"], model, ds)
        system = sparse.eye(s.size, format="csr") - dt * generator
        updated = np.asarray(spsolve(system, density_s), dtype=float)
        raw_minimum = float(np.min(updated))
        minimum_density = min(minimum_density, raw_minimum)
        if raw_minimum < -params.negative_tolerance:
            raise FloatingPointError(
                f"reduced density lost positivity ({raw_minimum:.3e})"
            )
        density_s = np.maximum(updated, 0.0)
        time = next_time
        profile = next_profile

    if not records["time"] or abs(records["time"][-1] - time) > 1.0e-12:
        append_record(time, load.value(time), profile)
    return {
        **{name: np.asarray(values, dtype=float) for name, values in records.items()},
        "s": s,
        "ds": np.asarray(ds),
        "density_s": density_s,
        "status": "experimental strain/configurational reduction; no crack model",
        "normal_domain": params.normal_domain,
        "effective_free_energy_symbol": "mathcal{F}_eff",
    }
