"""Deterministic, gauge-aware calibration of the analytic LJ--EAM hybrid.

The mapping is deliberately narrow: a uniform FCC (111) interface cell with
one atom per plane is projected onto the N=1 two-row energy surface.  It does
not turn the reduced model into a full three-dimensional aluminum potential.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable

import numpy as np
from scipy.optimize import least_squares, minimize_scalar
from scipy.special import zeta

from .analytic_lj_eam import (
    AnalyticLJEAM,
    SquareRootEmbedding,
    SquareRootLinearEmbedding,
)
from .exponential_density import (
    ExponentialDensityParams,
    exponential_cross_density,
    same_row_exponential_density,
)
from .model import ModelParams


EV_J = 1.602176634e-19
ANGSTROM_M = 1.0e-10
GPA_PA = 1.0e9


@dataclass(frozen=True)
class FCC111Geometry:
    """Microscopic FCC (111) geometry; unrelated to statistical area ``A_c``."""

    lattice_angstrom: float
    burgers_angstrom: float
    plane_spacing_angstrom: float
    atomic_cell_area_angstrom2: float


@dataclass(frozen=True)
class AluminumTargets:
    """Reduced targets and their declared model-discrepancy scales."""

    a0_over_b: float
    cohesive_energy_ev_per_atom: float
    normal_curvature_ev: float
    registry_barrier_ev_per_cell: float
    work_separation_ev_per_cell: float
    scales: np.ndarray
    names: tuple[str, ...] = (
        "a0_over_b",
        "cohesive_energy_ev_per_atom",
        "normal_curvature_ev",
        "registry_barrier_ev_per_cell",
        "work_separation_ev_per_cell",
    )

    def values(self) -> np.ndarray:
        return np.array(
            [
                self.a0_over_b,
                self.cohesive_energy_ev_per_atom,
                self.normal_curvature_ev,
                self.registry_barrier_ev_per_cell,
                self.work_separation_ev_per_cell,
            ],
            dtype=float,
        )


@dataclass(frozen=True)
class EffectiveHybridParameters:
    """Identifiable square-root parameters after fixing the density gauge.

    ``embedding_scale_ev`` is ``A`` after fixing ``rho_ref`` to the density at
    ideal FCC (111) spacing. ``density_decay_b`` is ``beta_rho*b/r_e``. Only
    these two combinations of the five raw density/embedding parameters enter
    the square-root embedding energy after this reference-density gauge.
    """

    epsilon_lj_ev: float
    sigma_lj_over_b: float
    density_decay_b: float
    embedding_scale_ev: float

    def as_array(self) -> np.ndarray:
        return np.array(
            [
                self.epsilon_lj_ev,
                self.sigma_lj_over_b,
                self.density_decay_b,
                self.embedding_scale_ev,
            ],
            dtype=float,
        )


@dataclass(frozen=True)
class ReducedObservables:
    a0_over_b: float
    cohesive_energy_ev_per_atom: float
    normal_curvature_ev: float
    registry_barrier_ev_per_cell: float
    work_separation_ev_per_cell: float
    min_hessian_eigenvalue_ev: float
    relaxed_kappa_ev: float

    def calibration_values(self) -> np.ndarray:
        return np.array(
            [
                self.a0_over_b,
                self.cohesive_energy_ev_per_atom,
                self.normal_curvature_ev,
                self.registry_barrier_ev_per_cell,
                self.work_separation_ev_per_cell,
            ],
            dtype=float,
        )


@dataclass(frozen=True)
class CalibrationFit:
    parameters: EffectiveHybridParameters | ExtendedHybridParameters
    observables: ReducedObservables
    normalized_residuals: np.ndarray
    sensitivity: np.ndarray
    singular_values: np.ndarray
    rank: int
    condition_number: float
    cost: float
    optimality: float
    nfev: int


@dataclass(frozen=True)
class ExtendedHybridParameters:
    epsilon_lj_ev: float
    sigma_lj_over_b: float
    density_decay_b: float
    embedding_scale_ev: float
    linear_ev: float

    def as_array(self) -> np.ndarray:
        return np.array(
            [self.epsilon_lj_ev, self.sigma_lj_over_b, self.density_decay_b,
             self.embedding_scale_ev, self.linear_ev],
            dtype=float,
        )


def fcc_111_geometry(lattice_angstrom: float) -> FCC111Geometry:
    """Return exact ideal-FCC geometry for one atom in a (111) primitive cell."""

    lattice = float(lattice_angstrom)
    if not np.isfinite(lattice) or lattice <= 0.0:
        raise ValueError("lattice parameter must be finite and positive")
    return FCC111Geometry(
        lattice_angstrom=lattice,
        burgers_angstrom=lattice / math.sqrt(2.0),
        plane_spacing_angstrom=lattice / math.sqrt(3.0),
        atomic_cell_area_angstrom2=math.sqrt(3.0) * lattice**2 / 4.0,
    )


def energy_density_to_cell_ev(value_j_m2: float, area_angstrom2: float) -> float:
    """Convert interfacial J/m^2 to eV per crystallographic interface cell."""

    return float(value_j_m2) * float(area_angstrom2) * ANGSTROM_M**2 / EV_J


def cell_ev_to_energy_density(value_ev: float, area_angstrom2: float) -> float:
    """Inverse of :func:`energy_density_to_cell_ev`."""

    return float(value_ev) * EV_J / (float(area_angstrom2) * ANGSTROM_M**2)


def longitudinal_modulus_111_gpa(c11: float, c12: float, c44: float) -> float:
    """Cubic-crystal longitudinal modulus for [111] strain at fixed lateral cell."""

    return (float(c11) + 2.0 * float(c12) + 4.0 * float(c44)) / 3.0


def normal_curvature_target_ev(
    geometry: FCC111Geometry, c111_gpa: float
) -> float:
    """Map fixed-lateral [111] elasticity to d^2 W/d(a/b)^2 in eV."""

    curvature_j_m2 = (
        float(c111_gpa)
        * GPA_PA
        * geometry.atomic_cell_area_angstrom2
        * ANGSTROM_M**2
        / (geometry.plane_spacing_angstrom * ANGSTROM_M)
    )
    curvature_ev_angstrom2 = curvature_j_m2 * ANGSTROM_M**2 / EV_J
    return curvature_ev_angstrom2 * geometry.burgers_angstrom**2


def mishin_lu_wei_targets() -> tuple[FCC111Geometry, AluminumTargets]:
    """Traceable mixed-reference target set used for identifiability assessment.

    Geometry, cohesion, and elastic constants use the 0 K Mishin-1999 Al EAM
    validation values.  The direct <110> GSF saddle uses Lu et al. DFT, and
    clean (111) surface energy uses Wei et al. DFT.  The mixture is explicit,
    not presented as a single internally consistent first-principles dataset.
    """

    geometry = fcc_111_geometry(4.05)
    c111 = longitudinal_modulus_111_gpa(114.0, 62.0, 32.0)
    registry = energy_density_to_cell_ev(0.250, geometry.atomic_cell_area_angstrom2)
    work_sep = energy_density_to_cell_ev(
        2.0 * 1.06, geometry.atomic_cell_area_angstrom2
    )
    values = np.array(
        [
            geometry.plane_spacing_angstrom / geometry.burgers_angstrom,
            3.36,
            normal_curvature_target_ev(geometry, c111),
            registry,
            work_sep,
        ]
    )
    # These are declared model-discrepancy/normalization scales, not claimed
    # measurement standard deviations: 0.5%, 2%, 5%, 10%, and 10%.
    scales = values * np.array([0.005, 0.02, 0.05, 0.10, 0.10])
    return geometry, AluminumTargets(*values, scales=scales)


def raw_to_effective_embedding(
    *, amplitude: float, rho0: float, beta_rho: float, r_e_over_b: float,
    rho_ref: float,
) -> tuple[float, float]:
    """Return the two identifiable square-root embedding combinations."""

    if min(rho0, beta_rho, r_e_over_b, rho_ref) <= 0.0 or amplitude < 0.0:
        raise ValueError("raw embedding scales must satisfy their positivity bounds")
    decay = beta_rho / r_e_over_b
    scale = amplitude * math.sqrt(rho0 * math.exp(beta_rho) / rho_ref)
    return float(decay), float(scale)


def model_from_effective_parameters(
    parameters: EffectiveHybridParameters,
    *,
    chi: float = 0.0,
    kT_ev: float = 0.02,
) -> AnalyticLJEAM:
    """Construct the gauge-fixed model with b=r_e=rho0=1 and x_ref=1."""

    values = parameters.as_array()
    if np.any(~np.isfinite(values)) or np.any(values <= 0.0):
        raise ValueError("all effective hybrid parameters must be finite and positive")
    beta = parameters.density_decay_b
    density_parameters = ExponentialDensityParams(
        rho0=1.0,
        beta_rho=beta,
        r_e=1.0,
        tol=2.0e-13,
        max_modes=128,
    )
    reference_a = math.sqrt(2.0 / 3.0)
    rho_ref = same_row_exponential_density(b=1.0, params=density_parameters)
    rho_ref += float(
        exponential_cross_density(
            reference_a, 0.0, b=1.0, params=density_parameters
        ).value
    )
    return AnalyticLJEAM(
        ModelParams(
            n_cells=1,
            b=1.0,
            epsilon=parameters.epsilon_lj_ev,
            sigma_lj=parameters.sigma_lj_over_b,
            chi_axial_projection=float(chi),
            kT=float(kT_ev),
            mobility_a=1.0,
            mobility_s=0.05,
            a_min=0.25,
            a_max=4.0,
        ),
        density_params=density_parameters,
        embedding=SquareRootEmbedding(
            amplitude=parameters.embedding_scale_ev, rho_ref=rho_ref
        ),
    )


def model_from_extended_parameters(
    parameters: ExtendedHybridParameters,
    *,
    chi: float = 0.0,
    kT_ev: float = 0.02,
) -> AnalyticLJEAM:
    """Construct the one-extra-parameter linear embedding experiment."""

    values = parameters.as_array()
    if np.any(~np.isfinite(values)) or np.any(values <= 0.0):
        raise ValueError("all extended hybrid parameters must be finite and positive")
    density_parameters = ExponentialDensityParams(
        rho0=1.0,
        beta_rho=parameters.density_decay_b,
        r_e=1.0,
        tol=2.0e-13,
        max_modes=128,
    )
    reference_a = math.sqrt(2.0 / 3.0)
    rho_ref = same_row_exponential_density(b=1.0, params=density_parameters)
    rho_ref += float(
        exponential_cross_density(
            reference_a, 0.0, b=1.0, params=density_parameters
        ).value
    )
    return AnalyticLJEAM(
        ModelParams(
            n_cells=1,
            b=1.0,
            epsilon=parameters.epsilon_lj_ev,
            sigma_lj=parameters.sigma_lj_over_b,
            chi_axial_projection=float(chi),
            kT=float(kT_ev),
            mobility_a=1.0,
            mobility_s=0.05,
            a_min=0.25,
            a_max=4.0,
        ),
        density_params=density_parameters,
        embedding=SquareRootLinearEmbedding(
            amplitude=parameters.embedding_scale_ev,
            linear=parameters.linear_ev,
            rho_ref=rho_ref,
        ),
    )


def _relaxed_energy(model: AnalyticLJEAM, s: float) -> tuple[float, float]:
    result = minimize_scalar(
        lambda a: model.local_energy(float(a), float(s)),
        bounds=(max(0.35, 0.45 * model.a0), 2.5),
        method="bounded",
        options={"xatol": 2.0e-11},
    )
    if not result.success:
        raise ValueError("normal branch minimization failed")
    curvature = float(model.local_hessian(float(result.x), float(s))[0, 0])
    if curvature <= 0.0:
        raise ValueError("normal branch is not locally stable")
    return float(result.fun), float(result.x)


def evaluate_reduced_observables(
    parameters: EffectiveHybridParameters,
) -> ReducedObservables:
    """Evaluate the five mapped static observables and stability diagnostics."""

    return evaluate_model_observables(model_from_effective_parameters(parameters))


def evaluate_model_observables(model: AnalyticLJEAM) -> ReducedObservables:
    """Evaluate mapped observables for either supported analytic embedding."""

    hessian = model.local_hessian(model.a0, 0.0)
    eigenvalues = np.linalg.eigvalsh(hessian)
    if np.min(eigenvalues) <= 0.0:
        raise ValueError("pristine hybrid Hessian is not positive definite")

    minimum_energy, _ = _relaxed_energy(model, 0.0)
    saddle_energy, _ = _relaxed_energy(model, 0.5 * model.p.b)
    registry_barrier = saddle_energy - minimum_energy

    parallel_density = same_row_exponential_density(
        b=model.p.b, params=model.density_params
    )
    separated_energy = float(2.0 * model.embedding.value(parallel_density))
    work_separation = separated_energy - minimum_energy

    same_row_pair_per_row = 4.0 * model.p.epsilon * (
        model.p.sigma_lj**12 * float(zeta(12.0, 1.0))
        - model.p.sigma_lj**6 * float(zeta(6.0, 1.0))
    )
    full_two_row_cell_energy = minimum_energy + 2.0 * same_row_pair_per_row
    atomized_cell_energy = float(2.0 * model.embedding.value(0.0))
    cohesive_per_atom = 0.5 * (atomized_cell_energy - full_two_row_cell_energy)

    if registry_barrier <= 0.0 or work_separation <= 0.0:
        raise ValueError("candidate lacks a positive registry or opening barrier")
    return ReducedObservables(
        a0_over_b=float(model.a0),
        cohesive_energy_ev_per_atom=float(cohesive_per_atom),
        normal_curvature_ev=float(hessian[0, 0]),
        registry_barrier_ev_per_cell=float(registry_barrier),
        work_separation_ev_per_cell=float(work_separation),
        min_hessian_eigenvalue_ev=float(np.min(eigenvalues)),
        relaxed_kappa_ev=float(model.sigma_over_E_force_scale()),
    )


def evaluate_extended_observables(
    parameters: ExtendedHybridParameters,
) -> ReducedObservables:
    return evaluate_model_observables(model_from_extended_parameters(parameters))


def _parameters_from_log(log_parameters: Iterable[float]) -> EffectiveHybridParameters:
    values = np.exp(np.asarray(tuple(log_parameters), dtype=float))
    return EffectiveHybridParameters(*map(float, values))


def _extended_from_log(log_parameters: Iterable[float]) -> ExtendedHybridParameters:
    values = np.exp(np.asarray(tuple(log_parameters), dtype=float))
    return ExtendedHybridParameters(*map(float, values))


def normalized_residuals(
    log_parameters: Iterable[float], targets: AluminumTargets
) -> np.ndarray:
    observables = evaluate_reduced_observables(_parameters_from_log(log_parameters))
    return (observables.calibration_values() - targets.values()) / targets.scales


def extended_normalized_residuals(
    log_parameters: Iterable[float], targets: AluminumTargets
) -> np.ndarray:
    observables = evaluate_extended_observables(_extended_from_log(log_parameters))
    return (observables.calibration_values() - targets.values()) / targets.scales


def _reject_infeasible_extended(
    log_parameters: Iterable[float], targets: AluminumTargets
) -> np.ndarray:
    """Reject unstable/no-basin optimizer proposals without clipping physics."""

    try:
        return extended_normalized_residuals(log_parameters, targets)
    except (ValueError, FloatingPointError):
        # The candidate is outside the feasible mechanical domain. No energy,
        # Hessian, or observable is clipped; the optimizer is simply told that
        # this trial is inadmissible. The returned final point is reevaluated
        # through the strict function before a CalibrationFit is constructed.
        return np.full(len(targets.names), 1.0e6, dtype=float)


def sensitivity_matrix(
    log_parameters: Iterable[float],
    targets: AluminumTargets,
    *,
    step: float = 2.0e-5,
) -> np.ndarray:
    """Central-difference normalized Jacobian with respect to log parameters."""

    center = np.asarray(tuple(log_parameters), dtype=float)
    matrix = np.empty((len(targets.names), center.size), dtype=float)
    for column in range(center.size):
        plus = center.copy()
        minus = center.copy()
        plus[column] += step
        minus[column] -= step
        matrix[:, column] = (
            normalized_residuals(plus, targets)
            - normalized_residuals(minus, targets)
        ) / (2.0 * step)
    return matrix


def extended_sensitivity_matrix(
    log_parameters: Iterable[float],
    targets: AluminumTargets,
    *,
    step: float = 2.0e-5,
) -> np.ndarray:
    center = np.asarray(tuple(log_parameters), dtype=float)
    matrix = np.empty((len(targets.names), center.size), dtype=float)
    for column in range(center.size):
        plus = center.copy()
        minus = center.copy()
        plus[column] += step
        minus[column] -= step
        matrix[:, column] = (
            extended_normalized_residuals(plus, targets)
            - extended_normalized_residuals(minus, targets)
        ) / (2.0 * step)
    return matrix


def calibrate_square_root_hybrid(
    targets: AluminumTargets | None = None,
    *,
    initial: EffectiveHybridParameters = EffectiveHybridParameters(
        0.20, 0.82, 3.0, 0.12
    ),
    max_nfev: int = 300,
) -> CalibrationFit:
    """Reproducible bounded least-squares calibration of the identifiable model."""

    if targets is None:
        _, targets = mishin_lu_wei_targets()
    lower = np.log([1.0e-5, 0.35, 0.20, 1.0e-5])
    upper = np.log([5.0, 1.35, 12.0, 5.0])
    solution = least_squares(
        lambda values: normalized_residuals(values, targets),
        np.log(initial.as_array()),
        bounds=(lower, upper),
        method="trf",
        x_scale="jac",
        ftol=2.0e-11,
        xtol=2.0e-11,
        gtol=2.0e-11,
        max_nfev=int(max_nfev),
    )
    if not solution.success:
        raise RuntimeError(f"deterministic calibration failed: {solution.message}")
    parameters = _parameters_from_log(solution.x)
    observables = evaluate_reduced_observables(parameters)
    residual = (observables.calibration_values() - targets.values()) / targets.scales
    jacobian = sensitivity_matrix(solution.x, targets)
    singular = np.linalg.svd(jacobian, compute_uv=False)
    tolerance = np.finfo(float).eps * max(jacobian.shape) * singular[0]
    rank = int(np.sum(singular > tolerance))
    condition = float(singular[0] / singular[-1])
    return CalibrationFit(
        parameters=parameters,
        observables=observables,
        normalized_residuals=residual,
        sensitivity=jacobian,
        singular_values=singular,
        rank=rank,
        condition_number=condition,
        cost=float(np.sum(residual**2)),
        optimality=float(solution.optimality),
        nfev=int(solution.nfev),
    )


def calibrate_linear_hybrid(
    targets: AluminumTargets | None = None,
    *,
    initial: ExtendedHybridParameters = ExtendedHybridParameters(
        0.15, 0.82, 3.0, 0.30, 0.05
    ),
    max_nfev: int = 500,
) -> CalibrationFit:
    """Fit the minimum fixed-LJ-gauge linear extension deterministically."""

    if targets is None:
        _, targets = mishin_lu_wei_targets()
    lower = np.log([1.0e-5, 0.35, 0.20, 1.0e-5, 1.0e-6])
    upper = np.log([5.0, 1.35, 12.0, 8.0, 8.0])
    solution = least_squares(
        lambda values: _reject_infeasible_extended(values, targets),
        np.log(initial.as_array()),
        bounds=(lower, upper),
        method="trf",
        x_scale="jac",
        ftol=2.0e-11,
        xtol=2.0e-11,
        gtol=2.0e-11,
        max_nfev=int(max_nfev),
    )
    if not solution.success:
        raise RuntimeError(f"extended deterministic calibration failed: {solution.message}")
    parameters = _extended_from_log(solution.x)
    observables = evaluate_extended_observables(parameters)
    residual = (observables.calibration_values() - targets.values()) / targets.scales
    jacobian = extended_sensitivity_matrix(solution.x, targets)
    singular = np.linalg.svd(jacobian, compute_uv=False)
    tolerance = np.finfo(float).eps * max(jacobian.shape) * singular[0]
    rank = int(np.sum(singular > tolerance))
    condition = float(singular[0] / singular[-1])
    return CalibrationFit(
        parameters=parameters,
        observables=observables,
        normalized_residuals=residual,
        sensitivity=jacobian,
        singular_values=singular,
        rank=rank,
        condition_number=condition,
        cost=float(np.sum(residual**2)),
        optimality=float(solution.optimality),
        nfev=int(solution.nfev),
    )


def barrier_ratio(
    configurational_barrier: float | None,
    opening_barrier: float | None,
) -> float:
    """Return Delta F_s / Delta G_open only for two positive barriers."""

    if configurational_barrier is None or opening_barrier is None:
        return math.nan
    if configurational_barrier <= 0.0 or opening_barrier <= 0.0:
        return math.nan
    return float(configurational_barrier / opening_barrier)
