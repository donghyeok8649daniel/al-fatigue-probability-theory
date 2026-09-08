"""Deterministic, gauge-fixed bulk calibration on the full FCC(111) stack.

Legacy normal/biaxial audit.  Its six residuals have at most FOUR independent
observables; use full_fcc_calibration_audit for the missing homogeneous shear.
Only bulk-representable targets enter this calibration.  GSF and cleavage are
held out for the localized active-interface model; homogeneous shear/dilation
are not silently relabeled as those interface observables.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable

import numpy as np
from scipy.optimize import least_squares, root_scalar

from .aluminum_calibration import EV_J, GPA_PA, mishin_lu_wei_targets
from .analytic_lj_eam import SquareRootEmbedding, SquareRootLinearEmbedding
from .exponential_density import ExponentialDensityParams
from .fcc111_full_energy import FullFCC111StackEnergy, StackSumConfig
from .model import ModelParams


CALIBRATION_STACK_CONFIG = StackSumConfig(tol=2e-9, max_layers=256)


@dataclass(frozen=True)
class FullFCCBulkTargets:
    a0_over_b: float
    cohesive_ev_atom: float
    strain_h_aa_ev: float
    strain_h_bb_ev: float
    strain_h_ab_ev: float
    residual_scales: np.ndarray
    names: tuple[str, ...] = (
        "normal_force_at_target",
        "lateral_force_at_target",
        "cohesive_energy_ev_atom",
        "strain_h_aa_ev",
        "strain_h_bb_ev",
        "strain_h_ab_ev",
    )

    def values(self) -> np.ndarray:
        return np.array(
            [0.0, 0.0, self.cohesive_ev_atom,
             self.strain_h_aa_ev, self.strain_h_bb_ev, self.strain_h_ab_ev]
        )


@dataclass(frozen=True)
class FullFCCParameters:
    epsilon_lj_ev: float
    sigma_lj_over_b: float
    density_decay_b: float
    embedding_sqrt_ev: float
    embedding_linear_ev: float = 0.0

    def as_array(self, *, extended: bool = True) -> np.ndarray:
        values = [self.epsilon_lj_ev, self.sigma_lj_over_b,
                  self.density_decay_b, self.embedding_sqrt_ev]
        if extended:
            values.append(self.embedding_linear_ev)
        return np.asarray(values, dtype=float)


@dataclass(frozen=True)
class FullFCCBulkObservables:
    normal_force: float
    lateral_force: float
    cohesive_ev_atom: float
    strain_h_aa_ev: float
    strain_h_bb_ev: float
    strain_h_ab_ev: float

    def values(self) -> np.ndarray:
        return np.asarray(tuple(self.__dict__.values()), dtype=float)


@dataclass(frozen=True)
class FullFCCCalibrationFit:
    parameters: FullFCCParameters
    observables: FullFCCBulkObservables
    normalized_residuals: np.ndarray
    jacobian: np.ndarray
    singular_values: np.ndarray
    rank: int
    condition_number: float
    cost: float
    optimality: float
    nfev: int
    stable: bool
    equilibrium_a0: float
    minimum_strain_hessian_eigenvalue: float
    start_index: int
    extended: bool


def full_fcc_bulk_targets() -> FullFCCBulkTargets:
    geometry, _ = mishin_lu_wei_targets()
    c11, c12, c44 = 114.0, 62.0, 32.0
    volume_m3 = (
        geometry.atomic_cell_area_angstrom2
        * geometry.plane_spacing_angstrom
        * 1.0e-30
    )
    convert = GPA_PA * volume_m3 / EV_J
    h_aa = (c11 + 2.0 * c12 + 4.0 * c44) / 3.0 * convert
    h_bb = 4.0 * (c11 + 2.0 * c12 + c44) / 3.0 * convert
    h_ab = 2.0 * (c11 + 2.0 * c12 - 2.0 * c44) / 3.0 * convert
    return FullFCCBulkTargets(
        a0_over_b=math.sqrt(2.0 / 3.0),
        cohesive_ev_atom=3.36,
        strain_h_aa_ev=h_aa,
        strain_h_bb_ev=h_bb,
        strain_h_ab_ev=h_ab,
        # Documented model-discrepancy scales: forces in eV per unit strain,
        # cohesion 2%, and each elastic curvature 5%.
        residual_scales=np.array(
            [0.10, 0.10, 0.0672, 0.05 * h_aa, 0.05 * h_bb, 0.05 * h_ab]
        ),
    )


def _density_reference(parameters: FullFCCParameters, target_a: float) -> float:
    density = ExponentialDensityParams(
        rho0=1.0, beta_rho=parameters.density_decay_b, r_e=1.0,
        tol=2e-11, max_modes=128,
    )
    probe = FullFCC111StackEnergy(
        ModelParams(n_cells=1, b=1.0, epsilon=parameters.epsilon_lj_ev,
                    sigma_lj=parameters.sigma_lj_over_b, a_min=0.25, a_max=4.0),
        density_params=density,
        embedding=SquareRootEmbedding(amplitude=0.0, rho_ref=1.0),
        stack_config=CALIBRATION_STACK_CONFIG,
        find_equilibrium=False,
    )
    return probe.full_environment_density(target_a, 0.0).value


def build_full_fcc_calibration_model(
    parameters: FullFCCParameters,
    *,
    b: float = 1.0,
    find_equilibrium: bool = False,
    rho_ref: float | None = None,
) -> FullFCC111StackEnergy:
    values = parameters.as_array(extended=True)
    if np.any(~np.isfinite(values)) or np.any(values[:4] <= 0.0):
        raise ValueError("positive finite LJ, decay, and square-root scales required")
    if parameters.embedding_linear_ev < 0.0:
        raise ValueError("linear embedding coefficient must be nonnegative")
    target_a = math.sqrt(2.0 / 3.0)
    density = ExponentialDensityParams(
        rho0=1.0, beta_rho=parameters.density_decay_b, r_e=1.0,
        tol=2e-11, max_modes=128,
    )
    if rho_ref is None:
        rho_ref = _density_reference(parameters, target_a)
    if parameters.embedding_linear_ev == 0.0:
        embedding = SquareRootEmbedding(
            amplitude=parameters.embedding_sqrt_ev, rho_ref=rho_ref
        )
    else:
        embedding = SquareRootLinearEmbedding(
            amplitude=parameters.embedding_sqrt_ev,
            linear=parameters.embedding_linear_ev,
            rho_ref=rho_ref,
        )
    return FullFCC111StackEnergy(
        ModelParams(
            n_cells=1, b=float(b), epsilon=parameters.epsilon_lj_ev,
            sigma_lj=parameters.sigma_lj_over_b,
            chi_axial_projection=0.0, kT=0.02, mobility_a=1.0,
            mobility_s=0.05, a_min=0.22, a_max=4.0,
        ),
        density_params=density,
        embedding=embedding,
        stack_config=CALIBRATION_STACK_CONFIG,
        find_equilibrium=find_equilibrium,
        require_stable_equilibrium=find_equilibrium,
    )


def evaluate_full_fcc_bulk_observables(
    parameters: FullFCCParameters,
    *,
    strain_step: float = 4.0e-4,
) -> FullFCCBulkObservables:
    targets = full_fcc_bulk_targets()
    a0 = targets.a0_over_b
    step = float(strain_step)
    cache: dict[tuple[int, int], float] = {}
    rho_ref = _density_reference(parameters, a0)

    def energy(i: int, j: int) -> float:
        key = (int(i), int(j))
        if key not in cache:
            alpha, beta = i * step, j * step
            model = build_full_fcc_calibration_model(
                parameters, b=1.0 + beta, find_equilibrium=False,
                rho_ref=rho_ref,
            )
            cache[key] = model.energy(a0 * (1.0 + alpha), 0.0)
        return cache[key]

    center = energy(0, 0)
    normal_force = (energy(1, 0) - energy(-1, 0)) / (2.0 * step)
    lateral_force = (energy(0, 1) - energy(0, -1)) / (2.0 * step)
    h_aa = (energy(1, 0) - 2.0 * center + energy(-1, 0)) / step**2
    h_bb = (energy(0, 1) - 2.0 * center + energy(0, -1)) / step**2
    h_ab = (
        energy(1, 1) - energy(1, -1) - energy(-1, 1) + energy(-1, -1)
    ) / (4.0 * step**2)
    atomized = -parameters.embedding_linear_ev
    return FullFCCBulkObservables(
        normal_force=float(normal_force), lateral_force=float(lateral_force),
        cohesive_ev_atom=float(atomized - center), strain_h_aa_ev=float(h_aa),
        strain_h_bb_ev=float(h_bb), strain_h_ab_ev=float(h_ab),
    )


def _from_log(values: Iterable[float], *, extended: bool) -> FullFCCParameters:
    physical = np.exp(np.asarray(tuple(values), dtype=float))
    if extended:
        return FullFCCParameters(*map(float, physical))
    return FullFCCParameters(*map(float, physical), embedding_linear_ev=0.0)


def normalized_residuals(
    log_parameters: Iterable[float], *, extended: bool
) -> np.ndarray:
    targets = full_fcc_bulk_targets()
    prediction = evaluate_full_fcc_bulk_observables(
        _from_log(log_parameters, extended=extended)
    )
    return (prediction.values() - targets.values()) / targets.residual_scales


def sensitivity_jacobian(
    parameters: FullFCCParameters,
    *,
    extended: bool,
    relative_step: float = 2.0e-4,
) -> np.ndarray:
    """Central log-parameter sensitivity of normalized bulk observables."""

    center = np.log(parameters.as_array(extended=extended))
    columns = []
    for index in range(center.size):
        offset = np.zeros_like(center)
        offset[index] = float(relative_step)
        columns.append(
            (normalized_residuals(center+offset, extended=extended)
             - normalized_residuals(center-offset, extended=extended))
            / (2.0*relative_step)
        )
    return np.column_stack(columns)


def full_fcc_equilibrium_near_target(
    parameters: FullFCCParameters,
) -> tuple[FullFCC111StackEnergy, float]:
    """Return the stable normal root on the target FCC branch at fixed b=1."""

    model = build_full_fcc_calibration_model(parameters, find_equilibrium=False)
    grid = np.linspace(0.65, 1.05, 33)
    derivative = np.array([model.local_deda(float(a), 0.0) for a in grid])
    roots = []
    for left, right, fleft, fright in zip(
        grid[:-1], grid[1:], derivative[:-1], derivative[1:]
    ):
        if fleft <= 0.0 <= fright:
            solved = root_scalar(
                lambda a: model.local_deda(float(a), 0.0),
                bracket=(float(left), float(right)), method="brentq",
            )
            if solved.converged:
                value = float(solved.root)
                if model.local_hessian(value, 0.0)[0, 0] > 0.0:
                    roots.append(value)
    if not roots:
        raise ValueError("no stable normal root near the target FCC branch")
    target = full_fcc_bulk_targets().a0_over_b
    equilibrium = min(roots, key=lambda value: abs(value-target))
    model.a0 = equilibrium
    model.equilibrium_error = None
    return model, equilibrium


def calibrate_full_fcc(
    *,
    extended: bool,
    starts: tuple[FullFCCParameters, ...],
    max_nfev: int = 220,
    fit_observation_indices: tuple[int, ...] | None = None,
) -> tuple[FullFCCCalibrationFit, ...]:
    fits = []
    for start_index, start in enumerate(starts):
        initial = np.log(start.as_array(extended=extended))
        lower_values = np.array([1e-4, 0.45, 0.35, 1e-4] + ([1e-4] if extended else []))
        upper_values = np.array([2.0, 1.35, 8.0, 12.0] + ([12.0] if extended else []))
        indices = (
            # Remove the duplicate cubic force/curvature constraints.  This
            # legacy mode remains underidentified for five parameters.
            np.asarray((0, 2, 3, 5), dtype=int)
            if fit_observation_indices is None
            else np.asarray(fit_observation_indices, dtype=int)
        )
        solved = least_squares(
            lambda value: normalized_residuals(value, extended=extended)[indices],
            initial, bounds=(np.log(lower_values), np.log(upper_values)),
            method="trf", x_scale="jac", ftol=2e-9, xtol=2e-9, gtol=2e-9,
            max_nfev=max_nfev,
        )
        parameters = _from_log(solved.x, extended=extended)
        observables = evaluate_full_fcc_bulk_observables(parameters)
        singular = np.linalg.svd(solved.jac, compute_uv=False)
        threshold = max(solved.jac.shape) * np.finfo(float).eps * singular[0]
        numerical_rank = int(np.sum(singular > threshold))
        # Force beta = 2 force alpha, Hbb = 2 Haa + Hab.  Do not count
        # finite-difference noise along these identities as physical information.
        structural_rows = np.array([
            [1, 0, 0, 0], [2, 0, 0, 0], [0, 1, 0, 0],
            [0, 0, 1, 0], [0, 0, 2, 1], [0, 0, 0, 1],
        ], dtype=float)
        rank = min(numerical_rank,
                   int(np.linalg.matrix_rank(structural_rows[indices])))
        condition = (float(singular[0] / singular[-1])
                     if rank == solved.x.size and singular[-1] > 0 else math.inf)
        stable = False
        equilibrium = math.nan
        minimum_eigenvalue = math.nan
        try:
            model, equilibrium = full_fcc_equilibrium_near_target(parameters)
            # Strain-space curvature, including lateral/cross components.
            h = np.array([[observables.strain_h_aa_ev, observables.strain_h_ab_ev],
                          [observables.strain_h_ab_ev, observables.strain_h_bb_ev]])
            minimum_eigenvalue = float(np.min(np.linalg.eigvalsh(h)))
            stable = minimum_eigenvalue > 0.0 and abs(equilibrium-full_fcc_bulk_targets().a0_over_b) < 0.08
        except (ValueError, RuntimeError):
            pass
        fits.append(
            FullFCCCalibrationFit(
                parameters=parameters, observables=observables,
                normalized_residuals=normalized_residuals(
                    solved.x, extended=extended
                ), jacobian=np.asarray(solved.jac),
                singular_values=singular, rank=rank, condition_number=condition,
                cost=float(solved.cost), optimality=float(solved.optimality),
                nfev=int(solved.nfev), stable=stable, equilibrium_a0=equilibrium,
                minimum_strain_hessian_eigenvalue=minimum_eigenvalue,
                start_index=start_index, extended=extended,
            )
        )
    return tuple(fits)
