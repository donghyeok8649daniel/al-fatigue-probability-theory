"""Fast/slow diagnostics for the N=1 Bessel-LJ probability model.

The routines in this module do not change the production probability model.
They separate three different reference problems:

* the linear response of the current finite-mobility two-coordinate dynamics,
* stable normal mechanical equilibrium at fixed registry ``s``, and
* finite-temperature conditional equilibrium in ``a`` at fixed ``s``.

The last of these is the actual adiabatic limit of the Smoluchowski equation as
``mobility_a / mobility_s`` tends to infinity.  The point equilibrium ``a*`` is
its low-temperature (Laplace) approximation, not an assumed replacement for
the conditional density.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import root, root_scalar

from .model import TwoRowLJ
from .probability_pde_2d import (
    Grid2D,
    energy_grid,
    observables,
    opening_intact_mask,
)


@dataclass(frozen=True)
class RelaxationSpectrum:
    """Pristine linear relaxation spectrum of ``dq/dt=-M H0 dq``."""

    hessian: np.ndarray
    mobility: np.ndarray
    symmetric_rate_matrix: np.ndarray
    eigenvalues: np.ndarray
    eigenvectors_symmetric: np.ndarray

    @property
    def lambda_slow(self) -> float:
        return float(self.eigenvalues[0])

    @property
    def lambda_fast(self) -> float:
        return float(self.eigenvalues[-1])

    @property
    def tau_slow(self) -> float:
        return 1.0 / self.lambda_slow

    @property
    def tau_fast(self) -> float:
        return 1.0 / self.lambda_fast


@dataclass(frozen=True)
class NormalEquilibrium:
    """Stable solution of ``partial G / partial a = 0`` at fixed ``s``."""

    a: float
    s: float
    force: float
    normal_stiffness: float
    residual: float


@dataclass(frozen=True)
class HarmonicResponse:
    """Complex fundamental response extracted using ``Re[z exp(i omega t)]``."""

    input_phasor: complex
    output_phasor: complex
    transfer: complex
    input_offset: float
    output_offset: float

    @property
    def magnitude(self) -> float:
        return float(abs(self.transfer))

    @property
    def phase_degrees(self) -> float:
        return float(np.degrees(np.angle(self.transfer)))


def pristine_relaxation_spectrum(model: TwoRowLJ) -> RelaxationSpectrum:
    """Return rates from the full symmetric matrix ``sqrt(M) H0 sqrt(M)``."""

    hessian = model.local_hessian(model.a0, 0.0)
    mobility = np.diag([model.p.mobility_a, model.p.mobility_s]).astype(float)
    if np.any(np.diag(mobility) <= 0.0):
        raise ValueError("coordinate mobilities must be positive")
    sqrt_mobility = np.diag(np.sqrt(np.diag(mobility)))
    rate_matrix = sqrt_mobility @ hessian @ sqrt_mobility
    eigenvalues, eigenvectors = np.linalg.eigh(rate_matrix)
    if not np.all(np.isfinite(eigenvalues)) or np.any(eigenvalues <= 0.0):
        raise FloatingPointError("pristine relaxation rates must be finite and positive")
    return RelaxationSpectrum(
        hessian=hessian,
        mobility=mobility,
        symmetric_rate_matrix=rate_matrix,
        eigenvalues=eigenvalues,
        eigenvectors_symmetric=eigenvectors,
    )


def linear_axial_transfer(model: TwoRowLJ, omega: float) -> complex:
    r"""Return ``G(omega)=epsilon_hat/(sigma_hat/E)`` for the full linear model.

    With ``c=(1,chi)^T`` and ``f_hat=kappa_axial sigma_hat/E``, the response is

    ``dq_hat=(i omega I + M H0)^(-1) M c f_hat``.

    This is a diagnostic of the current finite-mobility model.  It is not a
    physical-frequency calibration and is not an assertion that normal elastic
    relaxation should be the fatigue-memory mechanism.
    """

    if not np.isfinite(omega) or omega < 0.0:
        raise ValueError("omega must be finite and nonnegative")
    spectrum = pristine_relaxation_spectrum(model)
    c = np.array([1.0, model.p.chi_axial_projection], dtype=float)
    system = 1j * float(omega) * np.eye(2) + spectrum.mobility @ spectrum.hessian
    displacement_per_force = np.linalg.solve(system, spectrum.mobility @ c)
    return complex(
        model.sigma_over_E_force_scale()
        * (c @ displacement_per_force)
        / model.a0
    )


def model_frequency_diagnostics(model: TwoRowLJ, model_frequency: float) -> dict[str, float]:
    """Return local pristine diagnostics for cycles per model-time unit."""

    if not np.isfinite(model_frequency) or model_frequency <= 0.0:
        raise ValueError("model_frequency must be finite and positive")
    spectrum = pristine_relaxation_spectrum(model)
    period = 1.0 / float(model_frequency)
    omega = 2.0 * np.pi * float(model_frequency)
    transfer = linear_axial_transfer(model, omega)
    return {
        "model_frequency": float(model_frequency),
        "model_period": period,
        "omega": omega,
        "lambda_fast": spectrum.lambda_fast,
        "lambda_slow": spectrum.lambda_slow,
        "tau_fast": spectrum.tau_fast,
        "tau_slow": spectrum.tau_slow,
        "de_fast": omega * spectrum.tau_fast,
        "de_slow": omega * spectrum.tau_slow,
        "linear_transfer_magnitude": float(abs(transfer)),
        "linear_transfer_phase_degrees": float(np.degrees(np.angle(transfer))),
    }


def stable_normal_equilibrium(
    model: TwoRowLJ,
    s: float,
    force: float,
    *,
    previous_a: float | None = None,
    scan_points: int = 600,
) -> NormalEquilibrium | None:
    """Find the stable normal branch at fixed registry and applied force.

    All sign-changing stationary roots are checked and only roots with
    ``W_aa>0`` are eligible.  Continuation may be requested with ``previous_a``;
    otherwise the stable root nearest the pristine spacing is selected.  A
    return value of ``None`` denotes loss of a stable bound normal branch on the
    declared model domain.
    """

    if scan_points < 32:
        raise ValueError("scan_points must be at least 32")
    lower = float(model.p.a_min * (1.0 + 1.0e-10))
    upper = float(model.p.a_max)
    a_scan = np.linspace(lower, upper, scan_points)
    residual_scan = model.local_deda_array(a_scan, float(s)) - float(force)
    roots: list[NormalEquilibrium] = []
    for index in range(scan_points - 1):
        left_value = float(residual_scan[index])
        right_value = float(residual_scan[index + 1])
        if not np.isfinite(left_value) or not np.isfinite(right_value):
            continue
        if left_value == 0.0:
            candidate = float(a_scan[index])
        elif left_value * right_value > 0.0:
            continue
        else:
            solved = root_scalar(
                lambda a: model.local_deda(float(a), float(s)) - float(force),
                bracket=(float(a_scan[index]), float(a_scan[index + 1])),
                method="brentq",
                xtol=5.0e-15,
                rtol=4.0 * np.finfo(float).eps,
            )
            if not solved.converged:
                continue
            candidate = float(solved.root)
        stiffness = float(model.local_hessian(candidate, float(s))[0, 0])
        residual = abs(model.local_deda(candidate, float(s)) - float(force))
        if stiffness > 0.0 and np.isfinite(stiffness):
            roots.append(
                NormalEquilibrium(
                    a=candidate,
                    s=float(s),
                    force=float(force),
                    normal_stiffness=stiffness,
                    residual=float(residual),
                )
            )
    if not roots:
        return None
    target = model.a0 if previous_a is None else float(previous_a)
    return min(roots, key=lambda item: abs(item.a - target))


def stable_normal_equilibrium_curve(
    model: TwoRowLJ,
    s: float,
    forces: np.ndarray,
) -> tuple[NormalEquilibrium | None, ...]:
    """Follow one fixed-``s`` stable normal branch by force continuation."""

    values: list[NormalEquilibrium | None] = []
    previous_a: float | None = None
    for force in np.asarray(forces, dtype=float).reshape(-1):
        equilibrium = stable_normal_equilibrium(
            model, float(s), float(force), previous_a=previous_a
        )
        values.append(equilibrium)
        if equilibrium is not None:
            previous_a = equilibrium.a
    return tuple(values)


def stationary_total_equilibrium(
    model: TwoRowLJ,
    reduced_stress: float,
    *,
    initial_q: np.ndarray | None = None,
) -> dict[str, float]:
    """Solve the local two-coordinate stationary branch at a reduced stress."""

    force = float(model.force_from_sigma_over_E(reduced_stress))
    c = np.array([1.0, model.p.chi_axial_projection], dtype=float)
    if initial_q is None:
        initial_q = np.array([model.a0, 0.0]) + np.linalg.solve(
            model.local_hessian(model.a0, 0.0), force * c
        )
    solved = root(
        lambda q: np.array(
            [
                model.local_deda(float(q[0]), float(q[1])) - force,
                model.local_deds(float(q[0]), float(q[1]))
                - force * model.p.chi_axial_projection,
            ]
        ),
        np.asarray(initial_q, dtype=float),
        method="hybr",
        options={"xtol": 1.0e-11},
    )
    if not solved.success or np.linalg.norm(solved.fun, ord=np.inf) > 1.0e-9:
        raise RuntimeError("two-coordinate stationary equilibrium solve failed")
    a, s = (float(solved.x[0]), float(solved.x[1]))
    hessian = model.local_hessian(a, s)
    if np.min(np.linalg.eigvalsh(hessian)) <= 0.0:
        raise RuntimeError("stationary solution is not a stable local equilibrium")
    nwell = int(model.well_index(np.asarray(s)))
    xi = s - model.p.b * nwell
    normal = (a - model.a0) / model.a0
    intrawell = model.p.chi_axial_projection * xi / model.a0
    plastic = model.p.chi_axial_projection * model.p.b * nwell / model.a0
    return {
        "reduced_stress": float(reduced_stress),
        "force": force,
        "a": a,
        "s": s,
        "well_index": float(nwell),
        "normal_strain": float(normal),
        "intrawell_strain": float(intrawell),
        "plastic_strain": float(plastic),
        "strain": float(normal + intrawell + plastic),
        "stationarity_residual": float(np.linalg.norm(solved.fun, ord=np.inf)),
    }


def stationary_total_equilibrium_curve(
    model: TwoRowLJ,
    reduced_stresses: np.ndarray,
) -> dict[str, np.ndarray]:
    """Continue the stable two-coordinate equilibrium along a stress sequence."""

    rows: list[dict[str, float]] = []
    previous_q: np.ndarray | None = None
    for reduced_stress in np.asarray(reduced_stresses, dtype=float).reshape(-1):
        row = stationary_total_equilibrium(
            model, float(reduced_stress), initial_q=previous_q
        )
        rows.append(row)
        previous_q = np.array([row["a"], row["s"]], dtype=float)
    return {
        name: np.asarray([row[name] for row in rows], dtype=float)
        for name in rows[0]
    } if rows else {}


def conditional_fast_a_density(
    density: np.ndarray,
    model: TwoRowLJ,
    grid: Grid2D,
    force: float,
) -> tuple[np.ndarray, float]:
    r"""Project onto conditional Gibbs equilibrium in ``a`` at fixed ``P_s``.

    For the intact domain this constructs

    ``P_fast(a,s)=P_s(s) exp[-G(a,s)/kT] / Z_a(s)``.

    This is the finite-temperature adiabatic limit of the existing PDE.  It
    preserves every resolved ``s``-marginal column exactly whenever a bound
    normal cell exists.  Marginal mass for an ``s`` column with no intact
    normal cell is reported as ``unbound_mass`` and is not renormalized.
    """

    density = np.asarray(density, dtype=float)
    if density.shape != (grid.a.size, grid.s.size):
        raise ValueError("density shape does not match the two-dimensional grid")
    if model.p.kT <= 0.0:
        raise ValueError("conditional Gibbs equilibrium requires kT > 0")
    p_s = np.sum(density, axis=0) * grid.da
    energy = energy_grid(model, grid, float(force))
    intact = opening_intact_mask(model, grid, float(force))
    projected = np.zeros_like(density)
    unbound_mass = 0.0
    for j, marginal_density in enumerate(p_s):
        if marginal_density == 0.0:
            continue
        valid = intact[:, j]
        if not np.any(valid):
            unbound_mass += float(marginal_density * grid.ds)
            continue
        column_energy = energy[valid, j]
        exponent = -(column_energy - float(np.min(column_energy))) / model.p.kT
        weights = np.exp(np.clip(exponent, -745.0, 0.0))
        normalization = float(np.sum(weights) * grid.da)
        projected[valid, j] = float(marginal_density) * weights / normalization
    return projected, float(unbound_mass)


def point_fast_a_observables(
    density: np.ndarray,
    model: TwoRowLJ,
    grid: Grid2D,
    force: float,
) -> dict[str, float]:
    r"""Low-temperature ``a*(s,force)`` strain using the current ``P_s``.

    The configurational marginal is taken from ``density`` and is not evolved
    or re-equilibrated here.  This isolates the consequence of placing only the
    normal coordinate on its stable instantaneous mechanical branch.  It is a
    point-equilibrium reference; ``conditional_fast_a_observables`` is the
    finite-temperature adiabatic reference for the probability PDE.
    """

    density = np.asarray(density, dtype=float)
    if density.shape != (grid.a.size, grid.s.size):
        raise ValueError("density shape does not match the two-dimensional grid")
    p_s = np.sum(density, axis=0) * grid.da
    stable_a = np.full(grid.s.shape, np.nan, dtype=float)
    previous: float | None = None
    for j, s_value in enumerate(grid.s):
        equilibrium = stable_normal_equilibrium(
            model, float(s_value), float(force), previous_a=previous
        )
        if equilibrium is not None:
            stable_a[j] = equilibrium.a
            previous = equilibrium.a
    valid = np.isfinite(stable_a) & (p_s > 0.0)
    valid_mass = float(np.sum(p_s[valid]) * grid.ds)
    total_mass = float(np.sum(p_s) * grid.ds)
    if valid_mass <= 0.0:
        return {
            key: np.nan
            for key in (
                "normal_strain",
                "intrawell_strain",
                "plastic_strain",
                "strain",
                "mean_a",
            )
        } | {"unbound_marginal_mass": total_mass}
    weights = p_s[valid] * grid.ds / valid_mass
    s_valid = grid.s[valid]
    nwell = model.well_index(s_valid)
    xi = s_valid - model.p.b * nwell
    normal = float(np.sum(weights * (stable_a[valid] - model.a0) / model.a0))
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
        "normal_strain": normal,
        "intrawell_strain": intrawell,
        "plastic_strain": plastic,
        "strain": normal + intrawell + plastic,
        "mean_a": float(np.sum(weights * stable_a[valid])),
        "unbound_marginal_mass": total_mass - valid_mass,
    }


def conditional_fast_a_observables(
    density: np.ndarray,
    model: TwoRowLJ,
    grid: Grid2D,
    force: float,
) -> dict[str, float]:
    """Compare the current density with its exact conditional fast-``a`` limit."""

    projected, unbound_mass = conditional_fast_a_density(
        density, model, grid, float(force)
    )
    current = observables(density, model, grid, float(force))
    fast = observables(projected, model, grid, float(force))
    aa = grid.a[:, None]
    volume = grid.cell_volume
    current_mass = float(current["survival"])
    fast_mass = float(fast["survival"])
    current_mean_a = (
        float(np.sum(density * aa) * volume / current_mass)
        if current_mass > 0.0
        else np.nan
    )
    fast_mean_a = (
        float(np.sum(projected * aa) * volume / fast_mass)
        if fast_mass > 0.0
        else np.nan
    )
    return {
        "pde_mean_a": current_mean_a,
        "fast_a_mean_a": fast_mean_a,
        "mean_a_deviation": current_mean_a - fast_mean_a,
        "pde_normal_strain": float(current["normal_strain"]),
        "fast_a_normal_strain": float(fast["normal_strain"]),
        "pde_intrawell_strain": float(current["intrawell_strain"]),
        "fast_a_intrawell_strain": float(fast["intrawell_strain"]),
        "pde_plastic_strain": float(current["plastic_strain"]),
        "fast_a_plastic_strain": float(fast["plastic_strain"]),
        "pde_strain": float(current["strain"]),
        "fast_a_strain": float(fast["strain"]),
        "strain_deviation": float(current["strain"] - fast["strain"]),
        "unbound_marginal_mass": unbound_mass,
    }


def extract_harmonic_response(
    time: np.ndarray,
    input_signal: np.ndarray,
    output_signal: np.ndarray,
    omega: float,
    *,
    start_time: float = 0.0,
) -> HarmonicResponse:
    """Least-squares fundamental response after a caller-declared transient."""

    time = np.asarray(time, dtype=float)
    input_signal = np.asarray(input_signal, dtype=float)
    output_signal = np.asarray(output_signal, dtype=float)
    mask = time >= float(start_time)
    if np.count_nonzero(mask) < 5:
        raise ValueError("at least five samples are required for harmonic extraction")
    phase = float(omega) * time[mask]
    design = np.column_stack((np.ones(mask.sum()), np.cos(phase), np.sin(phase)))
    input_fit = np.linalg.lstsq(design, input_signal[mask], rcond=None)[0]
    output_fit = np.linalg.lstsq(design, output_signal[mask], rcond=None)[0]
    input_phasor = complex(input_fit[1], -input_fit[2])
    output_phasor = complex(output_fit[1], -output_fit[2])
    if abs(input_phasor) <= np.finfo(float).eps:
        raise ValueError("input fundamental amplitude is zero")
    return HarmonicResponse(
        input_phasor=input_phasor,
        output_phasor=output_phasor,
        transfer=output_phasor / input_phasor,
        input_offset=float(input_fit[0]),
        output_offset=float(output_fit[0]),
    )
