"""Deterministic probability-density reference solver for Theory Core v1.

This module solves the N=1, two-coordinate state q=(a,s) directly as a
Smoluchowski/Fokker--Planck PDE. It contains no random-number sampling.

The purpose is numerical validation before compressing the full N=3,
six-dimensional probability density with a sparse-grid or tensor method.

Discretization
--------------
A cell-centred conservative finite-volume method is used. Interior fluxes use
the Scharfetter--Gummel exponential fitting formula,

    J = D/h [ B(beta*DeltaG) P_L - B(-beta*DeltaG) P_R ],

where B(x)=x/(exp(x)-1), D=M*kT. This flux has the correct discrete Gibbs
stationary ratio on a fixed energy landscape and is much more robust than a
naive centred drift-diffusion stencil in steep LJ regions.

Crack first passage
-------------------
The mechanically defined local opening saddle from ``TwoRowLJ`` defines the
intact normal-opening basin. Probability outside that basin is absorbed and
is never renormalized into the surviving density. The cumulative initiated
probability is the mass explicitly removed by that opening operation. Physical
survival is one minus this absorbed mass; the direct intact-density integral is
retained independently so numerical mass residuals remain visible.

This N=1 module is a gold-standard/reference calculation, not the production
N=3 solver and not a calibrated pure-Al fatigue-life predictor.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

from .model import ModelParams, TwoRowLJ
from .configurational_plasticity import registry_rate_terms


@dataclass(frozen=True)
class Grid2DParams:
    """Cell-centred truncated domain for q=(a,s)."""

    n_a: int = 61
    n_s: int = 121
    a_lower_factor: float = 0.82
    a_upper: float = 1.80
    s_wells: int = 5


@dataclass(frozen=True)
class PDETimeParams:
    """Conservative integration controls.

    The validated production default remains explicit.  For ``integrator``
    equal to ``"explicit"``, ``max_dt`` is only an upper bound and the solver
    computes a positivity CFL bound from the Scharfetter--Gummel rates.  The
    optional backward-Euler path exists for stiff convergence diagnostics and
    uses the identical conservative spatial generator.
    """

    max_dt: float = 2.0e-3
    cfl: float = 0.45
    record_interval: float = 2.0e-2
    negative_tolerance: float = 2.0e-12
    integrator: str = "explicit"


@dataclass(frozen=True)
class CyclicLoad2D:
    force_min: float = 0.0
    force_max: float = 3.2
    period: float = 1.0
    cycles: float = 1.0
    phase_radians: float = 0.0
    value_function: Callable[[float], float] | None = None

    @property
    def duration(self) -> float:
        return float(self.period * self.cycles)

    def value(self, time: float) -> float:
        if self.value_function is not None:
            return float(self.value_function(time))
        midpoint = 0.5 * (self.force_max + self.force_min)
        amplitude = 0.5 * (self.force_max - self.force_min)
        return float(
            midpoint
            + amplitude
            * np.sin(2.0 * np.pi * time / self.period + self.phase_radians)
        )


def cyclic_load_from_sigma_over_E(
    model: TwoRowLJ,
    *,
    sigma_over_E_min: float,
    sigma_over_E_max: float,
    period: float = 1.0,
    cycles: float = 1.0,
    phase_radians: float = 0.0,
    value_function: Callable[[float], float] | None = None,
) -> CyclicLoad2D:
    r"""Build a force history from the signed reduced stress sigma/E.

    The conversion is not ``force = sigma/E``. The Bessel-LJ force coordinate
    has its own dimensionless tangent stiffness. We match the pristine relaxed
    total axial response to Young's law through

        f* = {a0 / [c^T H0^{-1} c]} (sigma/E),

    where ``c=(1,chi)^T`` and ``H0`` is the local ``(a,s)`` Hessian. The
    frozen-s value ``a0 W_aa`` remains a diagnostic rather than the canonical
    macroscopic mapping.

    This introduces no characteristic length, area, or volume. The tangent
    mapping only sets the dimensionless force coordinate; the PDE still uses
    the full nonlinear Bessel-LJ energy at every state and time.
    """

    scale = model.sigma_over_E_force_scale()
    mapped_function = None
    if value_function is not None:
        mapped_function = lambda time: float(
            model.force_from_sigma_over_E(value_function(time))
        )
    return CyclicLoad2D(
        force_min=float(scale * sigma_over_E_min),
        force_max=float(scale * sigma_over_E_max),
        period=period,
        cycles=cycles,
        phase_radians=phase_radians,
        value_function=mapped_function,
    )


@dataclass(frozen=True)
class Grid2D:
    a: np.ndarray
    s: np.ndarray
    da: float
    ds: float

    @property
    def cell_volume(self) -> float:
        return float(self.da * self.ds)


def _bernoulli(x: np.ndarray) -> np.ndarray:
    """Stable Bernoulli function B(x)=x/(exp(x)-1)."""

    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    small = np.abs(x) < 1.0e-5
    positive = x > 50.0
    negative = x < -50.0
    middle = ~(small | positive | negative)

    xs = x[small]
    out[small] = 1.0 - 0.5 * xs + xs * xs / 12.0 - xs**4 / 720.0

    xp = x[positive]
    out[positive] = xp * np.exp(-xp) / (1.0 - np.exp(-xp))

    xn = x[negative]
    out[negative] = -xn / (1.0 - np.exp(xn))

    xm = x[middle]
    out[middle] = xm / np.expm1(xm)
    return out


def build_grid(model: TwoRowLJ, params: Grid2DParams) -> Grid2D:
    if model.p.n_cells != 1:
        raise ValueError("probability_pde_2d requires ModelParams(n_cells=1)")
    if params.n_a < 5 or params.n_s < 5:
        raise ValueError("n_a and n_s must both be at least 5")
    if params.s_wells < 1 or params.s_wells % 2 == 0:
        raise ValueError("s_wells must be a positive odd integer")

    a_low = max(model.p.a_min * 1.001, params.a_lower_factor * model.a0)
    if params.a_upper <= a_low:
        raise ValueError("a_upper must exceed the lower a-domain bound")

    a_edges = np.linspace(a_low, params.a_upper, params.n_a + 1)
    s_half_span = 0.5 * params.s_wells * model.p.b
    s_edges = np.linspace(-s_half_span, s_half_span, params.n_s + 1)
    a = 0.5 * (a_edges[:-1] + a_edges[1:])
    s = 0.5 * (s_edges[:-1] + s_edges[1:])
    return Grid2D(
        a=a,
        s=s,
        da=float(a_edges[1] - a_edges[0]),
        ds=float(s_edges[1] - s_edges[0]),
    )


def _mesh(grid: Grid2D) -> tuple[np.ndarray, np.ndarray]:
    return np.meshgrid(grid.a, grid.s, indexing="ij")


def energy_grid(model: TwoRowLJ, grid: Grid2D, force: float) -> np.ndarray:
    """Evaluate the full N=1 effective energy G(a,s;force) at cell centres."""

    aa, ss = _mesh(grid)
    energy, _, _ = model.energy_gradient_batch(
        aa.reshape(-1, 1), ss.reshape(-1, 1), float(force)
    )
    return energy.reshape(aa.shape)


def configurational_interface_fluxes(
    density: np.ndarray,
    energy: np.ndarray,
    model: TwoRowLJ,
    grid: Grid2D,
) -> tuple[np.ndarray, np.ndarray]:
    """Return signed and gross probability rates at every ``s`` interface.

    Positive signed flux is toward increasing ``s``.  Gross flux is the sum of
    the two nonnegative SG one-way rates; unlike signed flux it detects
    balanced bidirectional thermal hopping.
    """

    density = np.asarray(density, dtype=float)
    energy = np.asarray(energy, dtype=float)
    if density.shape != energy.shape or density.shape != (grid.a.size, grid.s.size):
        raise ValueError("density/energy shape does not match the 2D grid")
    beta = 1.0 / model.p.kT
    diffusivity = model.p.mobility_s * model.p.kT
    psi = beta * (energy[:, 1:] - energy[:, :-1])
    forward = (diffusivity / grid.ds) * _bernoulli(psi) * density[:, :-1]
    reverse = (diffusivity / grid.ds) * _bernoulli(-psi) * density[:, 1:]
    signed = np.sum(forward - reverse, axis=0) * grid.da
    gross = np.sum(forward + reverse, axis=0) * grid.da
    return np.asarray(signed, dtype=float), np.asarray(gross, dtype=float)


def configurational_interface_one_way_fluxes(
    density: np.ndarray,
    energy: np.ndarray,
    model: TwoRowLJ,
    grid: Grid2D,
) -> tuple[np.ndarray, np.ndarray]:
    """Return nonnegative increasing-s and decreasing-s SG activities.

    These are reconstructed from the exact signed/gross SG decomposition,
    not from ``abs(net flux)``.
    """

    signed, gross = configurational_interface_fluxes(density, energy, model, grid)
    forward = 0.5 * (gross + signed)
    backward = 0.5 * (gross - signed)
    roundoff = 64.0 * np.finfo(float).eps * np.maximum(1.0, gross)
    if np.any(forward < -roundoff) or np.any(backward < -roundoff):
        raise FloatingPointError("invalid SG forward/backward flux decomposition")
    return np.maximum(forward, 0.0), np.maximum(backward, 0.0)


def configurational_well_structure(
    model: TwoRowLJ,
    grid: Grid2D,
) -> dict[str, np.ndarray]:
    """Map configurational wells to cells and their nearest FV interfaces."""

    cell_wells = np.asarray(model.well_index(grid.s), dtype=int)
    well_indices = np.arange(int(np.min(cell_wells)), int(np.max(cell_wells)) + 1)
    boundary_wells = well_indices[:-1]
    target_boundaries = (boundary_wells + 0.5) * model.p.b
    interfaces = 0.5 * (grid.s[:-1] + grid.s[1:])
    interface_indices = np.asarray(
        [int(np.argmin(np.abs(interfaces - boundary))) for boundary in target_boundaries],
        dtype=int,
    )
    represented_boundaries = interfaces[interface_indices]
    return {
        "cell_well_index": cell_wells,
        "well_indices": well_indices,
        "boundary_lower_well_index": boundary_wells,
        "boundary_target_s": np.asarray(target_boundaries, dtype=float),
        "boundary_interface_index": interface_indices,
        "boundary_interface_s": np.asarray(represented_boundaries, dtype=float),
        "boundary_alignment_error": np.asarray(
            represented_boundaries - target_boundaries, dtype=float
        ),
    }


def configurational_well_observables(
    density: np.ndarray,
    energy: np.ndarray,
    model: TwoRowLJ,
    grid: Grid2D,
) -> dict[str, np.ndarray]:
    """Return absolute well populations and nearest-boundary SG fluxes."""

    structure = configurational_well_structure(model, grid)
    populations = np.asarray(
        [
            np.sum(density[:, structure["cell_well_index"] == well])
            * grid.cell_volume
            for well in structure["well_indices"]
        ],
        dtype=float,
    )
    signed_all, gross_all = configurational_interface_fluxes(
        density, energy, model, grid
    )
    interface_indices = structure["boundary_interface_index"]
    return {
        **structure,
        "well_population": populations,
        "interwell_net_flux": signed_all[interface_indices],
        "interwell_gross_flux": gross_all[interface_indices],
        "interwell_forward_flux": 0.5
        * (gross_all[interface_indices] + signed_all[interface_indices]),
        "interwell_backward_flux": 0.5
        * (gross_all[interface_indices] - signed_all[interface_indices]),
    }


def opening_intact_mask(
    model: TwoRowLJ,
    grid: Grid2D,
    force: float,
) -> np.ndarray:
    """Return cells that remain inside the mechanically bound opening basin."""

    _, saddle, bound = model.opening_saddle_batch(grid.s[:, None], float(force))
    saddle = saddle[:, 0]
    bound = bound[:, 0]
    return bound[None, :] & (grid.a[:, None] < saddle[None, :])


def initial_gibbs_density(
    model: TwoRowLJ,
    grid: Grid2D,
    *,
    preload_force: float = 0.0,
    principal_well_only: bool = True,
) -> np.ndarray:
    """Conditional Gibbs density in the declared intact initial basin.

    The baseline uses the principal configurational well |s|<b/2. This is a
    conditional metastable initial ensemble, not an imposed Gaussian spacing
    law and not a product closure.
    """

    if model.p.kT <= 0.0:
        raise ValueError("finite-temperature Gibbs density requires kT > 0")

    g = energy_grid(model, grid, preload_force)
    mask = opening_intact_mask(model, grid, preload_force)
    if principal_well_only:
        mask &= np.abs(grid.s[None, :]) < 0.5 * model.p.b
    if not np.any(mask):
        raise ValueError("initial Gibbs basin contains no grid cells")

    g0 = float(np.min(g[mask]))
    exponent = -(g - g0) / model.p.kT
    density = np.zeros_like(g)
    density[mask] = np.exp(np.clip(exponent[mask], -745.0, 0.0))
    z = float(np.sum(density) * grid.cell_volume)
    if not np.isfinite(z) or z <= 0.0:
        raise FloatingPointError("failed to normalize the initial Gibbs density")
    return density / z


def _sg_rates_and_rhs(
    density: np.ndarray,
    energy: np.ndarray,
    model: TwoRowLJ,
    grid: Grid2D,
) -> tuple[np.ndarray, float]:
    """Return conservative SG semi-discrete RHS and maximum outgoing rate."""

    if model.p.kT <= 0.0:
        raise ValueError("Smoluchowski diffusion requires kT > 0")

    beta = 1.0 / model.p.kT
    d_a = model.p.mobility_a * model.p.kT
    d_s = model.p.mobility_s * model.p.kT
    rhs = np.zeros_like(density)
    outgoing = np.zeros_like(density)

    psi_a = beta * (energy[1:, :] - energy[:-1, :])
    bp_a = _bernoulli(psi_a)
    bm_a = _bernoulli(-psi_a)
    flux_a = (d_a / grid.da) * (
        bp_a * density[:-1, :] - bm_a * density[1:, :]
    )
    rhs[:-1, :] -= flux_a / grid.da
    rhs[1:, :] += flux_a / grid.da
    outgoing[:-1, :] += (d_a / grid.da**2) * bp_a
    outgoing[1:, :] += (d_a / grid.da**2) * bm_a

    psi_s = beta * (energy[:, 1:] - energy[:, :-1])
    bp_s = _bernoulli(psi_s)
    bm_s = _bernoulli(-psi_s)
    flux_s = (d_s / grid.ds) * (
        bp_s * density[:, :-1] - bm_s * density[:, 1:]
    )
    rhs[:, :-1] -= flux_s / grid.ds
    rhs[:, 1:] += flux_s / grid.ds
    outgoing[:, :-1] += (d_s / grid.ds**2) * bp_s
    outgoing[:, 1:] += (d_s / grid.ds**2) * bm_s

    return rhs, float(np.max(outgoing))


def _sg_generator_2d(
    energy: np.ndarray,
    model: TwoRowLJ,
    grid: Grid2D,
) -> sparse.csr_matrix:
    """Assemble the same SG semi-discretization as ``dp/dt=A p``."""

    if model.p.kT <= 0.0:
        raise ValueError("Smoluchowski diffusion requires kT > 0")
    shape = energy.shape
    index = np.arange(np.prod(shape), dtype=np.int64).reshape(shape)
    beta = 1.0 / model.p.kT
    rows: list[np.ndarray] = []
    columns: list[np.ndarray] = []
    values: list[np.ndarray] = []
    for axis, diffusivity, spacing in (
        (0, model.p.mobility_a * model.p.kT, grid.da),
        (1, model.p.mobility_s * model.p.kT, grid.ds),
    ):
        left = [slice(None)] * 2
        right = [slice(None)] * 2
        left[axis] = slice(0, -1)
        right[axis] = slice(1, None)
        left_t = tuple(left)
        right_t = tuple(right)
        left_index = index[left_t].ravel()
        right_index = index[right_t].ravel()
        psi = beta * (energy[right_t] - energy[left_t])
        rate_left_right = (diffusivity / spacing**2) * _bernoulli(psi)
        rate_right_left = (diffusivity / spacing**2) * _bernoulli(-psi)
        rate_left_right = rate_left_right.ravel()
        rate_right_left = rate_right_left.ravel()
        rows.extend((right_index, left_index, left_index, right_index))
        columns.extend((left_index, left_index, right_index, right_index))
        values.extend(
            (rate_left_right, -rate_left_right, rate_right_left, -rate_right_left)
        )
    return sparse.coo_matrix(
        (np.concatenate(values), (np.concatenate(rows), np.concatenate(columns))),
        shape=(index.size, index.size),
    ).tocsr()


def _implicit_step_2d(
    density: np.ndarray,
    energy: np.ndarray,
    model: TwoRowLJ,
    grid: Grid2D,
    dt: float,
    negative_tolerance: float,
) -> tuple[np.ndarray, float, float]:
    """Take one conservative backward-Euler step with the SG generator."""

    generator = _sg_generator_2d(energy, model, grid)
    system = sparse.eye(generator.shape[0], format="csr") - dt * generator
    updated = np.asarray(spsolve(system, density.ravel()), dtype=float).reshape(
        density.shape
    )
    minimum = float(np.min(updated))
    if minimum < -negative_tolerance:
        raise FloatingPointError(
            f"implicit 2D probability density lost positivity ({minimum:.3e})"
        )
    negative_mass_correction = float(
        np.sum(np.maximum(-updated, 0.0)) * grid.cell_volume
    )
    return np.maximum(updated, 0.0), negative_mass_correction, minimum


def _absorb_outside_opening_basin(
    density: np.ndarray,
    model: TwoRowLJ,
    grid: Grid2D,
    force: float,
) -> tuple[np.ndarray, float]:
    mask = opening_intact_mask(model, grid, force)
    removed = float(np.sum(density[~mask]) * grid.cell_volume)
    if removed == 0.0:
        return density, 0.0
    out = density.copy()
    out[~mask] = 0.0
    return out, removed


def observables(
    density: np.ndarray,
    model: TwoRowLJ,
    grid: Grid2D,
    force: float,
    *,
    first_passage_flux: float = 0.0,
) -> dict[str, float]:
    """Compute survivor-conditioned mechanical/configurational observables.

    The unwrapped configurational coordinate is split exactly as

        s = b*n + xi,  xi in [-b/2,b/2),

    so the axial strain bridge is reported as three additive parts:

        epsilon = epsilon_a + epsilon_xi + epsilon_p.

    ``plastic_strain`` is the signed mean well-index contribution. It is not
    the nonnegative activity diagnostic ``plastic_well_activity``.
    """

    volume = grid.cell_volume
    survival = float(np.sum(density) * volume)
    aa, ss = _mesh(grid)
    if survival > 0.0:
        conditional = density / survival
        nwell = model.well_index(ss)
        xi = ss - model.p.b * nwell

        normal_field = (aa - model.a0) / model.a0
        intrawell_field = model.p.chi_axial_projection * xi / model.a0
        plastic_field = (
            model.p.chi_axial_projection * model.p.b * nwell / model.a0
        )

        normal_strain = float(np.sum(conditional * normal_field) * volume)
        intrawell_strain = float(np.sum(conditional * intrawell_field) * volume)
        plastic_strain = float(np.sum(conditional * plastic_field) * volume)
        strain = normal_strain + intrawell_strain + plastic_strain
        strain_decomposition_residual = (
            strain - normal_strain - intrawell_strain - plastic_strain
        )
        mean_well_index = float(np.sum(conditional * nwell) * volume)
        well_activity = float(np.sum(conditional * np.abs(nwell)) * volume)
    else:
        strain = np.nan
        normal_strain = np.nan
        intrawell_strain = np.nan
        plastic_strain = np.nan
        strain_decomposition_residual = np.nan
        mean_well_index = np.nan
        well_activity = np.nan

    edge = max(1, min(2, grid.s.size // 4))
    s_boundary_mass = float(
        (np.sum(density[:, :edge]) + np.sum(density[:, -edge:])) * volume
    )
    a_lower_mass = float(np.sum(density[:edge, :]) * volume)
    a_upper_mass = float(np.sum(density[-edge:, :]) * volume)

    return {
        "force": float(force),
        "survival": survival,
        "initiation_probability": 1.0 - survival,
        "first_passage_flux": float(first_passage_flux),
        "strain": strain,
        "normal_strain": normal_strain,
        "intrawell_strain": intrawell_strain,
        "plastic_strain": plastic_strain,
        "strain_decomposition_residual": strain_decomposition_residual,
        "mean_well_index": mean_well_index,
        "plastic_well_activity": well_activity,
        "s_truncation_boundary_mass": s_boundary_mass,
        "a_lower_boundary_mass": a_lower_mass,
        "a_upper_boundary_mass": a_upper_mass,
    }


def run_probability_pde_2d(
    *,
    model_params: ModelParams | None = None,
    prepared_model: TwoRowLJ | None = None,
    grid_params: Grid2DParams = Grid2DParams(),
    time_params: PDETimeParams = PDETimeParams(),
    load: CyclicLoad2D = CyclicLoad2D(),
    preload_force: float = 0.0,
    record_callback: Callable[[dict[str, float]], None] | None = None,
    density_record_callback: (
        Callable[[float, float, np.ndarray, TwoRowLJ, Grid2D], None] | None
    ) = None,
    stop_requested: Callable[[], bool] | None = None,
) -> dict[str, np.ndarray | Grid2D | TwoRowLJ | str]:
    """Solve the deterministic N=1 probability PDE under a cyclic load."""

    if prepared_model is not None and model_params is not None:
        raise ValueError("pass either model_params or prepared_model, not both")
    p = prepared_model.p if prepared_model is not None else (
        model_params or ModelParams(n_cells=1)
    )
    if p.n_cells != 1:
        raise ValueError("N=1 reference solver requires n_cells=1")
    if time_params.max_dt <= 0.0 or not 0.0 < time_params.cfl < 1.0:
        raise ValueError("invalid PDE time controls")
    if time_params.record_interval <= 0.0:
        raise ValueError("record_interval must be positive")
    if time_params.integrator not in {"explicit", "implicit"}:
        raise ValueError("2D integrator must be 'explicit' or 'implicit'")
    if load.period <= 0.0 or load.cycles < 0.0:
        raise ValueError("invalid cyclic load duration")

    model = prepared_model if prepared_model is not None else TwoRowLJ(p)
    if not model._opening_table_ready:
        model._build_opening_table()
    grid = build_grid(model, grid_params)
    density = initial_gibbs_density(
        model, grid, preload_force=preload_force, principal_well_only=True
    )
    well_structure = configurational_well_structure(model, grid)
    well_indices = well_structure["well_indices"]
    boundary_interfaces = well_structure["boundary_interface_index"]
    cumulative_interwell_net = np.zeros(boundary_interfaces.size, dtype=float)
    cumulative_interwell_gross = np.zeros(boundary_interfaces.size, dtype=float)
    cumulative_opening_by_well = np.zeros(well_indices.size, dtype=float)
    interval_opening_by_well = np.zeros(well_indices.size, dtype=float)
    initial_well_population: np.ndarray | None = None
    initial_registry_moment: float | None = None
    well_population_records: list[np.ndarray] = []
    interwell_net_flux_records: list[np.ndarray] = []
    interwell_gross_flux_records: list[np.ndarray] = []
    cumulative_interwell_net_records: list[np.ndarray] = []
    cumulative_interwell_gross_records: list[np.ndarray] = []
    well_balance_residual_records: list[np.ndarray] = []
    cumulative_opening_by_well_records: list[np.ndarray] = []

    records: dict[str, list[float]] = {
        key: []
        for key in (
            "time",
            "force",
            "survival",
            "initiation_probability",
            "first_passage_flux",
            "strain",
            "normal_strain",
            "intrawell_strain",
            "plastic_strain",
            "strain_decomposition_residual",
            "mean_well_index",
            "plastic_well_activity",
            "s_truncation_boundary_mass",
            "a_lower_boundary_mass",
            "a_upper_boundary_mass",
            "intact_probability_mass",
            "cumulative_absorbed_mass",
            "absorbed_mass_increment",
            "initial_absorbed_mass",
            "integrated_first_passage_flux",
            "flux_consistency_residual",
            "mass_balance_residual",
            "negative_mass_correction",
            "cumulative_negative_mass_correction",
            "minimum_density",
            "unnormalized_registry_moment",
            "accumulated_net_registry_transfer",
            "absorbed_registry_moment",
            "registry_moment_balance_residual",
            "net_interwell_registry_rate",
            "gross_interwell_activity_rate",
            "net_plastic_flow_rate",
            "gross_configurational_slip_activity",
            "selective_opening_plastic_rate",
            "cumulative_forward_registry_activity",
            "cumulative_backward_registry_activity",
            "cumulative_gross_registry_activity",
        )
    }

    t = 0.0
    next_record = 0.0
    duration = load.duration
    cumulative_absorbed_mass = 0.0
    initial_absorbed_mass = 0.0
    interval_absorbed_mass = 0.0
    integrated_first_passage_flux = 0.0
    last_record_time = 0.0
    maximum_negative_mass_correction = 0.0
    cumulative_negative_mass_correction = 0.0
    minimum_density = float(np.min(density))

    def append_record(now: float, force: float) -> None:
        nonlocal interval_absorbed_mass, last_record_time
        nonlocal initial_well_population, initial_registry_moment
        elapsed = float(now - last_record_time)
        interval_flux = (
            interval_absorbed_mass / elapsed
            if elapsed > 10.0 * np.finfo(float).eps
            else 0.0
        )
        obs = observables(
            density, model, grid, force, first_passage_flux=interval_flux
        )
        well_obs = configurational_well_observables(
            density, energy_grid(model, grid, force), model, grid
        )
        population = well_obs["well_population"]
        if initial_well_population is None:
            initial_well_population = population.copy()
            initial_registry_moment = float(well_indices @ population)
        expected_population = initial_well_population - cumulative_opening_by_well
        if cumulative_interwell_net.size:
            expected_population[1:] += cumulative_interwell_net
            expected_population[:-1] -= cumulative_interwell_net
        well_population_records.append(population.copy())
        interwell_net_flux_records.append(well_obs["interwell_net_flux"].copy())
        interwell_gross_flux_records.append(well_obs["interwell_gross_flux"].copy())
        cumulative_interwell_net_records.append(cumulative_interwell_net.copy())
        cumulative_interwell_gross_records.append(cumulative_interwell_gross.copy())
        well_balance_residual_records.append(
            (population - expected_population).copy()
        )
        cumulative_opening_by_well_records.append(
            cumulative_opening_by_well.copy()
        )
        opening_rate_by_well = (
            interval_opening_by_well / elapsed
            if elapsed > 10.0 * np.finfo(float).eps
            else np.zeros_like(interval_opening_by_well)
        )
        registry_terms = registry_rate_terms(
            well_indices=well_indices,
            well_populations=population,
            interwell_net_flux=well_obs["interwell_net_flux"],
            interwell_gross_flux=well_obs["interwell_gross_flux"],
            opening_absorption_rate_by_well=opening_rate_by_well,
        )
        accumulated_registry_transfer = float(np.sum(cumulative_interwell_net))
        accumulated_gross_activity = float(np.sum(cumulative_interwell_gross))
        absorbed_registry_moment = float(well_indices @ cumulative_opening_by_well)
        registry_balance_residual = float(
            registry_terms.unnormalized_registry_moment
            - float(initial_registry_moment)
            - accumulated_registry_transfer
            + absorbed_registry_moment
        )
        plastic_factor = (
            model.p.chi_axial_projection * model.p.b / model.a0
        )
        snapshot = {
            "time": float(now),
            **{name: float(value) for name, value in obs.items()},
            # The physical first-passage fields are defined by probability
            # actually removed at the opening boundary.  The independently
            # integrated survivor mass is retained below as a numerical
            # diagnostic and may differ at roundoff/discretization level.
            "survival": float(1.0 - cumulative_absorbed_mass),
            "initiation_probability": float(cumulative_absorbed_mass),
            "intact_probability_mass": float(obs["survival"]),
            "cumulative_absorbed_mass": float(cumulative_absorbed_mass),
            "absorbed_mass_increment": float(interval_absorbed_mass),
            "initial_absorbed_mass": float(initial_absorbed_mass),
            "integrated_first_passage_flux": float(
                integrated_first_passage_flux
            ),
            "flux_consistency_residual": float(
                cumulative_absorbed_mass
                - initial_absorbed_mass
                - integrated_first_passage_flux
            ),
            "mass_balance_residual": float(
                obs["survival"] + cumulative_absorbed_mass - 1.0
            ),
            "negative_mass_correction": float(maximum_negative_mass_correction),
            "cumulative_negative_mass_correction": float(
                cumulative_negative_mass_correction
            ),
            "minimum_density": float(minimum_density),
            "unnormalized_registry_moment": float(
                registry_terms.unnormalized_registry_moment
            ),
            "accumulated_net_registry_transfer": accumulated_registry_transfer,
            "absorbed_registry_moment": absorbed_registry_moment,
            "registry_moment_balance_residual": registry_balance_residual,
            "net_interwell_registry_rate": float(
                registry_terms.net_registry_flow_rate
            ),
            "gross_interwell_activity_rate": float(
                registry_terms.gross_configurational_activity_rate
            ),
            "net_plastic_flow_rate": float(
                plastic_factor * registry_terms.conditional_net_registry_flow_rate
            ),
            "gross_configurational_slip_activity": float(
                abs(plastic_factor)
                * registry_terms.gross_configurational_activity_rate
                / max(float(np.sum(population)), np.finfo(float).tiny)
            ),
            "selective_opening_plastic_rate": float(
                plastic_factor * registry_terms.conditional_selective_opening_rate
            ),
            "cumulative_forward_registry_activity": 0.5
            * (accumulated_gross_activity + accumulated_registry_transfer),
            "cumulative_backward_registry_activity": 0.5
            * (accumulated_gross_activity - accumulated_registry_transfer),
            "cumulative_gross_registry_activity": accumulated_gross_activity,
        }
        for name in records:
            records[name].append(snapshot[name])
        if record_callback is not None:
            record_callback(dict(snapshot))
        if density_record_callback is not None:
            density_record_callback(
                float(now), float(force), density.copy(), model, grid
            )
        interval_absorbed_mass = 0.0
        interval_opening_by_well.fill(0.0)
        last_record_time = float(now)

    while True:
        force = load.value(t)
        density_before_instant_absorption = density
        density, instant_loss = _absorb_outside_opening_basin(
            density, model, grid, force
        )
        if instant_loss > 0.0:
            cumulative_absorbed_mass += instant_loss
            if t <= 10.0 * np.finfo(float).eps and not records["time"]:
                initial_absorbed_mass += instant_loss
            else:
                interval_absorbed_mass += instant_loss
                integrated_first_passage_flux += instant_loss
                removed_density = density_before_instant_absorption - density
                for index, well in enumerate(well_indices):
                    columns = well_structure["cell_well_index"] == well
                    cumulative_opening_by_well[index] += float(
                        np.sum(removed_density[:, columns]) * grid.cell_volume
                    )
                    interval_opening_by_well[index] += float(
                        np.sum(removed_density[:, columns]) * grid.cell_volume
                    )

        if t + 1.0e-14 >= next_record:
            append_record(t, force)
            next_record += time_params.record_interval

        if t >= duration - 1.0e-14 or np.sum(density) == 0.0:
            break
        if stop_requested is not None and stop_requested():
            break

        if time_params.integrator == "explicit":
            energy = energy_grid(model, grid, force)
            step_net_all, step_gross_all = configurational_interface_fluxes(
                density, energy, model, grid
            )
            rhs, max_outgoing_rate = _sg_rates_and_rhs(density, energy, model, grid)
            if max_outgoing_rate > 0.0:
                stable_dt = time_params.cfl / max_outgoing_rate
            else:
                stable_dt = time_params.max_dt
            dt = min(time_params.max_dt, stable_dt, duration - t)
            if not np.isfinite(dt) or dt <= 0.0:
                raise FloatingPointError("failed to obtain a positive stable PDE step")
            trial = density + dt * rhs
            min_value = float(np.min(trial))
            negative_mass_correction = float(
                np.sum(np.maximum(-trial, 0.0)) * grid.cell_volume
            )
            trial = np.maximum(trial, 0.0)
        else:
            dt = min(time_params.max_dt, duration - t)
            if not np.isfinite(dt) or dt <= 0.0:
                raise FloatingPointError("failed to obtain a positive implicit PDE step")
            implicit_force = load.value(t + dt)
            implicit_energy = energy_grid(model, grid, implicit_force)
            trial, negative_mass_correction, min_value = _implicit_step_2d(
                density,
                implicit_energy,
                model,
                grid,
                dt,
                time_params.negative_tolerance,
            )
            step_net_all, step_gross_all = configurational_interface_fluxes(
                trial, implicit_energy, model, grid
            )
        minimum_density = min(minimum_density, min_value)
        if min_value < -time_params.negative_tolerance:
            raise FloatingPointError(
                f"probability density became negative ({min_value:.3e}); refine dt/grid"
            )
        maximum_negative_mass_correction = max(
            maximum_negative_mass_correction, negative_mass_correction
        )
        cumulative_negative_mass_correction += negative_mass_correction
        cumulative_interwell_net += dt * step_net_all[boundary_interfaces]
        cumulative_interwell_gross += dt * step_gross_all[boundary_interfaces]

        next_t = t + dt
        next_force = load.value(next_t)
        before_absorb = float(np.sum(trial) * grid.cell_volume)
        trial_before_absorption = trial
        trial, removed = _absorb_outside_opening_basin(
            trial_before_absorption, model, grid, next_force
        )
        cumulative_absorbed_mass += removed
        interval_absorbed_mass += removed
        integrated_first_passage_flux += removed
        if removed > 0.0:
            removed_density = trial_before_absorption - trial
            for index, well in enumerate(well_indices):
                columns = well_structure["cell_well_index"] == well
                cumulative_opening_by_well[index] += float(
                    np.sum(removed_density[:, columns]) * grid.cell_volume
                )
                interval_opening_by_well[index] += float(
                    np.sum(removed_density[:, columns]) * grid.cell_volume
                )
        numerical_mass_error = abs(
            before_absorb - float(np.sum(density) * grid.cell_volume)
        )
        if numerical_mass_error > 5.0e-10:
            raise FloatingPointError(
                f"finite-volume mass conservation error {numerical_mass_error:.3e}"
            )

        density = trial
        t = next_t

    if not records["time"] or abs(records["time"][-1] - t) > 1.0e-12:
        append_record(t, load.value(t))

    return {
        **{name: np.asarray(values, dtype=float) for name, values in records.items()},
        "well_indices": well_indices.copy(),
        "well_populations": np.asarray(well_population_records, dtype=float),
        "interwell_boundary_lower_index": well_structure[
            "boundary_lower_well_index"
        ].copy(),
        "interwell_boundary_s": well_structure["boundary_target_s"].copy(),
        "interwell_boundary_grid_s": well_structure[
            "boundary_interface_s"
        ].copy(),
        "interwell_boundary_alignment_error": well_structure[
            "boundary_alignment_error"
        ].copy(),
        "interwell_net_flux": np.asarray(
            interwell_net_flux_records, dtype=float
        ),
        "interwell_gross_flux": np.asarray(
            interwell_gross_flux_records, dtype=float
        ),
        "interwell_forward_flux": 0.5
        * (
            np.asarray(interwell_gross_flux_records, dtype=float)
            + np.asarray(interwell_net_flux_records, dtype=float)
        ),
        "interwell_backward_flux": 0.5
        * (
            np.asarray(interwell_gross_flux_records, dtype=float)
            - np.asarray(interwell_net_flux_records, dtype=float)
        ),
        "cumulative_interwell_net_transfer": np.asarray(
            cumulative_interwell_net_records, dtype=float
        ),
        "cumulative_interwell_gross_transfer": np.asarray(
            cumulative_interwell_gross_records, dtype=float
        ),
        "cumulative_interwell_forward_transfer": 0.5
        * (
            np.asarray(cumulative_interwell_gross_records, dtype=float)
            + np.asarray(cumulative_interwell_net_records, dtype=float)
        ),
        "cumulative_interwell_backward_transfer": 0.5
        * (
            np.asarray(cumulative_interwell_gross_records, dtype=float)
            - np.asarray(cumulative_interwell_net_records, dtype=float)
        ),
        "cumulative_opening_absorption_by_well": cumulative_opening_by_well.copy(),
        "cumulative_opening_absorption_by_well_history": np.asarray(
            cumulative_opening_by_well_records, dtype=float
        ),
        "well_population_balance_residual": np.asarray(
            well_balance_residual_records, dtype=float
        ),
        "density": density,
        "grid": grid,
        "model": model,
        "first_passage_flux_definition": (
            "record-interval absorbed mass divided by record elapsed time"
        ),
        "absorbed_mass_definition": (
            "mass removed only by the opening absorbing-boundary mask"
        ),
    }
