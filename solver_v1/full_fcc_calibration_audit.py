"""Independent elastic information and deterministic full-FCC calibration.

Research only.  The production energy registry is not changed.  At fixed
density decay the atomistic observables are linear in the two LJ coefficients,
sqrt embedding amplitude, and linear embedding amplitude.  A deterministic
one-dimensional profile therefore exposes tradeoffs without random starts.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import math

import numpy as np
from scipy.optimize import least_squares, minimize_scalar, nnls, root_scalar

from .aluminum_calibration import EV_J, GPA_PA, mishin_lu_wei_targets
from .aluminum_full_fcc_calibration import (
    FullFCCParameters, build_full_fcc_calibration_model,
)


NAMES = ("normal_strain_force_ev", "cohesive_ev_atom",
         "hydrostatic_strain_curvature_ev", "normal_strain_curvature_ev",
         "simple_shear_strain_curvature_ev")
IDEAL_H = math.sqrt(2.0 / 3.0)


def atomic_volume_m3() -> float:
    geometry, _ = mishin_lu_wei_targets()
    return (geometry.atomic_cell_area_angstrom2
            * geometry.plane_spacing_angstrom * 1e-30)


def independent_bulk_targets() -> tuple[np.ndarray, np.ndarray]:
    """Zero-K Mishin EAM targets, not a claim of experimental exactness."""
    c11, c12, c44 = 114.0, 62.0, 32.0
    convert = atomic_volume_m3() * GPA_PA / EV_J
    values = np.array([
        0.0, 3.36,
        3.0 * (c11 + 2*c12) * convert,
        (c11 + 2*c12 + 4*c44) / 3.0 * convert,
        (c11 - c12 + c44) / 3.0 * convert,
    ])
    # eV/strain force, 2% cohesion, 5% each independent elastic measurement.
    scales = np.array([0.10, 0.0672, *(.05 * values[2:])])
    return values, scales


def cubic_constants_gpa(observables: np.ndarray, *, volume_scale: float = 1.) -> dict[str, float]:
    """Recover cubic constants at zero pressure from three independent modes.

    At nonzero pressure these are affine curvature diagnostics, not relaxed
    zero-pressure elastic constants; the caller must check the force residual.
    """
    y = np.asarray(observables, dtype=float)
    if not np.isfinite(volume_scale) or volume_scale <= 0:
        raise ValueError("volume_scale must be positive and finite")
    convert = EV_J / (atomic_volume_m3()*volume_scale) / GPA_PA
    x = y[2] / 3.0 * convert
    c44 = (3*y[3]*convert - x) / 4.0
    difference = 3*y[4]*convert - c44
    return {"C11_GPa": (x + 2*difference)/3,
            "C12_GPa": (x - difference)/3, "C44_GPa": c44,
            "bulk_modulus_GPa": x/3,
            "normal_stress_GPa": y[0]*convert}


@lru_cache(maxsize=1)
def _pair_columns() -> np.ndarray:
    p = FullFCCParameters(.1, 1., 3., 1.)
    model = build_full_fcc_calibration_model(p)
    columns = []
    for power, sign in ((6., 1.), (3., -1.)):
        result = model.full_power_sum(IDEAL_H, 0., power)
        # q=1+eta scales ALL distances.  r^(-2p) gives exact q derivatives.
        degree = 2*power
        columns.append(sign*np.array([
            IDEAL_H*result.d_da, -result.value,
            degree*(degree+1)*result.value,
            IDEAL_H**2*result.d2_daa, IDEAL_H**2*result.d2_dss,
        ]))
    return np.column_stack(columns)


@lru_cache(maxsize=512)
def bulk_basis(decay: float, strain_step: float = 8e-4) -> np.ndarray:
    """Return y = basis @ [4 eps sigma^12, 4 eps sigma^6, A, B].

    Normal/shear derivatives use the verified analytic reciprocal derivatives.
    Only isotropic environmental-density derivatives use a five-point stencil
    of the SAME infinite-lattice density.  These are checked by refinement.
    Density reference is fixed during strain, not renormalized at each point.
    """
    if not np.isfinite(decay) or decay <= 0:
        raise ValueError("density decay must be positive and finite")
    if not np.isfinite(strain_step) or strain_step <= 0:
        raise ValueError("strain step must be positive and finite")
    p = FullFCCParameters(.1, 1., float(decay), 1.)
    model = build_full_fcc_calibration_model(p)
    rho = model.full_environment_density(IDEAL_H, 0.)
    rho_ref = rho.value
    x_a = rho.d_da/rho_ref
    x_s = rho.d_ds/rho_ref
    x_aa = rho.d2_daa/rho_ref
    x_ss = rho.d2_dss/rho_ref

    xs = {0: 1.}
    for i in (-2, -1, 1, 2):
        stretch = 1 + i*strain_step
        shifted = build_full_fcc_calibration_model(
            p, b=stretch, rho_ref=rho_ref)
        xs[i] = shifted.full_environment_density(
            IDEAL_H*stretch, 0.).value/rho_ref

    def d2(values: dict[int, float]) -> float:
        return (-values[2]+16*values[1]-30*values[0]
                +16*values[-1]-values[-2])/(12*strain_step**2)

    sqrt_column = np.array([
        -.5*IDEAL_H*x_a, 1.,
        d2({i: -math.sqrt(x) for i, x in xs.items()}),
        IDEAL_H**2*(.25*x_a*x_a-.5*x_aa),
        IDEAL_H**2*(.25*x_s*x_s-.5*x_ss),
    ])
    linear_column = np.array([
        IDEAL_H*x_a, -1., d2(xs),
        IDEAL_H**2*x_aa, IDEAL_H**2*x_ss,
    ])
    return np.column_stack((_pair_columns(), sqrt_column, linear_column))


def parameter_coefficients(p: FullFCCParameters) -> np.ndarray:
    return np.array([4*p.epsilon_lj_ev*p.sigma_lj_over_b**12,
                     4*p.epsilon_lj_ev*p.sigma_lj_over_b**6,
                     p.embedding_sqrt_ev, p.embedding_linear_ev])


def independent_observables(
    p: FullFCCParameters, *, strain_step: float = 8e-4,
) -> np.ndarray:
    return bulk_basis(p.density_decay_b, strain_step) @ parameter_coefficients(p)


def independent_sensitivity(
    p: FullFCCParameters, *, extended: bool = True,
    log_step: float = 2e-4, strain_step: float = 8e-4,
) -> np.ndarray:
    """Central log sensitivities; scales are fixed physical target scales."""
    if log_step <= 0 or not np.isfinite(log_step):
        raise ValueError("log_step must be finite and positive")
    values = p.as_array(extended=extended)
    if np.any(values <= 0):
        raise ValueError("log sensitivities require positive parameters")
    _, scales = independent_bulk_targets()
    columns = []
    for i in range(len(values)):
        positive, negative = values.copy(), values.copy()
        positive[i] *= math.exp(log_step)
        negative[i] *= math.exp(-log_step)
        yp = independent_observables(FullFCCParameters(*positive),
                                    strain_step=strain_step)
        ym = independent_observables(FullFCCParameters(*negative),
                                    strain_step=strain_step)
        columns.append((yp-ym)/(2*log_step*scales))
    return np.column_stack(columns)


@dataclass(frozen=True)
class ProfilePoint:
    decay: float
    coefficients: np.ndarray
    prediction: np.ndarray
    residual: np.ndarray
    squared_loss: float
    parameters: FullFCCParameters | None
    admissible: bool
    reason: str


def profile_point(
    decay: float, *, extended: bool, observation_indices=(0, 1, 2, 3, 4),
    strain_step: float = 8e-4,
) -> ProfilePoint:
    """Constrained *linear* least squares at fixed decay; no penalty clipping."""
    indices = np.asarray(observation_indices, dtype=int)
    if indices.ndim != 1 or len(indices) == 0 or len(set(indices)) != len(indices):
        raise ValueError("observation indices must be a nonempty unique sequence")
    if np.any((indices < 0) | (indices >= len(NAMES))):
        raise ValueError("invalid observation index")
    target, scale = independent_bulk_targets()
    basis = bulk_basis(float(decay), strain_step)[:, :4 if extended else 3]
    coefficients, _ = nnls(basis[indices]/scale[indices,None],
                            target[indices]/scale[indices], maxiter=300)
    prediction = basis@coefficients
    residual = (prediction-target)/scale
    if not extended:
        coefficients = np.r_[coefficients, 0.]
    u, v, amplitude, linear = coefficients
    # Zero NNLS coefficients mean a boundary/infeasible LJ or embedding family,
    # not a tiny positive parameter manufactured by clipping.
    admissible = bool(u > 0 and v > 0 and amplitude > 0)
    parameters = None
    reason = "positive pair and embedding coefficients"
    if admissible:
        sigma = (u/v)**(1/6)
        epsilon = v*v/(4*u)
        parameters = FullFCCParameters(
            float(epsilon), float(sigma), float(decay),
            float(amplitude), float(linear))
    else:
        reason = "NNLS boundary: missing positive LJ/square-root coefficient"
    return ProfilePoint(float(decay), coefficients, prediction, residual,
                        float(residual[indices]@residual[indices]),
                        parameters, admissible, reason)


def deterministic_profile(
    *, extended: bool, decay_grid: np.ndarray,
    observation_indices=(0, 1, 2, 3, 4),
) -> tuple[ProfilePoint, list[ProfilePoint]]:
    """Evaluate every declared decay then refine each bracketed local minimum.

    The range is a declared research search range, not a physical prior.  End
    points and invalid coefficient boundaries remain visible in the record.
    """
    grid = np.asarray(decay_grid, dtype=float)
    if len(grid) < 3 or np.any(~np.isfinite(grid)) or np.any(grid <= 0):
        raise ValueError("need at least three finite positive decay points")
    if np.any(np.diff(grid) <= 0):
        raise ValueError("decay grid must be strictly increasing")
    rows = [profile_point(float(x), extended=extended,
                          observation_indices=observation_indices) for x in grid]
    candidates = rows.copy()
    for i in range(1, len(rows)-1):
        if (rows[i].squared_loss <= rows[i-1].squared_loss
                and rows[i].squared_loss <= rows[i+1].squared_loss):
            solved = minimize_scalar(
                lambda x: profile_point(float(x), extended=extended,
                    observation_indices=observation_indices).squared_loss,
                bounds=(grid[i-1], grid[i+1]), method="bounded",
                options={"xatol": 3e-7, "maxiter": 80})
            candidates.append(profile_point(
                float(solved.x), extended=extended,
                observation_indices=observation_indices))
    valid = [row for row in candidates if row.admissible]
    if not valid:
        raise ValueError("no admissible positive LJ/embedding profile point")
    return min(valid, key=lambda row: row.squared_loss), candidates


def bounded_independent_fit(
    start: FullFCCParameters, *, extended: bool = True, max_nfev: int = 180,
) -> tuple[FullFCCParameters, dict[str, object]]:
    """Fit independent modes within the PREVIOUS declared search box.

    These bounds are search restrictions, not independently established Al
    parameter uncertainties.  Compare to the unrestricted coefficient profile.
    """
    count = 5 if extended else 4
    lower = np.array([1e-4,.45,.35,1e-4,1e-4])[:count]
    upper = np.array([2.,1.35,8.,12.,12.])[:count]
    initial = start.as_array(extended=extended)
    if np.any(initial < lower) or np.any(initial > upper):
        raise ValueError("start must be inside the declared search box")
    target, scales = independent_bulk_targets()

    def residual(x):
        return (independent_observables(FullFCCParameters(*np.exp(x)))
                -target)/scales

    def jacobian(x):
        return independent_sensitivity(FullFCCParameters(*np.exp(x)),
                                       extended=extended)

    solved = least_squares(
        residual, np.log(initial), jac=jacobian,
        bounds=(np.log(lower),np.log(upper)), method="trf", x_scale="jac",
        max_nfev=max_nfev, ftol=1e-9, xtol=1e-9, gtol=1e-8)
    final = FullFCCParameters(*np.exp(solved.x))
    return final, {
        "start": initial.tolist(), "lower": lower.tolist(), "upper": upper.tolist(),
        "nfev": int(solved.nfev), "success": bool(solved.success),
        "message": solved.message, "active_mask": solved.active_mask.tolist(),
        "squared_loss": float(2*solved.cost),
        "normalized_residuals": solved.fun.tolist(),
        "parameters": final.as_array(extended=extended).tolist(),
    }


def isotropic_bulk_equilibrium(parameters: FullFCCParameters):
    """Relax FCC lattice constant, NOT just h at artificially fixed b.

    L0 remains the target b in metres.  b_model and h_model both change by
    stretch.  Density-reference gauge and dimensional potential are held fixed.
    """
    initial=build_full_fcc_calibration_model(parameters)
    rho_ref=initial.embedding.rho_ref

    def pressure(stretch):
        model=build_full_fcc_calibration_model(
            parameters,b=float(stretch),rho_ref=rho_ref)
        # Cubic symmetry: hydrostatic derivative = 3 h_target W_a.
        return 3*IDEAL_H*model.local_deda(IDEAL_H*stretch,0)

    scan=np.linspace(.90,1.10,21)
    forces=[pressure(x) for x in scan]
    roots=[]
    for l,r,fl,fr in zip(scan[:-1],scan[1:],forces[:-1],forces[1:]):
        if fl<=0<=fr:
            solved=root_scalar(pressure,bracket=(l,r),xtol=2e-12,method="brentq")
            if solved.converged:
                roots.append(float(solved.root))
    if not roots:
        raise ValueError("no stable isotropic FCC equilibrium near target")
    stretch=min(roots,key=lambda x:abs(x-1))
    model=build_full_fcc_calibration_model(parameters,b=stretch,rho_ref=rho_ref)
    model.a0=IDEAL_H*stretch
    model.equilibrium_error=None
    return model,stretch


def relaxed_bulk_observables(model, *, strain_step=8e-4):
    """Independent affine curvatures at the actual pressure-free FCC lattice."""
    from dataclasses import replace
    from .fcc111_full_energy import FullFCC111StackEnergy
    h=model.a0
    values={}
    for i in (-2,-1,0,1,2):
        stretch=1+i*strain_step
        shifted=FullFCC111StackEnergy(
            replace(model.p,b=model.p.b*stretch),
            density_params=model.density_params,embedding=model.embedding,
            stack_config=model.stack_config,find_equilibrium=False)
        values[i]=shifted.energy(h*stretch,0)
    hydro=(-values[2]+16*values[1]-30*values[0]+16*values[-1]-values[-2])/(12*strain_step**2)
    value=model.evaluate(h,0)
    atomized=float(model.embedding.value(0)) if model.embedding else 0.
    return np.array([h*value.d_da,atomized-value.energy,hydro,
                     h*h*value.d2_daa,h*h*value.d2_dss])
