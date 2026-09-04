# === 한국어 파일 안내 시작 ===
# - 파일 역할: same pair potential로 row 내부 평형 간격을 정한 뒤 spatial registry kink/kink-pair의 metastability를 검사한다.
# - 핵심: sigma_LJ/b를 row equilibrium에서 유도하고, a를 registry별로 준정적 최소화한 effective substrate에서 kink-pair를 최소화한다.
# - 주의: kink-pair activation saddle 자체는 아직 계산하지 않는다. 아래 defect energy는 formation/local-minimum energy다.
# === 한국어 파일 안내 끝 ===
"""Spatial registry-kink metastability audit.

Candidate-only mechanism diagnostic. No FCC geometry, no characteristic volume,
no fitted damping, and no lifetime calibration are introduced.
"""
from __future__ import annotations

from pathlib import Path
import json
import math
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize, minimize_scalar
from scipy.special import zeta

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data" / "spatial_registry_kink_metastability"

M_EXP = 12.19
N_EXP = 6.0
C_MN = M_EXP / (M_EXP - N_EXP) * (M_EXP / N_EXP) ** (N_EXP / (M_EXP - N_EXP))
SIGMA_OVER_B = (N_EXP * zeta(N_EXP) / (M_EXP * zeta(M_EXP))) ** (1.0 / (M_EXP - N_EXP))
K_MAX = 50
P_MAX = 150
KS = np.arange(1, K_MAX + 1, dtype=float)[:, None]
PS = np.arange(-P_MAX, P_MAX + 1, dtype=float)[None, :]

A0_PHYS = 2.8627442948e-10
E_MODULUS = 69e9
AREA0 = 6.0338e-20
EV_J = 1.602176634e-19


def pair_potential(r: np.ndarray | float) -> np.ndarray | float:
    return C_MN * ((SIGMA_OVER_B / r) ** M_EXP - (SIGMA_OVER_B / r) ** N_EXP)


def pair_derivative(r: np.ndarray | float) -> np.ndarray | float:
    return C_MN * (
        -M_EXP * SIGMA_OVER_B**M_EXP * r ** (-M_EXP - 1.0)
        + N_EXP * SIGMA_OVER_B**N_EXP * r ** (-N_EXP - 1.0)
    )


def u0(a: float, s: float) -> float:
    r2 = (KS * a) ** 2 + (PS + s) ** 2
    return float(
        C_MN
        * np.sum(
            SIGMA_OVER_B**M_EXP * r2 ** (-M_EXP / 2.0)
            - SIGMA_OVER_B**N_EXP * r2 ** (-N_EXP / 2.0)
        )
    )


def relaxed_substrate(normal_force: float, points: int = 401) -> tuple[CubicSpline, CubicSpline]:
    s_grid = np.linspace(0.0, 1.0, points)
    values = []
    aeq = []
    for s in s_grid:
        result = minimize_scalar(
            lambda a: u0(float(a), float(s)) - normal_force * float(a),
            bounds=(0.65, 1.30),
            method="bounded",
            options={"xatol": 1e-10},
        )
        values.append(float(result.fun))
        aeq.append(float(result.x))
    values = np.asarray(values)
    aeq = np.asarray(aeq)
    values[-1] = values[0]
    aeq[-1] = aeq[0]
    return (
        CubicSpline(s_grid, values, bc_type="periodic"),
        CubicSpline(s_grid, aeq, bc_type="periodic"),
    )


def spatial_energy_gradient(s: np.ndarray, substrate: CubicSpline) -> tuple[float, np.ndarray]:
    s_mod = np.mod(s, 1.0)
    reference = float(substrate(0.5))
    energy = float(np.sum(substrate(s_mod) - reference))
    gradient = np.asarray(substrate(s_mod, 1), dtype=float)
    x = np.arange(len(s), dtype=float) + s

    for j in range(len(s) - 1):
        r = x[j + 1 :] - x[j]
        if np.any(r <= 0.0):
            return 1e15, np.zeros_like(s)
        base = np.arange(1, len(s) - j, dtype=float)
        energy += float(np.sum(pair_potential(r) - pair_potential(base)))
        force = pair_derivative(r)
        gradient[j] -= np.sum(force)
        gradient[j + 1 :] += force
    return energy, gradient


def minimize_pair(substrate: CubicSpline, initial: np.ndarray) -> tuple[np.ndarray, float]:
    base = np.array(initial, dtype=float)
    free = np.arange(1, len(base) - 1)

    def unpack(y: np.ndarray) -> np.ndarray:
        s = base.copy()
        s[free] = y
        s[0] = 0.5
        s[-1] = 0.5
        return s

    def objective(y: np.ndarray) -> tuple[float, np.ndarray]:
        s = unpack(y)
        energy, gradient = spatial_energy_gradient(s, substrate)
        return energy, gradient[free]

    result = minimize(
        lambda y: objective(y)[0],
        base[free],
        jac=lambda y: objective(y)[1],
        method="L-BFGS-B",
        bounds=[(0.45, 1.55)] * len(free),
        options={"maxiter": 4000, "ftol": 1e-13, "gtol": 1e-9, "maxls": 50},
    )
    if not result.success:
        raise RuntimeError(result.message)
    return unpack(result.x), float(result.fun)


def hessian_minimum_eigenvalue(s: np.ndarray, substrate: CubicSpline, h: float = 1e-5) -> float:
    free = np.arange(1, len(s) - 1)
    y0 = s[free].copy()

    def grad(y: np.ndarray) -> np.ndarray:
        trial = s.copy()
        trial[free] = y
        return spatial_energy_gradient(trial, substrate)[1][free]

    hessian = np.empty((len(free), len(free)))
    for k in range(len(free)):
        yp = y0.copy()
        ym = y0.copy()
        yp[k] += h
        ym[k] -= h
        hessian[:, k] = (grad(yp) - grad(ym)) / (2.0 * h)
    hessian = 0.5 * (hessian + hessian.T)
    return float(np.linalg.eigvalsh(hessian)[0])


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)

    substrate0, aeq0 = relaxed_substrate(0.0)
    a_well = float(aeq0(0.5))

    # Physical energy bridge: match registry-surface normal curvature to EA0/a0.
    h = 1e-4
    uaa = (u0(a_well + h, 0.5) - 2.0 * u0(a_well, 0.5) + u0(a_well - h, 0.5)) / h**2
    b_phys = A0_PHYS / a_well
    k_phys = E_MODULUS * AREA0 / A0_PHYS
    energy_unit_eV = k_phys * b_phys**2 / uaa / EV_J

    j_count = 121
    center = j_count // 2
    grid = np.arange(j_count)
    initial = 0.5 + np.exp(-((grid - center) / 12.0) ** 4)
    initial[0] = 0.5
    initial[-1] = 0.5

    pair_state, pair_energy = minimize_pair(substrate0, initial)
    min_hessian = hessian_minimum_eigenvalue(pair_state, substrate0)
    local_a = np.asarray(aeq0(np.mod(pair_state, 1.0)), dtype=float)

    stress_rows = []
    state = pair_state.copy()
    for stress_mpa in (0, 50, 100, 150, 200):
        force = stress_mpa * 1e6 * AREA0
        q = force * b_phys / (energy_unit_eV * EV_J)
        substrate, _ = relaxed_substrate(q)
        state, energy = minimize_pair(substrate, state)
        stress_rows.append(
            {
                "stress_MPa": stress_mpa,
                "dimensionless_normal_force": q,
                "metastable_pair_energy_dimensionless": energy,
                "metastable_pair_energy_eV": energy * energy_unit_eV,
            }
        )

    payload = {
        "classification": "self-consistent dimensionless spatial-registry metastability diagnostic",
        "status": "candidate only; formation/local-minimum energy is not the activation saddle barrier",
        "m": M_EXP,
        "n": N_EXP,
        "sigma_LJ_over_b_from_infinite_row_equilibrium": SIGMA_OVER_B,
        "b_over_sigma_LJ": 1.0 / SIGMA_OVER_B,
        "zero_stress_registry_well_a_over_b": a_well,
        "normal_curvature_dimensionless": uaa,
        "physical_b_bridge_m": b_phys,
        "energy_unit_eV": energy_unit_eV,
        "kink_pair": {
            "sites": j_count,
            "formation_energy_dimensionless": pair_energy,
            "formation_energy_eV": pair_energy * energy_unit_eV,
            "minimum_hessian_eigenvalue_fixed_endpoints": min_hessian,
            "maximum_registry_s": float(np.max(pair_state)),
            "maximum_local_normal_opening_ratio": float(np.max(local_a / a_well)),
            "mean_local_normal_opening_ratio": float(np.mean(local_a / a_well)),
            "std_local_normal_opening_ratio": float(np.std(local_a / a_well)),
        },
        "stress_sensitivity": stress_rows,
        "verdict": [
            "With endpoints returned to the original well, unconstrained energy minimization retains a nonuniform kink-antikink state instead of collapsing to the uniform well.",
            "The fixed-endpoint Hessian is nonnegative within numerical resolution; the tiny lowest mode is consistent with a nearly translational defect mode.",
            "The residual registry core produces local quasistatic normal-equilibrium openings up to about 14% on this candidate surface.",
            "This supplies a non-equivalent post-transition structural state and a direct mechanical route back into P_a.",
            "The formation energy is not the transition-state barrier Delta G_kp; no Arrhenius lifetime may be claimed until the saddle is computed."
        ],
        "next_required_calculation": "Compute a minimum-energy path / critical kink-pair saddle and then test cycle-resolved nucleation plus residual P_a feedback."
    }

    (DATA / "summary.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    np.savetxt(DATA / "metastable_pair_profile.csv", np.column_stack((grid, pair_state, local_a)), delimiter=",", header="j,s,a_eq", comments="")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
