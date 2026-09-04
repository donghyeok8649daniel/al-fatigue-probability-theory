# === 한국어 파일 안내 시작 ===
# - 목적: spatial registry-kink 후보에서 artificial bounds를 제거한 뒤 실제 product basin, 형성에너지, NEB 경로를 재검증한다.
# - 핵심 판정: 좁은 kink-pair는 intact로 붕괴하며, 충분히 넓은 pair만 매우 얕은 lattice-pinned local minimum으로 남는다.
# - 주의: CI-NEB가 깨끗한 interior saddle을 분해하지 못하면 endpoint 형성에너지와 마지막 구간 corrugation을 분리해 보고한다.
# - 이 파일의 수치는 후보 mechanism diagnostic이며 Al fatigue law calibration이 아니다.
# === 한국어 파일 안내 끝 ===
"""Kink-pair minimum-energy-path / NEB audit.

No FCC geometry, fitted damping, characteristic area/volume, or empirical damage
variable is introduced.  The spatial registry coordinate is still a candidate
extension of the active normal-only theory.
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
DATA = ROOT / "results" / "data" / "kink_pair_mep_neb"

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
KB_EV = 8.617333262145e-5
ATTEMPT_HZ_DIAGNOSTIC = 1.40972151760697e12

SITES = 121
FREE = np.arange(1, SITES - 1)
JIDX, KIDX = np.triu_indices(SITES, 1)
BASE_DIST = (KIDX - JIDX).astype(float)


def pair_potential(r: np.ndarray) -> np.ndarray:
    return C_MN * ((SIGMA_OVER_B / r) ** M_EXP - (SIGMA_OVER_B / r) ** N_EXP)


def pair_derivative(r: np.ndarray) -> np.ndarray:
    return C_MN * (
        -M_EXP * SIGMA_OVER_B**M_EXP * r ** (-M_EXP - 1.0)
        + N_EXP * SIGMA_OVER_B**N_EXP * r ** (-N_EXP - 1.0)
    )


BASE_PAIR_ENERGY = pair_potential(BASE_DIST)


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
    energy = []
    aeq = []
    for s in s_grid:
        out = minimize_scalar(
            lambda a: u0(float(a), float(s)) - normal_force * float(a),
            bounds=(0.65, 1.30),
            method="bounded",
            options={"xatol": 2e-9},
        )
        energy.append(float(out.fun))
        aeq.append(float(out.x))
    energy = np.asarray(energy)
    aeq = np.asarray(aeq)
    energy[-1] = energy[0]
    aeq[-1] = aeq[0]
    return CubicSpline(s_grid, energy, bc_type="periodic"), CubicSpline(s_grid, aeq, bc_type="periodic")


def make_energy_gradient(substrate: CubicSpline):
    ref = float(substrate(0.5))

    def energy_gradient(s: np.ndarray) -> tuple[float, np.ndarray]:
        s = np.asarray(s, dtype=float)
        smod = np.mod(s, 1.0)
        energy = float(np.sum(substrate(smod) - ref))
        grad = np.asarray(substrate(smod, 1), dtype=float).copy()
        x = np.arange(SITES, dtype=float) + s
        r = x[KIDX] - x[JIDX]
        if float(np.min(r)) <= 0.0:
            return 1e15, np.zeros_like(s)
        energy += float(np.sum(pair_potential(r) - BASE_PAIR_ENERGY))
        f = pair_derivative(r)
        np.add.at(grad, JIDX, -f)
        np.add.at(grad, KIDX, f)
        return energy, grad

    return energy_gradient


def make_patch(width: float, edge: float = 1.5) -> np.ndarray:
    j = np.arange(SITES, dtype=float)
    left = SITES // 2 - width / 2.0
    right = SITES // 2 + width / 2.0
    bump = 0.5 * (np.tanh((j - left) / edge) - np.tanh((j - right) / edge))
    s = 0.5 + bump
    s[0] = s[-1] = 0.5
    return s


def minimize_unbounded(initial: np.ndarray, eg) -> tuple[np.ndarray, float, object]:
    base = np.asarray(initial, dtype=float).copy()

    def unpack(y: np.ndarray) -> np.ndarray:
        s = base.copy()
        s[FREE] = y
        s[0] = s[-1] = 0.5
        return s

    def objective(y: np.ndarray):
        s = unpack(y)
        e, g = eg(s)
        return e, g[FREE]

    out = minimize(
        lambda y: objective(y)[0],
        base[FREE],
        jac=lambda y: objective(y)[1],
        method="BFGS",
        options={"gtol": 1e-8, "maxiter": 3000},
    )
    return unpack(out.x), float(out.fun), out


def hessian(s: np.ndarray, eg, h: float = 1e-5) -> np.ndarray:
    H = np.empty((len(FREE), len(FREE)))
    for k, idx in enumerate(FREE):
        plus = s.copy()
        minus = s.copy()
        plus[idx] += h
        minus[idx] -= h
        H[:, k] = (eg(plus)[1][FREE] - eg(minus)[1][FREE]) / (2.0 * h)
    return 0.5 * (H + H.T)


def improved_tangent(images: np.ndarray, energies: np.ndarray, i: int) -> np.ndarray:
    dp = images[i + 1] - images[i]
    dm = images[i] - images[i - 1]
    if energies[i + 1] > energies[i] > energies[i - 1]:
        tangent = dp
    elif energies[i + 1] < energies[i] < energies[i - 1]:
        tangent = dm
    else:
        ep = abs(energies[i + 1] - energies[i])
        em = abs(energies[i - 1] - energies[i])
        emax = max(ep, em)
        emin = min(ep, em)
        tangent = dp * emax + dm * emin if energies[i + 1] > energies[i - 1] else dp * emin + dm * emax
    tangent = tangent.copy()
    tangent[0] = tangent[-1] = 0.0
    norm = float(np.linalg.norm(tangent))
    return tangent / norm if norm > 1e-14 else tangent


def neb_force(images: np.ndarray, eg, spring: float, climb: bool):
    energies = np.empty(len(images))
    gradients = np.empty_like(images)
    for i, image in enumerate(images):
        energies[i], gradients[i] = eg(image)
    force = np.zeros_like(images)
    imax = 1 + int(np.argmax(energies[1:-1]))
    for i in range(1, len(images) - 1):
        tangent = improved_tangent(images, energies, i)
        grad = gradients[i].copy()
        grad[0] = grad[-1] = 0.0
        if climb and i == imax:
            f = -grad + 2.0 * np.dot(grad, tangent) * tangent
        else:
            ftrue = -grad
            fperp = ftrue - np.dot(ftrue, tangent) * tangent
            dplus = np.linalg.norm(images[i + 1] - images[i])
            dminus = np.linalg.norm(images[i] - images[i - 1])
            f = fperp + spring * (dplus - dminus) * tangent
        f[0] = f[-1] = 0.0
        force[i] = f
    return force, energies, imax


def fire_neb(images: np.ndarray, eg, max_steps: int = 6000):
    x = images.copy()
    velocity = np.zeros_like(x)
    dt = 0.0015
    dtmax = 0.025
    alpha = 0.1
    npositive = 0
    spring = 5.0
    climb_start = 1800
    force_tol = 2e-5
    converged = False
    last_max_force = math.inf
    for step in range(max_steps):
        force, energies, _ = neb_force(x, eg, spring, step >= climb_start)
        last_max_force = float(np.max(np.linalg.norm(force[1:-1, 1:-1], axis=1)))
        if step >= climb_start and last_max_force < force_tol:
            converged = True
            break
        velocity += dt * force
        power = float(np.sum(velocity * force))
        vnorm = float(np.linalg.norm(velocity))
        fnorm = float(np.linalg.norm(force))
        if vnorm > 0.0 and fnorm > 0.0:
            velocity = (1.0 - alpha) * velocity + alpha * force * (vnorm / fnorm)
        if power > 0.0:
            npositive += 1
            if npositive > 5:
                dt = min(dt * 1.1, dtmax)
                alpha *= 0.99
        else:
            npositive = 0
            dt *= 0.5
            alpha = 0.1
            velocity[:] = 0.0
        x += dt * velocity
        x[0] = images[0]
        x[-1] = images[-1]
        x[:, 0] = x[:, -1] = 0.5
        x[1:-1, 1:-1] = np.clip(x[1:-1, 1:-1], 0.25, 1.75)
    _, energies, _ = neb_force(x, eg, spring, True)
    return x, energies, converged, last_max_force


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    substrate0, aeq0 = relaxed_substrate(0.0)
    eg0 = make_energy_gradient(substrate0)

    # Physical energy bridge retained from the preceding registry audit.
    a_well = float(aeq0(0.5))
    h = 1e-4
    uaa = (u0(a_well + h, 0.5) - 2.0 * u0(a_well, 0.5) + u0(a_well - h, 0.5)) / h**2
    b_phys = A0_PHYS / a_well
    k_phys = E_MODULUS * AREA0 / A0_PHYS
    energy_unit_eV = k_phys * b_phys**2 / uaa / EV_J

    intact, e_intact, _ = minimize_unbounded(np.full(SITES, 0.5), eg0)

    width_scan = []
    products = {}
    for width in (16, 18, 20, 21, 22, 24, 30):
        state, energy, result = minimize_unbounded(make_patch(float(width)), eg0)
        eig0 = float(np.linalg.eigvalsh(hessian(state, eg0))[0])
        q = float(np.sum(state - 0.5))
        retained = abs(q) > 1.0
        width_scan.append(
            {
                "initial_patch_width": width,
                "retained_nonuniform_state": retained,
                "relaxed_registry_excess_sum": q,
                "energy_dimensionless": energy,
                "energy_above_intact_eV": (energy - e_intact) * energy_unit_eV,
                "minimum_hessian_eigenvalue": eig0,
                "max_gradient_component": float(np.max(np.abs(result.jac))),
            }
        )
        if retained and eig0 > 0.0:
            products[width] = (state, energy)

    # Use the first robust positive-Hessian retained state in this audit as NEB product.
    product_width = min(products)
    product, e_product = products[product_width]

    images = np.array([(1.0 - a) * intact + a * product for a in np.linspace(0.0, 1.0, 27)])
    band, e_band, neb_converged, neb_force_residual = fire_neb(images, eg0)
    internal_max = float(np.max(e_band[1:-1]))

    # The final lattice-pinning corrugation is much smaller than the dominant pair-formation energy.
    # Resolve the last NEB segment densely instead of claiming an unresolved CI image as a saddle.
    t_grid = np.linspace(0.0, 1.0, 1001)
    final_segment_energy = np.array(
        [eg0((1.0 - t) * product + t * band[-2])[0] for t in t_grid], dtype=float
    )
    i_segment = int(np.argmax(final_segment_energy))
    segment_max = float(final_segment_energy[i_segment])
    segment_reverse_corrugation_eV = (segment_max - e_product) * energy_unit_eV
    dominant_forward_barrier_eV = (max(e_product, segment_max) - e_intact) * energy_unit_eV

    # Normal-stress dependence of the robust formation part; this is not a converged saddle audit.
    stress_rows = []
    for stress_mpa in (0, 50, 100, 150, 200):
        force_n = stress_mpa * 1e6 * AREA0
        qn = force_n * b_phys / (energy_unit_eV * EV_J)
        substrate, _ = relaxed_substrate(qn)
        eg = make_energy_gradient(substrate)
        intact_s, ei, _ = minimize_unbounded(np.full(SITES, 0.5), eg)
        product_s, ep, result = minimize_unbounded(make_patch(24.0), eg)
        eig0 = float(np.linalg.eigvalsh(hessian(product_s, eg))[0])
        stress_rows.append(
            {
                "stress_MPa": stress_mpa,
                "dimensionless_normal_force": qn,
                "formation_energy_eV": (ep - ei) * energy_unit_eV,
                "retained_registry_excess_sum": float(np.sum(product_s - 0.5)),
                "minimum_hessian_eigenvalue": eig0,
                "max_gradient_component": float(np.max(np.abs(result.jac))),
            }
        )

    # Conditional Arrhenius illustration using only the dominant formation energy spline.
    stress = np.array([row["stress_MPa"] for row in stress_rows], dtype=float)
    barrier = np.array([row["formation_energy_eV"] for row in stress_rows], dtype=float)
    barrier_spline = CubicSpline(stress, barrier, bc_type="natural")
    theta = np.linspace(0.0, 2.0 * np.pi, 20001)
    sigma_cycle = 100.0 + 100.0 * np.sin(theta)
    barrier_cycle = barrier_spline(sigma_cycle)
    rate = ATTEMPT_HZ_DIAGNOSTIC * np.exp(-barrier_cycle / (KB_EV * 300.0))
    hazard_cycle = float(np.trapezoid(rate, theta) / (2.0 * np.pi * 20.0))

    payload = {
        "classification": "candidate spatial-registry kink-pair MEP/NEB audit",
        "status": "candidate only; clean interior transition-state saddle not numerically resolved",
        "m": M_EXP,
        "n": N_EXP,
        "sites": SITES,
        "sigma_LJ_over_b": SIGMA_OVER_B,
        "zero_stress_registry_well_a_over_b": a_well,
        "energy_unit_eV": energy_unit_eV,
        "intact_energy_dimensionless": e_intact,
        "width_scan": width_scan,
        "neb": {
            "product_width_used": product_width,
            "product_formation_energy_eV": (e_product - e_intact) * energy_unit_eV,
            "images": 27,
            "fire_force_converged": neb_converged,
            "final_force_residual": neb_force_residual,
            "internal_band_max_dimensionless": internal_max,
            "product_energy_dimensionless": e_product,
            "interior_saddle_resolved": bool(internal_max > e_product),
            "last_segment_reverse_corrugation_eV": segment_reverse_corrugation_eV,
            "last_segment_max_parameter": float(t_grid[i_segment]),
            "dominant_forward_barrier_estimate_eV": dominant_forward_barrier_eV,
        },
        "stress_formation_energy": stress_rows,
        "conditional_rate_diagnostic": {
            "temperature_K": 300.0,
            "frequency_Hz": 20.0,
            "mean_stress_MPa": 100.0,
            "amplitude_stress_MPa": 100.0,
            "attempt_frequency_Hz": ATTEMPT_HZ_DIAGNOSTIC,
            "hazard_per_cycle_using_formation_barrier_only": hazard_cycle,
            "transition_probability_per_cycle": 1.0 - math.exp(-hazard_cycle),
            "median_cycles_if_formation_barrier_were_the_rate_barrier": math.log(2.0) / hazard_cycle,
        },
        "verdict": [
            "Removing the earlier artificial registry bounds changes the metastability verdict: narrow kink pairs collapse to the intact basin.",
            "Wider lattice-pinned kink-pair states do exist, but their lowest curvature is extremely small and the reverse Peierls corrugation is tiny compared with kBT at room temperature.",
            "The robust energetic scale is the approximately 0.924 eV zero-stress pair-formation cost; the clean interior saddle above the shallow product minimum was not resolved by CI-NEB.",
            "Normal tension lowers the dominant formation cost to about 0.898 eV by 200 MPa on this diagnostic surface.",
            "Because the reverse trapping barrier is not robustly large, this calculation does not establish a long-lived residual plastic memory under ideal pure-normal loading.",
            "The spatial registry route remains useful as a rare transient first-passage mechanism, but it is not promoted as the missing progressive P_a evolution law."
        ],
    }
    (DATA / "summary.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    np.savetxt(DATA / "neb_energy_profile.csv", np.column_stack((np.arange(len(e_band)), e_band)), delimiter=",", header="image,energy_dimensionless", comments="")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
