"""Audit registry stability under normal-only cyclic loading.

This script stays inside the active reduced row/layer U0(a,s).  It does not use
FCC geometry, q_s forcing, damping, random noise, Boltzmann statistics, or a
Fokker--Planck closure.

It performs three checks:
1. direct-sum roots of U_aa=0 and U_ss=0 at the symmetric registry well;
2. normal-only spatial-chain history of K_s(t)=U_ss(a_i(t),s0);
3. finite-time amplification of an infinitesimal registry perturbation
   mu_s*xi_ddot+K_s(t)*xi=0, reported without choosing a seed amplitude.

Because the undamped normal chain is not exactly period-one after ramping, the
last result is a finite-time variational amplification audit, not a Floquet
material prediction.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

M_EXP = 12.0
N_EXP = 6.0
B = 1.0
EPSILON = 1.0
SIGMA_LJ = 1.0
KMAX = 12
PMAX = 30
CELLS = 48
DT = 0.004
OMEGA = 0.35
NORMAL_MEAN = 0.2
NORMAL_AMPLITUDE = 0.2
RAMP_CYCLES = 1
A0 = 0.9919601753795769
S0 = 0.5

C_MN = M_EXP / (M_EXP - N_EXP) * (M_EXP / N_EXP) ** (N_EXP / (M_EXP - N_EXP))
_K = np.arange(1, KMAX + 1, dtype=float)[:, None, None]
_P = np.arange(-PMAX, PMAX + 1, dtype=float)[None, None, :]


def fields(a: np.ndarray, s: float = S0) -> tuple[np.ndarray, np.ndarray]:
    """Return dU/da and U_ss cellwise with the declared direct truncation."""
    aa = np.asarray(a, dtype=float)
    x = _K * aa[None, :, None]
    y = _P * B + s
    r2 = x*x + y*y
    radial = C_MN * EPSILON * (
        -M_EXP * SIGMA_LJ**M_EXP * r2 ** (-0.5 * (M_EXP + 2.0))
        + N_EXP * SIGMA_LJ**N_EXP * r2 ** (-0.5 * (N_EXP + 2.0))
    )
    d_u_da = np.sum(radial * (x * _K), axis=(0, 2))

    def hss(q: float) -> np.ndarray:
        return np.sum(
            -q * r2 ** (-0.5 * (q + 2.0))
            + q * (q + 2.0) * y*y * r2 ** (-0.5 * (q + 4.0)),
            axis=(0, 2),
        )

    u_ss = C_MN * EPSILON * (
        SIGMA_LJ**M_EXP * hss(M_EXP) - SIGMA_LJ**N_EXP * hss(N_EXP)
    )
    return d_u_da, u_ss


def direct_curvature(a: float, kind: str, kmax: int, pmax: int) -> float:
    k = np.arange(1, kmax + 1, dtype=float)[:, None]
    p = np.arange(-pmax, pmax + 1, dtype=float)[None, :]
    y = p * B + S0
    r2 = (k*a)**2 + y*y

    def term(q: float) -> float:
        if kind == "ss":
            value = (
                -q * r2 ** (-(q + 2.0)/2.0)
                + q*(q + 2.0)*y*y*r2 ** (-(q + 4.0)/2.0)
            )
        elif kind == "aa":
            value = (
                -q*k*k*r2 ** (-(q + 2.0)/2.0)
                + q*(q + 2.0)*k**4*a*a*r2 ** (-(q + 4.0)/2.0)
            )
        else:
            raise ValueError("kind must be 'aa' or 'ss'")
        return float(np.sum(value))

    return C_MN * EPSILON * (
        SIGMA_LJ**M_EXP * term(M_EXP) - SIGMA_LJ**N_EXP * term(N_EXP)
    )


def bisect(kind: str, lower: float, upper: float, kmax: int, pmax: int) -> float:
    lo, hi = lower, upper
    flo = direct_curvature(lo, kind, kmax, pmax)
    fhi = direct_curvature(hi, kind, kmax, pmax)
    if flo*fhi >= 0.0:
        raise ValueError("root bracket has no sign change")
    for _ in range(80):
        mid = 0.5*(lo + hi)
        fm = direct_curvature(mid, kind, kmax, pmax)
        if flo*fm <= 0.0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
    return 0.5*(lo + hi)


def normal_only_kss_history(cycles: int = 5) -> tuple[np.ndarray, int]:
    period = 2.0*math.pi/OMEGA
    steps = int(round(cycles*period/DT))
    steps_per_cycle = int(round(period/DT))
    x = np.arange(CELLS + 1, dtype=float) * A0
    velocity = np.zeros(CELLS + 1)

    def envelope(t: float) -> float:
        ramp_time = RAMP_CYCLES * period
        if t >= ramp_time:
            return 1.0
        return 0.5*(1.0 - math.cos(math.pi*t/ramp_time))

    def q_a(t: float) -> float:
        return envelope(t)*(NORMAL_MEAN + NORMAL_AMPLITUDE*math.sin(OMEGA*t))

    def acceleration(t: float) -> tuple[np.ndarray, np.ndarray]:
        spacing = np.diff(x)
        grad_a, kss = fields(spacing)
        acc = np.zeros(CELLS + 1)
        acc[1:CELLS] = grad_a[1:] - grad_a[:-1]
        acc[CELLS] = -grad_a[-1] + q_a(t)
        return acc, kss

    acc, kss = acceleration(0.0)
    history = np.empty((steps + 1, CELLS), dtype=np.float32)
    history[0] = kss
    for step in range(steps):
        t = step*DT
        velocity[1:] += 0.5*DT*acc[1:]
        x[1:] += DT*velocity[1:]
        new_acc, kss = acceleration(t + DT)
        velocity[1:] += 0.5*DT*new_acc[1:]
        acc = new_acc
        history[step + 1] = kss
    return history, steps_per_cycle


def finite_time_amplification(k_history: np.ndarray, mu: float) -> tuple[float, float, int]:
    """Worst cell amplification over the full history in an initial energy norm.

    The fundamental matrix maps [xi,xi_dot] at t0 to the final state.  The
    singular value is evaluated after scaling coordinates by
    [sqrt(K_s(t0)), sqrt(mu)].  The spectral radius is reported separately.
    No perturbation amplitude is chosen.
    """
    if mu <= 0.0:
        raise ValueError("mu must be positive")
    cells = k_history.shape[1]
    xi = np.zeros((cells, 2))
    vel = np.zeros((cells, 2))
    xi[:, 0] = 1.0
    vel[:, 1] = 1.0
    for step in range(k_history.shape[0] - 1):
        k0 = k_history[step, :, None].astype(float)
        k1 = k_history[step + 1, :, None].astype(float)
        vel -= 0.5*DT*(k0/mu)*xi
        xi += DT*vel
        vel -= 0.5*DT*(k1/mu)*xi

    k_initial = float(np.mean(k_history[0]))
    scale = np.diag([math.sqrt(k_initial), math.sqrt(mu)])
    inv_scale = np.diag([1.0/math.sqrt(k_initial), 1.0/math.sqrt(mu)])
    best_sv = 0.0
    best_rho = 0.0
    best_cell = 0
    for cell in range(cells):
        matrix = np.array([
            [xi[cell, 0], xi[cell, 1]],
            [vel[cell, 0], vel[cell, 1]],
        ])
        normalized = scale @ matrix @ inv_scale
        sv = float(np.linalg.svd(normalized, compute_uv=False)[0])
        rho = float(np.max(np.abs(np.linalg.eigvals(matrix))))
        if sv > best_sv:
            best_sv, best_cell = sv, cell
        best_rho = max(best_rho, rho)
    return best_sv, best_rho, best_cell


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    output = Path("results/data/registry_symmetry_stability")

    convergence = []
    for kmax, pmax in [(20, 50), (40, 100), (80, 200), (120, 300), (200, 500)]:
        normal_root = bisect("aa", 1.05, 1.20, kmax, pmax)
        registry_root = bisect("ss", 1.20, 1.50, kmax, pmax)
        convergence.append({
            "kmax": kmax,
            "pmax": pmax,
            "normal_curvature_zero_a": normal_root,
            "registry_curvature_zero_a": registry_root,
            "relative_gap": (registry_root - normal_root)/normal_root,
        })
    write_csv(output / "curvature_convergence.csv", convergence)

    k_history, steps_per_cycle = normal_only_kss_history(cycles=5)
    mismatch = []
    for cycle in range(1, 5):
        previous = k_history[(cycle - 1)*steps_per_cycle:cycle*steps_per_cycle]
        current = k_history[cycle*steps_per_cycle:(cycle + 1)*steps_per_cycle]
        rms = float(np.sqrt(np.mean((current - previous)**2)))
        relative = rms / float(np.sqrt(np.mean(previous**2)))
        mismatch.append({
            "cycle_pair": f"{cycle}->{cycle + 1}",
            "kss_rms_mismatch": rms,
            "relative_rms_mismatch": relative,
        })
    write_csv(output / "cycle_mismatch.csv", mismatch)

    scan = []
    for mu in [1, 10, 100, 500, 650, 700, 750, 800, 812.5, 850, 900, 1000, 1500, 2000]:
        sv, rho, cell = finite_time_amplification(k_history, float(mu))
        scan.append({
            "registry_inertia_mu": mu,
            "max_energy_norm_singular_amplification": sv,
            "max_transfer_spectral_radius": rho,
            "cell_at_max_singular_amplification": cell,
        })
    write_csv(output / "finite_time_amplification.csv", scan)
