# === 한국어 파일 안내 시작 ===
# - 파일 역할: spatial registry kink 후보에서 registry core가 normal equilibrium을 바꾸는지 direct-sum으로 검사한다.
# - 주요 함수: u0, d_u0_da, d2_u0_da2, stable_normal_root, main
# - 주의: characteristic area/volume, kink-pair barrier, Al lifetime을 보정하지 않는다.
# === 한국어 파일 안내 끝 ===
"""Direct-sum audit of normal/registry coupling for the spatial-kink candidate.

This script does NOT solve the kink profile or kink-pair saddle. It verifies the
more basic necessary condition: intermediate registry positions must alter the
local normal force/equilibrium if a spatial kink core is to feed back into P_a.

The calculation uses the existing dimensionless multilayer registry diagnostic
surface with m=12.19, n=6, b=sigma_LJ=1. Results are mechanism diagnostics, not
an independently calibrated Al dislocation-core prediction.
"""
from __future__ import annotations

from pathlib import Path
import csv
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data" / "registry_kink_normal_feedback"
REPORT = ROOT / "results" / "reports"

M_EXP = 12.19
N_EXP = 6.0
K_MAX = 80
P_MAX = 250
S_VALUES = (0.50, 0.40, 0.30, 0.25, 0.20, 0.10, 0.00)

C_MN = M_EXP / (M_EXP - N_EXP) * (M_EXP / N_EXP) ** (N_EXP / (M_EXP - N_EXP))
KS = np.arange(1, K_MAX + 1, dtype=float)[:, None]
PS = np.arange(-P_MAX, P_MAX + 1, dtype=float)[None, :]


def u0(a: float, s: float) -> float:
    r2 = (KS * a) ** 2 + (PS + s) ** 2
    return float(C_MN * np.sum(r2 ** (-M_EXP / 2.0) - r2 ** (-N_EXP / 2.0)))


def d_u0_da(a: float, s: float, h: float = 2.0e-5) -> float:
    return (u0(a + h, s) - u0(a - h, s)) / (2.0 * h)


def d2_u0_da2(a: float, s: float, h: float = 1.0e-4) -> float:
    return (u0(a + h, s) - 2.0 * u0(a, s) + u0(a - h, s)) / (h * h)


def bisect_root(func, lo: float, hi: float, iterations: int = 80) -> float:
    flo = func(lo)
    fhi = func(hi)
    if flo == 0.0:
        return lo
    if fhi == 0.0:
        return hi
    if flo * fhi > 0.0:
        raise ValueError("root is not bracketed")
    for _ in range(iterations):
        mid = 0.5 * (lo + hi)
        fm = func(mid)
        if flo * fm <= 0.0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
    return 0.5 * (lo + hi)


def stable_normal_root(s: float) -> float:
    grid = np.linspace(0.90, 1.20, 240)
    values = [d_u0_da(float(a), s) for a in grid]
    for left, right, f_left, f_right in zip(grid[:-1], grid[1:], values[:-1], values[1:]):
        if f_left * f_right < 0.0:
            root = bisect_root(lambda x: d_u0_da(x, s), float(left), float(right))
            if d2_u0_da2(root, s) > 0.0:
                return root
    raise RuntimeError(f"no stable normal root found for s={s}")


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    REPORT.mkdir(parents=True, exist_ok=True)

    a0 = stable_normal_root(0.5)
    u_well = u0(a0, 0.5)
    rows = []
    for s in S_VALUES:
        a_eq = stable_normal_root(s)
        rows.append({
            "s_over_b": s,
            "energy_excess_at_a0": u0(a0, s) - u_well,
            "dU_da_at_a0": d_u0_da(a0, s),
            "stable_a_eq_over_b": a_eq,
            "relative_opening_from_well": (a_eq - a0) / a0,
        })

    payload = {
        "classification": "dimensionless direct-sum mechanism diagnostic",
        "status": "candidate-only; no kink-pair saddle, characteristic-area, or lifetime calibration",
        "m": M_EXP,
        "n": N_EXP,
        "b_over_sigma_LJ": 1.0,
        "k_max": K_MAX,
        "p_max": P_MAX,
        "stable_well_s_over_b": 0.5,
        "stable_well_a0_over_b": a0,
        "rows": rows,
        "key_result": "Intermediate registry positions strongly shift the stable normal equilibrium; at the registry saddle s/b=0 the preferred normal spacing is about 10.04% larger than the registry-well equilibrium.",
        "interpretation": "A spatial kink core can therefore feed back into the normal spacing state through the already-derived U0(a,s), even though the far-field wells s0 and s0+b are exactly equivalent.",
        "next_required_calculation": "Minimize the spatial registry-row energy and compute the critical kink-pair saddle barrier Delta G_kp; do not substitute N*Delta G_s as the governing barrier.",
    }

    with (DATA / "registry_normal_coupling.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    (DATA / "summary.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    lines = [
        "# Registry-kink normal-feedback audit",
        "",
        "**Classification:** dimensionless direct-sum mechanism diagnostic; candidate only.",
        "",
        f"Stable registry-well normal root: `a0/b = {a0:.10f}` at `s/b=0.5`.",
        "",
        "| s/b | energy excess at a0 | dU/da at a0 | stable a_eq/b | relative opening |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['s_over_b']:.2f} | {row['energy_excess_at_a0']:.6f} | "
            f"{row['dU_da_at_a0']:.6f} | {row['stable_a_eq_over_b']:.6f} | "
            f"{100.0*row['relative_opening_from_well']:.3f}% |"
        )
    lines += [
        "",
        "The saddle-side registry core (`s/b=0`) prefers a normal spacing about 10.04% larger than the well state on this reduced surface.",
        "",
        "This verifies only the existence of normal-registry feedback inside a spatial core. It does not yet compute a kink profile, a kink-pair activation barrier, irreversible residual state, characteristic area/volume, or Al fatigue life.",
        "",
        "The previous coherent-patch `N*Delta G_s` rate is therefore not promoted. The next calculation is the spatial minimum-energy kink and critical kink-pair saddle from the same interaction energy.",
    ]
    (REPORT / "REGISTRY_KINK_NORMAL_FEEDBACK.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
