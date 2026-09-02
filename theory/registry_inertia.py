# === 한국어 파일 안내 시작 ===
# - 파일 역할: 현재 reduced row/repeat registry coordinate s의 관성척도를 kinetic energy에서 유도한다.
# - FCC/결정구조를 사용하지 않는다. 현재 W,U0의 'per upper atom/repeat' counting convention만 사용한다.
# - 핵심 함수: reduced_relative_mass, registry_inertia_ratio,
#   registry_frequency_from_curvature_ratio, principal_parametric_required_inertia_ratio.
# === 한국어 파일 안내 끝 ===
"""Registry inertia for the active reduced row/repeat coordinate.

The active row kernel W(d,s) is an interaction energy *per upper atom/repeat*.
The multilayer U0(a,s) keeps the same per-reference-repeat counting while the
same physical registry displacement s enters the interaction with each
background row.

A finite kinetic interpretation consistent with that counting is therefore a
reference repeat (or a finite coherent patch) translating by the physical
length s relative to its surrounding background.  For a reference repeat mass
m_r and an effective participating background mass M_b, the relative coordinate
s=y_r-y_b has reduced mass

    mu_s = m_r M_b / (m_r + M_b).

A frozen/heavy background is the limit M_b -> infinity, hence mu_s = m_r.
For two equal moving repeats, mu_s = m_r/2.  Thus the simplest local relative
coordinate gives 0 < mu_s/m_r <= 1; mu_s is not an arbitrary large inertia.

If a coherent patch contains N identical repeats, both its registry energy and
its mass scale by N.  Therefore the small-oscillation frequency is unchanged;
one must not multiply the inertia by N while keeping a one-repeat U0.

This module also gives a calibration-compatible frequency mapping.  Let the
active normalized U0 have equilibrium curvature ratio

    r_K = U_ss / U_aa.

If the normal stiffness is calibrated by K_a = E A0 / a0 and

    t0 = sqrt(m_r a0/(E A0)),

then for rho_mu = mu_s/m_r,

    omega_s t0 = sqrt(r_K / rho_mu),

without requiring a separate FCC geometry or an explicit value of b.
"""

from __future__ import annotations

import math


def reduced_relative_mass(reference_mass: float, background_mass: float) -> float:
    """Return the reduced mass for s=y_ref-y_bg.

    Both masses must be positive and finite.  A frozen-background limit is
    represented separately by ``fixed_background_registry_inertia``.
    """
    if not (math.isfinite(reference_mass) and reference_mass > 0.0):
        raise ValueError("reference_mass must be positive and finite")
    if not (math.isfinite(background_mass) and background_mass > 0.0):
        raise ValueError("background_mass must be positive and finite")
    return reference_mass * background_mass / (reference_mass + background_mass)


def fixed_background_registry_inertia(reference_mass: float) -> float:
    """Return mu_s=m_r when the surrounding background is fixed/heavy."""
    if not (math.isfinite(reference_mass) and reference_mass > 0.0):
        raise ValueError("reference_mass must be positive and finite")
    return float(reference_mass)


def registry_inertia_ratio(reference_mass: float, registry_inertia: float) -> float:
    """Return rho_mu=mu_s/m_r."""
    if not (math.isfinite(reference_mass) and reference_mass > 0.0):
        raise ValueError("reference_mass must be positive and finite")
    if not (math.isfinite(registry_inertia) and registry_inertia > 0.0):
        raise ValueError("registry_inertia must be positive and finite")
    return float(registry_inertia / reference_mass)


def registry_frequency_from_curvature_ratio(
    curvature_ratio: float,
    normal_atomic_time_scale_s: float,
    *,
    inertia_ratio: float = 1.0,
) -> float:
    """Return registry small-oscillation frequency in Hz.

    Uses

        omega_s t0 = sqrt((U_ss/U_aa)/(mu_s/m_r)).

    ``curvature_ratio`` must be positive and ``inertia_ratio`` is rho_mu.
    """
    if not (math.isfinite(curvature_ratio) and curvature_ratio > 0.0):
        raise ValueError("curvature_ratio must be positive and finite")
    if not (
        math.isfinite(normal_atomic_time_scale_s)
        and normal_atomic_time_scale_s > 0.0
    ):
        raise ValueError("normal_atomic_time_scale_s must be positive and finite")
    if not (math.isfinite(inertia_ratio) and inertia_ratio > 0.0):
        raise ValueError("inertia_ratio must be positive and finite")
    omega = math.sqrt(curvature_ratio / inertia_ratio) / normal_atomic_time_scale_s
    return omega / (2.0 * math.pi)


def principal_parametric_required_inertia_ratio(
    curvature_ratio: float,
    normal_atomic_time_scale_s: float,
    loading_frequency_hz: float,
) -> float:
    """Return rho_mu required for principal stiffness-parametric resonance.

    For the linearized equation mu_s xi_ddot + K_s(t) xi = 0, the principal
    small-modulation resonance condition is approximately

        omega_load = 2 omega_s.

    With omega_s t0=sqrt(r_K/rho_mu), this gives

        rho_mu = r_K / (pi f_load t0)^2.

    This function is a diagnostic of scale compatibility, not a claim that the
    coefficient is exactly Mathieu-periodic.
    """
    if not (math.isfinite(loading_frequency_hz) and loading_frequency_hz > 0.0):
        raise ValueError("loading_frequency_hz must be positive and finite")
    if not (math.isfinite(curvature_ratio) and curvature_ratio > 0.0):
        raise ValueError("curvature_ratio must be positive and finite")
    if not (
        math.isfinite(normal_atomic_time_scale_s)
        and normal_atomic_time_scale_s > 0.0
    ):
        raise ValueError("normal_atomic_time_scale_s must be positive and finite")
    denom = math.pi * loading_frequency_hz * normal_atomic_time_scale_s
    return float(curvature_ratio / (denom * denom))


def coherent_patch_frequency_scale(
    per_repeat_curvature: float,
    repeat_mass: float,
    repeats: int,
) -> float:
    """Return omega for a coherent N-repeat patch with extensive energy/mass.

    Since K_N=N K_1 and M_N=N m_1,

        omega_N=sqrt(K_N/M_N)=sqrt(K_1/m_1).

    The explicit ``repeats`` argument guards against the inconsistent operation
    of increasing mass without also increasing the per-repeat energy.
    """
    if not (math.isfinite(per_repeat_curvature) and per_repeat_curvature > 0.0):
        raise ValueError("per_repeat_curvature must be positive and finite")
    if not (math.isfinite(repeat_mass) and repeat_mass > 0.0):
        raise ValueError("repeat_mass must be positive and finite")
    if not isinstance(repeats, int) or repeats < 1:
        raise ValueError("repeats must be a positive integer")
    total_k = repeats * per_repeat_curvature
    total_m = repeats * repeat_mass
    return math.sqrt(total_k / total_m)
