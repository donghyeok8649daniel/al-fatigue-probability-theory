# === 한국어 파일 안내 시작 ===
# - 파일 역할: FCC {111}<110> slip-system geometry를 검증하고, 단축응력에서
#   slip-plane normal opening a와 in-plane registry s에 공액인 일반화힘을 계산한다.
# - 주요 함수: fcc_111_spacing, fcc_perfect_burgers_magnitude,
#   fcc_111_surface_primitive_area, validate_fcc_slip_system, schmid_factor,
#   uniaxial_plane_generalized_forces, two_rigid_plane_mass_metric.
# - 주의: 이 모듈은 기하/가상일만 다룬다. 현재 row-based U0(a,s)를 full FCC
#   stacking energy라고 가정하지 않는다.
# === 한국어 파일 안내 끝 ===
"""FCC slip-plane kinematics for the active spacing--registry coordinates.

This module makes one geometric interpretation explicit:

    a : relative opening normal to a declared slip plane n,
    s : relative translation along a declared in-plane slip direction d.

For a uniaxial Cauchy stress sigma * l⊗l acting on that plane, virtual work
through a patch of area A gives

    Q_a = A sigma (l·n)^2,
    Q_s = A sigma (l·n)(l·d).

The second factor is the signed Schmid factor.  Therefore the older pair

    Q_a = A sigma,
    Q_s = A M sigma

cannot both arise from one literal plane-opening/plane-slip displacement when
M != 0.  Q_a=A sigma is recovered only for n parallel to l, where every
in-plane d has l·d=0 and hence M=0.

This is a kinematic/virtual-work statement.  It does not by itself assert that
the current multiplicity-free row energy U0(a,s) is the exact three-dimensional
FCC energy of that plane patch.  FCC stacking and energy counting require a
separate consistency derivation.
"""

from __future__ import annotations

import math
import numpy as np


def _unit(values: tuple[int, int, int] | np.ndarray, name: str) -> np.ndarray:
    v = np.asarray(values, dtype=float)
    if v.shape != (3,) or not np.all(np.isfinite(v)):
        raise ValueError(f"{name} must contain three finite components")
    norm = float(np.linalg.norm(v))
    if norm == 0.0:
        raise ValueError(f"{name} cannot be the zero vector")
    return v / norm


def fcc_111_spacing(lattice_parameter_m: float) -> float:
    """Return d_111=a_lat/sqrt(3) for a cubic lattice."""
    if not math.isfinite(lattice_parameter_m) or lattice_parameter_m <= 0.0:
        raise ValueError("lattice_parameter_m must be positive and finite")
    return lattice_parameter_m / math.sqrt(3.0)


def fcc_perfect_burgers_magnitude(lattice_parameter_m: float) -> float:
    """Return |a_lat/2 <110>| = a_lat/sqrt(2)."""
    if not math.isfinite(lattice_parameter_m) or lattice_parameter_m <= 0.0:
        raise ValueError("lattice_parameter_m must be positive and finite")
    return lattice_parameter_m / math.sqrt(2.0)


def fcc_111_surface_primitive_area(lattice_parameter_m: float) -> float:
    """Return area of one triangular (111) surface primitive cell.

    The in-plane primitive translations have length a_lat/sqrt(2) and mutual
    angle 60 degrees, giving A_111=sqrt(3)*a_lat^2/4.
    """
    if not math.isfinite(lattice_parameter_m) or lattice_parameter_m <= 0.0:
        raise ValueError("lattice_parameter_m must be positive and finite")
    return math.sqrt(3.0) * lattice_parameter_m**2 / 4.0


def validate_fcc_slip_system(
    plane_hkl: tuple[int, int, int],
    direction_uvw: tuple[int, int, int],
) -> tuple[np.ndarray, np.ndarray]:
    """Validate one perfect-FCC {111}<110> slip system and return unit vectors."""
    hkl_abs = sorted(abs(int(x)) for x in plane_hkl)
    uvw_abs = sorted(abs(int(x)) for x in direction_uvw)
    if hkl_abs != [1, 1, 1]:
        raise ValueError("perfect FCC slip plane must belong to {111}")
    if uvw_abs != [0, 1, 1]:
        raise ValueError("perfect FCC slip direction must belong to <110>")
    n = _unit(plane_hkl, "plane_hkl")
    d = _unit(direction_uvw, "direction_uvw")
    if abs(float(n @ d)) > 1.0e-12:
        raise ValueError("slip direction must lie in the declared slip plane")
    return n, d


def schmid_factor(
    load_hkl: tuple[int, int, int],
    plane_hkl: tuple[int, int, int],
    direction_uvw: tuple[int, int, int],
) -> float:
    """Return signed M=(l·n)(l·d) for uniaxial loading."""
    l = _unit(load_hkl, "load_hkl")
    n, d = validate_fcc_slip_system(plane_hkl, direction_uvw)
    return float((l @ n) * (l @ d))


def uniaxial_plane_generalized_forces(
    stress_pa: float,
    area_m2: float,
    load_hkl: tuple[int, int, int],
    plane_hkl: tuple[int, int, int],
    direction_uvw: tuple[int, int, int],
) -> tuple[float, float]:
    """Return (Q_a,Q_s) from exact virtual work on the declared plane patch.

    The stress tensor is sigma*l⊗l.  The traction on plane n is t=sigma(l·n)l.
    For relative displacement dr=n da+d ds,

        A t·dr = Q_a da + Q_s ds.
    """
    if not math.isfinite(stress_pa) or not math.isfinite(area_m2) or area_m2 <= 0.0:
        raise ValueError("stress must be finite and area_m2 positive")
    l = _unit(load_hkl, "load_hkl")
    n, d = validate_fcc_slip_system(plane_hkl, direction_uvw)
    ln = float(l @ n)
    ld = float(l @ d)
    q_a = area_m2 * stress_pa * ln * ln
    q_s = area_m2 * stress_pa * ln * ld
    return q_a, q_s


def two_rigid_plane_mass_metric(
    upper_patch_mass_kg: float,
    lower_patch_mass_kg: float,
    plane_hkl: tuple[int, int, int],
    direction_uvw: tuple[int, int, int],
) -> np.ndarray:
    """Return finite reduced mass metric for relative plane displacement.

    The two rigid patches move with their center of mass fixed.  If their
    relative displacement is r=a*n+s*d, the reduced mass is

        mu = M1*M2/(M1+M2)

    and G_ij=mu e_i·e_j for e_a=n, e_s=d.  For a valid slip system n·d=0,
    so G=mu I exactly.  This metric is a kinematic candidate for a finite plane
    patch; pairing its mass with the current row-based U0 still requires an
    energy-area mapping.
    """
    masses = (upper_patch_mass_kg, lower_patch_mass_kg)
    if not all(math.isfinite(x) and x > 0.0 for x in masses):
        raise ValueError("patch masses must be positive and finite")
    n, d = validate_fcc_slip_system(plane_hkl, direction_uvw)
    mu = upper_patch_mass_kg * lower_patch_mass_kg / (
        upper_patch_mass_kg + lower_patch_mass_kg
    )
    basis = np.column_stack((n, d))
    metric = mu * (basis.T @ basis)
    return metric
