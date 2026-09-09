"""Exact FCC(111) geometry for the full infinite-plane-stack reference.

This module contains crystallography only.  ``atomic_cell_area`` is the area
of one primitive surface cell and is unrelated to the statistical specimen
correlation area ``A_c`` or to any FEM/FVM element area.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


DIRECT_110 = "direct_110"
SHOCKLEY_112 = "shockley_112"
REGISTRY_PATH_IDS = (DIRECT_110, SHOCKLEY_112)


@dataclass(frozen=True)
class RegistryPath:
    """One explicit scalar path through the two-dimensional registry plane."""

    path_id: str
    crystallographic_label: str
    direction_2d: tuple[float, float]
    direction_3d: tuple[float, float, float]
    period_over_b: float
    partial_increment_over_b: float | None
    calibration_status: str

    def direction(self) -> np.ndarray:
        return np.asarray(self.direction_2d, dtype=float)


@dataclass(frozen=True)
class FCC111Geometry:
    """Orthonormal basis and direct/reciprocal triangular primitive cells."""

    lattice_constant: float
    e1: np.ndarray
    e2: np.ndarray
    e3: np.ndarray
    b: float
    h111: float
    a1: np.ndarray
    a2: np.ndarray
    reciprocal_b1: np.ndarray
    reciprocal_b2: np.ndarray
    atomic_cell_area: float
    tau: np.ndarray

    def lattice_vector(self, m: int, n: int) -> np.ndarray:
        return int(m) * self.a1 + int(n) * self.a2

    def reciprocal_vector(self, h: int, k: int) -> np.ndarray:
        return int(h) * self.reciprocal_b1 + int(k) * self.reciprocal_b2

    def abc_shift(self, layer: int) -> np.ndarray:
        """Return A/B/C shift in a fixed primitive-cell gauge."""

        return int(layer) % 3 * self.tau

    def plane_basis_in_stacked_cubic_axes(self) -> np.ndarray:
        """Plane-frame axes expressed in the cubic axes of THIS +ABC stack.

        The stored e1,e2,e3 are the geometric/laboratory construction frame.
        With +tau=(a1+a2)/3 and +h, that frame is NOT the conventional cubic
        orientation often attached to the -tau FCC representative. A proper
        pi rotation about [111] supplies the cubic axes of the generated
        crystal: rows (-e1,-e2,e3). In these axes every neighbor is a_lat/2
        times an integer triplet of even sum. No atom/energy/path is changed.

        Use these rows for cubic elastic tensors and labeled cubic q paths.
        RegistryPath.direction_3d remains in the original geometric frame.
        """
        return np.stack((-self.e1, -self.e2, self.e3))


def fcc111_geometry(lattice_constant: float) -> FCC111Geometry:
    """Build the exact FCC(111) basis for a positive lattice constant."""

    alat = float(lattice_constant)
    if not np.isfinite(alat) or alat <= 0.0:
        raise ValueError("lattice_constant must be finite and positive")

    e1 = np.array([1.0, -1.0, 0.0]) / math.sqrt(2.0)
    e2 = np.array([1.0, 1.0, -2.0]) / math.sqrt(6.0)
    e3 = np.array([1.0, 1.0, 1.0]) / math.sqrt(3.0)
    b = alat / math.sqrt(2.0)
    h111 = alat / math.sqrt(3.0)

    # Components in the orthonormal (e1,e2) plane basis.
    a1 = np.array([b, 0.0])
    a2 = np.array([0.5 * b, 0.5 * math.sqrt(3.0) * b])
    area = 0.5 * math.sqrt(3.0) * b * b
    direct = np.column_stack((a1, a2))
    reciprocal = 2.0 * math.pi * np.linalg.inv(direct).T
    reciprocal_b1 = reciprocal[:, 0]
    reciprocal_b2 = reciprocal[:, 1]
    tau = (a1 + a2) / 3.0
    return FCC111Geometry(
        lattice_constant=alat,
        e1=e1,
        e2=e2,
        e3=e3,
        b=b,
        h111=h111,
        a1=a1,
        a2=a2,
        reciprocal_b1=reciprocal_b1,
        reciprocal_b2=reciprocal_b2,
        atomic_cell_area=area,
        tau=tau,
    )


def fcc111_geometry_from_b(b: float) -> FCC111Geometry:
    """Build the same geometry when the in-plane nearest spacing is primary."""

    return fcc111_geometry(float(b) * math.sqrt(2.0))


def registry_path(path_id: str) -> RegistryPath:
    """Return a named path without inferring an unspecified crystal orientation."""

    if path_id == DIRECT_110:
        return RegistryPath(
            path_id=DIRECT_110,
            crystallographic_label="(111) direct <110> full-Burgers path",
            direction_2d=(1.0, 0.0),
            direction_3d=(1.0 / math.sqrt(2.0), -1.0 / math.sqrt(2.0), 0.0),
            period_over_b=1.0,
            partial_increment_over_b=None,
            calibration_status=(
                "matches the repository's existing direct <110> scalar target"
            ),
        )
    if path_id == SHOCKLEY_112:
        # Unit vector parallel to a1+a2, equivalent to [2,-1,-1]/sqrt(6).
        return RegistryPath(
            path_id=SHOCKLEY_112,
            crystallographic_label="(111) <112> Shockley/minimum-GSF line",
            direction_2d=(math.sqrt(3.0) / 2.0, 0.5),
            direction_3d=(2.0 / math.sqrt(6.0), -1.0 / math.sqrt(6.0), -1.0 / math.sqrt(6.0)),
            period_over_b=math.sqrt(3.0),
            partial_increment_over_b=1.0 / math.sqrt(3.0),
            calibration_status=(
                "explicit held-out path; not the current scalar calibration path"
            ),
        )
    raise ValueError(f"unknown FCC(111) registry path: {path_id}")


def _nearest_lattice_representative(
    vector: np.ndarray, geometry: FCC111Geometry
) -> np.ndarray:
    """Choose a short lattice-equivalent representative of a 2D vector."""

    direct = np.column_stack((geometry.a1, geometry.a2))
    fractional = np.linalg.solve(direct, np.asarray(vector, dtype=float))
    center = np.rint(fractional).astype(int)
    candidates = []
    for dm in (-1, 0, 1):
        for dn in (-1, 0, 1):
            shift = direct @ (center + np.array([dm, dn]))
            candidate = np.asarray(vector, dtype=float) - shift
            candidates.append(candidate)
    return min(candidates, key=lambda value: float(value @ value)).copy()


def homogeneous_layer_registry(
    layer: int,
    s: float,
    *,
    geometry: FCC111Geometry,
    path_id: str = DIRECT_110,
    reduce_to_cell: bool = True,
) -> np.ndarray:
    r"""Return ``Delta_l = l(tau+s*e_s)`` modulo the triangular lattice.

    This is a homogeneous shear/registry increment: every adjacent (111)
    plane has the same additional displacement ``s*e_s``.  It is not a local
    single-interface fault construction.
    """

    path = registry_path(path_id)
    vector = int(layer) * (geometry.tau + float(s) * path.direction())
    if reduce_to_cell:
        return _nearest_lattice_representative(vector, geometry)
    return vector


def in_plane_lattice_norm_squared(m: int, n: int, b: float) -> float:
    """Exact ``|m*a1+n*a2|^2`` for the triangular primitive lattice."""

    return float(b) ** 2 * (int(m) ** 2 + int(m) * int(n) + int(n) ** 2)


def full_atomic_distance_squared(
    m: int,
    n: int,
    layer: int,
    a: float,
    s: float,
    *,
    geometry: FCC111Geometry,
    path_id: str = DIRECT_110,
) -> float:
    """Return the full reference-atom distance squared for ``(m,n,layer)``."""

    if int(m) == 0 and int(n) == 0 and int(layer) == 0:
        raise ValueError("the (0,0,0) self term is excluded")
    delta = homogeneous_layer_registry(
        int(layer), float(s), geometry=geometry, path_id=path_id
    )
    lateral = geometry.lattice_vector(int(m), int(n)) + delta
    return float((int(layer) * float(a)) ** 2 + lateral @ lateral)
