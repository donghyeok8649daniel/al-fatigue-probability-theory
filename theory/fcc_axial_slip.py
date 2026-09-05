"""FCC slip-system resolution for a single crystal under uniaxial stress.

The applied continuum stress remains purely axial.  Resolved shear is the
crystallographic driving stress on the twelve {111}<110> slip systems; it is
not an independently applied shear load and is not itself a fatigue law.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product

import numpy as np

from .cubic_normal_orientation import miller_unit_vector


@dataclass(frozen=True)
class FCCSlipSystem:
    plane_hkl: tuple[int, int, int]
    direction_uvw: tuple[int, int, int]
    signed_schmid_factor: float

    @property
    def schmid_factor(self) -> float:
        return abs(self.signed_schmid_factor)

    def resolved_shear_mpa(self, axial_stress_mpa: float) -> float:
        return self.signed_schmid_factor * axial_stress_mpa


def _canonical_sign(vector: tuple[int, int, int]) -> tuple[int, int, int]:
    for value in vector:
        if value < 0:
            return tuple(-item for item in vector)
        if value > 0:
            return vector
    raise ValueError("zero vector has no canonical sign")


def fcc_slip_system_indices() -> tuple[
    tuple[tuple[int, int, int], tuple[int, int, int]], ...
]:
    """Return the twelve unique undirected FCC {111}<110> systems."""
    planes = {
        _canonical_sign(tuple(signs))
        for signs in product((-1, 1), repeat=3)
    }
    directions: set[tuple[int, int, int]] = set()
    for zero_index in range(3):
        nonzero = [index for index in range(3) if index != zero_index]
        for signs in product((-1, 1), repeat=2):
            vector = [0, 0, 0]
            vector[nonzero[0]], vector[nonzero[1]] = signs
            directions.add(_canonical_sign(tuple(vector)))

    systems = []
    for plane in sorted(planes):
        for direction in sorted(directions):
            if sum(a * b for a, b in zip(plane, direction)) == 0:
                systems.append((plane, direction))
    if len(systems) != 12:
        raise RuntimeError(f"expected 12 FCC slip systems, generated {len(systems)}")
    return tuple(systems)


def fcc_axial_slip_systems(h: int, k: int, l: int) -> tuple[FCCSlipSystem, ...]:
    """Resolve uniaxial [h k l] stress onto all FCC slip systems."""
    loading = miller_unit_vector(h, k, l)
    resolved = []
    for plane, direction in fcc_slip_system_indices():
        normal = np.asarray(plane, dtype=float) / np.sqrt(3.0)
        slip = np.asarray(direction, dtype=float) / np.sqrt(2.0)
        signed = float(np.dot(loading, normal) * np.dot(loading, slip))
        resolved.append(FCCSlipSystem(plane, direction, signed))
    return tuple(sorted(resolved, key=lambda item: item.schmid_factor, reverse=True))


def maximum_schmid_factor(h: int, k: int, l: int) -> float:
    return fcc_axial_slip_systems(h, k, l)[0].schmid_factor
