"""Scalar normal-tension projection for a cubic single crystal.

This module supplies only the Young modulus along a declared crystallographic
loading direction.  It does not introduce shear fatigue, slip, a tensor damage
criterion, or the archived FCC pair-sum model.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np


@dataclass(frozen=True)
class CubicElasticConstants:
    c11_pa: float
    c12_pa: float
    c44_pa: float

    def validate(self) -> None:
        values = (self.c11_pa, self.c12_pa, self.c44_pa)
        if not all(math.isfinite(value) for value in values):
            raise ValueError("cubic elastic constants must be finite")
        # Born stability conditions for a cubic elastic solid.
        if not (self.c44_pa > 0 and self.c11_pa - self.c12_pa > 0
                and self.c11_pa + 2*self.c12_pa > 0):
            raise ValueError("cubic elastic constants violate mechanical stability")

    @property
    def compliances(self) -> tuple[float, float, float]:
        """Return S11, S12 and S44 in Pa^-1."""
        self.validate()
        denominator = (self.c11_pa-self.c12_pa)*(self.c11_pa+2*self.c12_pa)
        s11 = (self.c11_pa+self.c12_pa)/denominator
        s12 = -self.c12_pa/denominator
        s44 = 1/self.c44_pa
        return s11, s12, s44


def miller_unit_vector(h: int, k: int, l: int) -> np.ndarray:
    """Unit loading-axis vector parallel to the cubic direction [h k l]."""
    values = np.asarray([h, k, l], dtype=float)
    norm = float(np.linalg.norm(values))
    if norm == 0 or not np.all(np.isfinite(values)):
        raise ValueError("Miller loading direction cannot be [0 0 0]")
    return values/norm


def directional_young_modulus(
    constants: CubicElasticConstants,
    h: int,
    k: int,
    l: int,
) -> float:
    """Young modulus for uniaxial normal loading along cubic [h k l]."""
    direction = miller_unit_vector(h, k, l)
    l1, l2, l3 = direction
    s11, s12, s44 = constants.compliances
    orientation = l1*l1*l2*l2 + l2*l2*l3*l3 + l3*l3*l1*l1
    inverse_e = s11 - 2*(s11-s12-0.5*s44)*orientation
    if inverse_e <= 0 or not math.isfinite(inverse_e):
        raise ValueError("directional compliance is nonpositive")
    return 1/inverse_e
