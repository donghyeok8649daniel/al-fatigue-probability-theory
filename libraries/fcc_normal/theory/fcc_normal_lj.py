"""FCC generalized-Lennard-Jones normal-deformation lattice model.

Active normal-only model extension.

The model evaluates the homogeneous FCC lattice energy per atom

    U(F) = 1/2 sum_{R != 0} v(|F R|)

for a central generalized-LJ pair potential. For [001] normal loading the
homogeneous deformation gradient is diag(lambda_t, lambda_t, lambda_n), and
lambda_t is relaxed by minimizing U at fixed lambda_n. This enforces zero
homogeneous transverse nominal stress within the stated pair-potential model.

No damping, damage variable, slip coordinate, or empirical fatigue evolution
law is introduced here.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable

import numpy as np

EV_J = 1.602176634e-19


@dataclass(frozen=True)
class FCCNormalLJParameters:
    lattice_constant_m: float = 4.0495e-10
    repulsive_exponent: float = 12.19
    attractive_exponent: float = 6.0
    cutoff_lattice_constants: float = 12.0

    @property
    def atomic_volume_m3(self) -> float:
        return self.lattice_constant_m**3 / 4.0


def _validate_exponents(m: float, n: float) -> None:
    if not (m > n > 3.0):
        raise ValueError("for a convergent 3D infinite lattice sum require m > n > 3")


def generate_fcc_lattice_vectors(cutoff_lattice_constants: float) -> np.ndarray:
    """Return nonzero FCC Bravais vectors R/a_lat inside a spherical cutoff.

    FCC vectors are represented as (i,j,k)/2 with i+j+k even.
    """
    if cutoff_lattice_constants <= 0.0:
        raise ValueError("cutoff must be positive")
    lim = int(math.ceil(2.0 * cutoff_lattice_constants))
    vectors: list[tuple[float, float, float]] = []
    for i in range(-lim, lim + 1):
        for j in range(-lim, lim + 1):
            for k in range(-lim, lim + 1):
                if i == 0 and j == 0 and k == 0:
                    continue
                if (i + j + k) % 2 != 0:
                    continue
                x, y, z = 0.5 * i, 0.5 * j, 0.5 * k
                if x*x + y*y + z*z <= cutoff_lattice_constants**2:
                    vectors.append((x, y, z))
    return np.asarray(vectors, dtype=float)


def golden_section_minimize(
    function,
    lower: float,
    upper: float,
    *,
    tolerance: float = 1.0e-10,
    max_iterations: int = 200,
) -> float:
    """Dependency-free scalar minimization used for transverse relaxation."""
    if not lower < upper:
        raise ValueError("lower must be smaller than upper")
    ratio = (math.sqrt(5.0) - 1.0) / 2.0
    x1 = upper - ratio * (upper - lower)
    x2 = lower + ratio * (upper - lower)
    f1 = function(x1)
    f2 = function(x2)
    for _ in range(max_iterations):
        if abs(upper - lower) <= tolerance:
            break
        if f1 > f2:
            lower, x1, f1 = x1, x2, f2
            x2 = lower + ratio * (upper - lower)
            f2 = function(x2)
        else:
            upper, x2, f2 = x2, x1, f1
            x1 = upper - ratio * (upper - lower)
            f1 = function(x1)
    return 0.5 * (lower + upper)


class FCCNormalLJ:
    """Finite-cutoff approximation to the infinite FCC central-pair lattice."""

    def __init__(self, parameters: FCCNormalLJParameters = FCCNormalLJParameters()):
        self.parameters = parameters
        m = parameters.repulsive_exponent
        n = parameters.attractive_exponent
        _validate_exponents(m, n)
        self.vectors = generate_fcc_lattice_vectors(
            parameters.cutoff_lattice_constants
        )
        self._x, self._y, self._z = self.vectors.T
        radius = np.linalg.norm(self.vectors, axis=1)
        self._A_m = 0.5 * float(np.sum(radius ** (-m)))
        self._A_n = 0.5 * float(np.sum(radius ** (-n)))
        self.sigma_over_lattice_constant = (
            n * self._A_n / (m * self._A_m)
        ) ** (1.0 / (m - n))
        q = self.sigma_over_lattice_constant
        self._reference_energy_per_epsilon = (
            self._A_m * q**m - self._A_n * q**n
        )

    @property
    def m(self) -> float:
        return self.parameters.repulsive_exponent

    @property
    def n(self) -> float:
        return self.parameters.attractive_exponent

    @property
    def atomic_volume_m3(self) -> float:
        return self.parameters.atomic_volume_m3

    @property
    def sigma_lj_m(self) -> float:
        return self.sigma_over_lattice_constant * self.parameters.lattice_constant_m

    def cohesive_energy_j_per_atom(self, epsilon_lj_j: float) -> float:
        """Positive separation energy corresponding to U(infinity)-U(reference)."""
        return -epsilon_lj_j * self._reference_energy_per_epsilon

    def epsilon_for_cohesive_energy(self, cohesive_energy_j_per_atom: float) -> float:
        if cohesive_energy_j_per_atom <= 0.0:
            raise ValueError("cohesive energy must be positive")
        return -cohesive_energy_j_per_atom / self._reference_energy_per_epsilon

    def energy_per_atom_j(
        self,
        epsilon_lj_j: float,
        transverse_stretch: float,
        normal_stretch: float,
    ) -> float:
        if min(epsilon_lj_j, transverse_stretch, normal_stretch) <= 0.0:
            raise ValueError("energy scale and stretches must be positive")
        q = self.sigma_over_lattice_constant
        distance = np.sqrt(
            (transverse_stretch * self._x)**2
            + (transverse_stretch * self._y)**2
            + (normal_stretch * self._z)**2
        )
        return 0.5 * epsilon_lj_j * float(
            np.sum((q / distance)**self.m - (q / distance)**self.n)
        )

    def relaxed_transverse_stretch(
        self,
        epsilon_lj_j: float,
        normal_stretch: float,
        *,
        lower: float = 0.65,
        upper: float = 1.10,
    ) -> float:
        return golden_section_minimize(
            lambda transverse: self.energy_per_atom_j(
                epsilon_lj_j,
                transverse,
                normal_stretch,
            ),
            lower,
            upper,
        )

    def axial_engineering_stress_pa(
        self,
        epsilon_lj_j: float,
        normal_stretch: float,
        transverse_stretch: float | None = None,
        *,
        difference_step: float = 1.0e-5,
    ) -> float:
        if transverse_stretch is None:
            transverse_stretch = self.relaxed_transverse_stretch(
                epsilon_lj_j, normal_stretch
            )
        derivative = (
            self.energy_per_atom_j(
                epsilon_lj_j,
                transverse_stretch,
                normal_stretch + difference_step,
            )
            - self.energy_per_atom_j(
                epsilon_lj_j,
                transverse_stretch,
                normal_stretch - difference_step,
            )
        ) / (2.0 * difference_step)
        return derivative / self.atomic_volume_m3

    def small_strain_properties(
        self,
        epsilon_lj_j: float,
        *,
        strain_step: float = 5.0e-4,
    ) -> tuple[float, float]:
        plus = 1.0 + strain_step
        minus = 1.0 - strain_step
        transverse_plus = self.relaxed_transverse_stretch(epsilon_lj_j, plus)
        transverse_minus = self.relaxed_transverse_stretch(epsilon_lj_j, minus)
        stress_plus = self.axial_engineering_stress_pa(
            epsilon_lj_j, plus, transverse_plus
        )
        stress_minus = self.axial_engineering_stress_pa(
            epsilon_lj_j, minus, transverse_minus
        )
        young = (stress_plus - stress_minus) / (2.0 * strain_step)
        poisson = -(
            transverse_plus - transverse_minus
        ) / (2.0 * strain_step)
        return young, poisson

    def epsilon_for_youngs_modulus(
        self,
        target_youngs_modulus_pa: float,
        *,
        reference_epsilon_j: float = EV_J,
    ) -> float:
        if target_youngs_modulus_pa <= 0.0:
            raise ValueError("target modulus must be positive")
        reference_young, _ = self.small_strain_properties(reference_epsilon_j)
        return reference_epsilon_j * target_youngs_modulus_pa / reference_young

    def _energy_under_F(self, epsilon_lj_j: float, F: np.ndarray) -> float:
        deformed = self.vectors @ np.asarray(F, dtype=float).T
        distance = np.linalg.norm(deformed, axis=1)
        q = self.sigma_over_lattice_constant
        return 0.5 * epsilon_lj_j * float(
            np.sum((q / distance)**self.m - (q / distance)**self.n)
        )

    def cubic_elastic_constants_pa(
        self,
        epsilon_lj_j: float,
        *,
        strain_step: float = 2.0e-4,
    ) -> tuple[float, float, float]:
        """Return (C11,C12,C44) by finite differences at the zero-pressure state."""
        h = strain_step
        I = np.eye(3)
        U0 = self._energy_under_F(epsilon_lj_j, I)

        Fp = I.copy(); Fp[0, 0] += h
        Fm = I.copy(); Fm[0, 0] -= h
        C11 = (
            self._energy_under_F(epsilon_lj_j, Fp)
            - 2.0 * U0
            + self._energy_under_F(epsilon_lj_j, Fm)
        ) / h**2 / self.atomic_volume_m3

        def biaxial(ex: float, ey: float) -> float:
            F = I.copy()
            F[0, 0] += ex
            F[1, 1] += ey
            return self._energy_under_F(epsilon_lj_j, F)

        C12 = (
            biaxial(h, h)
            - biaxial(h, -h)
            - biaxial(-h, h)
            + biaxial(-h, -h)
        ) / (4.0 * h**2) / self.atomic_volume_m3

        Fp = I.copy(); Fp[0, 1] = h
        Fm = I.copy(); Fm[0, 1] = -h
        C44 = (
            self._energy_under_F(epsilon_lj_j, Fp)
            - 2.0 * U0
            + self._energy_under_F(epsilon_lj_j, Fm)
        ) / h**2 / self.atomic_volume_m3
        return C11, C12, C44

    def stress_strain_curve(
        self,
        epsilon_lj_j: float,
        normal_stretches: Iterable[float],
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        stretches = np.asarray(list(normal_stretches), dtype=float)
        transverse = np.empty_like(stretches)
        stress = np.empty_like(stretches)
        for i, normal in enumerate(stretches):
            transverse[i] = self.relaxed_transverse_stretch(epsilon_lj_j, normal)
            stress[i] = self.axial_engineering_stress_pa(
                epsilon_lj_j, normal, transverse[i]
            )
        return stretches, transverse, stress
