"""Analytic LJ--EAM hybrid for the uniform infinite two-row N=1 model.

The verified Poisson/Bessel LJ pair surface remains the base interaction.
This module adds an EAM-style environment density and differentiable
embedding energy; it is not a calibrated aluminum EAM parameterization.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import numpy as np

from .exponential_density import (
    ExponentialDensityParams,
    ExponentialDensityResult,
    exponential_cross_density,
    same_row_exponential_density,
)
from .model import ModelParams, TwoRowLJ


@runtime_checkable
class AnalyticEmbedding(Protocol):
    """Scalar differentiable embedding family."""

    def value(self, density): ...

    def first_derivative(self, density): ...

    def second_derivative(self, density): ...


@dataclass(frozen=True)
class SquareRootEmbedding:
    """Minimal analytic ``F(rho)=-A*sqrt(rho/rho_ref)`` family."""

    amplitude: float = 0.08
    rho_ref: float = 1.0

    def validate(self) -> None:
        if not np.isfinite(self.amplitude) or self.amplitude < 0.0:
            raise ValueError("embedding amplitude must be finite and nonnegative")
        if not np.isfinite(self.rho_ref) or self.rho_ref <= 0.0:
            raise ValueError("rho_ref must be finite and positive")

    def value(self, density):
        self.validate()
        rho = np.asarray(density, dtype=float)
        if np.any(rho < 0.0):
            raise ValueError("environment density cannot be negative")
        return -self.amplitude * np.sqrt(rho / self.rho_ref)

    def first_derivative(self, density):
        self.validate()
        rho = np.asarray(density, dtype=float)
        if np.any(rho <= 0.0):
            raise ValueError("square-root embedding derivative requires rho > 0")
        return -self.amplitude / (
            2.0 * np.sqrt(self.rho_ref) * np.sqrt(rho)
        )

    def second_derivative(self, density):
        self.validate()
        rho = np.asarray(density, dtype=float)
        if np.any(rho <= 0.0):
            raise ValueError("square-root embedding derivative requires rho > 0")
        return self.amplitude / (
            4.0 * np.sqrt(self.rho_ref) * rho ** 1.5
        )


@dataclass(frozen=True)
class SquareRootLinearEmbedding:
    r"""Minimal extension ``F(x)=-A*sqrt(x)+B*(x-1)``.

    Here ``x=rho/rho_ref`` and ``B>=0``.  In an unrestricted EAM decomposition
    this linear term is gauge-equivalent to a pair redistribution.  It is
    separately identifiable here only under the declared fixed-LJ pair gauge.
    This is an identifiability experiment, not an Al EAM claim.
    """

    amplitude: float
    linear: float
    rho_ref: float

    def validate(self) -> None:
        if not np.isfinite(self.amplitude) or self.amplitude < 0.0:
            raise ValueError("embedding amplitude must be finite and nonnegative")
        if not np.isfinite(self.linear) or self.linear < 0.0:
            raise ValueError("linear coefficient must be finite and nonnegative")
        if not np.isfinite(self.rho_ref) or self.rho_ref <= 0.0:
            raise ValueError("rho_ref must be finite and positive")

    def value(self, density):
        self.validate()
        rho = np.asarray(density, dtype=float)
        if np.any(rho < 0.0):
            raise ValueError("environment density cannot be negative")
        x = rho / self.rho_ref
        return -self.amplitude * np.sqrt(x) + self.linear * (x - 1.0)

    def first_derivative(self, density):
        self.validate()
        rho = np.asarray(density, dtype=float)
        if np.any(rho <= 0.0):
            raise ValueError("embedding derivative requires rho > 0")
        x = rho / self.rho_ref
        return (
            -self.amplitude / (2.0 * self.rho_ref * np.sqrt(x))
            + self.linear / self.rho_ref
        )

    def second_derivative(self, density):
        self.validate()
        rho = np.asarray(density, dtype=float)
        if np.any(rho <= 0.0):
            raise ValueError("embedding derivative requires rho > 0")
        x = rho / self.rho_ref
        return self.amplitude / (4.0 * self.rho_ref**2 * x ** 1.5)


@dataclass(frozen=True)
class HybridLocalTerms:
    density: np.ndarray
    density_cross: ExponentialDensityResult
    embedding_energy: np.ndarray
    embedding_d_da: np.ndarray
    embedding_d_ds: np.ndarray
    embedding_d2_daa: np.ndarray
    embedding_d2_das: np.ndarray
    embedding_d2_dss: np.ndarray


class AnalyticLJEAM(TwoRowLJ):
    """Separate analytic LJ pair plus EAM-style embedding energy surface.

    The reduced uniform two-row cell contains one upper and one lower atom.
    Their environment densities are identical by the even kernel and the
    integer substitution ``n -> -n-1``.  The cell embedding energy is thus
    exactly ``2*F(rho_parallel+rho_cross)``.  The variable LJ cross-row energy
    is the existing ``TwoRowLJ.local_energy`` and is not replaced or refitted.

    This first analytic hybrid is deliberately restricted to ``n_cells=1``;
    multi-cell environmental-density bookkeeping requires a separate
    derivation and is not silently approximated here.
    """

    def __init__(
        self,
        p: ModelParams,
        *,
        density_params: ExponentialDensityParams = ExponentialDensityParams(),
        embedding: AnalyticEmbedding = SquareRootEmbedding(),
    ):
        if p.n_cells != 1:
            raise ValueError("AnalyticLJEAM currently requires ModelParams(n_cells=1)")
        density_params.validate()
        if isinstance(embedding, (SquareRootEmbedding, SquareRootLinearEmbedding)):
            embedding.validate()
        if not isinstance(embedding, AnalyticEmbedding):
            raise TypeError("embedding must provide value and first/second derivatives")
        self.density_params = density_params
        self.embedding = embedding
        super().__init__(p)

    @property
    def embedding_is_zero(self) -> bool:
        if self.density_params.rho0 == 0.0:
            return True
        if isinstance(self.embedding, SquareRootLinearEmbedding):
            return self.embedding.amplitude == 0.0 and self.embedding.linear == 0.0
        amplitude = getattr(self.embedding, "amplitude", None)
        return amplitude == 0.0

    def local_environment_density(self, a, s) -> ExponentialDensityResult:
        cross = exponential_cross_density(
            a, s, b=self.p.b, params=self.density_params
        )
        parallel = same_row_exponential_density(
            b=self.p.b, params=self.density_params
        )
        return ExponentialDensityResult(
            value=cross.value + parallel,
            d_da=cross.d_da,
            d_ds=cross.d_ds,
            d2_daa=cross.d2_daa,
            d2_das=cross.d2_das,
            d2_dss=cross.d2_dss,
            modes_used=cross.modes_used,
            estimated_tail_absolute=cross.estimated_tail_absolute,
        )

    def embedding_terms(self, a, s) -> HybridLocalTerms:
        density = self.local_environment_density(a, s)
        if self.embedding_is_zero:
            zero = np.zeros_like(density.value)
            return HybridLocalTerms(
                density=density.value,
                density_cross=density,
                embedding_energy=zero,
                embedding_d_da=zero.copy(),
                embedding_d_ds=zero.copy(),
                embedding_d2_daa=zero.copy(),
                embedding_d2_das=zero.copy(),
                embedding_d2_dss=zero.copy(),
            )

        first = np.asarray(self.embedding.first_derivative(density.value))
        second = np.asarray(self.embedding.second_derivative(density.value))
        multiplicity = 2.0
        return HybridLocalTerms(
            density=density.value,
            density_cross=density,
            embedding_energy=multiplicity * np.asarray(
                self.embedding.value(density.value)
            ),
            embedding_d_da=multiplicity * first * density.d_da,
            embedding_d_ds=multiplicity * first * density.d_ds,
            embedding_d2_daa=multiplicity
            * (second * density.d_da**2 + first * density.d2_daa),
            embedding_d2_das=multiplicity
            * (second * density.d_da * density.d_ds + first * density.d2_das),
            embedding_d2_dss=multiplicity
            * (second * density.d_ds**2 + first * density.d2_dss),
        )

    def local_energy_gradient_array(self, a, s):
        pair_energy, pair_a, pair_s = super().local_energy_gradient_array(a, s)
        terms = self.embedding_terms(a, s)
        return (
            np.asarray(pair_energy) + terms.embedding_energy,
            np.asarray(pair_a) + terms.embedding_d_da,
            np.asarray(pair_s) + terms.embedding_d_ds,
        )

    def local_energy(self, a: float, s: float) -> float:
        energy, _, _ = self.local_energy_gradient_array(float(a), float(s))
        return float(energy)

    def local_energy_referenced(self, a, s):
        energy, _, _ = self.local_energy_gradient_array(a, s)
        reference, _, _ = self.local_energy_gradient_array(self.a0, 0.0)
        return np.asarray(energy) - float(reference)

    def local_deda(self, a: float, s: float) -> float:
        _, derivative, _ = self.local_energy_gradient_array(float(a), float(s))
        return float(derivative)

    def local_deds(self, a: float, s: float) -> float:
        _, _, derivative = self.local_energy_gradient_array(float(a), float(s))
        return float(derivative)

    def local_deda_array(self, a, s: float):
        aa = np.asarray(a, dtype=float)
        ss = np.full_like(aa, float(s), dtype=float)
        _, derivative, _ = self.local_energy_gradient_array(aa, ss)
        return np.asarray(derivative)

    def local_hessian(self, a: float, s: float) -> np.ndarray:
        pair = super().local_hessian(float(a), float(s))
        terms = self.embedding_terms(float(a), float(s))
        return pair + np.array(
            [
                [float(terms.embedding_d2_daa), float(terms.embedding_d2_das)],
                [float(terms.embedding_d2_das), float(terms.embedding_d2_dss)],
            ]
        )

    def normal_tangent_stiffness(self) -> float:
        stiffness = float(self.local_hessian(self.a0, 0.0)[0, 0])
        if not np.isfinite(stiffness) or stiffness <= 0.0:
            raise FloatingPointError("hybrid pristine normal stiffness must be positive")
        return stiffness

    def energy_gradient(self, a: np.ndarray, s: np.ndarray, force: float):
        energy, grad_a, grad_s = super().energy_gradient(a, s, force)
        terms = self.embedding_terms(a, s)
        return (
            float(energy + np.sum(terms.embedding_energy)),
            np.asarray(grad_a) + terms.embedding_d_da,
            np.asarray(grad_s) + terms.embedding_d_ds,
        )

    def energy_gradient_batch(self, a: np.ndarray, s: np.ndarray, force: float):
        energy, grad_a, grad_s = super().energy_gradient_batch(a, s, force)
        terms = self.embedding_terms(a, s)
        return (
            np.asarray(energy) + np.sum(terms.embedding_energy, axis=1),
            np.asarray(grad_a) + terms.embedding_d_da,
            np.asarray(grad_s) + terms.embedding_d_ds,
        )
