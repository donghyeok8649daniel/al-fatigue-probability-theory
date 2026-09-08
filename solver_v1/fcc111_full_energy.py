"""Full homogeneous FCC(111) infinite-plane-stack static reference energy.

The model is deliberately separate from ``TwoRowLJ`` and ``AnalyticLJEAM``.
It evaluates a per-atom infinite FCC plane stack.  The old reduced row energy
is neither replaced nor added to this energy.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import math

import numpy as np
from scipy.optimize import root_scalar
from scipy.special import zeta

from .analytic_lj_eam import AnalyticEmbedding, AnalyticLJEAM
from .exponential_density import ExponentialDensityParams
from .fcc111_geometry import (
    DIRECT_110,
    FCC111Geometry,
    fcc111_geometry_from_b,
    homogeneous_layer_registry,
    registry_path,
)
from .fcc111_lattice_sum import (
    PlaneKernelResult,
    ReciprocalSumConfig,
    plane_exponential_sum_direct,
    plane_exponential_sum_reciprocal,
    plane_power_sum_direct,
    plane_power_sum_reciprocal,
    same_plane_exponential_density_direct,
    same_plane_power_sum_direct,
    triangular_epstein_zeta,
)
from .model import ModelParams, TwoRowLJ


@dataclass(frozen=True)
class StackSumConfig:
    tol: float = 2.0e-12
    max_layers: int = 96
    consecutive_small_layers: int = 3
    reciprocal: ReciprocalSumConfig = ReciprocalSumConfig(
        tol=2.0e-13, max_shell_index=32, consecutive_small_shells=3
    )
    same_plane_density_radius: int = 32
    max_same_plane_density_radius: int = 256

    def validate(self) -> None:
        if not np.isfinite(self.tol) or self.tol <= 0.0:
            raise ValueError("stack tol must be finite and positive")
        if self.max_layers < 4 or self.consecutive_small_layers < 1:
            raise ValueError("invalid layer convergence controls")
        if self.same_plane_density_radius < 4:
            raise ValueError("same-plane density radius must be at least 4")
        if self.max_same_plane_density_radius < self.same_plane_density_radius:
            raise ValueError("maximum same-plane radius must not be smaller than initial radius")
        self.reciprocal.validate()


@dataclass(frozen=True)
class ScalarDerivatives:
    value: float
    d_da: float
    d_ds: float
    d2_daa: float
    d2_das: float
    d2_dss: float
    layers_used: int
    maximum_plane_shells: int
    estimated_tail_absolute: float


@dataclass(frozen=True)
class FullFCCEvaluation:
    pair: ScalarDerivatives
    density: ScalarDerivatives
    embedding_energy: float
    embedding_d_da: float
    embedding_d_ds: float
    embedding_d2_daa: float
    embedding_d2_das: float
    embedding_d2_dss: float
    energy: float
    d_da: float
    d_ds: float
    d2_daa: float
    d2_das: float
    d2_dss: float


@dataclass(frozen=True)
class DirectFullFCCEvaluation:
    pair_energy: float
    density: float
    embedding_energy: float
    energy: float
    d_da: float
    d_ds: float
    d2_daa: float
    d2_das: float
    d2_dss: float
    radial_index: int
    layers: int
    estimated_pair_tail: float
    estimated_density_tail: float


def _zero_mode_power_terms(
    d: float, p: float, area: float
) -> tuple[float, float, float]:
    alpha = 2.0 - 2.0 * float(p)
    value = math.pi / (float(area) * (float(p) - 1.0)) * float(d) ** alpha
    return value, alpha * value / d, alpha * (alpha - 1.0) * value / d**2


class FullFCC111StackEnergy:
    """Per-atom full FCC(111) infinite stack with optional EAM embedding.

    ``s`` is the additional homogeneous adjacent-plane displacement along an
    explicitly named in-plane path.  Plane ``l`` therefore has registry
    ``l*(tau+s*e_s)`` modulo the triangular lattice.  This is a homogeneous
    scalar shear path, not an isolated stacking-fault interface.
    """

    def __init__(
        self,
        p: ModelParams,
        *,
        path_id: str = DIRECT_110,
        density_params: ExponentialDensityParams | None = None,
        embedding: AnalyticEmbedding | None = None,
        stack_config: StackSumConfig = StackSumConfig(),
        require_stable_equilibrium: bool = True,
    ):
        if p.n_cells != 1:
            raise ValueError("full FCC homogeneous reference requires n_cells=1")
        if p.b <= 0.0 or p.epsilon <= 0.0 or p.sigma_lj <= 0.0:
            raise ValueError("b, epsilon, and sigma_lj must be positive")
        registry_path(path_id)
        stack_config.validate()
        if (density_params is None) != (embedding is None):
            raise ValueError("density_params and embedding must be supplied together")
        if density_params is not None:
            density_params.validate()
            if not isinstance(embedding, AnalyticEmbedding):
                raise TypeError("embedding must provide value and two derivatives")
        self.p = p
        self.path_id = path_id
        self.geometry = fcc111_geometry_from_b(p.b)
        self.density_params = density_params
        self.embedding = embedding
        self.stack_config = stack_config
        self.equilibrium_error: str | None = None
        try:
            self.a0 = self._find_reference_a()
        except ValueError as exc:
            if require_stable_equilibrium:
                raise
            self.a0 = float("nan")
            self.equilibrium_error = str(exc)

    @classmethod
    def from_reduced_model(
        cls,
        model: TwoRowLJ,
        *,
        path_id: str = DIRECT_110,
        stack_config: StackSumConfig = StackSumConfig(),
        require_stable_equilibrium: bool = True,
    ) -> "FullFCC111StackEnergy":
        """Evaluate exactly the same pair/density/embedding parameters in FCC."""

        if isinstance(model, AnalyticLJEAM):
            return cls(
                model.p,
                path_id=path_id,
                density_params=model.density_params,
                embedding=model.embedding,
                stack_config=stack_config,
                require_stable_equilibrium=require_stable_equilibrium,
            )
        return cls(
            model.p,
            path_id=path_id,
            stack_config=stack_config,
            require_stable_equilibrium=require_stable_equilibrium,
        )

    @property
    def has_embedding(self) -> bool:
        return self.density_params is not None and self.embedding is not None

    @property
    def embedding_is_zero(self) -> bool:
        if not self.has_embedding:
            return True
        assert self.density_params is not None and self.embedding is not None
        if self.density_params.rho0 == 0.0:
            return True
        amplitude = float(getattr(self.embedding, "amplitude", float("nan")))
        linear = float(getattr(self.embedding, "linear", 0.0))
        return amplitude == 0.0 and linear == 0.0

    def _layer_delta(self, layer: int, s: float) -> np.ndarray:
        return homogeneous_layer_registry(
            layer,
            s,
            geometry=self.geometry,
            path_id=self.path_id,
            reduce_to_cell=True,
        )

    def _transform_plane_result(
        self, plane: PlaneKernelResult, layer: int
    ) -> tuple[float, float, float, float, float]:
        direction = registry_path(self.path_id).direction()
        scale = float(layer)
        d_da = scale * plane.d_d
        d_ds = scale * float(direction @ plane.grad_delta)
        d2_daa = scale**2 * plane.d2_dd
        d2_das = scale**2 * float(direction @ plane.mixed_d_delta)
        d2_dss = scale**2 * float(direction @ plane.hess_delta @ direction)
        return d_da, d_ds, d2_daa, d2_das, d2_dss

    def full_power_sum(self, a: float, s: float, p: float) -> ScalarDerivatives:
        """Return the exact-mean/semi-analytic corrugation full power sum.

        The per-atom half-counting is ``0.5*same_plane + sum_(l>=1)``.
        Positive and negative layers are equal and cancel the original 1/2.
        """

        spacing = float(a)
        exponent = float(p)
        if spacing <= 0.0 or exponent <= 1.0:
            raise ValueError("require a>0 and p>1")
        same_plane = 0.5 * triangular_epstein_zeta(exponent, self.p.b)
        alpha = 2.0 - 2.0 * exponent
        mean = (
            math.pi
            / (self.geometry.atomic_cell_area * (exponent - 1.0))
            * spacing**alpha
            * float(zeta(2.0 * exponent - 2.0, 1.0))
        )
        value = same_plane + mean
        d_da = alpha * mean / spacing
        d2_daa = alpha * (alpha - 1.0) * mean / spacing**2
        d_ds = 0.0
        d2_das = 0.0
        d2_dss = 0.0
        small = 0
        layers_used = 0
        max_shells = 0
        recent: list[float] = []
        for layer in range(1, self.stack_config.max_layers + 1):
            distance = layer * spacing
            plane = plane_power_sum_reciprocal(
                distance,
                self._layer_delta(layer, s),
                p=exponent,
                geometry=self.geometry,
                config=self.stack_config.reciprocal,
            )
            zero, zero_d, zero_dd = _zero_mode_power_terms(
                distance, exponent, self.geometry.atomic_cell_area
            )
            transformed = self._transform_plane_result(plane, layer)
            terms = (
                plane.value - zero,
                transformed[0] - layer * zero_d,
                transformed[1],
                transformed[2] - layer**2 * zero_dd,
                transformed[3],
                transformed[4],
            )
            value += terms[0]
            d_da += terms[1]
            d_ds += terms[2]
            d2_daa += terms[3]
            d2_das += terms[4]
            d2_dss += terms[5]
            layers_used = layer
            max_shells = max(max_shells, plane.shells_used)
            bound = max(abs(term) for term in terms)
            recent.append(bound)
            recent = recent[-self.stack_config.consecutive_small_layers :]
            scale = max(
                1.0, abs(value), abs(d_da), abs(d_ds), abs(d2_daa),
                abs(d2_das), abs(d2_dss),
            )
            if bound <= self.stack_config.tol * scale:
                small += 1
                if small >= self.stack_config.consecutive_small_layers:
                    break
            else:
                small = 0
        else:
            raise RuntimeError("full FCC power corrugation did not meet layer tolerance")
        return ScalarDerivatives(
            value=float(value), d_da=float(d_da), d_ds=float(d_ds),
            d2_daa=float(d2_daa), d2_das=float(d2_das), d2_dss=float(d2_dss),
            layers_used=layers_used, maximum_plane_shells=max_shells,
            estimated_tail_absolute=float(max(recent, default=0.0)),
        )

    @lru_cache(maxsize=1)
    def _same_plane_density(self) -> tuple[float, int, float]:
        if self.density_params is None:
            return 0.0, 0, 0.0
        radius = self.stack_config.same_plane_density_radius
        while radius <= self.stack_config.max_same_plane_density_radius:
            result = same_plane_exponential_density_direct(
                kappa=self.density_params.kappa,
                amplitude=self.density_params.C_rho,
                geometry=self.geometry,
                radius_index=radius,
            )
            value, _, tail = result
            if tail <= self.stack_config.tol * max(1.0, abs(value)):
                return result
            radius *= 2
        raise RuntimeError("same-plane environment density did not meet tail tolerance")

    def full_environment_density(self, a: float, s: float) -> ScalarDerivatives:
        """Sum all neighbors first; no embedding function is applied here."""

        if self.density_params is None:
            return ScalarDerivatives(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0.0)
        spacing = float(a)
        if spacing <= 0.0:
            raise ValueError("a must be positive")
        value = self._same_plane_density()[0]
        d_da = d_ds = d2_daa = d2_das = d2_dss = 0.0
        small = 0
        layers_used = 0
        max_shells = 0
        recent: list[float] = []
        for layer in range(1, self.stack_config.max_layers + 1):
            plane = plane_exponential_sum_reciprocal(
                layer * spacing,
                self._layer_delta(layer, s),
                kappa=self.density_params.kappa,
                amplitude=self.density_params.C_rho,
                geometry=self.geometry,
                config=self.stack_config.reciprocal,
            )
            transformed = self._transform_plane_result(plane, layer)
            # Both +l and -l contribute to the environment density.
            terms = tuple(2.0 * term for term in (plane.value, *transformed))
            value += terms[0]
            d_da += terms[1]
            d_ds += terms[2]
            d2_daa += terms[3]
            d2_das += terms[4]
            d2_dss += terms[5]
            layers_used = layer
            max_shells = max(max_shells, plane.shells_used)
            bound = max(abs(term) for term in terms)
            recent.append(bound)
            recent = recent[-self.stack_config.consecutive_small_layers :]
            scale = max(
                1.0, abs(value), abs(d_da), abs(d_ds), abs(d2_daa),
                abs(d2_das), abs(d2_dss),
            )
            if bound <= self.stack_config.tol * scale:
                small += 1
                if small >= self.stack_config.consecutive_small_layers:
                    break
            else:
                small = 0
        else:
            raise RuntimeError("full FCC density did not meet layer tolerance")
        return ScalarDerivatives(
            value=float(value), d_da=float(d_da), d_ds=float(d_ds),
            d2_daa=float(d2_daa), d2_das=float(d2_das), d2_dss=float(d2_dss),
            layers_used=layers_used, maximum_plane_shells=max_shells,
            estimated_tail_absolute=float(max(recent, default=0.0)),
        )

    @lru_cache(maxsize=2048)
    def evaluate(self, a: float, s: float) -> FullFCCEvaluation:
        """Evaluate per-atom energy and analytic first/second derivatives."""

        spacing = float(a)
        registry = float(s)
        s6 = self.full_power_sum(spacing, registry, 6.0)
        s3 = self.full_power_sum(spacing, registry, 3.0)
        factor = 4.0 * self.p.epsilon

        def pair_component(name: str) -> float:
            return factor * (
                self.p.sigma_lj**12 * getattr(s6, name)
                - self.p.sigma_lj**6 * getattr(s3, name)
            )

        pair = ScalarDerivatives(
            value=pair_component("value"),
            d_da=pair_component("d_da"),
            d_ds=pair_component("d_ds"),
            d2_daa=pair_component("d2_daa"),
            d2_das=pair_component("d2_das"),
            d2_dss=pair_component("d2_dss"),
            layers_used=max(s6.layers_used, s3.layers_used),
            maximum_plane_shells=max(
                s6.maximum_plane_shells, s3.maximum_plane_shells
            ),
            estimated_tail_absolute=factor
            * (
                self.p.sigma_lj**12 * s6.estimated_tail_absolute
                + self.p.sigma_lj**6 * s3.estimated_tail_absolute
            ),
        )
        density = self.full_environment_density(spacing, registry)
        if self.embedding_is_zero:
            embedding_values = (0.0,) * 6
        else:
            assert self.embedding is not None
            first = float(self.embedding.first_derivative(density.value))
            second = float(self.embedding.second_derivative(density.value))
            embedding_values = (
                float(self.embedding.value(density.value)),
                first * density.d_da,
                first * density.d_ds,
                second * density.d_da**2 + first * density.d2_daa,
                second * density.d_da * density.d_ds + first * density.d2_das,
                second * density.d_ds**2 + first * density.d2_dss,
            )
        return FullFCCEvaluation(
            pair=pair,
            density=density,
            embedding_energy=embedding_values[0],
            embedding_d_da=embedding_values[1],
            embedding_d_ds=embedding_values[2],
            embedding_d2_daa=embedding_values[3],
            embedding_d2_das=embedding_values[4],
            embedding_d2_dss=embedding_values[5],
            energy=pair.value + embedding_values[0],
            d_da=pair.d_da + embedding_values[1],
            d_ds=pair.d_ds + embedding_values[2],
            d2_daa=pair.d2_daa + embedding_values[3],
            d2_das=pair.d2_das + embedding_values[4],
            d2_dss=pair.d2_dss + embedding_values[5],
        )

    def energy(self, a: float, s: float) -> float:
        return self.evaluate(float(a), float(s)).energy

    def local_energy(self, a: float, s: float) -> float:
        return self.energy(a, s)

    def grad(self, a: float, s: float) -> np.ndarray:
        result = self.evaluate(float(a), float(s))
        return np.array([result.d_da, result.d_ds])

    def local_deda(self, a: float, s: float) -> float:
        return float(self.evaluate(float(a), float(s)).d_da)

    def local_deds(self, a: float, s: float) -> float:
        return float(self.evaluate(float(a), float(s)).d_ds)

    def local_hessian(self, a: float, s: float) -> np.ndarray:
        result = self.evaluate(float(a), float(s))
        return np.array(
            [[result.d2_daa, result.d2_das], [result.d2_das, result.d2_dss]],
            dtype=float,
        )

    def hessian(self, a: float, s: float) -> np.ndarray:
        return self.local_hessian(a, s)

    def _find_reference_a(self) -> float:
        lower = max(0.15 * self.p.b, float(self.p.a_min))
        upper = min(float(self.p.a_max), 2.5 * self.p.b)
        grid = np.linspace(lower, upper, 96)
        derivative = np.array([self.local_deda(value, 0.0) for value in grid])
        roots = []
        for left, right, fleft, fright in zip(
            grid[:-1], grid[1:], derivative[:-1], derivative[1:]
        ):
            if not np.isfinite(fleft) or not np.isfinite(fright):
                continue
            if fleft <= 0.0 <= fright:
                solved = root_scalar(
                    lambda value: self.local_deda(float(value), 0.0),
                    bracket=(float(left), float(right)),
                    method="brentq",
                    xtol=8.0e-14,
                    rtol=8.0 * np.finfo(float).eps,
                )
                if solved.converged:
                    root = float(solved.root)
                    if self.local_hessian(root, 0.0)[0, 0] > 0.0:
                        roots.append(root)
        if not roots:
            raise ValueError("full FCC surface has no stable zero-load normal root")
        return min(roots, key=lambda value: self.local_energy(value, 0.0))

    def sigma_over_E_force_scale(self) -> float:
        hessian = self.local_hessian(self.a0, 0.0)
        eigenvalues = np.linalg.eigvalsh(hessian)
        if np.min(eigenvalues) <= 0.0:
            raise FloatingPointError("full FCC reference Hessian is not positive definite")
        c = np.array([1.0, self.p.chi_axial_projection])
        compliance = float(c @ np.linalg.solve(hessian, c))
        if compliance <= 0.0 or not np.isfinite(compliance):
            raise FloatingPointError("full FCC relaxed compliance is invalid")
        return float(self.a0 / compliance)

    def force_from_sigma_over_E(self, sigma_over_E):
        return self.sigma_over_E_force_scale() * np.asarray(sigma_over_E, dtype=float)

    def isolated_plane_energy(self) -> float:
        """Per-atom energy after homogeneous separation into isolated planes.

        This is a static uniform-plane-separation reference, not a finite
        single-cleavage-interface construction.
        """

        s6 = 0.5 * triangular_epstein_zeta(6.0, self.p.b)
        s3 = 0.5 * triangular_epstein_zeta(3.0, self.p.b)
        pair = 4.0 * self.p.epsilon * (
            self.p.sigma_lj**12 * s6 - self.p.sigma_lj**6 * s3
        )
        if self.embedding_is_zero:
            return float(pair)
        assert self.embedding is not None
        density = self._same_plane_density()[0]
        return float(pair + self.embedding.value(density))

    def direct_reference(
        self,
        a: float,
        s: float,
        *,
        radial_index: int = 80,
        layers: int = 80,
    ) -> DirectFullFCCEvaluation:
        """Independent truncated real-space full-stack reference evaluator."""

        spacing = float(a)
        registry = float(s)

        def direct_power(exponent: float) -> tuple[ScalarDerivatives, float]:
            same, _, same_tail = same_plane_power_sum_direct(
                p=exponent,
                geometry=self.geometry,
                radius_index=radial_index,
            )
            values = [0.5 * same, 0.0, 0.0, 0.0, 0.0, 0.0]
            radial_tail = 0.5 * same_tail
            for layer in range(1, int(layers) + 1):
                plane = plane_power_sum_direct(
                    layer * spacing,
                    self._layer_delta(layer, registry),
                    p=exponent,
                    geometry=self.geometry,
                    radius_index=radial_index,
                )
                transformed = self._transform_plane_result(plane, layer)
                terms = (plane.value, *transformed)
                values = [left + right for left, right in zip(values, terms)]
                radial_tail += plane.continuum_tail_estimate
            layer_tail = (
                math.pi
                / (self.geometry.atomic_cell_area * (exponent - 1.0))
                * spacing ** (2.0 - 2.0 * exponent)
                * layers ** (3.0 - 2.0 * exponent)
                / (2.0 * exponent - 3.0)
            )
            return ScalarDerivatives(
                value=values[0], d_da=values[1], d_ds=values[2],
                d2_daa=values[3], d2_das=values[4], d2_dss=values[5],
                layers_used=int(layers), maximum_plane_shells=0,
                estimated_tail_absolute=float(radial_tail + layer_tail),
            ), float(radial_tail + layer_tail)

        s6, tail6 = direct_power(6.0)
        s3, tail3 = direct_power(3.0)
        factor = 4.0 * self.p.epsilon
        pair_values = [
            factor
            * (self.p.sigma_lj**12 * getattr(s6, name) - self.p.sigma_lj**6 * getattr(s3, name))
            for name in ("value", "d_da", "d_ds", "d2_daa", "d2_das", "d2_dss")
        ]
        density_values = [0.0] * 6
        density_tail = 0.0
        if self.density_params is not None:
            same, _, same_tail = same_plane_exponential_density_direct(
                kappa=self.density_params.kappa,
                amplitude=self.density_params.C_rho,
                geometry=self.geometry,
                radius_index=radial_index,
            )
            density_values[0] = same
            density_tail = same_tail
            for layer in range(1, int(layers) + 1):
                plane = plane_exponential_sum_direct(
                    layer * spacing,
                    self._layer_delta(layer, registry),
                    kappa=self.density_params.kappa,
                    amplitude=self.density_params.C_rho,
                    geometry=self.geometry,
                    radius_index=radial_index,
                )
                transformed = self._transform_plane_result(plane, layer)
                terms = tuple(2.0 * value for value in (plane.value, *transformed))
                density_values = [
                    left + right for left, right in zip(density_values, terms)
                ]
                density_tail += 2.0 * plane.continuum_tail_estimate
            omitted = math.exp(-self.density_params.kappa * spacing * layers)
            density_tail += max(1.0, density_values[0]) * omitted
        if self.embedding_is_zero:
            embedding_values = [0.0] * 6
        else:
            assert self.embedding is not None
            rho = density_values[0]
            first = float(self.embedding.first_derivative(rho))
            second = float(self.embedding.second_derivative(rho))
            embedding_values = [
                float(self.embedding.value(rho)),
                first * density_values[1],
                first * density_values[2],
                second * density_values[1] ** 2 + first * density_values[3],
                second * density_values[1] * density_values[2] + first * density_values[4],
                second * density_values[2] ** 2 + first * density_values[5],
            ]
        total = [left + right for left, right in zip(pair_values, embedding_values)]
        return DirectFullFCCEvaluation(
            pair_energy=pair_values[0], density=density_values[0],
            embedding_energy=embedding_values[0], energy=total[0],
            d_da=total[1], d_ds=total[2], d2_daa=total[3],
            d2_das=total[4], d2_dss=total[5],
            radial_index=int(radial_index), layers=int(layers),
            estimated_pair_tail=float(
                factor * (self.p.sigma_lj**12 * tail6 + self.p.sigma_lj**6 * tail3)
            ),
            estimated_density_tail=float(density_tail),
        )
