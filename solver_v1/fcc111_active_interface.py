"""Static FCC(111) active-interface energy embedded in two bulk half-crystals.

The lower half has plane indices l<=0 and the upper half l>=1.  Opening and
rigid registry displacement affect only cross-interface interactions.  The
result is an energy difference per primitive interface cell, so infinite bulk
constants cancel analytically.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy.special import zeta

from .fcc111_full_energy import FullFCC111StackEnergy
from .fcc111_geometry import DIRECT_110, homogeneous_layer_registry, registry_path
from .fcc111_lattice_sum import (
    plane_exponential_sum_direct,
    plane_exponential_sum_reciprocal,
    plane_power_sum_direct,
    plane_power_sum_reciprocal,
)


@dataclass(frozen=True)
class ActiveInterfaceEvaluation:
    energy: float
    d_da: float
    d_ds: float
    d2_daa: float
    d2_das: float
    d2_dss: float
    pair_energy: float
    embedding_energy: float
    layers_used: int
    maximum_plane_shells: int
    estimated_tail_absolute: float


def _power_mean_difference(
    *, exponent: float, h: float, delta_a: float, area: float
) -> tuple[float, float, float]:
    """Exact ``sum k[(kh+da)^-r-(kh)^-r]`` and opening derivatives."""

    r = 2.0 * float(exponent) - 2.0
    x = float(delta_a) / float(h)
    if x <= -1.0:
        raise ValueError("active interface spacing must be positive")
    prefactor = math.pi / (float(area) * (float(exponent) - 1.0))
    value_sum = (
        zeta(r - 1.0, 1.0 + x)
        - x * zeta(r, 1.0 + x)
        - zeta(r - 1.0, 1.0)
    )
    first_sum = -r * (
        zeta(r, 1.0 + x) - x * zeta(r + 1.0, 1.0 + x)
    )
    second_sum = r * (r + 1.0) * (
        zeta(r + 1.0, 1.0 + x) - x * zeta(r + 2.0, 1.0 + x)
    )
    return (
        float(prefactor * h**(-r) * value_sum),
        float(prefactor * h**(-r - 1.0) * first_sum),
        float(prefactor * h**(-r - 2.0) * second_sum),
    )


class FCC111ActiveInterface:
    """One locally opening/slipping interface between infinite bulk halves."""

    def __init__(
        self,
        bulk_model: FullFCC111StackEnergy,
        *,
        path_id: str = DIRECT_110,
        tolerance: float = 2.0e-9,
        max_layers: int = 256,
        consecutive_small_layers: int = 4,
    ):
        if not np.isfinite(bulk_model.a0):
            raise ValueError("active interface requires a stable bulk reference")
        registry_path(path_id)
        if tolerance <= 0.0 or max_layers < 8:
            raise ValueError("invalid active-interface convergence controls")
        self.bulk = bulk_model
        self.path_id = path_id
        self.geometry = bulk_model.geometry
        self.h = float(bulk_model.a0)
        self.tolerance = float(tolerance)
        self.max_layers = int(max_layers)
        self.consecutive_small_layers = int(consecutive_small_layers)
        self.direction = registry_path(path_id).direction()
        self.rho_bulk = (
            bulk_model.full_environment_density(self.h, 0.0).value
            if bulk_model.has_embedding else 0.0
        )

    @property
    def period(self) -> float:
        return registry_path(self.path_id).period_over_b * self.bulk.p.b

    def _baseline_delta(self, layer: int) -> np.ndarray:
        return homogeneous_layer_registry(
            layer, 0.0, geometry=self.geometry, path_id=self.path_id,
            reduce_to_cell=True,
        )

    def _active_delta(self, layer: int, s: float) -> np.ndarray:
        return self._baseline_delta(layer) + float(s) * self.direction

    def _power_cross_difference(
        self, a: float, s: float, exponent: float
    ) -> tuple[np.ndarray, int, int, float]:
        delta_a = float(a) - self.h
        mean = _power_mean_difference(
            exponent=exponent, h=self.h, delta_a=delta_a,
            area=self.geometry.atomic_cell_area,
        )
        total = np.array([mean[0], mean[1], 0.0, mean[2], 0.0, 0.0])
        small = 0
        recent = []
        maximum_shells = 0
        for layer in range(1, self.max_layers + 1):
            base_d = layer * self.h
            active_d = base_d + delta_a
            base = plane_power_sum_reciprocal(
                base_d, self._baseline_delta(layer), p=exponent,
                geometry=self.geometry, config=self.bulk.stack_config.reciprocal,
            )
            active = plane_power_sum_reciprocal(
                active_d, self._active_delta(layer, s), p=exponent,
                geometry=self.geometry, config=self.bulk.stack_config.reciprocal,
            )
            base_zero = math.pi / (
                self.geometry.atomic_cell_area * (exponent - 1.0)
            ) * base_d ** (2.0 - 2.0 * exponent)
            active_zero = math.pi / (
                self.geometry.atomic_cell_area * (exponent - 1.0)
            ) * active_d ** (2.0 - 2.0 * exponent)
            multiplicity = float(layer)
            term = multiplicity * np.array(
                [
                    (active.value-active_zero) - (base.value-base_zero),
                    active.d_d - (2.0-2.0*exponent)*active_zero/active_d,
                    float(self.direction @ active.grad_delta),
                    active.d2_dd - ((2.0-2.0*exponent)*(1.0-2.0*exponent)
                                    * active_zero/active_d**2),
                    float(self.direction @ active.mixed_d_delta),
                    float(self.direction @ active.hess_delta @ self.direction),
                ]
            )
            total += term
            maximum_shells = max(maximum_shells, base.shells_used, active.shells_used)
            bound = float(np.max(np.abs(term)))
            recent.append(bound)
            recent = recent[-self.consecutive_small_layers:]
            if bound <= self.tolerance * max(1.0, float(np.max(np.abs(total)))):
                small += 1
                if small >= self.consecutive_small_layers:
                    return total, layer, maximum_shells, max(recent)
            else:
                small = 0
        raise RuntimeError("active-interface pair corrugation did not converge")

    def _density_changes(
        self, a: float, s: float
    ) -> tuple[list[np.ndarray], int, int, float]:
        if not self.bulk.has_embedding:
            return [], 0, 0, 0.0
        assert self.bulk.density_params is not None
        delta_a = float(a) - self.h
        plane_terms = []
        small = 0
        recent = []
        maximum_shells = 0
        for layer in range(1, self.max_layers + 1):
            base = plane_exponential_sum_reciprocal(
                layer*self.h, self._baseline_delta(layer),
                kappa=self.bulk.density_params.kappa,
                amplitude=self.bulk.density_params.C_rho,
                geometry=self.geometry, config=self.bulk.stack_config.reciprocal,
            )
            active = plane_exponential_sum_reciprocal(
                layer*self.h+delta_a, self._active_delta(layer, s),
                kappa=self.bulk.density_params.kappa,
                amplitude=self.bulk.density_params.C_rho,
                geometry=self.geometry, config=self.bulk.stack_config.reciprocal,
            )
            term = np.array(
                [active.value-base.value, active.d_d,
                 float(self.direction @ active.grad_delta), active.d2_dd,
                 float(self.direction @ active.mixed_d_delta),
                 float(self.direction @ active.hess_delta @ self.direction)]
            )
            plane_terms.append(term)
            maximum_shells = max(maximum_shells, base.shells_used, active.shells_used)
            bound = float(np.max(np.abs(term[1:]))) + abs(float(term[0]))
            recent.append(bound)
            recent = recent[-self.consecutive_small_layers:]
            if bound <= self.tolerance * max(1.0, self.rho_bulk):
                small += 1
                if small >= self.consecutive_small_layers:
                    break
            else:
                small = 0
        else:
            raise RuntimeError("active-interface density difference did not converge")
        # A plane at depth r sees every opposite-half plane with k>=r+1.
        cumulative = np.cumsum(np.asarray(plane_terms)[::-1], axis=0)[::-1]
        return [row for row in cumulative], layer, maximum_shells, max(recent)

    def evaluate(self, a: float, s: float) -> ActiveInterfaceEvaluation:
        if float(a) <= 0.0:
            raise ValueError("active interface spacing must be positive")
        s6, layers6, shells6, tail6 = self._power_cross_difference(a, s, 6.0)
        s3, layers3, shells3, tail3 = self._power_cross_difference(a, s, 3.0)
        factor = 4.0 * self.bulk.p.epsilon
        pair = factor * (
            self.bulk.p.sigma_lj**12 * s6 - self.bulk.p.sigma_lj**6 * s3
        )
        embedding = np.zeros(6)
        density_layers, layers_rho, shells_rho, tail_rho = self._density_changes(a, s)
        if density_layers:
            assert self.bulk.embedding is not None
            for change in density_layers:
                density = self.rho_bulk + change[0]
                if density <= 0.0:
                    raise FloatingPointError("active-interface site density is nonpositive")
                first = float(self.bulk.embedding.first_derivative(density))
                second = float(self.bulk.embedding.second_derivative(density))
                # Two symmetry-related atoms occur at every distance from the cut.
                embedding += 2.0 * np.array(
                    [
                        float(self.bulk.embedding.value(density)
                              - self.bulk.embedding.value(self.rho_bulk)),
                        first*change[1], first*change[2],
                        second*change[1]**2 + first*change[3],
                        second*change[1]*change[2] + first*change[4],
                        second*change[2]**2 + first*change[5],
                    ]
                )
        total = pair + embedding
        return ActiveInterfaceEvaluation(
            energy=float(total[0]), d_da=float(total[1]), d_ds=float(total[2]),
            d2_daa=float(total[3]), d2_das=float(total[4]), d2_dss=float(total[5]),
            pair_energy=float(pair[0]), embedding_energy=float(embedding[0]),
            layers_used=max(layers6, layers3, layers_rho),
            maximum_plane_shells=max(shells6, shells3, shells_rho),
            estimated_tail_absolute=float(
                factor*(self.bulk.p.sigma_lj**12*tail6+self.bulk.p.sigma_lj**6*tail3)
                + tail_rho
            ),
        )

    def energy(self, a: float, s: float) -> float:
        return self.evaluate(a, s).energy

    def grad(self, a: float, s: float) -> np.ndarray:
        value = self.evaluate(a, s)
        return np.array([value.d_da, value.d_ds])

    def hessian(self, a: float, s: float) -> np.ndarray:
        value = self.evaluate(a, s)
        return np.array([[value.d2_daa, value.d2_das],
                         [value.d2_das, value.d2_dss]])

    def separated_limit(self) -> float:
        """Exact-mean, tolerance-controlled work of separating the two halves."""

        power_limits = {}
        for exponent in (3.0, 6.0):
            r = 2.0*exponent-2.0
            prefactor = math.pi/(
                self.geometry.atomic_cell_area*(exponent-1.0)
            )
            value = -prefactor*self.h**(-r)*float(zeta(r-1.0, 1.0))
            small = 0
            for layer in range(1, self.max_layers+1):
                base_d = layer*self.h
                base = plane_power_sum_reciprocal(
                    base_d, self._baseline_delta(layer), p=exponent,
                    geometry=self.geometry, config=self.bulk.stack_config.reciprocal,
                )
                zero = prefactor*base_d**(-r)
                term = -float(layer)*(base.value-zero)
                value += term
                if abs(term) <= self.tolerance*max(1.0,abs(value)):
                    small += 1
                    if small >= self.consecutive_small_layers:
                        break
                else:
                    small = 0
            else:
                raise RuntimeError("separated pair corrugation did not converge")
            power_limits[exponent] = value
        pair = 4.0*self.bulk.p.epsilon*(
            self.bulk.p.sigma_lj**12*power_limits[6.0]
            - self.bulk.p.sigma_lj**6*power_limits[3.0]
        )
        if not self.bulk.has_embedding:
            return float(pair)
        assert self.bulk.density_params is not None and self.bulk.embedding is not None
        plane_values = []
        small = 0
        for layer in range(1, self.max_layers+1):
            base = plane_exponential_sum_reciprocal(
                layer*self.h, self._baseline_delta(layer),
                kappa=self.bulk.density_params.kappa,
                amplitude=self.bulk.density_params.C_rho,
                geometry=self.geometry, config=self.bulk.stack_config.reciprocal,
            )
            plane_values.append(base.value)
            if abs(base.value) <= self.tolerance*max(1.0,self.rho_bulk):
                small += 1
                if small >= self.consecutive_small_layers:
                    break
            else:
                small = 0
        else:
            raise RuntimeError("separated density neighborhood did not converge")
        lost = np.cumsum(np.asarray(plane_values)[::-1])[::-1]
        embedding = 0.0
        for value in lost:
            density = self.rho_bulk-value
            embedding += 2.0*float(
                self.bulk.embedding.value(density)-self.bulk.embedding.value(self.rho_bulk)
            )
        return float(pair+embedding)

    def direct_reference(
        self, a: float, s: float, *, radial_index: int = 60, layers: int = 40
    ) -> ActiveInterfaceEvaluation:
        """Independent finite real-space/layer reference for validation only."""

        delta_a = float(a) - self.h
        pair_by_p = {}
        pair_tail = 0.0
        max_points = 0
        for exponent in (3.0, 6.0):
            total = np.zeros(6)
            last = 0.0
            for layer in range(1, int(layers) + 1):
                base = plane_power_sum_direct(
                    layer*self.h, self._baseline_delta(layer), p=exponent,
                    geometry=self.geometry, radius_index=radial_index,
                )
                active = plane_power_sum_direct(
                    layer*self.h+delta_a, self._active_delta(layer, s), p=exponent,
                    geometry=self.geometry, radius_index=radial_index,
                )
                multiplicity = float(layer)
                term = multiplicity * np.array(
                    [active.value-base.value, active.d_d,
                     float(self.direction @ active.grad_delta), active.d2_dd,
                     float(self.direction @ active.mixed_d_delta),
                     float(self.direction @ active.hess_delta @ self.direction)]
                )
                total += term
                last = float(np.max(np.abs(term)))
                max_points = max(max_points, active.points_used, base.points_used)
            pair_by_p[exponent] = total
            pair_tail += last
        factor = 4.0*self.bulk.p.epsilon
        pair = factor*(self.bulk.p.sigma_lj**12*pair_by_p[6.0]
                       - self.bulk.p.sigma_lj**6*pair_by_p[3.0])

        embedding = np.zeros(6)
        density_tail = 0.0
        if self.bulk.has_embedding:
            assert self.bulk.density_params is not None and self.bulk.embedding is not None
            bulk_direct = self.bulk.direct_reference(
                self.h, 0.0, radial_index=radial_index, layers=layers
            ).density
            plane_terms = []
            for layer in range(1, int(layers) + 1):
                base = plane_exponential_sum_direct(
                    layer*self.h, self._baseline_delta(layer),
                    kappa=self.bulk.density_params.kappa,
                    amplitude=self.bulk.density_params.C_rho,
                    geometry=self.geometry, radius_index=radial_index,
                )
                active = plane_exponential_sum_direct(
                    layer*self.h+delta_a, self._active_delta(layer, s),
                    kappa=self.bulk.density_params.kappa,
                    amplitude=self.bulk.density_params.C_rho,
                    geometry=self.geometry, radius_index=radial_index,
                )
                plane_terms.append(np.array(
                    [active.value-base.value, active.d_d,
                     float(self.direction @ active.grad_delta), active.d2_dd,
                     float(self.direction @ active.mixed_d_delta),
                     float(self.direction @ active.hess_delta @ self.direction)]
                ))
                density_tail += active.continuum_tail_estimate + base.continuum_tail_estimate
            cumulative = np.cumsum(np.asarray(plane_terms)[::-1], axis=0)[::-1]
            for change in cumulative:
                density = bulk_direct + change[0]
                first = float(self.bulk.embedding.first_derivative(density))
                second = float(self.bulk.embedding.second_derivative(density))
                embedding += 2.0*np.array(
                    [float(self.bulk.embedding.value(density)
                           - self.bulk.embedding.value(bulk_direct)),
                     first*change[1], first*change[2],
                     second*change[1]**2+first*change[3],
                     second*change[1]*change[2]+first*change[4],
                     second*change[2]**2+first*change[5]]
                )
        total = pair+embedding
        return ActiveInterfaceEvaluation(
            energy=float(total[0]), d_da=float(total[1]), d_ds=float(total[2]),
            d2_daa=float(total[3]), d2_das=float(total[4]), d2_dss=float(total[5]),
            pair_energy=float(pair[0]), embedding_energy=float(embedding[0]),
            layers_used=int(layers), maximum_plane_shells=max_points,
            estimated_tail_absolute=float(abs(pair_tail)+density_tail),
        )
