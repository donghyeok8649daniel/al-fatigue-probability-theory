"""One-axis contract adapter only. No PDE, material closure, or spatial operator.

All q coordinates and probability masses are dimensionless. This module never
imports the production solver. Scalars must be supplied explicitly; tensors are
not projected onto a scalar load. See C_V1_1_CROSS_REVIEW_KO.txt.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import math
from typing import Callable


class ContractError(ValueError):
    pass


def scalar(name: str, value: object, *, positive: bool = False) -> float:
    if type(value) not in (int, float) or not math.isfinite(value):
        raise ContractError(f"{name}: explicit finite scalar required")
    if positive and value <= 0:
        raise ContractError(f"{name}: positive value required")
    return float(value)


def fingerprint(value: object) -> str:
    return hashlib.sha256(
        json.dumps(asdict(value), sort_keys=True, allow_nan=False).encode("utf-8")
    ).hexdigest()


@dataclass(frozen=True)
class Node:
    a_star: float
    s_star: float
    well_index: int
    xi_star: float


@dataclass(frozen=True)
class EnergyBundleReference:
    # All three identifiers are caller declarations, not independently verified.
    energy_model_id: str | None = None
    parameter_fingerprint: str | None = None
    energy_counting_source: str | None = None


@dataclass(frozen=True)
class ModelDeclaration:
    a0_star: float
    b_star: float
    chi: float
    kappa: float
    source: str
    n_original_cells: int = 1
    bundle: EnergyBundleReference = EnergyBundleReference()


@dataclass(frozen=True)
class ScaleDeclaration:
    # Unknown values stay None. No default material calibration is supplied.
    E0_J_per_representative: float | None = None
    L0_m: float | None = None
    same_L0_for_a_and_s: bool | None = None
    source: str | None = None
    bundle: EnergyBundleReference = EnergyBundleReference()


@dataclass(frozen=True)
class CellDeclaration:
    # An energy multiplicity does not establish event independence.
    representative_count: int | None = None
    actual_geometric_volume_m3: float | None = None
    actual_interface_area_m2: float | None = None
    actual_loaded_area_m2: float | None = None
    actual_axis_gauge_length_m: float | None = None
    E_x_Pa: float | None = None
    geometry_source: str | None = None
    E_x_direction_temperature_constraint_source: str | None = None
    axis_l: tuple[float, float, float] | None = None
    interface_normal_n: tuple[float, float, float] | None = None
    slip_direction_m: tuple[float, float, float] | None = None
    single_interface_supplies_axial_strain: bool | None = None
    loaded_prism_geometry_declared: bool = False
    reference_small_strain_work_pair_declared: bool = False
    caller_work_audit_accepted: bool = False
    energy_basis: str = "E0_per_original_representative"
    bundle: EnergyBundleReference = EnergyBundleReference()
    probability_basis: str = "one_original_representative_unmapped_to_physical_cell"


@dataclass(frozen=True)
class DynamicsDeclaration:
    # theta is the thermal energy ratio, not the time coordinate or load scenario.
    theta_kBT_over_E0: float | None = None
    mobility_a_internal: float | None = None
    mobility_s_internal: float | None = None
    internal_time_definition: str | None = None
    seconds_per_internal_time: float | None = None
    source: str | None = None


@dataclass(frozen=True)
class Declaration:
    model: ModelDeclaration
    scales: ScaleDeclaration
    cell: CellDeclaration
    dynamics: DynamicsDeclaration
    nodes: tuple[Node, ...]
    support_source: str
    # p entries are probability masses, including quadrature weights if needed.
    probability_measure: str = "node_probability_masses"
    evolution_mode: str = "interface_only"
    tolerance: float = 1.0e-12


@dataclass(frozen=True)
class ProbabilityState:
    time_internal: float
    alive_mass: tuple[float, ...]
    opening_absorbed_mass: float


@dataclass(frozen=True)
class LoadKnot:
    time_internal: float
    sigma_x_over_E_x: float


@dataclass(frozen=True)
class LoadHistory:
    knots: tuple[LoadKnot, ...]
    interpolation: str
    axis_source: str
    stress_input_kind: str = "explicit_signed_one_axis_ratio"


@dataclass(frozen=True)
class TrialRequest:
    declaration: Declaration
    base_revision: int
    base_fingerprint: str
    base_state: ProbabilityState
    load_history: LoadHistory
    force_star_at_knots: tuple[float, ...]
    purpose: str


@dataclass(frozen=True)
class Transition:
    alive_mass: tuple[float, ...]
    opening_absorbed_mass: float
    backend_report: str
    purpose: str


@dataclass(frozen=True)
class Trial:
    token: int
    request: TrialRequest
    endpoint: ProbabilityState
    backend_report: str


ProbabilityBackend = Callable[[TrialRequest], Transition]


def bundle_provenance_status(d: Declaration) -> dict:
    bundles = (d.model.bundle, d.scales.bundle, d.cell.bundle)
    fields = ("energy_model_id", "parameter_fingerprint", "energy_counting_source")
    missing, mismatches = [], []
    for name in fields:
        values = [getattr(bundle, name) for bundle in bundles]
        if any(value is None for value in values):
            missing.append(name)
        if len({value for value in values if value is not None}) > 1:
            mismatches.append(name)
    return {
        "basis": "caller_declaration_only",
        "identifiers_match_and_complete": not missing and not mismatches,
        "independently_verified": False,
        "original_parameter_fingerprint_recomputed": False,
        "missing": missing,
        "mismatches": mismatches,
        "declared_bundles": [asdict(bundle) for bundle in bundles],
    }


def validate_declaration(d: Declaration) -> None:
    m = d.model
    if type(m.n_original_cells) is not int or m.n_original_cells != 1:
        raise ContractError("adapter is restricted to original N=1, not spatial cell count")
    for name in ("a0_star", "b_star", "kappa"):
        scalar(name, getattr(m, name), positive=True)
    scalar("chi", m.chi)
    scalar("tolerance", d.tolerance, positive=True)
    if d.tolerance > 1.0e-8:
        raise ContractError("tolerance too large for this interface audit")
    if not m.source or not d.support_source:
        raise ContractError("original-model and support provenance required")
    for bundle in (d.model.bundle, d.scales.bundle, d.cell.bundle):
        if not isinstance(bundle, EnergyBundleReference):
            raise ContractError("explicit immutable EnergyBundleReference required")
        for name, value in asdict(bundle).items():
            if value is not None and (type(value) is not str or not value.strip()):
                raise ContractError(name + ": nonempty caller identifier or None required")
    provenance = bundle_provenance_status(d)
    if provenance["mismatches"]:
        raise ContractError("mixed energy bundle: " + ", ".join(provenance["mismatches"]))
    if d.cell.energy_basis != "E0_per_original_representative":
        raise ContractError("aggregate E0 and representative count must not be double counted")
    if d.probability_measure != "node_probability_masses":
        raise ContractError("densities require explicitly weighted masses before entry")
    if d.evolution_mode not in ("interface_only", "declared_dimensionless_backend"):
        raise ContractError("unsupported evolution mode; no physical/spatial solve mode exists")
    if type(d.nodes) is not tuple or not d.nodes:
        raise ContractError("immutable nonempty support tuple required")
    for q in d.nodes:
        for name in ("a_star", "s_star", "xi_star"):
            scalar(name, getattr(q, name))
        if type(q.well_index) is not int:
            raise ContractError("integer well memory required")
        if not (-m.b_star / 2 <= q.xi_star < m.b_star / 2):
            raise ContractError("xi must use declared half-open [-b/2,b/2) convention")
        if not math.isclose(q.s_star, m.b_star * q.well_index + q.xi_star,
                            rel_tol=d.tolerance, abs_tol=d.tolerance):
            raise ContractError("s != b*n + xi; do not discard well memory")
    if len(set(d.nodes)) != len(d.nodes):
        raise ContractError("duplicate support nodes require an explicit measure decision")
    for obj, names in (
        (d.scales, ("E0_J_per_representative", "L0_m")),
        (d.cell, ("actual_geometric_volume_m3", "actual_interface_area_m2",
                  "actual_loaded_area_m2", "actual_axis_gauge_length_m", "E_x_Pa")),
        (d.dynamics, ("mobility_a_internal", "mobility_s_internal", "seconds_per_internal_time")),
    ):
        for name in names:
            value = getattr(obj, name)
            if value is not None:
                scalar(name, value, positive=True)
    count = d.cell.representative_count
    if count is not None and (type(count) is not int or count <= 0):
        raise ContractError("representative_count must be an explicitly declared positive integer")
    if d.scales.same_L0_for_a_and_s is not None and type(d.scales.same_L0_for_a_and_s) is not bool:
        raise ContractError("same_L0_for_a_and_s must be bool or None")
    theta = d.dynamics.theta_kBT_over_E0
    if theta is not None and scalar("theta_kBT_over_E0", theta) < 0:
        raise ContractError("thermal energy ratio cannot be negative")
    for name in ("loaded_prism_geometry_declared", "reference_small_strain_work_pair_declared",
                 "caller_work_audit_accepted"):
        if type(getattr(d.cell, name)) is not bool:
            raise ContractError(name + " must be bool")
    if (d.cell.single_interface_supplies_axial_strain is not None
            and type(d.cell.single_interface_supplies_axial_strain) is not bool):
        raise ContractError("single_interface_supplies_axial_strain must be bool or None")
    for name in ("axis_l", "interface_normal_n", "slip_direction_m"):
        vector = getattr(d.cell, name)
        if vector is not None:
            if type(vector) is not tuple or len(vector) != 3:
                raise ContractError(name + ": immutable 3-vector required")
            for value in vector:
                scalar(name, value)
            if not math.isclose(math.fsum(x*x for x in vector), 1.0,
                                rel_tol=d.tolerance, abs_tol=d.tolerance):
                raise ContractError(name + ": declared unit direction required")
    n, slip = d.cell.interface_normal_n, d.cell.slip_direction_m
    if n is not None and slip is not None:
        if abs(math.fsum(x*y for x, y in zip(n, slip))) > d.tolerance:
            raise ContractError("declared slip direction must be tangent to interface")
    if d.cell.probability_basis != "one_original_representative_unmapped_to_physical_cell":
        raise ContractError("physical-cell event aggregation is outside this adapter")


def validate_state(d: Declaration, p: ProbabilityState) -> None:
    scalar("time_internal", p.time_internal)
    if type(p.alive_mass) is not tuple or len(p.alive_mass) != len(d.nodes):
        raise ContractError("immutable masses must match fixed support")
    values = p.alive_mass + (p.opening_absorbed_mass,)
    for value in values:
        if scalar("probability mass", value) < 0:
            raise ContractError("negative mass is rejected, not clipped")
    if abs(math.fsum(values) - 1.0) > d.tolerance:
        raise ContractError("alive + cumulative opening absorption must equal 1")


def force_from_ratio(d: Declaration, ratio: float) -> float:
    """The only stress input accepted here is a declared signed scalar sigma_x/E_x."""
    return scalar("f_star", d.model.kappa * scalar("sigma_x/E_x", ratio))


def dimensionless_evolution_missing(d: Declaration) -> list[str]:
    missing = []
    for name in ("theta_kBT_over_E0", "mobility_a_internal", "mobility_s_internal",
                 "internal_time_definition", "source"):
        if getattr(d.dynamics, name) in (None, ""):
            missing.append("dynamics." + name)
    for name in bundle_provenance_status(d)["missing"]:
        missing.append("energy_bundle." + name)
    return missing


def normal_slip_geometry_status(d: Declaration) -> dict:
    """A's single-interface example only; no geometry is generated or fitted."""
    validate_declaration(d)
    c, s, m = d.cell, d.scales, d.model
    missing = []
    for name in ("actual_geometric_volume_m3", "actual_interface_area_m2", "axis_l",
                 "interface_normal_n", "slip_direction_m", "geometry_source"):
        if getattr(c, name) in (None, ""):
            missing.append("cell." + name)
    if s.L0_m is None or s.same_L0_for_a_and_s is not True:
        missing.append("common_L0")
    if c.single_interface_supplies_axial_strain is not True:
        missing.append("single_interface_supplies_axial_strain=True")
    normal_scaled = slip_scaled = None
    normal_consistent = slip_consistent = False
    if not missing:
        alpha = math.fsum(x*y for x, y in zip(c.axis_l, c.interface_normal_n))
        beta = math.fsum(x*y for x, y in zip(c.axis_l, c.slip_direction_m))
        prefactor = scalar("Aint*L0*a0/Vgeom",
                           c.actual_interface_area_m2 * s.L0_m * m.a0_star
                           / c.actual_geometric_volume_m3, positive=True)
        normal_scaled = scalar("normal scaled coefficient", prefactor * alpha * alpha)
        slip_scaled = scalar("slip scaled coefficient", prefactor * alpha * beta)
        normal_consistent = math.isclose(normal_scaled, 1.0, rel_tol=d.tolerance, abs_tol=d.tolerance)
        slip_consistent = math.isclose(slip_scaled, m.chi, rel_tol=d.tolerance, abs_tol=d.tolerance)
    return {
        "scope": "declared single-interface jump n*da+m*ds supplies adopted axial strain",
        "missing": missing,
        "normal_coefficient_checked": not missing,
        "slip_coefficient_checked": not missing,
        "normal_scaled_coefficient_target_1": normal_scaled,
        "slip_scaled_coefficient_target_chi": slip_scaled,
        "normal_coefficient_consistent": normal_consistent,
        "slip_coefficient_consistent": slip_consistent,
        "normal_and_slip_consistent_on_caller_geometry": normal_consistent and slip_consistent,
        "actual_normal_slip_geometry_independently_verified": False,
    }


def work_normalization_status(d: Declaration) -> dict:
    """Conditional algebraic audit; candidate scales are not measured geometry."""
    validate_declaration(d)
    s, c, m = d.scales, d.cell, d.model
    provenance = bundle_provenance_status(d)
    missing = []
    for prefix, obj, names in (
        ("scales", s, ("E0_J_per_representative", "L0_m", "source")),
        ("cell", c, ("representative_count", "actual_geometric_volume_m3",
                     "actual_interface_area_m2", "actual_loaded_area_m2",
                     "actual_axis_gauge_length_m", "E_x_Pa", "geometry_source",
                     "E_x_direction_temperature_constraint_source")),
    ):
        for name in names:
            if getattr(obj, name) in (None, ""):
                missing.append(prefix + "." + name)
    if s.same_L0_for_a_and_s is not True:
        missing.append("scales.same_L0_for_a_and_s=True")
    if not c.caller_work_audit_accepted:
        missing.append("cell.caller_work_audit_accepted")
    if not c.reference_small_strain_work_pair_declared:
        missing.append("cell.reference_small_strain_work_pair_declared")
    missing.extend("energy_bundle." + name for name in provenance["missing"])
    geometry = normal_slip_geometry_status(d)
    missing.extend(geometry["missing"])
    candidate_volume = candidate_area = None
    coefficient_ratio = None
    volume_residual = None
    loaded_prism_ratio = None
    if (provenance["identifiers_match_and_complete"]
            and all(value is not None for value in
                    (s.E0_J_per_representative, c.representative_count, c.E_x_Pa))):
        candidate_volume = scalar("conditional work-effective volume",
                                  c.representative_count * s.E0_J_per_representative
                                  * m.kappa * m.a0_star / c.E_x_Pa, positive=True)
        if s.L0_m is not None and s.same_L0_for_a_and_s is True:
            candidate_area = scalar("conditional work-effective area",
                                    candidate_volume / (s.L0_m * m.a0_star), positive=True)
        if c.actual_geometric_volume_m3 is not None:
            coefficient_ratio = scalar("Veff/Vgeom", candidate_volume / c.actual_geometric_volume_m3)
            volume_residual = scalar("Ex*Vgeom-r*E0*kappa*a0",
                                     c.E_x_Pa * c.actual_geometric_volume_m3
                                     - c.representative_count * s.E0_J_per_representative
                                     * m.kappa * m.a0_star)
    if all(value is not None for value in
           (c.actual_geometric_volume_m3, c.actual_loaded_area_m2, c.actual_axis_gauge_length_m)):
        loaded_prism_ratio = scalar("Aloaded*Hactual/Vgeom",
                                    c.actual_loaded_area_m2 * c.actual_axis_gauge_length_m
                                    / c.actual_geometric_volume_m3, positive=True)
    mismatches = []
    checks = [("work_effective_volume_vs_actual_volume", coefficient_ratio)]
    if c.loaded_prism_geometry_declared:
        checks.append(("declared_loaded_prism_Aloaded_Hactual_vs_Vgeom", loaded_prism_ratio))
    for name, value in checks:
        if value is not None and not math.isclose(value, 1.0, rel_tol=d.tolerance,
                                                  abs_tol=d.tolerance):
            mismatches.append(name)
    if geometry["normal_coefficient_checked"] and not geometry["normal_coefficient_consistent"]:
        mismatches.append("single_interface_normal_coefficient")
    if geometry["slip_coefficient_checked"] and not geometry["slip_coefficient_consistent"]:
        mismatches.append("single_interface_slip_coefficient")
    return {
        "conditional_identity": "Veff = r*E0*kappa*a0_star/E_x",
        "conditional_work_effective_area_identity": "Aeff = Veff/(L0*a0_star)",
        "work_effective_volume_m3_not_actual_geometry": candidate_volume,
        "work_effective_area_m2_not_actual_loaded_or_interface_area": candidate_area,
        "actual_geometric_volume_m3_declared": c.actual_geometric_volume_m3,
        "actual_interface_area_m2_declared": c.actual_interface_area_m2,
        "actual_loaded_area_m2_declared": c.actual_loaded_area_m2,
        "actual_axis_gauge_length_m_declared": c.actual_axis_gauge_length_m,
        "model_reference_spacing_m_not_axis_gauge_length": (
            s.L0_m * m.a0_star if s.L0_m is not None and s.same_L0_for_a_and_s is True else None),
        "gauge_length_inferred_from_model_spacing": False,
        "work_effective_to_actual_volume_ratio": coefficient_ratio,
        "actual_volume_residual_Pa_m3": volume_residual,
        "declared_loaded_prism_Aloaded_Hactual_over_Vgeom": loaded_prism_ratio,
        "loaded_prism_identity_checked": c.loaded_prism_geometry_declared and loaded_prism_ratio is not None,
        "bundle_provenance": provenance,
        "normal_slip_geometry_audit": geometry,
        "missing": sorted(set(missing)),
        "mismatches": mismatches,
        "one_axis_declared_checks_consistent": not missing and not mismatches,
        "physical_work_ready": False,
        "physical_ready_reason": "caller-declared algebra is not independent geometry/material calibration",
        "spatial_solver_ready": False,
        "physical_cell_event_mapping_ready": False,
    }


def require_physical_work_ready(d: Declaration) -> None:
    result = work_normalization_status(d)
    raise ContractError("physical work not certified by this adapter: "
                        + ", ".join(result["missing"] + result["mismatches"]
                                    + [result["physical_ready_reason"]]))


def require_one_axis_declared_checks_consistent(d: Declaration) -> None:
    result = work_normalization_status(d)
    if not result["one_axis_declared_checks_consistent"]:
        raise ContractError("declared axis checks blocked: "
                            + ", ".join(result["missing"] + result["mismatches"]))


def observe(d: Declaration, p: ProbabilityState, ratio: float) -> dict:
    validate_declaration(d)
    validate_state(d, p)
    m = d.model
    f_star = force_from_ratio(d, ratio)
    survival = math.fsum(p.alive_mass)
    terms = {
        "opening": tuple((q.a_star - m.a0_star) / m.a0_star for q in d.nodes),
        "slip_well_memory": tuple(m.chi * m.b_star * q.well_index / m.a0_star for q in d.nodes),
        "slip_wrapped": tuple(m.chi * q.xi_star / m.a0_star for q in d.nodes),
    }
    moments = {key: math.fsum(w * v for w, v in zip(p.alive_mass, values))
               for key, values in terms.items()}
    moments["total"] = math.fsum(moments.values())
    conditional = ({key: value / survival for key, value in moments.items()}
                   if survival > d.tolerance else None)
    mean_h = m.a0_star * moments["total"]
    return {
        "f_star": f_star,
        "survival_mass": survival,
        "opening_absorbed_mass_cumulative": p.opening_absorbed_mass,
        "eps_x_alive_unnormalized_moments": moments,
        "eps_x_survival_conditional": conditional,
        "conditional_mean_status": "defined" if conditional is not None else "insufficient_survival_mass",
        "endpoint_external_potential_star_alive_unnormalized": -f_star * mean_h,
        "endpoint_external_potential_star_survival_conditional": (
            -f_star * mean_h / survival if conditional is not None else None),
        "metadata": {
            "schema_version": "one-axis-1.1",
            "declaration_fingerprint": fingerprint(d),
            "state_fingerprint": fingerprint(p),
            "original_energy_bundle_provenance": bundle_provenance_status(d),
            "time_internal": p.time_internal,
            "force_definition": "f_star=kappa*(sigma_x/E_x)",
            "strain_definition": "[(a_star-a0_star)+chi*s_star]/a0_star",
            "energy_definition": "G_star=W_star-f_star*h_star",
            "external_potential_is_not_accumulated_path_work": True,
            "absorbed_state_strain_not_assigned": True,
            "dimensionless_evolution_missing": dimensionless_evolution_missing(d),
            "seconds_per_internal_time": d.dynamics.seconds_per_internal_time,
            "physical_clock_declared_not_validated": bool(d.dynamics.seconds_per_internal_time
                                                          and d.dynamics.source
                                                          and not dimensionless_evolution_missing(d)),
            "physical_evolution_validated": False,
            "work_normalization": work_normalization_status(d),
            "no_bulk_tangent_or_spatial_operator": True,
        },
    }


class TrialCommitInterface:
    """Each candidate starts from one immutable P_n. Commit invalidates siblings.

    A supplied callable owns any probability integration. This adapter checks
    only the contract and mass accounting, not the accuracy of that integrator.
    """

    def __init__(self, declaration: Declaration, initial_state: ProbabilityState):
        validate_declaration(declaration)
        validate_state(declaration, initial_state)
        self._declaration = declaration
        self._state = initial_state
        self._revision = 0
        self._next_token = 0
        self._issued: dict[int, Trial] = {}

    @property
    def declaration(self) -> Declaration:
        return self._declaration

    @property
    def state(self) -> ProbabilityState:
        return self._state

    @property
    def revision(self) -> int:
        return self._revision

    def trial(self, history: LoadHistory, backend: ProbabilityBackend, *, purpose: str) -> Trial:
        d = self.declaration
        if history.stress_input_kind != "explicit_signed_one_axis_ratio":
            raise ContractError("raw weak stress and tensors require an external audited interface")
        if not history.axis_source or history.interpolation not in (
                "piecewise_linear", "piecewise_constant_left"):
            raise ContractError("explicit axis source and load interpolation required")
        if type(history.knots) is not tuple or len(history.knots) < 2:
            raise ContractError("immutable load history needs at least two knots")
        times = tuple(scalar("knot time", k.time_internal) for k in history.knots)
        if times[0] != self._state.time_internal or any(b <= a for a, b in zip(times, times[1:])):
            raise ContractError("history must start exactly at P_n time and strictly advance")
        forces = tuple(force_from_ratio(d, k.sigma_x_over_E_x) for k in history.knots)
        if purpose == "declared_dimensionless_backend":
            missing = dimensionless_evolution_missing(d)
            if d.evolution_mode != purpose or missing:
                raise ContractError("probability evolution blocked: " + ", ".join(missing or ["evolution_mode"]))
        elif purpose != "interface_logical_fixture" or d.evolution_mode != "interface_only":
            raise ContractError("purpose/mode mismatch; no physical/spatial solve is provided")
        request = TrialRequest(d, self._revision, fingerprint(self._state), self._state,
                               history, forces, purpose)
        # The callable receives the same immutable state for every corrector.
        transition = backend(request)
        if not isinstance(transition, Transition) or transition.purpose != purpose or not transition.backend_report:
            raise ContractError("backend must declare its purpose and report provenance")
        endpoint = ProbabilityState(times[-1], transition.alive_mass, transition.opening_absorbed_mass)
        validate_state(d, endpoint)
        if endpoint.opening_absorbed_mass < self._state.opening_absorbed_mass:
            raise ContractError("cumulative opening absorption cannot decrease")
        trial = Trial(self._next_token, request, endpoint, transition.backend_report)
        self._next_token += 1
        self._issued[trial.token] = trial
        return trial

    def commit(self, trial: Trial, *, coupled_residuals_accepted: bool) -> dict:
        # Caller must explicitly assert mechanical/probability/work acceptance.
        # This is an assertion record, not a residual convergence computation.
        if coupled_residuals_accepted is not True:
            raise ContractError("commit requires explicit external residual acceptance")
        if self._issued.get(trial.token) is not trial:
            raise ContractError("trial is foreign, stale, or already committed")
        if (trial.request.base_revision != self._revision
                or trial.request.base_fingerprint != fingerprint(self._state)):
            raise ContractError("trial does not start at current P_n")
        old_absorbed = self._state.opening_absorbed_mass
        self._state = trial.endpoint
        self._revision += 1
        self._issued.clear()
        return {
            "committed_revision": self._revision,
            "base_fingerprint": trial.request.base_fingerprint,
            "endpoint_fingerprint": fingerprint(self._state),
            "opening_absorbed_mass_increment": self._state.opening_absorbed_mass - old_absorbed,
            "coupled_residuals_accepted_by_caller": True,
            "residuals_computed_by_adapter": False,
            "purpose": trial.request.purpose,
        }
