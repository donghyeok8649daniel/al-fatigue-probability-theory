"""Isolated SI unit contract. No production imports, material defaults or solver.

Numbers must be Decimal, strings or integers. Energy representative labels are
preserved: converting eV to J never converts per-cell energy to per-atom energy.
"""
from dataclasses import dataclass
from decimal import Decimal


class UnitContractError(ValueError):
    pass


def number(value):
    if isinstance(value, bool) or not isinstance(value, (str, int, Decimal)):
        raise UnitContractError("Use a decimal string, integer or Decimal")
    try:
        result = Decimal(value)
    except Exception as exc:
        raise UnitContractError("Invalid decimal value") from exc
    if not result.is_finite():
        raise UnitContractError("A finite value is required")
    return result


# Conversion constants already declared in solver_v1/aluminum_calibration.py.
EV_J = Decimal("1.602176634e-19")
ANGSTROM_M = Decimal("1e-10")
_UNITS = {
    "m": ("length", Decimal(1)),
    "angstrom": ("length", ANGSTROM_M),
    "nm": ("length", Decimal("1e-9")),
    "m^2": ("area", Decimal(1)),
    "angstrom^2": ("area", Decimal("1e-20")),
    "nm^2": ("area", Decimal("1e-18")),
    "m^3": ("volume", Decimal(1)),
    "angstrom^3": ("volume", Decimal("1e-30")),
    "nm^3": ("volume", Decimal("1e-27")),
    "J": ("energy", Decimal(1)),
    "eV": ("energy", EV_J),
    "Pa": ("stress", Decimal(1)),
    "MPa": ("stress", Decimal("1e6")),
    "GPa": ("stress", Decimal("1e9")),
}


@dataclass(frozen=True)
class SIQuantity:
    value: Decimal
    dimension: str
    representative: str | None = None


def to_si(value, unit, *, representative=None):
    normalized = unit.replace("Å", "angstrom").replace("²", "^2").replace("³", "^3")
    if normalized not in _UNITS:
        raise UnitContractError("Unsupported unit; keep energy basis separate from its unit")
    dimension, multiplier = _UNITS[normalized]
    if dimension == "energy" and (not isinstance(representative, str) or not representative.strip()):
        raise UnitContractError("Energy needs an explicit representative label")
    return SIQuantity(number(value) * multiplier, dimension, representative)


def assert_compatible(left, right):
    if left.dimension != right.dimension:
        raise UnitContractError("Physical dimensions differ")
    if left.representative != right.representative:
        raise UnitContractError("Representative bases differ; no implicit count conversion")


@dataclass(frozen=True)
class ReducedUnits:
    model_id: str
    parameter_reference: str
    coordinate_reference: str
    energy_representative: str
    length_scale_m: Decimal | str | int | None
    energy_scale_J: Decimal | str | int | None

    def __post_init__(self):
        for name in ("model_id", "parameter_reference", "coordinate_reference", "energy_representative"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise UnitContractError("All provenance labels are required")
        for name in ("length_scale_m", "energy_scale_J"):
            value = getattr(self, name)
            if value is not None:
                parsed = number(value)
                if parsed <= 0:
                    raise UnitContractError("Declared scales must be positive")
                object.__setattr__(self, name, parsed)

    @property
    def basis(self):
        return repr((self.model_id, self.parameter_reference, self.coordinate_reference, self.energy_representative))

    def _require(self, name):
        value = getattr(self, name)
        if value is None:
            raise UnitContractError("Unresolved scale: " + name)
        return value

    def energy(self, w_star):
        return SIQuantity(number(w_star) * self._require("energy_scale_J"), "energy", self.basis)

    def energy_gradient(self, derivative_star):
        # Positive derivative of energy, not automatically the negative force.
        value = number(derivative_star) * self._require("energy_scale_J") / self._require("length_scale_m")
        return SIQuantity(value, "energy/length", self.basis)

    def energy_curvature(self, curvature_star):
        value = number(curvature_star) * self._require("energy_scale_J") / self._require("length_scale_m")**2
        return SIQuantity(value, "energy/length^2", self.basis)


# No automatic actual volume, actual interface area, temperature, mobility,
# physical clock, strain tensor, energy multiplicity or production activation.
