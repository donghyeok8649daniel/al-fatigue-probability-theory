"""Material identity and wafer setup, separate from the Al energy registry.

Wafer doping is user-supplied composition metadata. It does not select a
validated silicon potential, carrier density, resistivity, or physical clock.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from math import isfinite
from numbers import Real


ALUMINUM = "aluminum"
SILICON_WAFER = "silicon_wafer"
MATERIAL_IDS = (ALUMINUM, SILICON_WAFER)
DOPANT_SPECIES = ("B", "P", "As", "Sb")


class MaterialInputError(ValueError):
    def __init__(self, key: str):
        super().__init__(key)
        self.ui_error_data = {"key": key, "values": {}}


@dataclass(frozen=True)
class MaterialSelection:
    material_id: str = ALUMINUM
    doping_enabled: bool = False
    dopant_species: str | None = None
    dopant_concentration_cm3: float | None = None

    def validate(self) -> None:
        if self.material_id not in MATERIAL_IDS or type(self.doping_enabled) is not bool:
            raise MaterialInputError("error.material_selection")
        if self.material_id == ALUMINUM and self.doping_enabled:
            raise MaterialInputError("error.material_selection")
        if not self.doping_enabled:
            if self.dopant_species is not None or self.dopant_concentration_cm3 is not None:
                raise MaterialInputError("error.material_inactive_doping")
            return
        if self.dopant_species not in DOPANT_SPECIES:
            raise MaterialInputError("error.dopant_species")
        value = self.dopant_concentration_cm3
        if (isinstance(value, bool) or not isinstance(value, Real)
                or not isfinite(value) or value <= 0):
            raise MaterialInputError("error.dopant_concentration")

    def require_backend(self) -> None:
        self.validate()
        if self.material_id == SILICON_WAFER:
            raise MaterialInputError("error.silicon_backend_unavailable")

    def metadata(self) -> dict[str, object]:
        self.validate()
        return {**asdict(self), "material_display_key": f"material.{self.material_id}",
                "dopant_concentration_basis": (
                    "dopant_atoms_per_cm3" if self.doping_enabled else None)}


def default_material_draft() -> dict[str, object]:
    return {"material_id": ALUMINUM, "doping_enabled": False,
            "dopant_species": "B", "dopant_concentration_cm3": ""}


def validate_material_draft(value) -> dict[str, object]:
    """Validate project structure while allowing unfinished numerical input.

    The concentration text is preserved verbatim, like other project entry fields.
    Only material_from_draft turns active inputs into a validated configuration.
    """
    if (not isinstance(value, dict) or set(value) != set(default_material_draft())
            or value["material_id"] not in MATERIAL_IDS
            or type(value["doping_enabled"]) is not bool
            or value["dopant_species"] not in DOPANT_SPECIES
            or not isinstance(value["dopant_concentration_cm3"], str)):
        raise MaterialInputError("error.material_selection")
    return dict(value)


def material_from_draft(value) -> MaterialSelection:
    draft = validate_material_draft(value)
    active = draft["material_id"] == SILICON_WAFER and draft["doping_enabled"]
    concentration = None
    if active:
        try:
            concentration = float(draft["dopant_concentration_cm3"])
        except (ValueError, OverflowError):
            raise MaterialInputError("error.dopant_concentration") from None
    result = MaterialSelection(draft["material_id"], bool(active),
                               draft["dopant_species"] if active else None, concentration)
    result.validate()
    return result
