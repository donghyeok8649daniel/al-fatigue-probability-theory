from dataclasses import replace

import pytest

from app.materials import (
    ALUMINUM, SILICON_WAFER, DOPANT_SPECIES, MaterialInputError,
    MaterialSelection, default_material_draft, material_from_draft,
)
from app.solver_adapter import UIAnalysisConfig, physical_load_conversion, run_ui_analysis
from solver_v1.energy_model_registry import ENERGY_MODEL_IDS


@pytest.mark.parametrize('species', DOPANT_SPECIES)
def test_dopant_density_and_inactive_drafts(species):
    draft = dict(material_id=SILICON_WAFER, doping_enabled=True,
                 dopant_species=species, dopant_concentration_cm3=' 2.5E16 ')
    selected = material_from_draft(draft)
    assert selected == MaterialSelection(SILICON_WAFER, True, species, 2.5e16)
    assert selected.metadata()['dopant_concentration_basis'] == 'dopant_atoms_per_cm3'
    # Hidden entries are preserved in the draft but never become active Al data.
    draft['material_id'] = ALUMINUM
    draft['dopant_concentration_cm3'] = 'unfinished'
    assert material_from_draft(draft) == MaterialSelection()
    draft.update(material_id=SILICON_WAFER, doping_enabled=False)
    assert material_from_draft(draft) == MaterialSelection(SILICON_WAFER)
    assert material_from_draft(default_material_draft()) == MaterialSelection()


@pytest.mark.parametrize('value', ['', 'NaN', 'Inf', '-inf', '0', '-1', '1e999', '1e-999', '1e16 cm-3'])
def test_invalid_active_concentration_is_not_coerced(value):
    draft = dict(material_id=SILICON_WAFER, doping_enabled=True,
                 dopant_species='B', dopant_concentration_cm3=value)
    with pytest.raises(MaterialInputError, match='error.dopant_concentration'):
        material_from_draft(draft)


@pytest.mark.parametrize('energy_model', ENERGY_MODEL_IDS)
@pytest.mark.parametrize('doped', [False, True])
def test_silicon_cannot_enter_al_mechanics_pde_or_clock(energy_model, doped, monkeypatch):
    from app.axial_specimen import run_axial_probability
    from app.solid_mechanics import run_solid_probability
    from app.convergence_check import run_convergence_check

    def forbidden(*_a, **_k):
        pytest.fail('Si reached Al model construction, meshing, or PDE')

    monkeypatch.setattr('app.solver_adapter.build_time_basis_model', forbidden)
    monkeypatch.setattr('app.solver_adapter.run_probability_pde_2d', forbidden)
    monkeypatch.setattr('solver_v1.kinetic_calibration_workflow.build_time_basis_model', forbidden)
    monkeypatch.setattr('app.solid_mechanics.tetrahedralize', forbidden)
    config = UIAnalysisConfig(material_id=SILICON_WAFER, energy_model=energy_model,
        doping_enabled=doped, dopant_species='P' if doped else None,
        dopant_concentration_cm3=1e17 if doped else None)
    config.validate()  # A valid material configuration is distinct from an available backend.
    for run in (
        lambda: physical_load_conversion(config),
        lambda: run_ui_analysis(config),
        lambda: run_ui_analysis(config, _prepared_model=object()),
        lambda: run_axial_probability(config, None),
        lambda: run_solid_probability(config, None, [], None,
            poisson=.33, target_mm=3., direction=[1., 0., 0.]),
        lambda: run_convergence_check(config, reference_result={}),
        lambda: replace(config, time_basis='physical').validate_time_basis(),
    ):
        with pytest.raises(MaterialInputError, match='error.silicon_backend_unavailable'):
            run()


@pytest.mark.parametrize('selection', [
    MaterialSelection('unknown'), MaterialSelection(ALUMINUM, True, 'B', 1e16),
    MaterialSelection(SILICON_WAFER, False, 'P', 1e16),
    MaterialSelection(SILICON_WAFER, True, 'Ga', 1e16),
    MaterialSelection(SILICON_WAFER, True, 'B', True),
])
def test_invalid_material_combinations_are_rejected(selection):
    with pytest.raises(MaterialInputError):
        selection.validate()
