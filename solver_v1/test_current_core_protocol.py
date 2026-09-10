"""Static-core observations do not manufacture time, pinning or activation."""
import argparse

import numpy as np
import pytest

from .current_core_diagnostics import (
    phase_winding_cells, adjacent_registry_profile, inner_field_difference, translated_core_seed,
    twofold_core_partner,
    reference_configuration_seed,
)
from .current_material_rows import CurrentMaterialScrewCore
from .run_current_material_core import load_current_material, run
from .test_current_material_rows import model, reference_far


def test_phase_circulation_and_row_period_invariance(model):
    for sign in (-1, 0, 1):
        core = CurrentMaterialScrewCore(model, reference_far(model), free_radius=1.3, ring=1, burgers_sign=sign)
        q = core.initial
        first = phase_winding_cells(core, q)
        assert first['total_winding'] == pytest.approx(sign, abs=1e-12)
        shifted = q.copy(); shifted[0, 0] += core.rows.b
        second = phase_winding_cells(core, shifted)
        assert second['total_winding'] == pytest.approx(first['total_winding'], abs=1e-12)
        assert second['cells'] == first['cells']
        compare = inner_field_difference(core, q, shifted, radius=1.3)
        assert compare['maximum_vector_change_over_L0'] < 1e-14


def test_profile_keeps_real_values_and_affine_subtraction(model):
    core = CurrentMaterialScrewCore(model, reference_far(model), free_radius=1.3, ring=1,
                                   burgers_sign=0, shear_traction=.017)
    records = adjacent_registry_profile(core, core.initial)
    assert any(row['both_rows_free'] for row in records)
    assert max(abs(row['raw_delta_x']) for row in records) > 0
    assert max(abs(row['affine_subtracted_delta_x']) for row in records) < 1e-14
    assert all(row['raw_delta_y'] == row['raw_delta_z'] == 0 for row in records)


def test_glide_initialization_moves_defect_not_affine_or_boundary(model):
    far = reference_far(model)
    affine = CurrentMaterialScrewCore(model, far, free_radius=3., ring=1,
                                     burgers_sign=0, shear_traction=.01)
    translated = translated_core_seed(affine, affine.initial, 1)
    np.testing.assert_allclose(translated, affine.initial, atol=1e-15, rtol=1e-15)
    core = CurrentMaterialScrewCore(model, far, free_radius=3., ring=1)
    q = core.initial
    first = phase_winding_cells(core, q)
    translated = translated_core_seed(core, q, 1)
    second = phase_winding_cells(core, translated)
    assert len(first['cells']) == len(second['cells']) == 1
    assert second['cells'][0]['y_reference']-first['cells'][0]['y_reference'] == pytest.approx(core.rows.d)
    assert second['total_winding'] == pytest.approx(1., abs=1e-12)
    outside = np.ones(len(core.indices), bool); outside[core.free_ids] = False
    np.testing.assert_array_equal(core.full_field(translated)[outside], core.boundary[outside])
    with pytest.raises(ValueError):
        translated_core_seed(core, q, 10)


def test_exact_twofold_symmetry_of_complete_current_energy(model):
    core = CurrentMaterialScrewCore(model,reference_far(model),free_radius=1.3,ring=2)
    q = core.initial+.012*np.sin(np.arange(core.initial.size)*.7).reshape(core.initial.shape)
    partner = twofold_core_partner(core,q)
    np.testing.assert_allclose(twofold_core_partner(core,partner),q,atol=1e-15,rtol=1e-15)
    assert core.evaluate(q)['energy'] == pytest.approx(core.evaluate(partner)['energy'],abs=2e-12)


def test_reference_multistart_preserves_own_exterior_and_units(model):
    core=CurrentMaterialScrewCore(model,reference_far(model),free_radius=1.3,ring=1)
    q=core.initial+.001*np.sin(np.arange(core.initial.size)).reshape(core.initial.shape)
    saved={tuple(i):u for i,u in zip(core.indices[core.free_ids],q)}
    meta=dict(shear_traction_MPa=0.,burgers_sign=1,length_scale_m=2e-10,
              b_reduced=core.rows.b,h_reduced=core.rows.h,center_over_L0=core.far_field.center)
    before=core.boundary.copy()
    result=reference_configuration_seed(core,meta,saved,length_scale_m=2e-10)
    np.testing.assert_array_equal(result,q)
    np.testing.assert_array_equal(core.boundary,before)
    with pytest.raises(ValueError):reference_configuration_seed(core,meta,saved,length_scale_m=3e-10)
    with pytest.raises(ValueError):reference_configuration_seed(core,meta,{},length_scale_m=2e-10)


def test_reference_multistart_cannot_change_load_or_topology(model):
    core=CurrentMaterialScrewCore(model,reference_far(model),free_radius=1.3,ring=1)
    saved={tuple(i):u for i,u in zip(core.indices[core.free_ids],core.initial)}
    meta=dict(shear_traction_MPa=0.,burgers_sign=1,length_scale_m=2e-10,
              b_reduced=core.rows.b,h_reduced=core.rows.h,center_over_L0=core.far_field.center)
    for key,value in [('shear_traction_MPa',1.),('burgers_sign',-1),('h_reduced',1.)]:
        with pytest.raises(ValueError):
            reference_configuration_seed(core,dict(meta,**{key:value}),saved,length_scale_m=2e-10)


def test_loader_preserves_current_parameter_binding_and_uncalibrated_status(monkeypatch):
    from .run_current_material_core import ROOT, DEFAULT_MATERIAL
    monkeypatch.chdir(ROOT)
    _, tensor, metadata = load_current_material(DEFAULT_MATERIAL.relative_to(ROOT))
    assert tensor.shape == (3, 3, 3, 3)
    assert metadata['parameter_sha256'] == '8736cb9d28f1430e3991b1046ef42544cefaa3c9a4f5d743e0d02931329f5c2e'
    for key in ('energy_terms_omitted', 'material_accepted', 'actual_yield_calibrated', 'physical_time', 'physical_Hz'):
        assert metadata[key] is False


@pytest.mark.parametrize('key,value', [('radius',0),('max_seconds',float('nan')),
    ('tolerance',-1),('escape_amplitude',0),('iterations',0),('ring',0),('shear_mpa',float('inf'))])
def test_runner_invalid_controls_refused_before_creating_output(tmp_path, key, value):
    values = dict(out=tmp_path/'untouched', radius=2., max_seconds=10., tolerance=2e-12,
        force_tolerance=1e-6, escape_amplitude=.02, iterations=2, ring=2, shear_mpa=0., center_y=0.)
    values[key] = value
    with pytest.raises(ValueError):
        run(argparse.Namespace(**values))
    assert not values['out'].exists()
