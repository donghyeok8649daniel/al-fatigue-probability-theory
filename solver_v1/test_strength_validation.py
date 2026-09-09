from dataclasses import replace
import numpy as np
import pytest

from .strength_validation import StrengthConditions, compare_strength, stress_at_plastic_strain


def conditions():
    return StrengthConditions('stress_at_plastic_strain','resolved_shear',.002,
        'specified_test_Al','hypothetical_test_microstructure','test_temperature','test_load_protocol')


def test_uniform_fold_cannot_be_reported_as_measured_yield():
    ref=conditions(); ideal=replace(ref,observable='uniform_ideal_fold',plastic_strain_level=None)
    r=compare_strength(ref,ideal,4.,3000.)
    assert not r['comparable'] and r['relative_error'] is None


def test_axial_and_resolved_shear_cannot_be_equated():
    ref=conditions()
    assert not compare_strength(ref,replace(ref,stress_component='axial'),4.,4.)['comparable']


def test_two_missing_microstructures_do_not_make_a_match():
    ref=replace(conditions(),microstructure_binding=None)
    assert not compare_strength(ref,ref,4.,4.)['comparable']


def test_matched_comparison_returns_actual_error():
    ref=conditions(); r=compare_strength(ref,ref,5.,6.)
    assert r['comparable'] and r['relative_error']==pytest.approx(.2)


def test_stress_criterion_is_not_first_nonzero_plastic_strain():
    r=stress_at_plastic_strain([0.,1.,5.,7.],[0.,1e-12,.001,.003],.002,maximum_bracket_width=.0021)
    assert r['resolved'] and r['stress_mpa']==pytest.approx(6.)


def test_insufficient_sampling_refuses_yield_interpolation():
    r=stress_at_plastic_strain([0.,10.],[0.,.02],.002,maximum_bracket_width=.001)
    assert not r['resolved'] and r['stress_mpa'] is None


def test_unreached_criterion_is_not_zero_strength():
    r=stress_at_plastic_strain([0.,10.],[0.,.0001],.002,maximum_bracket_width=.001)
    assert not r['resolved'] and r['stress_mpa'] is None


def test_nonmonotonic_hold_is_not_a_loading_branch():
    with pytest.raises(ValueError):
        stress_at_plastic_strain([0,10,0],[0,.002,0],.002,maximum_bracket_width=.001)


def test_nonfinite_resolution_is_rejected():
    with pytest.raises(ValueError):
        stress_at_plastic_strain([0,10],[0,.002],.001,maximum_bracket_width=np.nan)


def test_digitization_rejects_occlusion_and_ignores_isolated_artifact():
    from .run_strength_benchmark import trace_pixel_y
    rgb=np.full((40,5,3),255)
    rgb[20:25,:,0:2]=0; rgb[20:25,:,2]=120
    rgb[2,0,:]=[0,0,120]
    y, count, status=trace_pixel_y(rgb,offset=0)
    assert y==22 and count==25
    rgb[20:25]=[0,255,255]
    assert trace_pixel_y(rgb,offset=0)[0] is None
    rgb[20:25]=[0,0,120]; rgb[8:13]=[0,0,120]
    assert trace_pixel_y(rgb,offset=0)[0] is None
