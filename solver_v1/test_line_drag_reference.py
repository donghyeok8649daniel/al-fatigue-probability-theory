import numpy as np
import pytest

from .line_drag_reference import (
    load_line_drag_references,source_drag_pa_s,source_line_velocity_over_b,
    registry_profile_metric,conditional_registry_metric,
)


def test_published_reference_units_and_temperature_ranges():
    rows=load_line_drag_references()
    assert len(rows)==6
    assert source_drag_pa_s('olmsted2005_edge_015',300.)==pytest.approx(1.11e-5)
    assert source_drag_pa_s('olmsted2005_screw_030',300.)==pytest.approx(2.25e-5)
    for r in rows.values():
        if 'value_original' in r:
            # 1 dyn/cm² = .1 Pa, so retain .1 in the friction conversion.
            assert r['value_original']*.1==pytest.approx(r['value_SI'])
    with pytest.raises(ValueError):source_drag_pa_s('gorman1969_123K',300.)
    with pytest.raises(ValueError):source_drag_pa_s('olmsted2005_edge_015',50.)


def test_source_velocity_range_and_sign_are_not_a_solver_clock():
    result=source_line_velocity_over_b('olmsted2005_edge_015',shear_mpa=[-10.,0.,10.],temperature_K=300.)
    np.testing.assert_allclose(result,np.array([-10.,0.,10.])*1e6/1.11e-5)
    with pytest.raises(ValueError):
        source_line_velocity_over_b('olmsted2005_edge_015',shear_mpa=100.,temperature_K=300.)


def test_profile_metric_converges_to_analytic_arctangent_and_scales_with_length():
    # Synthetic mathematical test, NOT a calibrated core width.
    b=2.86e-10;width=3*b;exact=b*b/(2*np.pi*width)
    errors=[]
    for cells in (401,801,1601):
        x=np.linspace(-40*width,40*width,cells)
        s=b*(.5+np.arctan(x/width)/np.pi)
        measured=registry_profile_metric(x,s)
        errors.append(abs(measured-exact))
        assert registry_profile_metric(2*x,2*s)==pytest.approx(2*measured)
        assert registry_profile_metric(x,s+b)==pytest.approx(measured)
    assert errors[2]<errors[1]<errors[0]
    assert errors[-1]/exact<2e-4


def test_conditional_drag_projection_preserves_dissipation_not_global_time():
    B=2e-5;I=4e-11;area=7e-20
    d=conditional_registry_metric(line_drag_pa_s=B,profile_metric_m=I,atomic_cell_area_m2=area)
    eta=d['slip_friction_per_area_J_s_m4'];M=d['conditional_cell_mobility_m2_J_s']
    assert eta*I==pytest.approx(B)
    assert M*eta*area==pytest.approx(1.)
    assert d['production_calibration_available'] is False
    assert d['normal_mobility_available'] is False
    assert d['t0_seconds'] is None
    with pytest.raises(ValueError):
        conditional_registry_metric(line_drag_pa_s=B,profile_metric_m=0.,atomic_cell_area_m2=area)
