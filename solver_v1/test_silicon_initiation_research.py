import numpy as np
import pytest

from solver_v1.silicon_initiation_research import (
    intact_prism, prescribed_grips, grip_observables, connectivity_diagnostics)
from solver_v1.silicon_oxide_research import (
    axial_composite, rectangular_shell, inverse_equivalent_flaw_nm, precipitate_radius_nm)


def test_reference_prism_has_no_internal_seed():
    atoms,g=intact_prism()
    assert len(atoms)==360
    assert g['free'].sum()==216
    assert not atoms.pbc.any()
    checks=connectivity_diagnostics(atoms.positions,lower=g['lower'],upper=g['upper'])
    assert all(r['components']==1 and r['grip_connected'] for r in checks)
    assert all(r['pair_count']==592 for r in checks)
    assert np.isclose(np.linalg.det(g['directions_cubic']),2)


def test_grips_translate_without_changing_internal_distances():
    atoms,g=intact_prism()
    r=atoms.positions
    s=prescribed_grips(r,lower=g['lower'],upper=g['upper'],extension_A=3.7)
    for key in ('lower','upper'):
        ids=np.flatnonzero(g[key])
        assert np.allclose(s[ids]-s[ids[0]],r[ids]-r[ids[0]],atol=1e-14)
    assert np.array_equal(s[g['free']],r[g['free']])
    assert np.isclose((s[g['upper'],2].mean()-s[g['lower'],2].mean())-g['gauge_length_A'],3.7)


@pytest.mark.parametrize('delta',[-.1,0.,.4])
def test_grip_conjugacy_independent_two_spring_energy(delta):
    r=np.array([[0.,0.,-delta/2],[0.,0.,1.],[0.,0.,2+delta/2]])
    lower,upper=np.array([True,False,False]),np.array([False,False,True])
    k=7.; bonds=np.diff(r[:,2])-1
    f=np.zeros_like(r);f[0,2]=k*bonds[0];f[1,2]=k*(bonds[1]-bonds[0]);f[2,2]=-k*bonds[1]
    out=grip_observables(r,f,lower=lower,upper=upper,area_A2=2.)
    def energy(d):
        p=prescribed_grips(np.array([[0.,0.,0.],[0.,0.,1.],[0.,0.,2.]]),lower=lower,upper=upper,extension_A=d)
        return .5*k*np.sum((np.diff(p[:,2])-1)**2)
    h=1e-5
    assert out['conjugate_force_eV_A']==pytest.approx((energy(delta+h)-energy(delta-h))/(2*h),abs=1e-10)
    assert out['nominal_stress_GPa']==pytest.approx(k*delta/4*160.2176634,abs=1e-10)
    assert out['reaction_force_residual_eV_A']<1e-14
    assert out['total_internal_torque_eV']<1e-14


def test_connectivity_exposes_threshold_dependence():
    r=np.array([[0,0,0],[0,0,2.4],[0,0,5.4]])
    rows=connectivity_diagnostics(r,lower=[True,False,False],upper=[False,False,True],cutoffs_A=[2.8,3.2])
    assert not rows[0]['grip_connected']
    assert rows[1]['grip_connected']


def test_cross_section_and_composite_energy_force_consistency():
    area=rectangular_shell(outer_width=4.11,outer_height=5.11,thickness_per_face=.2)
    assert area==pytest.approx([17.4741,3.528])
    e=np.array([168.9,70]);e0=np.array([0.,.46/70])
    out=axial_composite(areas=area,young_moduli=e,eigenstrains=e0,nominal_stress=3.27)
    eps=out['strain'];h=1e-7
    energy=lambda x: .5*np.sum(area*e*(x-e0)**2)
    assert (energy(eps+h)-energy(eps-h))/(2*h)==pytest.approx(3.27*sum(area),rel=1e-10)
    assert abs(out['force_residual'])<2e-14
    assert out['effective_modulus']==pytest.approx(152.28646135386464)


def test_zero_external_force_balances_oxide_residual_tension():
    out=axial_composite(areas=[19.,1.],young_moduli=[169.,70.],eigenstrains=[0.,.46/70],nominal_stress=0.)
    assert out['phase_stress'][0]>0 and out['phase_stress'][1]<0
    assert abs(out['total_force'])<1e-15


def test_phase_subdivision_and_area_scaling_do_not_change_stress():
    a=axial_composite(areas=[19.,1.],young_moduli=[169.,70.],eigenstrains=[0.,.004],nominal_stress=3.)
    b=axial_composite(areas=[190.,4.,6.],young_moduli=[169.,70.,70.],eigenstrains=[0.,.004,.004],nominal_stress=3.)
    assert a['strain']==pytest.approx(b['strain'])
    assert a['effective_modulus']==pytest.approx(b['effective_modulus'])
    assert a['phase_stress']==pytest.approx(b['phase_stress'][:2])


def test_zero_oxide_has_pure_silicon_response():
    out=axial_composite(areas=[20.,0.],young_moduli=[168.9,70.],eigenstrains=[0.,.1],nominal_stress=4.09)
    assert out['effective_modulus']==pytest.approx(168.9)
    assert out['phase_stress'][0]==pytest.approx(4.09)


def test_paper_units_reconstruct_nanometres():
    assert inverse_equivalent_flaw_nm(4.09)==pytest.approx(19.02845,abs=1e-5)
    assert inverse_equivalent_flaw_nm(2.55)==pytest.approx(48.95192,abs=1e-5)
    assert precipitate_radius_nm(5.5)==pytest.approx(6.77594,abs=1e-5)
    assert precipitate_radius_nm(90)==pytest.approx(27.41002,abs=1e-5)
    assert precipitate_radius_nm(0)==0


@pytest.mark.parametrize('call',[
    lambda: intact_prism(repeats=(1,1,1)),
    lambda: rectangular_shell(outer_width=1,outer_height=1,thickness_per_face=.5),
    lambda: axial_composite(areas=[-1,2],young_moduli=[1,2],eigenstrains=[0,0],nominal_stress=0),
    lambda: inverse_equivalent_flaw_nm(0),
    lambda: precipitate_radius_nm(-1),
])
def test_invalid_physical_inputs_fail(call):
    with pytest.raises(ValueError):call()
