"""New density product forces and preservation of zero-angular bulk anchors."""
import numpy as np
import pytest
from results.silicon_wafer_feasibility.run_static_probe import source_parameters
from .silicon_sw_material_fit import force_features
from .silicon_sw_density_screening import (cubic_density_reference,screened_site,ScreenedDiamondCell,
                                          ScreenedSW111BesselInterface)
from .silicon_environment_research import DiamondCell,environment_jet


def test_zero_exponent_recovers_original_exact_force_channels():
    p,_=source_parameters();cell=np.eye(3)*20.;x=np.array([[10.,10.,10.],[12.,10.,10.],[10.,12.,10.]])
    plain=force_features(x,cell,np.zeros(3,bool),p,strain_derivative=True)
    screened=force_features(x,cell,np.zeros(3,bool),p,strain_derivative=True,density_power=0.,density_reference=.8)
    for a,b in zip(plain,screened):np.testing.assert_array_equal(a,b)


@pytest.mark.parametrize('density_gamma',[None,.2])
def test_density_gradient_in_conservative_force_and_stress(density_gamma):
    p,_=source_parameters();cell=np.eye(3)*20.;x=np.array([[10.,10.,10.],[12.3,10.2,10.1],[10.1,12.1,10.3],[9.7,9.9,12.4]])
    pbc=np.zeros(3,bool);reference=cubic_density_reference(p,5.431,density_gamma=density_gamma);power=2.;h=2e-5
    options=dict(density_power=power,density_reference=reference,density_gamma=density_gamma)
    e,f,s=force_features(x,cell,pbc,p,strain_derivative=True,**options)
    direct=sum(screened_site(np.delete(x,i,axis=0)-x[i],p,power,reference,density_gamma=density_gamma) for i in range(len(x)))
    np.testing.assert_allclose(e.sum(),direct,atol=2e-12)
    for atom,axis in [(0,0),(1,2),(2,1)]:
        d=np.zeros_like(x);d[atom,axis]=h
        ep,em=[force_features(q,cell,pbc,p,**options)[0] for q in [x+d,x-d]]
        np.testing.assert_allclose(-(ep-em)/(2*h),f[atom,axis],atol=2e-8,rtol=2e-7)
    mode=np.diag([1.,-.3,.2])
    ep,em=[force_features(x@(np.eye(3)+sign*h*mode).T,cell@(np.eye(3)+sign*h*mode).T,pbc,p,
              **options)[0] for sign in [1,-1]]
    np.testing.assert_allclose((ep-em)/(2*h),np.einsum('ab,abk->k',mode,s),atol=3e-7,rtol=3e-7)
    np.testing.assert_allclose(f.sum(axis=0),0.,atol=2e-13)
    for k in range(3):np.testing.assert_allclose(np.cross(x-x.mean(axis=0),f[:,:,k]).sum(axis=0),0.,atol=2e-12)


def test_perfect_tetrahedral_bulk_full_hessian_is_preserved():
    p,_=source_parameters();lattice=5.431;reference=cubic_density_reference(p,lattice)
    base=DiamondCell(p,lattice).evaluate(np.zeros(9),'direct')
    for power in [0.,1.,2.,4.]:
        value=ScreenedDiamondCell(p,lattice,density_power=power,density_reference=reference).evaluate(np.zeros(9))
        np.testing.assert_allclose(value.value,base.value,atol=2e-12)
        np.testing.assert_allclose(value.gradient,base.gradient,atol=2e-10)
        np.testing.assert_allclose(value.hessian,base.hessian,atol=3e-9)


def test_cutoff_underflow_and_nonzero_screened_bessel_derivatives():
    p,_=source_parameters();rc=p['sigma']*p['a'];cell=np.eye(3)*30
    x=np.array([[10.,10.,10.],[10.+rc-1e-6,10.,10.],[10.,10.+rc-1e-6,10.]])
    e,f,s=force_features(x,cell,np.zeros(3,bool),p,strain_derivative=True,
                         density_power=2.,density_reference=cubic_density_reference(p,5.431))
    assert np.all(np.isfinite(e)) and np.all(np.isfinite(f)) and np.all(np.isfinite(s))
    model=ScreenedSW111BesselInterface(p,5.431,density_power=2.,
                                     density_reference=cubic_density_reference(p,5.431))
    b,d=[model.evaluate([.25,.15,-.07],method=m) for m in ['bessel','direct']]
    np.testing.assert_allclose(b.value,d.value,atol=2e-8)
    np.testing.assert_allclose(b.gradient,d.gradient,atol=2e-7)
    np.testing.assert_allclose(b.hessian,d.hessian,atol=3e-6)


def test_separate_coordination_density_kernel_has_same_bessel_derivatives():
    p,_=source_parameters();reference=cubic_density_reference(p,5.431,density_gamma=.2)
    model=ScreenedSW111BesselInterface(p,5.431,density_power=2.,density_reference=reference,density_gamma=.2)
    b,d=[model.evaluate([.25,.15,-.07],method=m) for m in ['bessel','direct']]
    np.testing.assert_allclose(b.value,d.value,atol=2e-8)
    np.testing.assert_allclose(b.gradient,d.gradient,atol=2e-7)
    np.testing.assert_allclose(b.hessian,d.hessian,atol=3e-6)
