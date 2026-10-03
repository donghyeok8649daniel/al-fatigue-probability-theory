"""Independent direct triples, conservative forces, bulk and J0..J4 checks."""
import numpy as np
import pytest
from results.silicon_wafer_feasibility.run_static_probe import source_parameters
from .silicon_environment_research import DiamondCell,Jet,_state,_finish
from .silicon_bessel_reference import SW111BesselInterface
from .silicon_sw_material_fit import force_features
from .silicon_sw_anharmonic import (beta_contract,angular_site,AngularDiamondCell,
                                   AngularSW111BesselInterface,MONOMIALS,MULTIPLICITIES)


def test_zero_beta_force_channels_and_angular_shape_bounds():
    p,_=source_parameters();x=np.array([[10.,10.,10.],[12.,10.,10.],[10.,12.,10.]])
    original=force_features(x,np.eye(3)*20,np.zeros(3,bool),p,strain_derivative=True)
    control=force_features(x,np.eye(3)*20,np.zeros(3,bool),p,strain_derivative=True,angle_beta=0.)
    for a,b in zip(original,control):np.testing.assert_array_equal(a,b)
    c=p['costheta0'];theta=np.linspace(-1-c,1-c,2001)
    for beta in [-.375,0.,.375]:
        beta_contract(beta,c)
        assert np.all(theta**2*(1-beta*theta)**2>=0)
        assert np.min(2*(1-beta*theta)*(1-2*beta*theta))>=-1e-15
    for beta in [np.inf,np.nan,.376,-.8]:
        with pytest.raises(ValueError):beta_contract(beta,c)


def test_rank_four_identity_keeps_tetrahedral_rank_three_and_self_subtraction():
    p,_=source_parameters();directions=np.array([[1.,1.,1.],[1.,-1.,-1.],[-1.,1.,-1.],[-1.,-1.,1.]])/np.sqrt(3)
    assert len(MONOMIALS)==35
    assert np.sum(MULTIPLICITIES[[sum(a)==3 for a in MONOMIALS]])==27
    # T_xyz is nonzero in diamond; bulk rank-three terms cannot be omitted.
    assert abs(np.sum(np.prod(directions,axis=1)))>.7
    vectors=2.35*directions+np.array([[.12,-.05,.03],[.03,.1,-.07],[-.1,.08,.02],[.05,-.11,.04]])
    coordinates=_state(vectors.ravel(),12,True);jets=[coordinates[i:i+3] for i in range(0,12,3)]
    for beta in [-.375,0.,.375]:
        d,m=[_finish(angular_site(jets,p,beta,r),12,True) for r in ['direct','moments']]
        np.testing.assert_allclose(m.value,d.value,atol=2e-12)
        np.testing.assert_allclose(m.gradient,d.gradient,atol=2e-11)
        np.testing.assert_allclose(m.hessian,d.hessian,atol=2e-10)
    one=[jets[0]]
    angular=angular_site(one,p,.375,'moments')-angular_site(one,dict(p,**{'lambda':0.}),.375,'direct')
    np.testing.assert_allclose(angular.value,0.,atol=2e-12)
    np.testing.assert_allclose(angular.hessian,0.,atol=1e-10)


@pytest.mark.parametrize('beta',[-.375,.375])
def test_angular_conservative_force_and_static_stress(beta):
    p,_=source_parameters();cell=np.eye(3)*20;x=np.array([[10.,10.,10.],[12.3,10.2,10.1],[10.1,12.1,10.3],[9.7,9.9,12.4]])
    pbc=np.zeros(3,bool);h=2e-5;e,f,s=force_features(x,cell,pbc,p,strain_derivative=True,angle_beta=beta)
    independent=sum(angular_site(np.delete(x,i,axis=0)-x[i],p,beta) for i in range(len(x)))
    np.testing.assert_allclose(e.sum(),independent,atol=3e-12)
    for atom,axis in [(0,0),(1,2),(2,1)]:
        d=np.zeros_like(x);d[atom,axis]=h
        ep,em=[force_features(q,cell,pbc,p,angle_beta=beta)[0] for q in [x+d,x-d]]
        np.testing.assert_allclose(-(ep-em)/(2*h),f[atom,axis],atol=3e-8,rtol=3e-7)
    mode=np.array([[1.,.2,-.1],[.2,-.3,.13],[-.1,.13,.2]])
    ep,em=[force_features(x@(np.eye(3)+sign*h*mode).T,cell@(np.eye(3)+sign*h*mode).T,pbc,p,angle_beta=beta)[0] for sign in [1,-1]]
    np.testing.assert_allclose((ep-em)/(2*h),np.einsum('ab,abk->k',mode,s),atol=3e-7,rtol=3e-7)
    np.testing.assert_allclose(f.sum(axis=0),0.,atol=3e-13)
    for k in range(3):np.testing.assert_allclose(np.cross(x-x.mean(axis=0),f[:,:,k]).sum(axis=0),0.,atol=3e-12)


def test_tetrahedral_full_nine_dimensional_hessian_and_cubic_eos_preserved():
    p,_=source_parameters();q=np.zeros(9)
    for lattice in [5.05,5.431,5.8]:
        original=DiamondCell(p,lattice).evaluate(q,'direct')
        for beta in [-.375,.375]:
            value=AngularDiamondCell(p,lattice,angle_beta=beta).evaluate(q)
            np.testing.assert_allclose(value.value,original.value,atol=2e-12)
            np.testing.assert_allclose(value.gradient,original.gradient,atol=3e-10)
            np.testing.assert_allclose(value.hessian,original.hessian,atol=3e-8)


def test_nonzero_angular_bessel_energy_gradient_full_hessian_and_sw_control():
    p,_=source_parameters();state=[.25,.15,-.07]
    model=AngularSW111BesselInterface(p,5.431,angle_beta=.375)
    b,d=[model.evaluate(state,method=m) for m in ['bessel','direct']]
    np.testing.assert_allclose(b.value,d.value,atol=2e-8)
    np.testing.assert_allclose(b.gradient,d.gradient,atol=2e-7)
    np.testing.assert_allclose(b.hessian,d.hessian,atol=3e-6)
    model.beta=0.
    control=model.evaluate(state);plain=SW111BesselInterface(p,5.431).evaluate(state)
    np.testing.assert_allclose(control.value,plain.value,atol=2e-10)
    np.testing.assert_allclose(control.gradient,plain.gradient,atol=2e-9)
    np.testing.assert_allclose(control.hessian,plain.hessian,atol=3e-8)
