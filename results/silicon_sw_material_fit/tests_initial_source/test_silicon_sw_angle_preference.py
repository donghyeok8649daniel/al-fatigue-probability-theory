"""Nonzero bulk angular background: independent finite differences and sums."""
import numpy as np
from results.silicon_wafer_feasibility.run_static_probe import source_parameters
from .silicon_environment_research import _state,_finish
from .silicon_sw_anharmonic import angular_site,AngularDiamondCell,AngularSW111BesselInterface


def test_non_tetrahedral_preference_keeps_full_tensor_background():
    p,_=source_parameters();p=dict(p,costheta0=-.25)
    vectors=2.35*np.array([[1.,1.,1.],[1.,-1.,-1.],[-1.,1.,-1.],[-1.,-1.,1.]])/np.sqrt(3)
    q=_state(vectors.ravel(),12,True);jets=[q[i:i+3] for i in range(0,12,3)]
    direct,moments=[_finish(angular_site(jets,p,0.,r),12,True) for r in ['direct','moments']]
    np.testing.assert_allclose(moments.value,direct.value,atol=2e-12)
    np.testing.assert_allclose(moments.gradient,direct.gradient,atol=2e-11)
    np.testing.assert_allclose(moments.hessian,direct.hessian,atol=2e-10)
    # The angular energy no longer vanishes at exact tetrahedral directions.
    no_angle=_finish(angular_site(jets,dict(p,**{'lambda':0.}),0.,'direct'),12,True)
    assert direct.value-no_angle.value>1e-3
    assert np.max(abs(direct.hessian-no_angle.hessian))>1e-3


def test_preferred_cosine_bulk_hessian_by_energy_and_gradient_differences():
    p,_=source_parameters();model=AngularDiamondCell(dict(p,costheta0=-.25),5.461,angle_beta=0.)
    q=np.array([.005,-.002,.001,.003,.001,-.004,.008,-.007,.006]);jet=model.evaluate(q)
    step=1e-5
    for axis in [0,3,6]:
        delta=step*np.eye(9)[axis];plus,minus=[model.evaluate(state) for state in [q+delta,q-delta]]
        np.testing.assert_allclose((plus.value-minus.value)/(2*step),jet.gradient[axis],atol=2e-7)
        np.testing.assert_allclose((plus.gradient-minus.gradient)/(2*step),jet.hessian[:,axis],atol=2e-6)


def test_nonzero_background_bessel_interface_full_hessian():
    p,_=source_parameters();model=AngularSW111BesselInterface(dict(p,costheta0=-.25),5.461,
        angle_beta=0.,cut_kind='glide')
    b,d=[model.evaluate([.25,.15,-.07],method=method) for method in ['bessel','direct']]
    np.testing.assert_allclose(b.value,d.value,atol=2e-8)
    np.testing.assert_allclose(b.gradient,d.gradient,atol=2e-7)
    np.testing.assert_allclose(b.hessian,d.hessian,atol=3e-6)
