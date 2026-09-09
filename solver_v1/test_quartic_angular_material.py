import numpy as np
import pytest

from solver_v1.quartic_angular_material import QuarticSymmetryInterface, per_site_quartic_jet, QuarticSymmetryObservationCache
from solver_v1.symmetry_resolved_material import SymmetryResolvedInterface


def test_quartic_is_per_site_not_squared_total_and_rotation_invariant():
    values=np.zeros((2,10,3));values[0,0]=[1.,0.,0.];values[1,0]=[0.,2.,0.]
    assert per_site_quartic_jet(values)[0]==17.
    assert per_site_quartic_jet(values)[0]!=25.
    Q=np.array([[0.,-1.,0.],[1.,0.,0.],[0.,0.,1.]])
    np.testing.assert_allclose(per_site_quartic_jet(values@Q),per_site_quartic_jet(values))


def test_quartic_zero_limit_and_pristine_hessian_are_exact():
    c=np.array([.3,.5,2.,1.,.2,3.,-.1,2.,.03,1200.]);decays=(3.,6.,3.5)
    model=QuarticSymmetryInterface(decays,c);base=SymmetryResolvedInterface(decays,c[:9])
    q=(model.h,0.,0.)
    a,b=model.evaluate(q),base.evaluate(q)
    np.testing.assert_array_equal(a.hessian,b.hessian)
    np.testing.assert_array_equal(a.gradient,b.gradient)
    c[-1]=0.;zero=QuarticSymmetryInterface(decays,c)
    q=(1.12*model.h,.13,.09)
    assert zero.evaluate(q).energy==base.evaluate(q).energy


def test_quartic_surface_analytic_force_and_hessian_and_coefficient_basis():
    from solver_v1.run_vector_material_calibration import source_and_targets
    c=np.array([.3,.5,2.,1.,.2,3.,-.1,2.,.03,1200.]);decays=(3.,6.,3.5)
    model=QuarticSymmetryInterface(decays,c)
    q=np.array([1.12*model.h,.13,.09]);value=model.evaluate(q)
    for axis in range(3):
        h=np.eye(3)[axis]*2e-5;plus=model.evaluate(q+h);minus=model.evaluate(q-h)
        assert (plus.energy-minus.energy)/4e-5==pytest.approx(value.gradient[axis],rel=3e-6,abs=3e-7)
        np.testing.assert_allclose((plus.gradient-minus.gradient)/4e-5,value.hessian[:,axis],rtol=3e-6,atol=3e-6)
    _,obs,_=source_and_targets();cache=QuarticSymmetryObservationCache(obs)
    mat=cache.matrix(decays)
    np.testing.assert_array_equal(mat[:5,-1],np.zeros(5))
    # Source and candidate compute h through different floating-point paths;
    # this fourth-order term's Hessian is quadratic in that tiny displacement.
    np.testing.assert_allclose(mat[5:7,-1],np.zeros(2),atol=np.finfo(float).eps**2,rtol=0.)
    for i,o in enumerate(obs):
        if o.bulk_index is None:
            v=model.evaluate(o.state)
            jet=np.r_[v.energy,v.gradient,v.hessian[0],v.hessian[1,1:],v.hessian[2,2]]
            assert mat[i]@c==pytest.approx(np.asarray(o.jet_weights)@jet,rel=2e-9,abs=2e-9)


def test_rational_extension_analytic_jet_and_quartic_limit():
    c=np.array([.3,.5,2.,1.,.2,3.,-.1,2.,.03,1200.]);decays=(3.,6.,3.5)
    model=QuarticSymmetryInterface(decays,c,saturation=5000.)
    q=np.array([1.15*model.h,.16,-.03]);v=model.evaluate(q)
    for axis in range(3):
        h=np.eye(3)[axis]*1e-5;p,m=model.evaluate(q+h),model.evaluate(q-h)
        assert (p.energy-m.energy)/2e-5==pytest.approx(v.gradient[axis],rel=1e-6,abs=1e-7)
        np.testing.assert_allclose((p.gradient-m.gradient)/2e-5,v.hessian[:,axis],rtol=1e-6,atol=1e-6)
    unbounded=QuarticSymmetryInterface(decays,c)
    assert 0<model.quartic.jet(tuple(q))[0][0]<unbounded.quartic.jet(tuple(q))[0][0]
    values=np.zeros((2,10,3));values[:,0]=[[1.,2.,3.],[-.4,.2,.1]]
    np.testing.assert_allclose(per_site_quartic_jet(values,1e-10),per_site_quartic_jet(values),rtol=2e-9)
