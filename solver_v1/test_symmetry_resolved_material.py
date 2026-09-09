import numpy as np
import pytest

from solver_v1.symmetry_resolved_material import (
    SymmetryResolvedInterface, SymmetryTailBulkBasis, quadrupole_channel_gauge,
    SymmetryObservationCache,
)
from solver_v1.range_resolved_material import build_range_surface
from solver_v1.tail_constrained_material import normalized_amplitude
from solver_v1.angular_environment_reference import traceless_second
from solver_v1.yield_elastic_metric import cubic_to_mode_matrix


@pytest.fixture(scope='module')
def model():
    return SymmetryResolvedInterface((3.,6.,3.5),[.3,.5,2.,1.,.2,3.,-.1,2.,.03])


def test_symmetry_channel_is_independent_and_has_no_T2g_response():
    gauge=quadrupole_channel_gauge(6.,3.5)
    assert gauge['normal_response_residual']<1e-10
    cubic=np.linalg.solve(cubic_to_mode_matrix(),gauge['bulk_column'][2:5])
    assert cubic[0]>0 and cubic[1]<0
    assert abs(cubic[2])<1e-14
    assert abs(cubic[0]+2*cubic[1])<1e-13
    with pytest.raises(ValueError): quadrupole_channel_gauge(4.,4.)


def test_zero_extension_recovers_entire_verified_energy_jet(model):
    c=model.coefficients.copy();c[-1]=0
    new=SymmetryResolvedInterface(model.decays,c)
    old=build_range_surface(*model.decays,c[:8])
    q=(1.1*model.h,.17,-.09)
    n,o=new.evaluate(q),old.evaluate(q)
    assert n.energy==o.energy
    np.testing.assert_array_equal(n.gradient,o.gradient)
    np.testing.assert_array_equal(n.hessian,o.hessian)


def test_extra_energy_analytic_gradient_hessian_periodicity(model):
    for q in (np.array([model.h,0.,0.]),np.array([1.08*model.h,.18,-.11])):
        analytic=model.evaluate(q); step=2e-5
        for axis in range(3):
            delta=np.eye(3)[axis]*step
            plus,minus=model.evaluate(q+delta),model.evaluate(q-delta)
            assert (plus.energy-minus.energy)/(2*step)==pytest.approx(analytic.gradient[axis],rel=2e-6,abs=2e-7)
            np.testing.assert_allclose((plus.gradient-minus.gradient)/(2*step),analytic.hessian[:,axis],rtol=2e-6,atol=2e-6)
        np.testing.assert_allclose(analytic.hessian,analytic.hessian.T,atol=1e-14)
        shifted=model.evaluate(q+np.array([0.,1.,0.]))
        assert shifted.energy==pytest.approx(analytic.energy,abs=2e-10)
    at_bulk=model.evaluate([model.h,0.,0.])
    assert at_bulk.energy==pytest.approx(0.,abs=1e-12)


def test_combination_is_per_atom_after_neighbor_sum(model):
    q=np.array([1.08*model.h,.19,-.08]); face=model.face
    m,n=np.meshgrid(np.arange(-18,19),np.arange(-18,19),indexing='ij')
    R=m.ravel()[:,None]*model.geometry.a1+n.ravel()[:,None]*model.geometry.a2
    gauge=model.moment.gauge
    def plane(d,shift,k):
        v=np.column_stack([R+shift,np.full(len(R),d)])
        w=normalized_amplitude(k)*np.exp(-k*np.linalg.norm(v,axis=1))
        return traceless_second(np.einsum('n,ni,nj->ij',w,v,v))
    changes=[]
    for layer in range(1,23):
        shift=model.geometry.abc_shift(layer)
        parts=[]
        for k in model.decays[1:]:
            parts.append(plane(q[0]+(layer-1)*model.h,shift+q[1:],k)-plane(layer*model.h,shift,k))
        changes.append((parts[0]-gauge['eta']*parts[1])/gauge['normalization'])
    depths=np.cumsum(np.array(changes)[::-1],axis=0)[::-1]
    direct=2*np.sum(depths**2)
    assert model.extra_jet(q)[0][0]==pytest.approx(direct,rel=2e-10,abs=2e-11)
    # Squaring each plane before summing is a different, forbidden energy.
    assert abs(direct-2*np.sum(np.array(changes)**2))>1e-5


def test_extra_finite_q_hessian_matches_actual_sinusoidal_energy(model):
    basis=SymmetryTailBulkBasis(model.decays,radius=9.)
    q_cubic=np.array([.3,.17,.24]);q=basis.wavevector(q_cubic)
    polarization=np.array([1.,-2.,.4]);polarization/=np.linalg.norm(polarization)
    R=basis.R; phase=R@q;gauge=basis.gauge
    def energy(eps,theta):
        v=R+eps*(np.cos(phase+theta)-np.cos(theta))[:,None]*polarization
        rad=np.linalg.norm(v,axis=1);parts=[]
        for k in model.decays[1:]:
            w=normalized_amplitude(k)*np.exp(-k*rad)
            parts.append(traceless_second(np.einsum('n,ni,nj->ij',w,v,v)))
        Q=(parts[0]-gauge['eta']*parts[1])/gauge['normalization']
        return np.sum(Q**2)
    step=2e-5
    numeric=4*(np.mean([energy(step,t) for t in (0,np.pi/2,np.pi,3*np.pi/2)])-energy(0,0))/step**2
    columns,tail=basis.evaluate(q_cubic)
    assert numeric==pytest.approx(polarization@columns[-1]@polarization,rel=2e-6,abs=2e-7)
    fine=SymmetryTailBulkBasis(model.decays,radius=12.)
    cf,_=fine.evaluate(q_cubic)
    assert np.linalg.norm(cf[-1]-columns[-1],2)<=tail[-1]+1e-10


def test_cubic_embedding_isolated_atom_and_zero_reference_jet():
    from solver_v1.symmetry_resolved_material import CubicDensityBasis
    basis=CubicDensityBasis(2.)
    assert basis.value(0.)==-1.
    assert basis.value(2.)==basis.first_derivative(2.)==basis.second_derivative(2.)==0.
    for x in (.3,1.4):
        h=1e-5
        assert (basis.value(x+h)-basis.value(x-h))/(2*h)==pytest.approx(basis.first_derivative(x),rel=1e-8)
        assert (basis.first_derivative(x+h)-basis.first_derivative(x-h))/(2*h)==pytest.approx(basis.second_derivative(x),rel=1e-8)


def test_cubic_embedding_complete_surface_and_coefficient_jet(model):
    from solver_v1.symmetry_resolved_material import CubicSymmetryInterface,CubicSymmetryObservationCache
    from solver_v1.run_vector_material_calibration import source_and_targets
    _,obs,_=source_and_targets()
    c=np.r_[model.coefficients,.7];new=CubicSymmetryInterface(model.decays,c)
    q=np.array([1.25*model.h,.1,.04]);v=new.evaluate(q)
    for axis in range(3):
        h=np.eye(3)[axis]*2e-5
        plus,minus=new.evaluate(q+h),new.evaluate(q-h)
        assert (plus.energy-minus.energy)/4e-5==pytest.approx(v.gradient[axis],rel=3e-6,abs=2e-7)
        np.testing.assert_allclose((plus.gradient-minus.gradient)/4e-5,v.hessian[:,axis],rtol=3e-6,atol=2e-6)
    cache=CubicSymmetryObservationCache(obs);mat=cache.matrix(model.decays)
    assert mat[1,-1]==-1.  # isolated atom reference, NOT zero!
    np.testing.assert_array_equal(mat[[0,2,3,4],-1],np.zeros(4))
    for i,o in enumerate(obs):
        if o.bulk_index is None:
            v=new.evaluate(o.state)
            jet=np.r_[v.energy,v.gradient,v.hessian[0],v.hessian[1,1:],v.hessian[2,2]]
            assert mat[i]@c==pytest.approx(np.asarray(o.jet_weights)@jet,rel=1e-9,abs=1e-9)
