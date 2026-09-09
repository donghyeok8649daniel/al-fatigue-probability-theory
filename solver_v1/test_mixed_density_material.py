import numpy as np
import pytest

from .mixed_density_material import (
    combine_density_jets, mixed_bulk_embedding_columns, MixedDensityQuarticInterface,
    MixedDensityObservationCache, MixedDensityTailBulkBasis,
)
from .quartic_angular_material import QuarticSymmetryInterface, QuarticSymmetryObservationCache
from .run_vector_material_calibration import source_and_targets
from .vector_interface_reference import _embedding_jet, HESSIAN_INDICES
from .joint_fcc_interface_calibration import MinimalConvexEmbedding


SHAPES=(3.1,4.8,6.2,.35)
COEFFICIENTS=np.array([.3,1.4,4.,3.,.2,1.,.2,.4,.25,100.])


def test_environment_combines_before_embedding_not_after():
    first=np.zeros((2,10));second=np.zeros((3,10))
    first[:,0]=[-.3,-.1];second[:,0]=[-.7,-.2,-.05]
    first[:,1]=[.2,.1];second[:,1]=[.3,.05,.01]
    w=.4;emb=MinimalConvexEmbedding(2.,0.,0.,1.)
    changes=combine_density_jets(first,second,w)
    value=_embedding_jet(changes,1.,emb)
    assert value[0]==pytest.approx(2*np.sum(emb.value(1+changes[:,0])-emb.value(1)))
    wrong=(1-w)*_embedding_jet(first,1.,emb)+w*_embedding_jet(second,1.,emb)
    assert abs(value[0]-wrong[0])>.01
    with pytest.raises(ValueError):combine_density_jets(first,second,1.1)


def test_zero_mixture_preserves_all_old_jets_exactly():
    old=QuarticSymmetryInterface(SHAPES[:3],COEFFICIENTS)
    new=MixedDensityQuarticInterface((*SHAPES[:3],0.),COEFFICIENTS)
    for state in ((old.h,0.,0.),(.92,.23,-.04),(1.4,.17,.11)):
        a,b=old.evaluate(state),new.evaluate(state)
        assert a.energy==b.energy
        np.testing.assert_array_equal(a.gradient,b.gradient)
        np.testing.assert_array_equal(a.hessian,b.hessian)


def test_unit_weight_and_coalescent_density_gauge():
    q=(.95,.23,.06)
    new=MixedDensityQuarticInterface((*SHAPES[:3],1.),COEFFICIENTS)
    old=QuarticSymmetryInterface((SHAPES[1],SHAPES[1],SHAPES[2]),COEFFICIENTS)
    a,b=new.evaluate(q),old.evaluate(q)
    assert a.energy==pytest.approx(b.energy,rel=1e-12,abs=1e-12)
    np.testing.assert_allclose(a.hessian,b.hessian,rtol=1e-12,atol=1e-12)
    a=MixedDensityQuarticInterface((3.1,3.1,6.2,.1),COEFFICIENTS).evaluate(q)
    b=MixedDensityQuarticInterface((3.1,3.1,6.2,.8),COEFFICIENTS).evaluate(q)
    assert a.energy==pytest.approx(b.energy,rel=1e-13)
    np.testing.assert_allclose(a.hessian,b.hessian,rtol=1e-13,atol=1e-12)


def test_analytic_mixed_energy_gradient_hessian_and_periodicity():
    model=MixedDensityQuarticInterface(SHAPES,COEFFICIENTS)
    for q in (np.array([.89,.23,-.05]),np.array([1.08,.31,.12])):
        value=model.evaluate(q);step=2e-5
        for axis in range(3):
            d=np.eye(3)[axis]*step
            left,right=model.evaluate(q-d),model.evaluate(q+d)
            assert value.gradient[axis]==pytest.approx((right.energy-left.energy)/(2*step),rel=3e-7,abs=2e-8)
            np.testing.assert_allclose(value.hessian[:,axis],(right.gradient-left.gradient)/(2*step),rtol=5e-7,atol=2e-6)
        shifted=model.evaluate(q+np.array([0.,1.,0.]))
        assert value.energy==pytest.approx(shifted.energy,rel=1e-12,abs=1e-12)
        np.testing.assert_array_equal(value.hessian,value.hessian.T)


def test_perfect_density_gauge_and_bulk_force_not_renormalized_with_state():
    model=MixedDensityQuarticInterface(SHAPES,COEFFICIENTS)
    changes,_=model.density_changes((model.h,0.,0.))
    np.testing.assert_allclose(changes[:,0],0.,atol=1e-14)
    bulk=mixed_bulk_embedding_columns(SHAPES[0],SHAPES[1],SHAPES[3])
    emb=_embedding_jet(changes,1.,model.normalized_embedding)
    assert emb[0]==pytest.approx(0.,abs=1e-13)
    assert model.h*emb[1]==pytest.approx(bulk[0]@COEFFICIENTS[2:5],rel=2e-10,abs=2e-11)
    opened,_=model.density_changes((1.5,0.,0.))
    assert np.all(opened[:,0]<=0.)
    assert np.min(opened[:,0])<-.1


def test_hydrostatic_refinement_and_atomization_reference():
    a=mixed_bulk_embedding_columns(SHAPES[0],SHAPES[1],SHAPES[3],strain_step=8e-4)
    b=mixed_bulk_embedding_columns(SHAPES[0],SHAPES[1],SHAPES[3],strain_step=4e-4)
    np.testing.assert_allclose(a,b,rtol=1e-8,atol=3e-8)
    np.testing.assert_array_equal(a[1],[1.,-1.,1.])


def test_observation_matrix_matches_actual_energy_not_mixture_surrogate():
    _,obs,_=source_and_targets()
    new=MixedDensityObservationCache(obs)
    raw=new.matrix(SHAPES)
    model=MixedDensityQuarticInterface(SHAPES,COEFFICIENTS)
    for index in (5,8,15):
        o=obs[index];v=model.evaluate(o.state)
        jet=np.r_[v.energy,v.gradient,v.hessian[0],v.hessian[1,1:],v.hessian[2,2]]
        assert raw[index]@COEFFICIENTS==pytest.approx(np.asarray(o.jet_weights)@jet,rel=5e-10,abs=2e-9)
    old=QuarticSymmetryObservationCache(obs).matrix(SHAPES[:3])
    np.testing.assert_array_equal(new.matrix((*SHAPES[:3],0.)),old)


def test_mixed_scalar_tail_bounds_independent_radius_change_and_acoustic_zero():
    coarse=MixedDensityTailBulkBasis(SHAPES,radius=7.)
    fine=MixedDensityTailBulkBasis(SHAPES,radius=11.)
    for q in ((.01,.01,.0),(.375,.375,0.),(.7,.2,.1)):
        columns,errors=coarse.evaluate(q);ref,_=fine.evaluate(q)
        difference=np.linalg.norm(columns-ref,ord=2,axis=(-2,-1))
        assert np.all(difference<=errors+2e-12)
    columns,errors=coarse.evaluate((0.,0.,0.))
    np.testing.assert_array_equal(columns,np.zeros_like(columns))
    np.testing.assert_array_equal(errors,np.zeros_like(errors))
