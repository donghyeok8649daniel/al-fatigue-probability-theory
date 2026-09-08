"""Analytic, symmetry, and independent-direct checks for research I3 only."""
import math
import numpy as np
import pytest

from solver_v1.angular_environment_reference import (
    AngularInterfaceInvariant, AngularInterfaceResearchSurface, traceless_third,
)
from solver_v1.fcc111_active_interface import FCC111ActiveInterface
from solver_v1.fcc111_geometry import SHOCKLEY_112
from solver_v1.interface_static_scenarios import packed_evaluation
from solver_v1.joint_fcc_interface_calibration import IDEAL_H, JointFit, build_joint_model


@pytest.fixture(scope="module")
def invariant():
    probe=JointFit("unfitted_test_only",3.,np.array([.1,.15,4.,2.,0.]),0.,None,None)
    bulk=build_joint_model(probe); bulk.a0=IDEAL_H
    face=FCC111ActiveInterface(bulk,path_id=SHOCKLEY_112,tolerance=2e-11)
    return AngularInterfaceInvariant(face,tolerance=2e-11)


def test_traceless_third_moment_rotation_and_inversion():
    vectors=np.array([[.3,.2,.7],[-.1,.6,.2],[.8,-.3,.4]])
    weights=np.array([.3,.5,.2])
    def tensor(r):
        return traceless_third(np.einsum("n,ni,nj,nk->ijk",weights,r,r,r))
    original=tensor(vectors)
    np.testing.assert_allclose(np.einsum("ijj->i",original),0.,atol=1e-16)
    np.testing.assert_allclose(tensor(-vectors),-original,atol=0.)
    angle=.71
    rotation=np.array([[math.cos(angle),-math.sin(angle),0],
                       [math.sin(angle),math.cos(angle),0],[0,0,1]])
    assert np.sum(tensor(vectors@rotation.T)**2)==pytest.approx(np.sum(original**2),rel=3e-15)
    centrosymmetric=np.r_[vectors,-vectors]
    tensor0=np.einsum("ni,nj,nk->ijk",centrosymmetric,centrosymmetric,centrosymmetric)
    np.testing.assert_allclose(tensor0,0,atol=1e-16)


def test_angular_reciprocal_matches_independent_direct_and_refinement(invariant):
    a,s=1.03*IDEAL_H,.19
    reciprocal=invariant.evaluate(a,s)[0][0]
    coarse=invariant.direct_value(a,s,radius=10,layers=10)
    fine=invariant.direct_value(a,s,radius=20,layers=20)
    assert abs(fine-reciprocal)<abs(coarse-reciprocal)
    assert abs(fine-reciprocal)<2e-12


def test_angular_gradient_and_hessian_are_analytic(invariant):
    a,s=1.03*IDEAL_H,.19; step=1e-5
    v=invariant.evaluate(a,s)[0]
    ap=invariant.evaluate(a+step,s)[0]; am=invariant.evaluate(a-step,s)[0]
    sp=invariant.evaluate(a,s+step)[0]; sm=invariant.evaluate(a,s-step)[0]
    np.testing.assert_allclose(v[1:3],[(ap[0]-am[0])/(2*step),(sp[0]-sm[0])/(2*step)],atol=1e-10)
    np.testing.assert_allclose([v[3],v[4]],(ap[1:3]-am[1:3])/(2*step),atol=3e-9)
    np.testing.assert_allclose([v[4],v[5]],(sp[1:3]-sm[1:3])/(2*step),atol=3e-9)


def test_perfect_bulk_limit_and_partial_registry_are_distinct(invariant):
    pristine=invariant.evaluate(IDEAL_H,0)[0]
    assert pristine[0]<1e-24
    np.testing.assert_allclose(pristine[1:3],0.,atol=1e-15)
    fault=invariant.evaluate(IDEAL_H,1/math.sqrt(3))[0][0]
    assert fault>1e-5
    assert invariant.evaluate(IDEAL_H,invariant.interface.period)[0][0]<1e-24


def test_zero_angular_amplitude_preserves_original_energy(invariant):
    scalar=invariant.interface
    extended=AngularInterfaceResearchSurface(scalar,0.)
    np.testing.assert_array_equal(extended.packed(1.03*IDEAL_H,.19),
                                  packed_evaluation(scalar,1.03*IDEAL_H,.19))


def test_first_moment_direct_and_analytic_derivatives(invariant):
    vector=AngularInterfaceInvariant(invariant.interface,rank=1,tolerance=2e-11)
    a,s=1.03*IDEAL_H,.19; step=1e-5
    value=vector.evaluate(a,s)[0]
    assert value[0]==pytest.approx(vector.direct_value(a,s,radius=24,layers=24),abs=2e-12)
    ap=vector.evaluate(a+step,s)[0]; am=vector.evaluate(a-step,s)[0]
    sp=vector.evaluate(a,s+step)[0]; sm=vector.evaluate(a,s-step)[0]
    np.testing.assert_allclose(value[1:3],[(ap[0]-am[0])/(2*step),(sp[0]-sm[0])/(2*step)],atol=1e-10)
    np.testing.assert_allclose([value[3],value[4]],(ap[1:3]-am[1:3])/(2*step),atol=2e-9)
    np.testing.assert_allclose([value[4],value[5]],(sp[1:3]-sm[1:3])/(2*step),atol=2e-9)


def test_signed_vector_term_is_separate_and_no_bulk_energy_is_added(invariant):
    scalar=invariant.interface; a,s=1.03*IDEAL_H,.19
    extended=AngularInterfaceResearchSurface(scalar,3.,vector_amplitude_ev=-.2)
    expected=(packed_evaluation(scalar,a,s)+3*extended.angular.evaluate(a,s)[0]
              -.2*extended.vector.evaluate(a,s)[0])
    np.testing.assert_array_equal(extended.packed(a,s),expected)
    assert abs(extended.packed(IDEAL_H,0)[0])<1e-20
    # Whether a signed coefficient yields physical TOTAL stability must still
    # be checked. No empirical plasticity or positivity clamp is inserted.


def test_local_interface_curvature_basis_matches_the_scalar_surface(invariant):
    from solver_v1.odd_moment_calibration import scalar_curvature_row
    scalar=invariant.interface
    predicted=scalar_curvature_row(3.)@np.array([.1,.15,4.,2.])
    assert predicted==pytest.approx(scalar.hessian(IDEAL_H,0)[0,0],abs=2e-10)
