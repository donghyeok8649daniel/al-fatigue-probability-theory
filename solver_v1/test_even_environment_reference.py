"""Independent checks for the nested, non-production even environment term."""
import numpy as np
import pytest

from solver_v1.angular_environment_reference import AngularInterfaceInvariant, traceless_second
from solver_v1.even_moment_calibration import quadrupole_bulk_curvatures
from solver_v1.fcc111_active_interface import FCC111ActiveInterface
from solver_v1.joint_fcc_interface_calibration import IDEAL_H, JointFit, build_joint_model
from solver_v1.static_bulk_stability import neighbors_in_plane_frame


@pytest.fixture(scope="module")
def invariant():
    bulk=build_joint_model(JointFit("unfitted_test_only",3.,np.array([.1,.15,4.,2.,0.]),0.,None,None))
    bulk.a0=IDEAL_H
    return AngularInterfaceInvariant(FCC111ActiveInterface(bulk,tolerance=2e-11),rank=2,tolerance=2e-11)


def test_quadrupole_trace_rotation_and_even_parity():
    r=np.array([[.3,.2,.7],[-.1,.6,.2],[.8,-.3,.4]])
    def tensor(v):
        return traceless_second(v.T@v)
    q=tensor(r)
    np.testing.assert_allclose(np.trace(q),0.,atol=5e-16)
    np.testing.assert_array_equal(tensor(-r),q)
    rotation=np.array([[.6,-.8,0],[.8,.6,0],[0,0,1.]])
    assert np.sum(tensor(r@rotation.T)**2)==pytest.approx(np.sum(q*q),rel=2e-15)


def test_even_interface_direct_derivatives_and_perfect_limit(invariant):
    a,s=1.07*IDEAL_H,.23; step=1e-5
    v=invariant.evaluate(a,s)[0]
    assert v[0]==pytest.approx(invariant.direct_value(a,s,radius=24,layers=24),abs=2e-12)
    ap=invariant.evaluate(a+step,s)[0]; am=invariant.evaluate(a-step,s)[0]
    sp=invariant.evaluate(a,s+step)[0]; sm=invariant.evaluate(a,s-step)[0]
    np.testing.assert_allclose(v[1:3],[(ap[0]-am[0])/(2*step),(sp[0]-sm[0])/(2*step)],atol=1e-10)
    np.testing.assert_allclose([v[3],v[4]],(ap[1:3]-am[1:3])/(2*step),atol=3e-9)
    np.testing.assert_allclose([v[4],v[5]],(sp[1:3]-sm[1:3])/(2*step),atol=3e-9)
    # a+(k-1)h and kh differ by rounding in some layers. Do not require
    # sub-roundoff exact subtraction or insert a force-zero clamp in physics.
    np.testing.assert_allclose(invariant.evaluate(IDEAL_H,0.)[0][:3],0.,atol=32*np.finfo(float).eps)
    assert invariant.evaluate(IDEAL_H,invariant.interface.period)[0][0]<1e-25


def test_bulk_even_elasticity_matches_independent_direct_affine_derivative(invariant):
    r=neighbors_in_plane_frame(invariant.interface.bulk,14.)
    length=np.linalg.norm(r,axis=1)
    w=invariant.amplitude*np.exp(-invariant.kappa*length)
    raw=np.einsum("ni,nj->nij",r,r)
    base=traceless_second(raw)
    np.testing.assert_allclose(np.einsum("n,nij->ij",w,base),0.,atol=3e-12)
    actual=quadrupole_bulk_curvatures(invariant)
    expected=[]
    for mode in (np.eye(3),np.diag([0.,0.,1.]),np.outer([1.,0.,0.],[0.,0.,1.])):
        er=r@mode.T
        derivative=traceless_second(np.einsum("ni,nj->nij",er,r)+np.einsum("ni,nj->nij",r,er))
        derivative-=invariant.kappa*(np.sum(r*er,axis=1)/length)[:,None,None]*base
        qder=np.einsum("n,nij->ij",w,derivative)
        expected.append(2*np.sum(qder*qder))
    np.testing.assert_allclose(actual,expected,atol=2e-10,rtol=1e-9)


def test_exact_constraint_solver_does_not_certify_material_from_small_loss():
    from solver_v1.even_moment_calibration import even_basis
    from solver_v1.run_constrained_odd_calibration import constrained_matrix_fit
    matrix=even_basis(3.,5.)
    # Synthetic self-consistency data, explicitly NOT a physical Al fit.
    pair=np.linalg.solve(matrix[:2,:2],[0.,3.36])
    target=matrix[:,:2]@pair
    result=constrained_matrix_fit(matrix,target,np.ones(11),
        coefficient_order=["u","v","A","B","D3","D1","D2"],
        lower_rest=[0.,0.,0.,-np.inf,0.])
    np.testing.assert_allclose(result["predictions"],target,atol=1e-11)
    assert result["strictly_positive_LJ_resolved"]
    assert not result["physical_candidate_accepted"]


def test_even_interface_refuses_nonzero_unimplemented_bulk_background():
    model=build_joint_model(JointFit("noncubic_test",3.,np.array([.1,.15,4.,2.,0.]),0.,None,None))
    model.a0=1.01*IDEAL_H
    with pytest.raises(ValueError,match="cubic FCC bulk"):
        AngularInterfaceInvariant(FCC111ActiveInterface(model),rank=2)


def test_opening_inequality_columns_match_actual_analytic_forces():
    from solver_v1.angular_environment_reference import AngularInterfaceResearchSurface
    from solver_v1.monotone_opening_calibration import opening_force_basis
    coefficients=np.array([.1,.15,4.,2.,.3,2.,-.3,.1])
    bulk=build_joint_model(JointFit("force_basis_test",3.,coefficients[:5],0.,None,None))
    bulk.a0=IDEAL_H
    surface=AngularInterfaceResearchSurface(FCC111ActiveInterface(bulk,tolerance=2e-11),
        2.,angular_decay=3.,vector_amplitude_ev=-.3,quadrupole_amplitude_ev=.1,tolerance=2e-11)
    ratios=(1.1,2.7)
    np.testing.assert_allclose(opening_force_basis(3.,3.,ratios)@coefficients,
        [surface.packed(ratio*IDEAL_H,0.)[1] for ratio in ratios],atol=5e-11,rtol=2e-11)


def test_sampled_positive_traction_is_not_a_monotonicity_certificate():
    from solver_v1.monotone_opening_calibration import OPENING_RATIOS
    from solver_v1.run_monotone_opening_refinement import opening_force_extrema
    # Synthetic polynomial ONLY tests the stationary-force detector. This is
    # never an empirical constitutive energy or an Al calibration input.
    force=np.polynomial.Polynomial.fromroots([1.,3.3,3.7])
    energy=force.integ(); curvature=force.deriv()
    class TestSurface:
        h=1.
        def packed(self,a,s):
            return np.array([energy(a),force(a),0.,curvature(a),0.,1.])
    assert min(force(np.asarray(OPENING_RATIOS)))>0
    coarse=min(opening_force_extrema(TestSurface(),80),key=lambda r:r["force_ev_coordinate"])
    refined=min(opening_force_extrema(TestSurface(),160),key=lambda r:r["force_ev_coordinate"])
    assert coarse["force_ev_coordinate"]<-.09
    assert coarse["a_over_h"]==pytest.approx(refined["a_over_h"],abs=2e-11)
