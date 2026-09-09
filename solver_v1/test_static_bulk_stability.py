"""Independent second-variation tests; no inertial mass or physical time."""
import numpy as np
import pytest

from solver_v1.angular_environment_reference import traceless_third, traceless_second
from solver_v1.joint_fcc_interface_calibration import JointFit, build_joint_model, IDEAL_H
from solver_v1.static_bulk_stability import StaticBulkHessian


@pytest.fixture(scope="module")
def bulk():
    model=build_joint_model(JointFit("unfitted_test",3.,np.array([.1,.15,4.,2.,0.]),0.,None,None))
    model.a0=IDEAL_H
    return model


def test_translational_invariance_parity_and_positive_I3_increment(bulk):
    base=StaticBulkHessian(bulk,cutoff=6.)
    extended=StaticBulkHessian(bulk,cutoff=6.,D3=2.)
    np.testing.assert_array_equal(base.evaluate([0,0,0])["matrix"],np.zeros((3,3)))
    q=base.crystallographic_wavevector([.3,.4,.2])
    v=extended.evaluate(q)
    np.testing.assert_allclose(v["matrix"],v["matrix"].T,atol=1e-14)
    np.testing.assert_allclose(v["matrix"],extended.evaluate(-q)["matrix"],atol=1e-14)
    assert np.linalg.eigvalsh(v["matrix"]-base.evaluate(q)["matrix"])[0]>=-1e-14


def test_cubic_wavevector_labels_match_positive_abc_crystal(bulk):
    evaluator=StaticBulkHessian(bulk,cutoff=3.)
    g=bulk.geometry
    fraction=np.array([.2,.4,.1])
    q=evaluator.crystallographic_wavevector(fraction)
    qcubic=fraction*2*np.pi/g.lattice_constant
    np.testing.assert_allclose(evaluator.R@q,
        (evaluator.R@g.plane_basis_in_stacked_cubic_axes())@qcubic,atol=1e-15)


@pytest.mark.parametrize("D2",[0.,.3])
def test_force_constant_matrix_matches_direct_sinusoidal_energy_variation(bulk,D2):
    evaluator=StaticBulkHessian(bulk,cutoff=10.,D3=2.,D1=-.15,D2=D2)
    q=evaluator.crystallographic_wavevector([.23,.37,.17])
    polarization=np.array([1.,2.,-.5]); polarization/=np.linalg.norm(polarization)
    R=evaluator.R; phase=R@q; p=bulk.p; density=bulk.density_params
    def energy(displacement,angle):
        shifted=R+displacement*(np.cos(phase+angle)-np.cos(angle))[:,None]*polarization
        r=np.linalg.norm(shifted,axis=1)
        pair=2*p.epsilon*np.sum((p.sigma_lj/r)**12-(p.sigma_lj/r)**6)
        rho=density.C_rho*np.exp(-density.kappa*r)
        scalar=float(bulk.embedding.value(np.sum(rho)))
        w=rho/bulk.embedding.rho_ref
        vector=np.einsum("n,ni->i",w,shifted)
        tensor=traceless_third(np.einsum("n,ni,nj,nk->ijk",w,shifted,shifted,shifted))
        quadrupole=traceless_second(np.einsum("n,ni,nj->ij",w,shifted,shifted))
        return pair+scalar+2*np.sum(tensor*tensor)-.15*np.sum(vector*vector)+D2*np.sum(quadrupole*quadrupole)
    step=1e-4
    baseline=energy(0,0)
    average=np.mean([energy(step,angle) for angle in (0,np.pi/2,np.pi,3*np.pi/2)])
    measured=4*(average-baseline)/step**2
    exact=polarization@evaluator.evaluate(q)["matrix"]@polarization
    assert measured==pytest.approx(exact,rel=4e-6,abs=2e-5)
