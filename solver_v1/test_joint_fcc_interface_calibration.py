"""Matched-source, units, analytic extension, and scenario regressions."""
import math
import numpy as np
import pytest

from solver_v1.reference_eam_targets import DEFAULT_CACHE, MishinRigidFCCReference
from solver_v1.joint_fcc_interface_calibration import (
    IDEAL_H, INTERFACE_STATES, JointFit, MinimalConvexEmbedding,
    build_joint_model, fixed_decay_fit, joint_basis, matched_targets,
)
from solver_v1.fcc111_active_interface import FCC111ActiveInterface
from solver_v1.fcc111_geometry import DIRECT_110, SHOCKLEY_112, registry_path
from solver_v1.full_fcc_calibration_audit import relaxed_bulk_observables


@pytest.fixture(scope="module")
def reference():
    if not DEFAULT_CACHE.exists():
        pytest.skip("optional verified NIST reference file not downloaded; see source ledger")
    return MishinRigidFCCReference()


def test_external_reference_reproduces_published_cohesion(reference):
    assert reference.cohesion()==pytest.approx(3.36,abs=3e-7)
    assert reference._rho_bulk==pytest.approx(1.,abs=2e-5)
    assert reference.interface_energy(reference.h,0)==0.


@pytest.mark.parametrize("path",(DIRECT_110,SHOCKLEY_112))
def test_external_reference_periodic_registry_and_surface_normalization(reference,path):
    period=registry_path(path).period_over_b*reference.geometry.b
    left=reference.interface_energy(1.1*reference.h,.12,path_id=path)
    right=reference.interface_energy(1.1*reference.h,.12+period,path_id=path)
    assert left==pytest.approx(right,abs=1e-12)
    work=reference.work_separation()*16.02176634/reference.geometry.atomic_cell_area
    assert work==pytest.approx(1.7412853,abs=1e-6)
    # This is source data, NOT a probability-PDE-selectable energy interface.
    assert not hasattr(reference,"grad")


def test_source_interpolation_error_is_small_and_fault_sign_resolved(reference):
    linear=MishinRigidFCCReference(interpolation="linear")
    slip=reference.geometry.b/math.sqrt(3)
    cubic=reference.interface_energy(reference.h,slip,path_id=SHOCKLEY_112)
    alternative=linear.interface_energy(linear.h,slip,path_id=SHOCKLEY_112)
    assert abs(cubic-alternative)<1e-7
    assert cubic*16.02176634/reference.geometry.atomic_cell_area==pytest.approx(.15663043,abs=1e-7)


def test_quadratic_embedding_derivatives_and_atomized_reference():
    embedding=MinimalConvexEmbedding(3.,2.,.4,1.7)
    assert embedding.value(0)==pytest.approx(-1.6)
    rho=1.35; h=1e-5
    first=(embedding.value(rho+h)-embedding.value(rho-h))/(2*h)
    second=(embedding.first_derivative(rho+h)-embedding.first_derivative(rho-h))/(2*h)
    assert embedding.first_derivative(rho)==pytest.approx(first,abs=1e-9)
    assert embedding.second_derivative(rho)==pytest.approx(second,abs=1e-9)
    with pytest.raises(ValueError):
        MinimalConvexEmbedding(3.,2.,-.1,1.)


@pytest.mark.parametrize("quadratic",(0.,.3))
def test_joint_linear_basis_matches_actual_bulk_and_interface(quadratic):
    coefficients=np.array([.1,.15,4.,2.,quadratic])
    fit=JointFit("probe",3.,coefficients,0.,np.zeros(10),np.zeros(10))
    bulk=build_joint_model(fit); bulk.a0=IDEAL_H
    expected=joint_basis(3.)@coefficients
    np.testing.assert_allclose(expected[:5],relaxed_bulk_observables(bulk),atol=2e-7,rtol=1e-8)
    for i,(_,opening,slip,path) in enumerate(INTERFACE_STATES):
        face=FCC111ActiveInterface(bulk,path_id=path,tolerance=2e-11)
        value=face.separated_limit() if math.isinf(opening) else face.energy(opening*IDEAL_H,slip)
        assert expected[5+i]==pytest.approx(value,abs=2e-10)


def test_matched_targets_are_physical_cell_energies_not_unit_area_numbers(reference):
    targets,scales=matched_targets(reference)
    assert targets[5]==pytest.approx(.06943467959,abs=1e-10)
    assert len(targets)==10 and np.all(scales>0)
    assert targets[-1]==pytest.approx(reference.work_separation())
    assert targets[1]==3.36


def test_joint_coefficient_calibration_is_deterministic(reference):
    target,scales=matched_targets(reference)
    first=fixed_decay_fit(2.8,target,scales,family="linear")
    second=fixed_decay_fit(2.8,target,scales,family="linear")
    np.testing.assert_array_equal(first.coefficients,second.coefficients)
    assert first.squared_loss==pytest.approx(float(first.residuals@first.residuals))
    # Improvements in opening must not hide the failed intrinsic-fault target.
    assert first.predictions[-1]>0
    assert abs(first.residuals[5])>9
