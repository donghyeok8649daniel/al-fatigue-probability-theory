"""Unit, continuation, reference-derivative, and density-reduction checks."""
import math
import numpy as np
import pytest

from solver_v1.analytic_density_feasibility import (
    intrinsic_fault_upper_bound, mixture_basis, squared_envelope_basis,
)
from solver_v1.interface_static_scenarios import (
    InterfaceUnits, continue_static_branch, resolved_uniaxial_tractions,
)
from solver_v1.joint_fcc_interface_calibration import joint_basis, matched_targets
from solver_v1.reference_eam_targets import DEFAULT_CACHE, MishinRigidFCCReference


def test_traction_work_and_energy_units_round_trip():
    units=InterfaceUnits(2.8e-10,7.1e-20)
    load=np.array([-.9,.05,.15,1.,5.])
    np.testing.assert_allclose(units.force_to_traction(units.traction_to_force(load)),load,rtol=2e-16)
    force=units.traction_to_force(.15)
    assert force*.1==pytest.approx(.15e9*7.1e-20*2.8e-11/1.602176634e-19)
    assert units.energy_to_surface(1.)==pytest.approx(1.602176634e-19/7.1e-20)


def test_orientation_is_explicit_and_not_an_adjustable_chi():
    normal=[0,0,1]; slip=[1,0,0]
    np.testing.assert_allclose(resolved_uniaxial_tractions(1.,normal,normal,slip),[1,0])
    loading=np.array([1.,0,1])/math.sqrt(2)
    np.testing.assert_allclose(resolved_uniaxial_tractions(2.,loading,normal,slip),[1,1])
    with pytest.raises(ValueError):
        resolved_uniaxial_tractions(1.,[1,0,1],normal,slip)


def test_harmonic_static_load_and_unload_is_reversible():
    h=.8; stiffness=np.array([[3.,.2],[.2,2.]])
    def evaluate(a,s):
        q=np.array([a-h,s]); force=stiffness@q
        return np.array([.5*q@force,*force,3.,.2,2.])
    units=InterfaceUnits(1e-10,1e-20)
    loads=np.array([[0,0],[-1,0],[0,0],[1,.5],[0,0]])
    rows=continue_static_branch(evaluate,h=h,period=1.,units=units,loads_gpa=loads)
    for row,load in zip(rows,loads):
        assert row["status"]=="stable_local_equilibrium"
        expected=np.linalg.solve(stiffness,units.traction_to_force(load))
        np.testing.assert_allclose([row["a"]-h,row["s"]],expected,atol=1e-12)
    assert rows[-1]["s"]==pytest.approx(0,abs=1e-12)


def test_branch_failure_is_not_claimed_as_a_spinodal_or_plastic_jump():
    def unbound(a,s):
        return np.array([a,1.,0.,-1.,0.,1.])
    rows=continue_static_branch(unbound,h=1.,period=1.,units=InterfaceUnits(1e-10,1e-20),
                                loads_gpa=[[0,0],[1,1]])
    assert rows[0]["status"]!="stable_local_equilibrium"
    assert rows[0]["a"] is None
    assert rows[1]["status"]=="not_evaluated_after_branch_loss"


def test_positive_density_mixture_recovers_single_exponential_gauge():
    for weight in (0.,.4,1.):
        np.testing.assert_allclose(mixture_basis(2.8,2.8,weight),joint_basis(2.8),atol=1e-6,rtol=2e-7)
    np.testing.assert_allclose(mixture_basis(2.,4.,.3),mixture_basis(4.,2.,.7),atol=3e-14)


def test_squared_density_envelope_is_nonnegative_and_recovers_baseline():
    q=np.exp(-2.1*(np.linspace(.01,8,200)-1))
    shape=1.3
    np.testing.assert_allclose(q*(1-shape*q)**2,q-2*shape*q*q+shape*shape*q**3,atol=2e-13)
    assert np.all(q*(1-shape*q)**2>=0)
    np.testing.assert_allclose(squared_envelope_basis(2.1,0.),joint_basis(2.1),atol=1e-6,rtol=2e-7)


@pytest.mark.skipif(not DEFAULT_CACHE.exists(),reason="optional NIST source not downloaded")
def test_source_interface_derivatives_against_independent_finite_differences():
    reference=MishinRigidFCCReference(); a=1.04*reference.h; s=.19
    x=reference.interface_derivatives(a,s); step=2e-5
    ea=(reference.interface_energy(a+step,s)-reference.interface_energy(a-step,s))/(2*step)
    es=(reference.interface_energy(a,s+step)-reference.interface_energy(a,s-step))/(2*step)
    np.testing.assert_allclose(x[1:3],[ea,es],atol=2e-8,rtol=2e-7)
    ha=(reference.interface_derivatives(a+step,s)[1:3]-reference.interface_derivatives(a-step,s)[1:3])/(2*step)
    hs=(reference.interface_derivatives(a,s+step)[1:3]-reference.interface_derivatives(a,s-step)[1:3])/(2*step)
    np.testing.assert_allclose([x[3],x[4]],ha,atol=1e-7,rtol=2e-6)
    np.testing.assert_allclose([x[4],x[5]],hs,atol=1e-7,rtol=2e-6)


@pytest.mark.skipif(not DEFAULT_CACHE.exists(),reason="optional NIST source not downloaded")
def test_fixed_decay_observable_incompatibility_not_hidden_by_optimization():
    target,scales=matched_targets()
    bound=intrinsic_fault_upper_bound(2.8,target,scales,allowance=3.)
    assert bound["success"]
    assert bound["intrinsic_upper_ev_cell"]<0<target[5]
    assert bound["status"]==0
