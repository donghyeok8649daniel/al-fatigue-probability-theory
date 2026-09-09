"""Exact coefficient/target bookkeeping; tests do not certify material strength."""
from dataclasses import replace
import json

import numpy as np
import pytest

from .run_low_stress_cyclic_diagnostic import FIT, build_surface
from .reference_eam_targets import MishinRigidFCCReference
from .vector_interface_reference import FullRegistryInterface, MishinVectorInterfaceReference
from .vector_material_calibration import (
    IDEAL_H, LENGTH_M, UNITS, VectorCoefficientBasis, matched_observations,
    observation_matrix, exact_constraint_tangent, fit_coefficients,
    build_coefficient_surface, coefficient_identifiability,
)


@pytest.fixture(scope='module')
def setup():
    old = json.loads(FIT.read_bytes())['angular_monotone_opening']
    c = np.array(old['coefficients'])
    basis = VectorCoefficientBasis(old['scalar_decay'],old['angular_decay'])
    source = MishinVectorInterfaceReference(MishinRigidFCCReference(),LENGTH_M/1e-10)
    obs, states = matched_observations(source)
    matrix = observation_matrix(basis,obs)
    return basis,c,source,obs,states,matrix


@pytest.mark.parametrize('q', [(IDEAL_H,0.,0.),(1.08*IDEAL_H,.237,-.113), (1.7*IDEAL_H,.3,.19)])
def test_linear_basis_recovers_unmodified_candidate(setup,q):
    basis,c,*_ = setup
    direct = FullRegistryInterface(build_surface(tolerance=2e-11)[0]).evaluate(q)
    value = basis.evaluate(q,c)
    np.testing.assert_allclose(value.energy,direct.energy,rtol=0,atol=2e-11)
    np.testing.assert_allclose(value.gradient,direct.gradient,rtol=0,atol=3e-10)
    np.testing.assert_allclose(value.hessian,direct.hessian,rtol=0,atol=6e-10)


def test_coefficient_derivatives_are_exact(setup):
    basis,c,*_=setup; q=(1.05*IDEAL_H,.21,.121)
    matrix=basis.jet(q)
    for i in range(8):
        step=1e-4*max(1,abs(c[i])); p=c.copy(); m=c.copy(); p[i]+=step; m[i]-=step
        plus=basis.evaluate(q,p); minus=basis.evaluate(q,m)
        np.testing.assert_allclose((plus.gradient-minus.gradient)/(2*step),matrix[1:4,i],rtol=1e-7,atol=1e-8)


def test_reference_density_gauge_is_fixed_during_displacement(setup):
    basis=setup[0]; face=basis.model.face; rho=face.bulk.embedding.rho_ref
    for q in ((IDEAL_H,0.,0.), (1.3*IDEAL_H,.3,.1)):
        basis.jet(q)
        assert face.bulk.embedding.rho_ref == rho


def test_exact_bulk_constraints_and_tangent(setup):
    _,_,_,obs,_,matrix=setup
    target=np.array([o.target for o in obs]); offset,T=exact_constraint_tangent(matrix,target)
    np.testing.assert_allclose(matrix[:2]@offset,target[:2],atol=2e-14)
    np.testing.assert_allclose(matrix[:2]@T,0,atol=3e-14)
    np.testing.assert_array_equal(T[2:],np.eye(6))


def test_source_stationary_states_and_dependent_reverse_barrier(setup):
    _,_,source,obs,states,_=setup
    for name in ('saddle','fault'):
        assert max(abs(source.evaluate(states[name]).gradient)) < 2e-9
    assert not any('reverse' in o.name for o in obs)
    assert not any('yield' in o.name or 'fatigue' in o.name for o in obs)
    assert sum(o.role=='heldout' for o in obs)==9


def test_heldout_never_enters_loss_and_fit_is_deterministic(setup):
    _,_,_,obs,_,matrix=setup
    result=fit_coefficients(matrix,obs)
    perturbed=[replace(o,target=o.target+100.) if o.role=='heldout' else o for o in obs]
    again=fit_coefficients(matrix,perturbed)
    np.testing.assert_array_equal(result['coefficients'],again['coefficients'])
    assert result['squared_loss']==again['squared_loss']
    assert not result['physical_candidate_accepted']
    assert all(obs[i].role!='heldout' for i in result['selected_rows'])


def test_opening_inequalities_do_not_clip_forces(setup):
    basis,_,_,obs,_,matrix=setup
    rows=np.array([basis.jet((a*IDEAL_H,0.,0.))[1] for a in (1.1,2.,4.,8.)])
    result=fit_coefficients(matrix,obs,opening_inequalities=rows)
    assert result['admissible']
    c=result['coefficients']; direct=build_coefficient_surface(basis.scalar_decay,basis.angular_decay,c)
    actual=np.array([direct.evaluate((a*IDEAL_H,0.,0.)).gradient[0] for a in (1.1,2.,4.,8.)])
    np.testing.assert_allclose(actual,rows@c,atol=4e-11,rtol=0)


def test_identifiability_excludes_fixed_constraints(setup):
    _,c,_,obs,_,matrix=setup
    info=coefficient_identifiability(matrix,obs,c)
    assert info['jacobian'].shape==(sum(o.role=='fit' for o in obs),6)
    assert len(info['singular_values'])==6
    assert not info['statistical_confidence_claimed']
    assert not info['full_radial_identifiability_checked']


def test_physical_traction_units_are_not_specimen_aggregation():
    stress=np.array([-150.,0.,25.,150.])
    np.testing.assert_allclose(UNITS.force_to_traction_mpa(UNITS.traction_mpa_to_force(stress)),stress,rtol=2e-15)
    assert UNITS.atomic_cell_area_m2==pytest.approx(np.sqrt(3)/2*LENGTH_M**2)
    import inspect
    from . import vector_material_calibration as module
    text=inspect.getsource(module)
    assert 'A_c' not in text and 'specimen_aggregation' not in text


@pytest.mark.parametrize('q', [(0.,0.,0.), (1.,float('nan'),0.)])
def test_invalid_basis_state_rejected(setup,q):
    with pytest.raises(ValueError):
        setup[0].jet(q)


def test_extremum_scan_detects_hidden_negative_traction():
    from .vector_interface_reference import VectorInterfaceEvaluation
    from .vector_material_calibration import opening_traction_extrema
    class Example:
        h=1.
        def evaluate(self,q):
            a=q[0]; H=np.zeros((3,3)); H[0,0]=2*(a-2.8)
            return VectorInterfaceEvaluation(0.,np.array([(a-2.8)**2-.01,0,0]),H,{}, {})
    m=Example()
    assert all(m.evaluate([a,0,0]).gradient[0]>0 for a in (2.5,3.))
    for n in (31,61):
        roots=opening_traction_extrema(m,n)
        assert len(roots)==1
        assert roots[0]['a_over_h']==pytest.approx(2.8)
        assert roots[0]['force_eV_L0']==pytest.approx(-.01)


def test_nonfinite_density_decay_rejected():
    with pytest.raises(ValueError):
        VectorCoefficientBasis(np.nan,4.)


def test_source_no_neighbors_is_exact_empty_sum_not_force_clipping(setup,monkeypatch):
    source=setup[2]
    def forbidden(*args,**kwargs):
        raise AssertionError('empty SOURCE sum must not depend on an einsum reduction')
    monkeypatch.setattr(np,'einsum',forbidden)
    pair,rho=source._plane(2*source.reference.cutoff,np.array([.37,.18]))
    np.testing.assert_array_equal(pair,np.zeros(10))
    np.testing.assert_array_equal(rho,np.zeros(10))


def test_separated_source_plane_has_zero_force_for_multiple_allocations(setup):
    source=setup[2]
    for _ in range(30):
        value=source.evaluate((40*source.h,.31,-.27))
        np.testing.assert_array_equal(value.gradient,np.zeros(3))
        np.testing.assert_array_equal(value.hessian,np.zeros((3,3)))
