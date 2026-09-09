"""The actually fitted research candidate: scoped calibration, not Al yield."""
from dataclasses import replace
from pathlib import Path
import hashlib
import json

import numpy as np
import pytest

from solver_v1.validate_tail_calibration import load_material
from solver_v1.run_vector_material_calibration import source_and_targets
from solver_v1.tail_constrained_material import convex_profile
from solver_v1.quartic_angular_material import QuarticSymmetryInterface
from solver_v1.vector_material_calibration import UNITS
from solver_v1.vector_registry_audit import stationary_state
from solver_v1.yield_elastic_metric import cubic_metric_problem,cubic_to_mode_matrix


ROOT=Path(__file__).resolve().parents[1]
RESULT=ROOT/'results/fcc111_active_interface/tail_calibration_v14/quartic_exact_bulk'


@pytest.fixture(scope='module')
def fitted():
    return load_material(RESULT)


def test_actual_calibration_reproducible_with_exact_bulk_and_heldout_excluded(fitted):
    data,definition,_,cache_type,_=fitted
    _,observations,_=source_and_targets()
    best=data['best'];raw=cache_type(observations).matrix(best['decays'])
    matrix,obs=cubic_metric_problem(raw[:,:8],observations)
    extra=raw[:,8:].copy();extra[2:5]=np.linalg.solve(cubic_to_mode_matrix(),extra[2:5])
    matrix=np.column_stack([matrix,extra])
    obs=[replace(o,role='exact') if i<5 else o for i,o in enumerate(obs)]
    fit=convex_profile(matrix,obs,nonnegative=(0,1,2,4,5,8,9))
    np.testing.assert_allclose(fit['coefficients'],best['coefficients'],rtol=2e-7,atol=2e-7)
    np.testing.assert_allclose(fit['predictions'][:5],[0.,3.36,114.,62.,32.],atol=2e-10)
    assert fit['kkt_residual']<1e-6
    assert all(obs[i].role=='fit' for i in fit['selected_rows'])
    assert not data['material_accepted']
    assert not definition['physical_kinetics_calibrated']


def test_fitted_equilibria_are_actual_roots_not_source_point_snapshots(fitted):
    _,_,model,_,_=fitted
    _,_,states=source_and_targets()
    energies=[]
    for label,morse in [('perfect',0),('fault',0),('saddle',1)]:
        root=stationary_state(model,states[label],expected_index=morse)
        assert root['valid'] and root['morse_index']==morse
        assert np.linalg.norm(root['evaluation'].gradient)<1e-9
        energies.append(float(UNITS.energy_to_surface(root['evaluation'].energy)))
    np.testing.assert_allclose(energies,[0.,.1492292528432122,.17912665414840465],atol=2e-8)


def test_fitted_large_nonlinear_amplitude_has_converged_analytic_jets(fitted):
    data,_,model,_,_=fitted
    b=data['best'];fine=QuarticSymmetryInterface(b['decays'],b['coefficients'],tolerance=2e-13)
    for q in ([model.h,0.,0.],[1.04*model.h,.34,.20],[1.5*model.h,.1,-.03]):
        q=np.asarray(q);value=model.evaluate(q);ref=fine.evaluate(q)
        assert value.energy==pytest.approx(ref.energy,abs=2e-8)
        np.testing.assert_allclose(value.hessian,ref.hessian,rtol=1e-8,atol=2e-7)
        for axis in range(3):
            d=np.eye(3)[axis]*1e-5
            plus,minus=model.evaluate(q+d),model.evaluate(q-d)
            assert (plus.energy-minus.energy)/2e-5==pytest.approx(value.gradient[axis],rel=3e-6,abs=2e-7)
            np.testing.assert_allclose((plus.gradient-minus.gradient)/2e-5,value.hessian[:,axis],rtol=3e-6,atol=4e-6)
    np.testing.assert_array_equal(model.coefficients,b['coefficients'])


def test_fitted_finite_q_refinement_above_analytic_tail_without_whole_zone_claim(fitted):
    data,_,_,_,bulk_type=fitted
    b=data['best'];c=np.asarray(b['coefficients'])
    coarse=bulk_type(b['decays'],radius=12.);fine=bulk_type(b['decays'],radius=16.)
    for q in ((.01,0.,0.),(.375,.375,0.),(.8,.5,.1)):
        H,e=coarse.evaluate(q);Hf,ef=fine.evaluate(q)
        actual=np.einsum('c,cij->ij',c,H)
        reference=np.einsum('c,cij->ij',c,Hf)
        assert np.linalg.eigvalsh(actual)[0]>e@abs(c)
        assert np.linalg.eigvalsh(reference)[0]>ef@abs(c)
        assert np.linalg.norm(reference-actual,2)<=e@abs(c)+ef@abs(c)
    decision=json.loads((RESULT/'validation/decision.json').read_bytes())
    assert not decision['whole_zone_proof'] and not decision['full_Al_material_accepted']


def test_historical_calibration_and_physical_time_status_unchanged():
    original=ROOT/'results/fcc111_active_interface/matched_v3/monotone_opening/fitted_candidates.json'
    assert hashlib.sha256(original.read_bytes()).hexdigest()=='9d00fbf54831c134fe9961e17f0603defdc7958bc3590743b2094264a31b6655'
    decision=json.loads((RESULT/'validation/decision.json').read_bytes())
    assert not decision['experimental_yield_validated']
    assert not decision['physical_seconds'] and not decision['physical_Hz']
