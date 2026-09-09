from dataclasses import replace
import json

import numpy as np
import pytest

from .full_fcc_calibration_audit import independent_bulk_targets, cubic_constants_gpa
from .run_low_stress_cyclic_diagnostic import FIT
from .run_vector_material_calibration import source_and_targets
from .vector_material_calibration import VectorCoefficientBasis, observation_matrix, fit_coefficients
from .yield_elastic_metric import cubic_to_mode_matrix, cubic_metric_problem, transformed_mode_covariance, fit_exact_bulk


def test_exact_elastic_conversion_and_zero_pressure():
    H=cubic_to_mode_matrix(); values,_=independent_bulk_targets()
    np.testing.assert_allclose(H@[114,62,32],values[2:5],atol=2e-14)
    for c in ([114,62,32],[87.816,64.833,49.271],[110,65,38]):
        observed=np.r_[0,3.36,H@c]
        recovered=cubic_constants_gpa(observed)
        np.testing.assert_allclose([recovered[f'C{i}_GPa'] for i in (11,12,44)],c,atol=3e-14)


def test_covariance_transform_preserves_old_loss_not_diagonal_approximation():
    target,scales=independent_bulk_targets(); H=cubic_to_mode_matrix()
    covariance=transformed_mode_covariance(scales[2:5])
    delta=np.array([-28.74,2.24,20.54])
    old=np.sum((H@delta/scales[2:5])**2)
    same=delta@np.linalg.solve(covariance,delta)
    assert same==pytest.approx(old,rel=2e-14)
    assert np.max(abs(covariance-np.diag(np.diag(covariance))))>1
    # Changing to independently scaled C is not an invariant re-labeling.
    different=np.sum((delta/(.05*np.array([114,62,32])))**2)
    assert abs(different-old)>50


@pytest.fixture(scope='module')
def actual():
    _,obs,_=source_and_targets()
    basis=VectorCoefficientBasis(2.906199643752374,5.308553273298447)
    mat=observation_matrix(basis,obs)
    return basis,mat,obs


def test_metric_changes_only_three_rows_and_keeps_heldout(actual):
    _,mat,obs=actual; original=mat.copy()
    changed,new=cubic_metric_problem(mat,obs)
    np.testing.assert_array_equal(mat,original)
    np.testing.assert_array_equal(changed[:2],mat[:2])
    np.testing.assert_array_equal(changed[5:],mat[5:])
    assert new[5:]==obs[5:]
    assert [o.units for o in new[2:5]]==['GPa']*3
    assert not any('yield' in o.name for o in new)


def test_five_exact_bulk_constraints_and_fit_reproducibility(actual):
    _,mat,obs=actual; cm,co=cubic_metric_problem(mat,obs)
    fit=fit_exact_bulk(cm,co)
    assert fit['admissible'] and fit['strictly_positive_LJ_resolved']
    np.testing.assert_allclose(cm[2:5]@fit['coefficients'],[114,62,32],atol=4e-11)
    assert fit['exact_bulk_rank']==5 and fit['reduced_jacobian'].shape[1]==3
    changed=[replace(o,target=o.target+1000) if o.role=='heldout' else o for o in co]
    again=fit_exact_bulk(cm,changed)
    np.testing.assert_array_equal(fit['coefficients'],again['coefficients'])
    assert fit['squared_loss']==again['squared_loss']
    assert not fit['physical_candidate_accepted']


def test_old_parameter_file_remains_unchanged():
    import hashlib
    assert hashlib.sha256(FIT.read_bytes()).hexdigest()=='9d00fbf54831c134fe9961e17f0603defdc7958bc3590743b2094264a31b6655'


def test_no_strength_or_kinetic_target_and_no_area_aggregation():
    import inspect
    from . import yield_elastic_metric as module
    text=inspect.getsource(module)
    assert 'A_c' not in text and 'M_a' not in text and 'M_s' not in text


@pytest.mark.parametrize('scale',[0,np.nan,-1])
def test_invalid_metric_scale_refused(actual,scale):
    _,mat,obs=actual
    with pytest.raises(ValueError): cubic_metric_problem(mat,obs,relative_scale=scale)
