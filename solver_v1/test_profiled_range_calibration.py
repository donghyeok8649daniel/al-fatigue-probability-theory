import numpy as np
import pytest

from .profiled_range_calibration import fit_profiled_coefficients
from .range_resolved_material import RangeObservationCache
from .run_vector_material_calibration import source_and_targets
from .vector_material_calibration import fit_coefficients
from .yield_elastic_metric import cubic_metric_problem


@pytest.fixture(scope='module')
def problem():
    obs=source_and_targets()[1]
    return cubic_metric_problem(RangeObservationCache(obs).matrix((2.579101374767287,5.996534839345317,4.962329189737055)),obs)


def test_profile_agrees_with_previous_constrained_fit_without_new_physics(problem):
    matrix,obs=problem; actual=fit_profiled_coefficients(matrix,obs)
    old=fit_coefficients(matrix,obs)
    assert actual['squared_loss']<=old['squared_loss']+2e-9
    np.testing.assert_allclose(actual['predictions'][old['selected_rows']],old['predictions'],atol=3e-5,rtol=2e-6)
    assert actual['equality_residual']<1e-10
    assert actual['kkt_free_stationarity']<1e-8
    assert actual['strictly_positive_LJ_resolved']
    assert not actual['physical_candidate_accepted']


def test_heldout_does_not_enter_profile_and_repeat_is_deterministic(problem):
    matrix,obs=problem; expected=fit_profiled_coefficients(matrix,obs)
    alternate=matrix.copy()
    alternate[[i for i,o in enumerate(obs) if o.role=='heldout']]*=3.1
    altered=fit_profiled_coefficients(alternate,obs)
    np.testing.assert_array_equal(expected['coefficients'],altered['coefficients'])
    np.testing.assert_array_equal(expected['coefficients'],fit_profiled_coefficients(matrix,obs)['coefficients'])


def test_original_parameter_file_stays_bound_to_historical_result():
    import hashlib
    from .run_low_stress_cyclic_diagnostic import FIT
    assert hashlib.sha256(FIT.read_bytes()).hexdigest()=='9d00fbf54831c134fe9961e17f0603defdc7958bc3590743b2094264a31b6655'
