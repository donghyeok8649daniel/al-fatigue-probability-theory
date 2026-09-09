import numpy as np
import pytest

from .spectral_material_constraints import BulkHessianCoefficientBasis,stability_halfspace,fit_spectral_faces
from .range_resolved_material import RangeBulkValidation,build_range_surface,RangeObservationCache
from .profiled_range_calibration import fit_profiled_coefficients
from .run_vector_material_calibration import source_and_targets
from .yield_elastic_metric import cubic_metric_problem


def test_finite_q_exact_coefficient_linearity_with_negative_angular_terms():
    decays=(2.3,5.2,3.1); coefficients=[.15,.6,3.,2.,1.,20.,-2.,-8.]
    basis=BulkHessianCoefficientBasis(decays,radius=4.)
    direct=RangeBulkValidation(build_range_surface(*decays,coefficients),cutoff=4.)
    q=[.375,.375,0.]; cols=basis.matrices(q)
    expected=direct.evaluate(direct.crystallographic_wavevector(q))['matrix']
    np.testing.assert_allclose(np.einsum('c,cij->ij',coefficients,cols),expected,rtol=2e-12,atol=3e-11)
    probe=stability_halfspace(cols,coefficients)
    assert np.dot(probe['halfspace'],coefficients)==pytest.approx(np.linalg.eigvalsh(expected)[0],abs=3e-11)


def test_no_spectral_constraints_matches_same_convex_profile():
    obs=source_and_targets()[1]
    matrix,obs=cubic_metric_problem(RangeObservationCache(obs).matrix((2.58,6.,4.96)),obs)
    original=fit_profiled_coefficients(matrix,obs)
    test=fit_spectral_faces(matrix,obs,[])
    np.testing.assert_array_equal(original['coefficients'],test['coefficients'])
    assert original['squared_loss']==test['squared_loss']
    # A necessary test halfspace is not added to the physical energy.
    lower=np.zeros(8); lower[6]=1.
    bounded=fit_spectral_faces(matrix,obs,[lower])
    assert bounded['coefficients'][6]>=-1e-10
    assert bounded['kkt_stationarity']<1e-8
    assert not bounded['physical_candidate_accepted']
    altered=matrix.copy()
    altered[[i for i,o in enumerate(obs) if o.role=='heldout']]*=2.3
    np.testing.assert_array_equal(bounded['coefficients'],fit_spectral_faces(altered,obs,[lower])['coefficients'])
