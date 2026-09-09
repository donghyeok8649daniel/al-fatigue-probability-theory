import numpy as np
import pytest
from scipy.integrate import quad

from solver_v1.tail_constrained_material import (
    FCCTailEnvelope, TailBulkCoefficientBasis, convex_profile, spectral_profile,
    declared_wavepoints,
)
from solver_v1.spectral_material_constraints import BulkHessianCoefficientBasis
from solver_v1.profiled_range_calibration import fit_profiled_coefficients
from solver_v1.range_resolved_material import RangeObservationCache
from solver_v1.run_vector_material_calibration import source_and_targets
from solver_v1.yield_elastic_metric import cubic_metric_problem


def test_analytic_tail_integrals_and_monotonicity():
    tail=FCCTailEnvelope(8.)
    for p in (6,8,12,14):
        numeric=tail.factor*quad(lambda x:(x+tail.cover)**2*x**-p,tail.lower,np.inf,
                                epsabs=1e-16)[0]
        assert tail.power(p)==pytest.approx(numeric,rel=1e-10,abs=1e-17)
        assert FCCTailEnvelope(10.).power(p)<tail.power(p)
    for k in (2.,4.,6.):
        for m in range(6):
            numeric=tail.factor*quad(lambda x:(x+tail.cover)**2*x**m*np.exp(-k*x),
                                    tail.lower,np.inf,epsabs=1e-20)[0]
            assert tail.exponential(k,m)==pytest.approx(numeric,rel=1e-9,abs=1e-18)
    with pytest.raises(ValueError): tail.exponential(.1,3)


def test_exact_coefficient_operator_matches_independent_legacy():
    decays=(3.1,4.8,3.7)
    new=TailBulkCoefficientBasis(decays,radius=6.)
    old=BulkHessianCoefficientBasis(decays,radius=6.)
    for q in ((.23,.37,.17),(.375,.375,0.)):
        columns,_=new.evaluate(q)
        np.testing.assert_allclose(columns,old.matrices(q),rtol=2e-10,atol=3e-12)
    c,t=new.evaluate((0.,0.,0.))
    np.testing.assert_array_equal(c,np.zeros_like(c))
    np.testing.assert_array_equal(t,np.zeros_like(t))


def test_analytic_tail_bounds_independent_radius_changes():
    decays=(3.,4.,5.)
    coarse=TailBulkCoefficientBasis(decays,radius=7.)
    fine=TailBulkCoefficientBasis(decays,radius=11.)
    for q in ((.01,.01,0.),(.375,.375,0.),(.7,.1,.3)):
        columns,bound=coarse.evaluate(q); reference,_=fine.evaluate(q)
        actual=np.linalg.norm(reference-columns,ord=2,axis=(-2,-1))
        assert np.all(actual<=bound+2e-12)
        assert np.all(bound>0)


def test_whitened_profile_matches_independent_active_face_and_excludes_holdout():
    from dataclasses import replace
    _,obs,_=source_and_targets()
    matrix,obs=cubic_metric_problem(RangeObservationCache(obs).matrix((3.,5.,4.)),obs)
    first=convex_profile(matrix,obs); reference=fit_profiled_coefficients(matrix,obs)
    np.testing.assert_allclose(first['coefficients'],reference['coefficients'],rtol=1e-8,atol=2e-8)
    modified=[replace(o,target=o.target+123.) if o.role=='heldout' else o for o in obs]
    second=convex_profile(matrix,modified)
    np.testing.assert_array_equal(first['coefficients'],second['coefficients'])
    assert first['exact_residual']<1e-9


def test_full_polarization_with_tail_is_not_single_zero_eigenvalue():
    from solver_v1.vector_material_calibration import MaterialObservation
    obs=[MaterialObservation('a',1.,1.,'test','exact'),
         MaterialObservation('b',1.,1.,'test','exact')]
    obs += [MaterialObservation(str(i),(-2. if i==7 else 1.),1.,'test','fit') for i in range(2,8)]
    matrix=np.eye(8)
    columns=np.zeros((1,8,3,3));columns[0,0]=np.eye(3);columns[0,7,0,0]=1.
    tails=np.full((1,8),.001)
    fit=spectral_profile(matrix,obs,columns,tails)
    c=fit['coefficients']; eig=np.linalg.eigvalsh(np.einsum('c,cij->ij',c,columns[0]))
    assert eig[0]>=tails[0]@abs(c)-1e-10
    assert eig[0]>0
    assert not fit['whole_zone_proved']


def test_wavepoint_geometry_contains_previous_instability():
    points=declared_wavepoints()
    assert any(np.array_equal(q,[.375,.375,0.]) for q in points)
    assert np.all(points[:,0]>=points[:,1]) and np.all(points[:,1]>=points[:,2])
    assert np.max(points.sum(axis=1))<=1.5+1e-12


def test_zero_dual_near_active_constraint_does_not_change_verified_optimum():
    """Regression for the actual v15 active-face failure, in a solvable QP.

    The large whitened-coordinate norm makes the positive third coefficient
    'near active' under the screening tolerance. Its exact dual is ZERO:
    screening is not authority to force that coefficient to zero. The old
    physical-face polish did that and turned a verified optimum into KKT=1.
    No material tolerance or physical sign constraint is relaxed here.
    """
    from solver_v1.vector_material_calibration import MaterialObservation
    targets=(1.,1.,1.,1e8)
    obs=[MaterialObservation(str(i),y,1.,'test',
         'exact' if i<2 else 'fit') for i,y in enumerate(targets)]
    matrix=np.diag([1.,1.,1.,1e-6])
    fit=convex_profile(matrix,obs,nonnegative=(0,1,2,3))
    np.testing.assert_allclose(fit['coefficients'],[1.,1.,1.,1e14],rtol=2e-14)
    assert fit['kkt_residual']<1e-6
    assert fit['squared_loss']<1e-12
    assert fit['scaled_primal_violation']==0.


def test_coefficient_unit_rescaling_preserves_physical_profile_and_active_zero():
    from solver_v1.vector_material_calibration import MaterialObservation
    obs=[MaterialObservation(str(i),y,1.,'test','exact' if i<2 else 'fit')
         for i,y in enumerate((1.,1.,-2.,3.))]
    reference=convex_profile(np.eye(4),obs,nonnegative=(0,1,2,3))
    for unit_scale in (1e-5,1.,1e5):
        matrix=np.diag([1.,1.,unit_scale,1/unit_scale])
        fit=convex_profile(matrix,obs,nonnegative=(0,1,2,3))
        np.testing.assert_allclose(fit['predictions'],reference['predictions'],atol=2e-12)
        assert fit['coefficients'][2]==0.
        assert fit['squared_loss']==pytest.approx(4.,abs=1e-12)
        assert fit['kkt_residual']<1e-6


def test_nearly_parallel_affine_face_does_not_square_condition_number():
    from solver_v1.tail_constrained_material import project_affine_face
    A=np.array([[1.,0.,0.],[1.,1e-8,0.]])
    h=np.array([0.,1e-8]);start=np.array([-1.,-3.,4.])
    actual=project_affine_face(A,h,start)
    np.testing.assert_allclose(actual,[0.,1.,4.],atol=1e-14)
    assert np.max(abs(A@actual-h))<1e-22
    # Independent well-known failure mode: Gram normal equations discard
    # the 1e-16 singular direction although A itself still resolves 1e-8.
    gram=start+A.T@np.linalg.lstsq(A@A.T,h-A@start,rcond=None)[0]
    assert abs(gram[1]-actual[1])>1


def test_counterexample_wavepoint_refinement_preserves_base_and_rejects_bad_coordinates():
    from solver_v1.tail_constrained_material import augmented_wavepoints
    point=[.14463521202111632]*3
    base=declared_wavepoints(.1);actual=augmented_wavepoints(.1,[point,point])
    assert len(actual)==len(base)+1
    assert all(any(np.array_equal(q,r) for r in actual) for q in base)
    assert any(np.array_equal(q,point) for q in actual)
    for invalid in ([[0.,0.,0.]],[[.1,.2,0.]],[[.8,.8,.8]],[[1.,0.,np.nan]]):
        with pytest.raises(ValueError):augmented_wavepoints(additional=invalid)
