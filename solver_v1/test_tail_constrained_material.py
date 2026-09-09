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
