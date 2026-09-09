import hashlib
import json

import numpy as np
import pytest

from .range_resolved_material import (
    RangeResolvedBasis, RangeObservationCache, build_range_surface,
    range_observation_matrix,RangeBulkValidation,
)
from .run_low_stress_cyclic_diagnostic import FIT,ROOT
from .run_vector_material_calibration import source_and_targets
from .vector_material_calibration import (
    VectorCoefficientBasis, build_coefficient_surface, observation_matrix,
)


COEFF = [.09, .24, 3.9, 3.1, 1.3, 100., -3., 50.]


@pytest.fixture(scope='module')
def observations():
    return source_and_targets()[1]


def test_equal_ranges_recover_existing_entire_analytic_jet():
    old = build_coefficient_surface(1.3, 4.9, COEFF)
    new = build_range_surface(1.3, 4.9, 4.9, COEFF)
    for q in [(old.h, 0., 0.), (1.1*old.h, .23, -.12), (1.5*old.h, .5, .2)]:
        a, b = old.evaluate(q), new.evaluate(q)
        assert a.energy == b.energy
        np.testing.assert_array_equal(a.gradient, b.gradient)
        np.testing.assert_array_equal(a.hessian, b.hessian)


def test_observable_column_separation_and_exact_cache(observations):
    old = observation_matrix(VectorCoefficientBasis(1.3, 4.9), observations)
    new = range_observation_matrix(RangeResolvedBasis(1.3, 4.9, 3.2), observations)
    np.testing.assert_array_equal(old[:, :7], new[:, :7])
    assert np.max(abs(old[:, 7]-new[:, 7])) > .001
    np.testing.assert_array_equal(new[:2, 7], [0., 0.])
    cached = RangeObservationCache(observations).matrix((1.3, 4.9, 3.2))
    np.testing.assert_allclose(cached, new, rtol=2e-13, atol=2e-13)


def test_range_derivatives_and_periodic_registry():
    model = build_range_surface(1.3, 4.9, 3.2, COEFF)
    q = np.array([1.11*model.h, .23, -.12]); value = model.evaluate(q)
    errors = []
    for step in (4e-5, 2e-5):
        grad = []; hessian = []
        for d in np.eye(3)*step:
            plus, minus = model.evaluate(q+d), model.evaluate(q-d)
            grad.append((plus.energy-minus.energy)/(2*step))
            hessian.append((plus.gradient-minus.gradient)/(2*step))
        errors.append([np.max(abs(grad-value.gradient)),
                       np.max(abs(np.array(hessian).T-value.hessian))])
    assert np.all(np.array(errors[1]) < .3*np.array(errors[0]))
    for delta in (model.geometry.a1, model.geometry.a2):
        equivalent = model.evaluate(q+np.r_[0., delta])
        np.testing.assert_allclose(value.hessian, equivalent.hessian, atol=2e-12)
    np.testing.assert_array_equal(value.hessian, value.hessian.T)


def test_rank2_direct_sum_independent_and_not_global_embedding():
    model = build_range_surface(1.3, 4.9, 3.2, COEFF)
    obj = model.surface.quadrupole
    a, s = 1.14*model.h, .27
    exact = obj.evaluate(a, s)[0][0]
    errors = [abs(obj.direct_value(a, s, radius=n, layers=n)-exact) for n in (8, 16, 24)]
    assert errors[-1] < 2e-13
    assert errors[1] < errors[0]
    assert model.moments[-1][1].invariant is obj


@pytest.mark.parametrize('decays', [(1., 2., 0.), (np.nan, 3., 4.), (1., -1., 4.)])
def test_invalid_ranges_refused(decays):
    with pytest.raises(ValueError):
        build_range_surface(*decays, COEFF)


def test_historical_parameters_not_overwritten():
    assert hashlib.sha256(FIT.read_bytes()).hexdigest() == '9d00fbf54831c134fe9961e17f0603defdc7958bc3590743b2094264a31b6655'


def test_exact_cubic_fit_cannot_hide_its_finite_q_instability():
    study=json.loads((ROOT/'results/fcc111_active_interface/range_core_v12/material/calibration.json').read_bytes())
    fit=study['best']['bulk_exact']
    np.testing.assert_allclose(fit['cubic_GPa'],[114,62,32],atol=2e-10)
    model=build_range_surface(fit['scalar_decay'],fit['odd_decay'],fit['quadrupole_decay'],fit['coefficients'])
    values=[]
    for radius in (6.,8.):
        check=RangeBulkValidation(model,cutoff=radius)
        q=check.crystallographic_wavevector([.375,.375,0.])
        values.append(check.evaluate(q)['eigenvalues'][0])
    assert max(values)<-25.
    assert abs(values[1]-values[0])<.05
    assert not study['material_accepted']


def test_full_local_range_rank_does_not_imply_material_calibration():
    data=json.loads((ROOT/'results/fcc111_active_interface/range_core_v12/material/identifiability.json').read_bytes())
    info=data['cubic_5pct']; J=np.asarray(info['full_jacobian'])
    np.testing.assert_allclose(np.linalg.svd(J,compute_uv=False),info['full_singular_values'],rtol=2e-13)
    assert info['full_rank']==9
    assert not info['on_common_range_submanifold']
    assert not info['confidence_intervals_available']
