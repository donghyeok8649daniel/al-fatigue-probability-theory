import numpy as np
import pytest

from .coordination_screening import CoordinationScreenedInterface
from .current_material_rows import CurrentMaterialSiteLaw
from .quartic_jet_ablation import RationalQuarticStaticInterface, RationalQuarticJetColumn
from .vector_material_calibration import MaterialObservation as O


SHAPE = (3., 6., 8., 2., 5., -.2)
COEFFICIENTS = np.array([.005, .012, 2., 1., 2., 10., 2., 2., 1., 100.])


@pytest.fixture(scope='module')
def base():
    return CoordinationScreenedInterface(SHAPE, COEFFICIENTS, law='power')


def test_zero_ablation_recovers_full_existing_jet(base):
    trial = RationalQuarticStaticInterface(SHAPE, COEFFICIENTS, quartic_saturation=0.)
    for q in [(base.h, 0., 0.), (1.12*base.h, .17, -.04)]:
        old, new = base.evaluate(q), trial.evaluate(q)
        assert new.energy == pytest.approx(old.energy, abs=1e-14)
        np.testing.assert_allclose(new.gradient, old.gradient, atol=1e-13)
        np.testing.assert_allclose(new.hessian, old.hessian, atol=1e-12)


def test_nonzero_ablation_keeps_initial_jet_but_changes_finite_state(base):
    trial = RationalQuarticStaticInterface(SHAPE, COEFFICIENTS, quartic_saturation=1000.)
    old, new = base.evaluate((base.h, 0., 0.)), trial.evaluate((base.h, 0., 0.))
    np.testing.assert_allclose(new.hessian, old.hessian, atol=1e-12)
    q = (1.12*base.h, .17, -.04)
    assert abs(base.evaluate(q).energy-trial.evaluate(q).energy) > 1e-6
    with pytest.raises(ValueError, match='not implemented'):
        CurrentMaterialSiteLaw(trial)


def test_cap_column_is_same_full_model_difference_and_fd(base):
    q = np.array([1.12*base.h, .17, -.04])
    observations = [O('E', 1., 1., 'eV', 'fit', tuple(q), tuple(np.eye(10)[0])),
                    O('Haa', 1., 1., 'eV/L0^2', 'fit', tuple(q), tuple(np.eye(10)[4]))]
    cache = RationalQuarticJetColumn(base, observations)
    trial = RationalQuarticStaticInterface(SHAPE, COEFFICIENTS, quartic_saturation=1000.)
    difference = COEFFICIENTS[9]*(cache.column(1000.)-cache.column(0.))
    old, new = base.evaluate(q), trial.evaluate(q)
    np.testing.assert_allclose(difference, [new.energy-old.energy, new.hessian[0, 0]-old.hessian[0, 0]],
                               atol=2e-12, rtol=1e-10)
    step = 2e-6
    numerical = np.column_stack([(trial.evaluate(q+step*e).gradient-trial.evaluate(q-step*e).gradient)/(2*step)
                                 for e in np.eye(3)])
    np.testing.assert_allclose(numerical, new.hessian, atol=4e-7, rtol=2e-7)
    assert cache.reference_invariant(q) > 0
