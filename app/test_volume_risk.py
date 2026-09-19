import numpy as np
import pytest

from .volume_risk import volume_risk


@pytest.mark.parametrize('p', [0., 1e-16, .01, .7, 1.])
def test_hazard_density_and_combined_risk_are_subdivision_invariant(p):
    h, risk, unit = volume_risk(p, 8., 2.)
    h_split, split, unit_split = volume_risk(p, np.array([1., 2., 5.]), 2.)
    np.testing.assert_array_equal(h_split, np.full(3, h))
    np.testing.assert_array_equal(unit_split, np.full(3, unit))
    with np.errstate(divide='ignore'):
        combined = -np.expm1(np.log1p(-split).sum())
    assert combined == pytest.approx(risk, rel=1e-13, abs=1e-30)
    assert risk == pytest.approx(-np.expm1(4*np.log1p(-p)) if p < 1 else 1)


def test_cell_risk_differs_from_density_and_no_area_enters():
    h, p, unit = volume_risk(.1, [1., 10.], .5)
    np.testing.assert_allclose(h, -np.log(.9)/.5)
    assert p[1] > p[0]
    assert unit[0] == pytest.approx(unit[1])
    h, p, _ = volume_risk(1e-300, 1e300, 1e300)
    assert p == pytest.approx(1e-300, rel=1e-12, abs=0)
    assert h == 0  # underflow in per-mm^3 display does not erase the cell signal


@pytest.mark.parametrize('vc', [0., -1., np.nan, np.inf])
def test_no_implicit_correlation_volume(vc):
    with pytest.raises(ValueError, match='risk.need_volume'):
        volume_risk(.1, 1., vc)
