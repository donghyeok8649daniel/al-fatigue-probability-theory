"""Reject nonfinite scalar comparisons and invalid independent-control records."""
import math
import pytest
from results.silicon_wafer_feasibility.replay_mpa_qe_v13 import scalar_difference,check_controls


def test_nan_can_hide_from_legacy_max_but_is_rejected():
    # Reproduce the specific max(known_error, NaN) masking behavior.
    assert max(0.,abs(float('nan')-1.))==0.
    with pytest.raises(ValueError,match='nonfinite'):scalar_difference(float('nan'),1.)


@pytest.mark.parametrize('actual,expected',[(1.,float('nan')),(float('inf'),1.),(1.,-float('inf'))])
def test_reject_nonfinite_either_side(actual,expected):
    with pytest.raises(ValueError,match='nonfinite'):scalar_difference(actual,expected)


def test_preserves_small_finite_discrepancy_and_control_sign():
    assert math.isclose(scalar_difference(2.,2.+1e-12),1e-12,rel_tol=.001)
    check_controls([dict(frame=7,energy_difference_eV=-1e-10,force_difference_eV_A=1e-10)],[7])


@pytest.mark.parametrize('kind',['duplicate','missing','nan_energy','nan_force','negative_norm','large_error'])
def test_invalid_control_records(kind):
    rows=[dict(frame=i,energy_difference_eV=0.,force_difference_eV_A=0.) for i in [2,7,11]]
    if kind=='duplicate':rows[1]['frame']=2
    elif kind=='missing':rows.pop()
    elif kind=='nan_energy':rows[1]['energy_difference_eV']=float('nan')
    elif kind=='nan_force':rows[1]['force_difference_eV_A']=float('nan')
    elif kind=='negative_norm':rows[1]['force_difference_eV_A']=-1e-12
    else:rows[1]['force_difference_eV_A']=1e-4
    with pytest.raises(ValueError):check_controls(rows,[2,7,11])
