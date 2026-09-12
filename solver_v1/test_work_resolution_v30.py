import inspect
import numpy as np
import pytest
from scipy.integrate import trapezoid
from .run_reference_thermostat_md import work_trapezoid_from_right_sum,run
from .run_work_resolution_v30 import analyze


@pytest.mark.parametrize('stride',[1,2,10,40])
def test_block_mean_accumulation_retains_all_internal_powers(stride):
    dt=.001; t=np.arange(401)*dt
    power=3*np.cos(29*t)+.4*t
    right=0.
    for k in range(stride,len(t),stride):
        right+=stride*dt*power[k-stride+1:k+1].mean()
        actual=work_trapezoid_from_right_sum(right,dt,power[0],power[k])
        assert actual==pytest.approx(trapezoid(power[:k+1],t[:k+1]),abs=2e-15)


def test_work_measurement_is_opt_in_and_unresolved_result_preserved(tmp_path):
    assert inspect.signature(run).parameters['timestep_work'].default is False
    out=tmp_path/'user';out.mkdir()
    with pytest.raises(FileExistsError):analyze(tmp_path,out)
    with pytest.raises(ValueError):work_trapezoid_from_right_sum(0,0,0,0)


def test_signed_work_not_clipped():
    assert work_trapezoid_from_right_sum(-3,.1,2,4)==pytest.approx(-3.1)
