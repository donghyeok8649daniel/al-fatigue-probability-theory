import math
from dataclasses import replace

import numpy as np
import pytest

from solver_v1.aluminum_calibration import cell_ev_to_energy_density
from solver_v1.aluminum_full_fcc_calibration import (
    FullFCCParameters,
    build_full_fcc_calibration_model,
)
from solver_v1.fcc111_active_interface import FCC111ActiveInterface
from solver_v1.fcc111_geometry import DIRECT_110, SHOCKLEY_112
from solver_v1.fcc111_full_energy import FullFCC111StackEnergy


CANDIDATE = FullFCCParameters(
    0.051463163185591194, 0.9278241602685452, 4.6853176139362835,
    5.679412614403284, 2.760290609784449,
)
A0 = 0.816496172881018


@pytest.fixture(scope="module")
def interface():
    bulk = build_full_fcc_calibration_model(CANDIDATE, find_equilibrium=False)
    bulk.a0 = A0
    return FCC111ActiveInterface(bulk, tolerance=2e-8)


def test_perfect_bulk_limit_has_zero_energy_and_force(interface):
    value = interface.evaluate(interface.h, 0.0)
    assert value.energy == pytest.approx(0.0, abs=2e-14)
    assert value.pair_energy == pytest.approx(0.0, abs=2e-14)
    assert value.embedding_energy == pytest.approx(0.0, abs=2e-14)
    assert abs(value.d_da) < 5e-6
    assert abs(value.d_ds) < 1e-10
    assert np.min(np.linalg.eigvalsh(interface.hessian(interface.h, 0.0))) > 0


@pytest.mark.parametrize("path_id", (DIRECT_110, SHOCKLEY_112))
def test_active_interface_registry_periodicity(path_id):
    bulk = build_full_fcc_calibration_model(CANDIDATE, find_equilibrium=False)
    bulk.a0 = A0
    model = FCC111ActiveInterface(bulk, path_id=path_id, tolerance=2e-8)
    left = model.evaluate(model.h*1.03, 0.117)
    right = model.evaluate(model.h*1.03, 0.117+model.period)
    assert left.energy == pytest.approx(right.energy, abs=2e-10)
    assert left.d_ds == pytest.approx(right.d_ds, abs=2e-9)


def test_active_interface_gradient_and_hessian_are_analytic(interface):
    a, s, step = 1.05*interface.h, 0.17, 8e-5
    value = interface.evaluate(a, s)
    energy = interface.energy
    da = (energy(a+step,s)-energy(a-step,s))/(2*step)
    ds = (energy(a,s+step)-energy(a,s-step))/(2*step)
    daa = (energy(a+step,s)-2*energy(a,s)+energy(a-step,s))/step**2
    dss = (energy(a,s+step)-2*energy(a,s)+energy(a,s-step))/step**2
    das = (energy(a+step,s+step)-energy(a+step,s-step)
           -energy(a-step,s+step)+energy(a-step,s-step))/(4*step**2)
    np.testing.assert_allclose([value.d_da,value.d_ds],[da,ds],rtol=2e-6,atol=2e-6)
    np.testing.assert_allclose(
        [value.d2_daa,value.d2_das,value.d2_dss], [daa,das,dss],
        rtol=2e-5, atol=2e-5,
    )


def test_active_interface_reciprocal_matches_direct_reference(interface):
    a, s = 1.05*interface.h, 0.17
    value = interface.evaluate(a, s)
    direct = interface.direct_reference(a, s, radial_index=45, layers=32)
    assert value.energy == pytest.approx(direct.energy, abs=2e-6)
    assert value.d_da == pytest.approx(direct.d_da, abs=3e-5)
    assert value.d_ds == pytest.approx(direct.d_ds, abs=2e-7)
    assert value.d2_daa == pytest.approx(direct.d2_daa, abs=1e-5)


def test_opening_limit_and_gsf_conversion_use_atomic_area_not_Ac(interface):
    work_sep = interface.separated_limit()
    area_angstrom2 = 7.102490842787127
    work_j_m2 = cell_ev_to_energy_density(work_sep, area_angstrom2)
    assert 0.45 < work_sep < 0.60
    assert 1.0 < work_j_m2 < 1.3
    direct_usf = interface.energy(interface.h, 0.5*interface.period)
    assert direct_usf > 0.0
    assert "A_c" not in FCC111ActiveInterface.__init__.__annotations__


def test_interface_layer_shell_and_neighborhood_refinement(interface):
    bulk = interface.bulk
    config = replace(bulk.stack_config, tol=2e-12,
                     reciprocal=replace(bulk.stack_config.reciprocal,tol=2e-14))
    fine_bulk = FullFCC111StackEnergy(
        bulk.p, density_params=bulk.density_params, embedding=bulk.embedding,
        stack_config=config, find_equilibrium=False)
    fine_bulk.a0 = bulk.a0
    fine = FCC111ActiveInterface(fine_bulk,tolerance=2e-12,
                                 consecutive_small_layers=6)
    fields = ("energy","d_da","d_ds","d2_daa","d2_das","d2_dss")
    for a,s in ((interface.h*1.05,.17),(interface.h*2.,.33)):
        first, second = interface.evaluate(a,s),fine.evaluate(a,s)
        np.testing.assert_allclose(
            [getattr(first,k) for k in fields],
            [getattr(second,k) for k in fields],rtol=1e-7,atol=2e-8)
    assert interface.separated_limit() == pytest.approx(
        fine.separated_limit(),abs=1e-8)


def test_real_space_radius_and_layer_refinement_independently(interface):
    a,s=1.05*interface.h,.17
    reference=interface.evaluate(a,s)
    coarse=interface.direct_reference(a,s,radial_index=20,layers=20)
    finer_r=interface.direct_reference(a,s,radial_index=40,layers=20)
    finer_l=interface.direct_reference(a,s,radial_index=20,layers=40)
    combined=interface.direct_reference(a,s,radial_index=60,layers=60)
    # Derivative tail is nonnegative here; both omitted-domain refinements
    # independently improve its value.
    error=lambda x:abs(x.d_da-reference.d_da)
    assert error(finer_r)<error(coarse)
    assert error(finer_l)<error(coarse)
    assert error(combined)<error(finer_r)
    assert error(combined)<error(finer_l)


def test_opening_asymptote_is_separation_work_not_a_finite_saddle(interface):
    limit=interface.separated_limit()
    near=interface.energy(10*interface.h,0)
    far=interface.energy(40*interface.h,0)
    assert abs(far-limit)<abs(near-limit)
    assert interface.grad(40*interface.h,0)[0]>0
