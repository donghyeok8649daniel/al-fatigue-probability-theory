"""Physical unit, independent-mode, and deterministic fitting regressions."""
import math

import numpy as np
import pytest

from solver_v1.aluminum_calibration import EV_J
from solver_v1.aluminum_full_fcc_calibration import (
    FullFCCParameters, build_full_fcc_calibration_model,
    evaluate_full_fcc_bulk_observables,
)
from solver_v1.full_fcc_calibration_audit import (
    IDEAL_H, atomic_volume_m3, bulk_basis, cubic_constants_gpa,
    deterministic_profile, independent_bulk_targets, independent_observables,
    independent_sensitivity, parameter_coefficients, profile_point,
    bounded_independent_fit,
    isotropic_bulk_equilibrium, relaxed_bulk_observables,
)


LEGACY = FullFCCParameters(.051463163185591194, .9278241602685452,
                          4.6853176139362835, 5.679412614403284,
                          2.760290609784449)


def test_three_independent_elastic_modes_and_unit_conversion():
    target, _ = independent_bulk_targets()
    recovered = cubic_constants_gpa(target)
    assert recovered["C11_GPa"] == pytest.approx(114)
    assert recovered["C12_GPa"] == pytest.approx(62)
    assert recovered["C44_GPa"] == pytest.approx(32)
    assert atomic_volume_m3() == pytest.approx(4.05**3/4*1e-30, rel=2e-14)
    assert target[4] == pytest.approx(28e9*atomic_volume_m3()/EV_J)
    # (X,C44,D) -> (hydrostatic,normal,simple shear) is invertible.
    assert np.linalg.matrix_rank([[3, 0, 0], [1/3, 4/3, 0],
                                  [0, 1/3, 1/3]]) == 3


def test_legacy_fit_does_not_determine_all_cubic_constants():
    c = cubic_constants_gpa(independent_observables(LEGACY))
    assert c["C44_GPa"] == pytest.approx(32, abs=.001)
    assert abs(c["C11_GPa"]-114) > 18
    assert abs(c["C12_GPa"]-62) > 9


def test_basis_matches_original_reciprocal_energy_derivatives():
    y = independent_observables(LEGACY)
    model = build_full_fcc_calibration_model(LEGACY)
    center = model.evaluate(IDEAL_H, 0)
    np.testing.assert_allclose(
        y[[0,1,3,4]],
        [IDEAL_H*center.d_da, -LEGACY.embedding_linear_ev-center.energy,
         IDEAL_H**2*center.d2_daa, IDEAL_H**2*center.d2_dss],
        rtol=3e-10, atol=2e-11)
    old = evaluate_full_fcc_bulk_observables(LEGACY, strain_step=1e-4)
    assert y[2] == pytest.approx(
        old.strain_h_aa_ev+2*old.strain_h_ab_ev+old.strain_h_bb_ev,
        abs=2e-5)


def test_simple_shear_is_homogeneous_strain_not_local_gsf():
    model = build_full_fcc_calibration_model(LEGACY)
    dg = 2e-4
    w = model.energy(IDEAL_H, 0)
    numerical = (model.energy(IDEAL_H, dg*IDEAL_H)-2*w
                 +model.energy(IDEAL_H,-dg*IDEAL_H))/dg**2
    assert numerical == pytest.approx(independent_observables(LEGACY)[4],
                                      rel=2e-6)


def test_density_gauge_is_fixed_under_physical_strain():
    model = build_full_fcc_calibration_model(LEGACY)
    ref = model.embedding.rho_ref
    shifted = build_full_fcc_calibration_model(LEGACY, b=1.01, rho_ref=ref)
    assert shifted.embedding.rho_ref == ref
    assert shifted.full_environment_density(IDEAL_H*1.01,0).value/ref < 1


def test_hydrostatic_stencil_refines_without_changing_the_energy():
    values = [independent_observables(LEGACY,strain_step=h)
              for h in (1.6e-3,8e-4,4e-4)]
    np.testing.assert_allclose(values[0],values[2],rtol=2e-8,atol=2e-7)
    np.testing.assert_allclose(values[1],values[2],rtol=2e-8,atol=2e-7)


def test_log_sensitivity_predicts_small_parameter_change():
    y0 = independent_observables(LEGACY)
    jac = independent_sensitivity(LEGACY)
    _, scales = independent_bulk_targets()
    direction = np.array([.3,-.1,.2,.1,-.2])
    delta = 1e-5
    p1 = FullFCCParameters(*(LEGACY.as_array()*np.exp(delta*direction)))
    difference = (independent_observables(p1)-y0)/scales
    np.testing.assert_allclose(difference,delta*jac@direction,
                               rtol=8e-4,atol=2e-7)


def test_profile_is_reproducible_and_coefficient_mapping_is_exact():
    first = profile_point(5.2,extended=True)
    second = profile_point(5.2,extended=True)
    assert first.admissible and first.parameters is not None
    np.testing.assert_array_equal(first.coefficients,second.coefficients)
    np.testing.assert_allclose(parameter_coefficients(first.parameters),
                               first.coefficients,rtol=3e-15)
    np.testing.assert_allclose(independent_observables(first.parameters),
                               first.prediction,rtol=1e-12,atol=3e-13)


def test_profile_uses_only_declared_observations():
    # Held-out data do not exist in the bulk basis or loss.
    point = profile_point(4.7,extended=True,observation_indices=(0,1))
    assert point.squared_loss == pytest.approx(
        float(point.residual[:2]@point.residual[:2]),abs=1e-25)
    assert bulk_basis(4.7).shape == (5,4)


def test_deterministic_profile_refines_declared_interval():
    best, rows = deterministic_profile(
        extended=True,decay_grid=np.array([4.8,5.,5.2,5.4,5.6]))
    assert best.squared_loss < 1e-7
    assert 5 < best.decay < 5.4
    assert len(rows) >= 5
    assert best.parameters.epsilon_lj_ev > 2  # flag, not hide large cancellation


def test_bounded_optimizer_reports_limits_and_actual_stopping_status():
    p, report=bounded_independent_fit(LEGACY,max_nfev=2)
    assert report["nfev"]<=2
    assert report["squared_loss"]>=0
    assert "success" in report and "message" in report
    assert np.all(p.as_array()>=report["lower"])
    assert np.all(p.as_array()<=report["upper"])
    with pytest.raises(ValueError):
        bounded_independent_fit(FullFCCParameters(20.,.7,5.,4.,55.))


def test_bulk_reference_relaxes_isotropically_and_preserves_density_gauge():
    # An imperfect target fit must not become an axially strained "FCC bulk".
    p=FullFCCParameters(.8466141,.8069657,1.36765,10.1664764,12.)
    original=build_full_fcc_calibration_model(p)
    model,stretch=isotropic_bulk_equilibrium(p)
    assert model.a0/model.p.b == pytest.approx(IDEAL_H,rel=1e-14)
    assert model.embedding.rho_ref == original.embedding.rho_ref
    assert abs(model.grad(model.a0,0)[0])<1e-8
    assert abs(stretch-1)>1e-5
    c=cubic_constants_gpa(relaxed_bulk_observables(model),
                          volume_scale=stretch**3)
    assert abs(c["normal_stress_GPa"])<1e-7
    assert c["C11_GPa"]>c["C12_GPa"]>0 and c["C44_GPa"]>0


def test_positive_local_hessian_does_not_certify_positive_cleavage_work():
    from solver_v1.fcc111_active_interface import FCC111ActiveInterface
    p=FullFCCParameters(.8466141,.8069657,1.36765,10.1664764,12.)
    bulk,_=isotropic_bulk_equilibrium(p)
    interface=FCC111ActiveInterface(bulk,tolerance=2e-10)
    assert np.min(np.linalg.eigvalsh(interface.hessian(interface.h,0)))>0
    # This is a model-failure regression, not permission to use this candidate.
    assert interface.separated_limit()<0


def test_opening_diagnostic_does_not_miss_a_narrow_metastable_peak():
    from solver_v1.fcc111_active_interface import FCC111ActiveInterface
    from solver_v1.run_full_fcc_calibration_audit import fixed_registry_opening_barriers
    p=FullFCCParameters(.8466141,.8069657,1.36765,10.1664764,12.)
    bulk,_=isotropic_bulk_equilibrium(p)
    interface=FCC111ActiveInterface(bulk,tolerance=2e-10)
    # Despite negative *global* separation work, the positive intact curvature
    # produces a small local metastable barrier. A global unimodal search misses it.
    rows=fixed_registry_opening_barriers("failure_probe",interface,1.,1.)
    assert len(rows)==4
    for row in rows:
        assert row["fixed_s_spinodal_force_eV_per_reduced_a"]>0
        assert row["opening_barrier_eV_cell"]>0
        assert row["bound_a_reduced_L0"]<row["saddle_a_reduced_L0"]
        assert abs(interface.hessian(row["fixed_s_spinodal_a_reduced_L0"],0)[0,0])<1e-6


@pytest.mark.parametrize("decay", (0.,-1.,math.nan))
def test_invalid_decay_rejected(decay):
    with pytest.raises(ValueError):
        bulk_basis(decay)
