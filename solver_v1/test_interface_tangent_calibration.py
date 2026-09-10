"""The v20 constrained-calibration math is not an Al material acceptance test."""
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from .interface_tangent_calibration import (
    TangentCalibrationProblem, impose_tangents, equality_sign_feasibility,
    make_residual_table, tangent_summary, shape_selection_status,
)
from .vector_material_calibration import MaterialObservation
from .tail_constrained_material import convex_profile


def test_imposed_tangents_preserve_reference_and_reject_implicit_controls():
    original = [MaterialObservation('perfect_Haa', 20., 2., 'eV/L0^2', 'fit'),
                MaterialObservation('perfect_Hxx', 4., .4, 'eV/L0^2', 'fit')]
    changed = impose_tangents(original, {'perfect_Haa': 25.})
    assert original[0].target == 20. and original[0].role == 'fit'
    assert changed[0].target == 25. and changed[0].role == 'exact'
    assert changed[1] is original[1]
    for bad in ({'unexplained': 1.}, {'perfect_Haa': -1.}, {'perfect_Hxx': float('nan')}):
        with pytest.raises(ValueError):
            impose_tangents(original, bad)


def test_lp_separates_sign_infeasibility_from_spectral_claim():
    obs = [MaterialObservation('exact', 1., 1., 'test', 'exact'),
           MaterialObservation('other', 0., 1., 'test', 'fit')]
    M = np.array([[1., 0.], [0., 1.]])
    feasible = equality_sign_feasibility(M, obs, nonnegative=(0, 1))
    assert feasible['success'] and feasible['exact_rank'] == 1
    assert not feasible['spectral_constraints_included']
    bad = [replace(obs[0], target=-1), obs[1]]
    impossible = equality_sign_feasibility(M, bad, nonnegative=(0, 1))
    assert not impossible['success'] and impossible['status'] == 2
    assert equality_sign_feasibility(M, bad, nonnegative=(1,))['success']


def test_imposed_and_source_residuals_are_distinct_in_exports():
    source = [MaterialObservation('perfect_Haa', 2., .1, 'eV/L0^2', 'fit'),
              MaterialObservation('other', 3., 1., 'eV', 'fit')]
    obs = impose_tangents(source, {'perfect_Haa': 4.})
    fit = convex_profile(np.eye(2), obs, nonnegative=(0, 1))
    rows = make_residual_table(source, obs, fit)
    assert rows[0]['imposed_target'] == 4. and rows[0]['source_target'] == 2.
    assert abs(rows[0]['imposed_normalized_residual']) < 1e-12
    assert rows[0]['source_normalized_residual'] == pytest.approx(20.)
    assert rows[0]['units'] == 'eV/L0^2'


def test_convex_control_value_function_and_heldout_exclusion():
    # E fixes c0; h fixes c1. Other loss favors c1=3. Distinct imposed
    # tangents yield the known convex quadratic, not an altered target scale.
    M = np.array([[1., 0., 0.], [0., 1., 0.], [0., 1., 1.], [0., 0., 1.], [0., 1., 2.]])
    source = [MaterialObservation('bulk', 1., 1., 'test', 'exact'),
              MaterialObservation('perfect_Haa', 2., 1., 'test', 'fit'),
              MaterialObservation('fit', 3., 1., 'test', 'fit'),
              MaterialObservation('fit2', 0., 1., 'test', 'fit'),
              MaterialObservation('heldout', 99., 1., 'test', 'heldout')]
    losses = []
    for t in (1., 2., 3.):
        obs = impose_tangents(source, {'perfect_Haa': t})
        result = convex_profile(M, obs, nonnegative=(0, 1))
        losses.append(result['squared_loss'])
        altered = obs[:-1]+[replace(obs[-1], target=-1e8)]
        other = convex_profile(M, altered, nonnegative=(0, 1))
        np.testing.assert_array_equal(result['coefficients'], other['coefficients'])
    np.testing.assert_allclose(losses, [2., .5, 0.], atol=1e-12)
    assert losses[1] <= (losses[0]+losses[2])/2


@pytest.fixture(scope='module')
def actual_problem():
    root = Path(__file__).resolve().parents[1]
    return TangentCalibrationProblem(root/'results/fcc111_active_interface/coordination_screening_v19/joint_refinement')


def test_actual_source_units_and_predeclared_validation(actual_problem):
    p = actual_problem
    assert len(p.development) == 104 and len(p.other_rows) == 102
    assert all(o.units == 'GPa' for o in p.observations[2:5])
    np.testing.assert_allclose([o.target for o in p.observations[2:5]], [114., 62., 32.])
    fresh = [o for o in p.observations if o.name.startswith('v20_new_')]
    assert len(fresh) == 48 and all(o.role == 'heldout' for o in fresh)
    assert not {o.state for o in fresh} & {o.state for o in p.raw_observations[:-48]}


def test_actual_old_candidate_replay_and_no_energy_extension(actual_problem):
    p = actual_problem; parent = p.parents['old_family']; shape = tuple(parent['shape'])
    prediction = p.matrix(shape)@np.array(parent['coefficients'])
    residual = (prediction-np.array([o.target for o in p.observations]))/np.array([o.scale for o in p.observations])
    assert residual[p.development]@residual[p.development] == pytest.approx(parent['source_loss'], rel=1e-12)
    assert max(abs(residual[:5])) < 1e-10
    assert not p.definition['energy_family_changed']


def test_unverified_profile_has_no_fake_success_values():
    row = tangent_summary('bad', 'toy', dict(completed=False, classification='spectral_profile_not_verified'), None)
    assert row['source_loss_104'] is None and row['Haa'] is None
    assert not row['material_accepted']


def test_actual_two_tangent_constraints_are_verified_not_material_adoption(actual_problem):
    p = actual_problem
    result = p.fit(p.parents['old_family']['shape'], p.source_tangents)
    assert result['completed'] and result['feasibility']['success']
    assert result['feasibility']['exact_rank'] == 7
    assert result['feasibility']['null_dimension'] == 3
    for name, target in p.source_tangents.items():
        assert result['predictions'][p.tangent_rows[name]] == pytest.approx(target, abs=2e-10)
    assert len(result['selected_rows']) == 102
    assert result['exact_residual'] < 1e-8 and result['kkt_residual'] < 1e-6
    assert result['minimum_robust_margin'] > 0
    assert not result['strictly_positive_LJ'] and not result['material_accepted']
    assert result['source_loss_104'] > p.parents['old_family']['source_loss']


def test_saved_fixed_shape_controls_keep_convex_tradeoff_and_source_residuals():
    import csv
    root = Path(__file__).resolve().parents[1]
    with (root/'results/fcc111_active_interface/tangent_calibration_v20/controls/summary.csv').open(
            encoding='utf-8', newline='') as stream:
        rows = list(csv.DictReader(stream))
    for family in ('old_family', 'power_family'):
        selected = {r['case']: r for r in rows if r['family'] == family}
        names = ['normal_scan_0', 'normal_scan_0.25', 'normal_scan_0.5', 'normal_scan_0.75', 'both_exact']
        losses = np.array([float(selected[n]['other_loss_102']) for n in names])
        assert np.all(np.diff(losses, n=2) > 0)
        assert selected['baseline']['positive_LJ'] == 'True'
        assert selected['both_exact']['positive_LJ'] == 'False'
        assert all(selected[n]['material_accepted'] == 'False' for n in names)


def test_positive_trial_is_not_mislabeled_converged_optimizer_endpoint():
    best = dict(shape=[2., 3., 4., 5., 6., 0.], strictly_positive_LJ=True)
    data = dict(best=best, optimizer=dict(success=True,
        last_accepted_coordinates=np.log([2., 3., 4., 5., 7.])))
    state = shape_selection_status(data)
    assert state['positive_profile_found']
    assert not state['selected_profile_is_optimizer_endpoint']
    assert not state['positive_material_shape_optimizer_converged']
    data['optimizer']['last_accepted_coordinates'] = np.log(best['shape'][:5])
    assert shape_selection_status(data)['positive_material_shape_optimizer_converged']
    data['optimizer']['success'] = False
    assert not shape_selection_status(data)['positive_material_shape_optimizer_converged']
