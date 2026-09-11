import numpy as np
import pytest

from .core_interface_compatibility import minimax_compatibility, interval_profile
from .vector_material_calibration import MaterialObservation as O


def test_exact_minimax_dual_for_conflicting_curvatures():
    # c0=1, y1=c1, y2=c1; targets -1,+1 cannot both be matched.
    M = np.array([[1., 0.], [0., 1.], [0., 1.]])
    obs = [O('anchor', 1., 1., 'x', 'exact'),
           O('saddle', -1., .1, 'curvature', 'fit'),
           O('opening', 1., .1, 'curvature', 'heldout')]
    r = minimax_compatibility(M, obs, [1, 2], nonnegative=[0])
    assert r['completed'] and not r['target_box_feasible']
    assert r['minimax_normalized_error'] == pytest.approx(10.)
    assert r['dual_lower_bound'] == pytest.approx(10.)
    assert r['duality_gap'] == pytest.approx(0., abs=1e-11)
    assert obs[2].role == 'heldout'  # diagnosis does not relabel input records


def test_sign_constraint_is_distinct_from_general_linear_family():
    M = np.eye(2)
    obs = [O('anchor', 2., .1, 'x', 'exact'), O('target', -3., .2, 'x', 'fit')]
    signed = minimax_compatibility(M, obs, [1], nonnegative=[0, 1])
    unrestricted = minimax_compatibility(M, obs, [1])
    assert signed['minimax_normalized_error'] == pytest.approx(15.)
    assert unrestricted['minimax_normalized_error'] < 1e-10


def test_units_rescale_without_changing_compatibility():
    M = np.array([[1., 0.], [0., 1.], [0., 1.]])
    obs = [O('eq', 1., 1., 'a', 'exact'), O('j1', -2., .2, 'b', 'fit'),
           O('j2', 4., .4, 'b', 'fit')]
    old = minimax_compatibility(M, obs, [1, 2])
    factors = np.array([1e9, 1e-19, 1e-19])
    scaled = [O(o.name, o.target*f, o.scale*f, o.units, o.role) for o, f in zip(obs, factors)]
    new = minimax_compatibility(M*factors[:, None], scaled, [1, 2])
    assert new['minimax_normalized_error'] == pytest.approx(old['minimax_normalized_error'])
    np.testing.assert_allclose(new['coefficients'], old['coefficients'], atol=2e-13)


def test_inconsistent_anchors_are_not_a_material_certificate():
    M = np.array([[1.], [1.], [1.]])
    obs = [O('eq1', 1., 1., 'x', 'exact'), O('eq2', 2., 1., 'x', 'exact'),
           O('test', 0., 1., 'x', 'fit')]
    r = minimax_compatibility(M, obs, [2])
    assert not r['completed'] and not r['certified_lower_bound_available']


@pytest.mark.parametrize('scale', [1., 1e-5, 1e5])
def test_affine_interval_profile_exact_boundary_and_objective(scale):
    M = np.array([[1., 0.], [0., scale], [0., scale]])
    obs = [O('eq', 2., 1., 'x', 'exact'), O('fit', 0., scale, 'x', 'fit'),
           O('box', 3*scale, scale, 'x', 'heldout')]
    fit = interval_profile(M, obs, [2], 1., nonnegative=(0, 1))
    np.testing.assert_allclose(fit['coefficients'], [2., 2.], atol=1e-10)
    assert fit['squared_loss'] == pytest.approx(4.)
    assert fit['bookkeeping_constant'] == pytest.approx(1.)
    assert fit['selected_rows'] == [1]
    assert not fit['energy_term_added'] and obs[2].role == 'heldout'


def test_affine_interval_and_sampled_stability_are_both_retained():
    M = np.array([[1., 0.], [0., 1.], [0., 1.]])
    obs = [O('eq', 2., 1., 'x', 'exact'), O('fit', -2., 1., 'x', 'fit'),
           O('box', 1., 1., 'x', 'heldout')]
    ops = np.array([[np.eye(3), np.eye(3)]])
    fit = interval_profile(M, obs, [2], 1., nonnegative=(0,),
                           operators=ops, tails=np.zeros((1, 2)))
    np.testing.assert_allclose(fit['coefficients'], [2., 0.], atol=1e-10)
    assert fit['minimum_robust_margin'] == pytest.approx(2.)
    assert not fit['whole_zone_proved']


def test_invalid_rows_and_widths_rejected():
    obs = [O('eq', 1., 1., 'x', 'exact'), O('fit', 1., 1., 'x', 'fit')]
    with pytest.raises(ValueError):
        minimax_compatibility(np.eye(2), obs, [1, 1])
    with pytest.raises(ValueError):
        interval_profile(np.eye(2), obs, [1], -1., nonnegative=(0, 1))


def test_interval_scale_is_independent_of_objective_block_normalization():
    M = np.array([[1., 0.], [0., 1.], [0., 1.]])
    obs = [O('eq', 2., 1., 'x', 'exact'), O('fit', 0., 50., 'x', 'fit'),
           O('box', 3., 50., 'x', 'fit')]
    fit = interval_profile(M, obs, [2], 1., nonnegative=(0, 1), interval_scales=[.1])
    assert fit['predictions'][2] == pytest.approx(2.9)
    assert fit['squared_loss'] == pytest.approx((2.9**2+.1**2)/50.**2)


def test_recorded_physical_design_recomputes_two_jet_feasibility_and_full_conflict():
    import csv
    from pathlib import Path
    from .current_core_coefficient_basis import NAMES
    from .interface_tangent_calibration import NONNEGATIVE
    path = Path(__file__).resolve().parents[1]/'results/core_interface_compatibility_v23/fixed_shape_audit/candidate_1_design.csv'
    with path.open(encoding='utf8', newline='') as stream:
        rows = list(csv.DictReader(stream))
    M = np.array([[float(r[k]) for k in NAMES] for r in rows])
    obs = [O(r['observable'], float(r['target']), float(r['scale']), r['units'], r['role']) for r in rows]
    pair = [next(i for i, o in enumerate(obs) if o.name == name)
            for name in ('saddle_Hxx', 'v19_new_4_Haa')]
    all_jets = [i for i, o in enumerate(obs) if o.units == 'eV/L0^2' and o.role != 'exact']
    feasible = minimax_compatibility(M, obs, pair, nonnegative=NONNEGATIVE)
    conflict = minimax_compatibility(M, obs, all_jets, nonnegative=NONNEGATIVE)
    assert feasible['target_box_feasible'] and feasible['minimax_normalized_error'] < 1e-9
    assert conflict['exact_rank'] == 7 and len(all_jets) == 115
    assert conflict['minimax_normalized_error'] == pytest.approx(18.4644467070752, rel=2e-10)
    assert not conflict['target_box_feasible'] and not conflict['whole_family_impossibility_proved']


def test_selected_bessel_jets_keep_full_bulk_unit_conversion():
    import csv
    import json
    from pathlib import Path
    from types import SimpleNamespace
    from .core_interface_compatibility import ExistingJetSubset
    from .current_core_coefficient_basis import NAMES
    root = Path(__file__).resolve().parents[1]
    definition = json.loads((root/'results/fcc111_active_interface/coordination_screening_v19/joint_refinement/definition.json').read_bytes())
    raw = [O(**dict(o, state=tuple(o['state']) if o['state'] is not None else None,
                   jet_weights=tuple(o['jet_weights']) if o['jet_weights'] is not None else None))
           for o in definition['observations']]
    critical = [next(i for i, o in enumerate(raw) if o.name == name)
                for name in ('perfect_Haa', 'saddle_Hxx', 'v19_new_4_Haa')]
    subset = ExistingJetSubset(SimpleNamespace(raw_observations=raw, observations=raw), critical)
    snapshot = json.loads((root/'results/current_material_core_v22/wider_probe_validation/research_candidate_snapshot.json').read_bytes())
    actual = subset.matrix(tuple(snapshot['best']['shape']))
    with (root/'results/core_interface_compatibility_v23/fixed_shape_audit/candidate_1_design.csv').open(
            encoding='utf8', newline='') as stream:
        by_name = {r['observable']: r for r in csv.DictReader(stream)}
    names = [raw[i].name if i not in (2, 3, 4) else ('C11_GPa', 'C12_GPa', 'C44_GPa')[i-2]
             for i in subset.rows]
    expected = np.array([[float(by_name[name][k]) for k in NAMES] for name in names])
    np.testing.assert_allclose(actual, expected, atol=2e-10, rtol=2e-12)
