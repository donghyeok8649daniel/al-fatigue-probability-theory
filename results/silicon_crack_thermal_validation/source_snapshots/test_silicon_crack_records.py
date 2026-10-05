import numpy as np
import pytest
from scipy.linalg import expm

from solver_v1.silicon_initiation_probability import initiation_collectors, evolve_first_passage
from solver_v1.silicon_crack_records import joint_crack_records, evolve_joint_records


def generator(size, edges):
    matrix = np.zeros((size, size))
    for source, destination, rate in edges:
        matrix[destination, source] += rate
        matrix[source, source] -= rate
    return matrix


def test_joint_histories_preserve_the_unmodified_healing_dynamics():
    # A intact, C open crack, D closed crack, B independently qualified crack.
    g = generator(4, [(0,1,2), (1,0,3), (1,2,4), (2,1,1),
                      (2,0,2), (2,3,.7), (3,0,5)])
    model = joint_crack_records(g, formation_states=[1,2,3], qualified_states=[3])
    times = [0., .02, .3, 1., 4.]
    result = evolve_joint_records(model, [1,0,0,0], times)
    for index, time in enumerate(times):
        np.testing.assert_allclose(result['physical_mass'][index],
                                   expm(time*g)@np.array([1,0,0,0]), atol=2e-14)
    assert np.max(abs(result['mass_residual'])) < 2e-14
    assert result['minimum_mass'] >= 0
    assert np.all(np.diff(result['formation_cumulative']) >= 0)
    assert np.all(np.diff(result['qualified_cumulative']) >= 0)
    assert np.all(result['qualified_cumulative'] <= result['formation_cumulative'])
    # Qualification is a history, even if the underlying crack later heals.
    assert result['qualified_cumulative'][-1] > result['physical_mass'][-1,3]
    assert result['formation_cumulative'][-1] > result['present_crack_probability'][-1]


def test_each_marginal_matches_a_separate_first_passage_problem():
    g = generator(3, [(0,1,2), (1,0,3), (1,2,.4), (2,1,9)])
    result = evolve_joint_records(joint_crack_records(g, formation_states=[1,2],
                                  qualified_states=[2]), [1,0,0], [0,.1,1,10])
    formed = evolve_first_passage(initiation_collectors(g, {'formed':[1,2]}), [1], [0,.1,1,10])
    qualified = evolve_first_passage(initiation_collectors(g, {'qualified':[2]}), [1,0], [0,.1,1,10])
    np.testing.assert_allclose(result['formation_cumulative'], formed['cause_cumulative'][:,0], atol=1e-14)
    np.testing.assert_allclose(result['qualified_cumulative'], qualified['cause_cumulative'][:,0], atol=1e-14)
    np.testing.assert_allclose(result['formed_not_yet_qualified'],
                               result['formation_cumulative']-result['qualified_cumulative'], atol=1e-14)


def test_closed_crack_is_formation_without_requiring_an_open_gap():
    g = generator(3, [(0,1,1), (1,0,2), (1,2,.3)])
    # State 1 explicitly denotes a closed but unhealed crack, not intact material.
    result = evolve_joint_records(joint_crack_records(g, formation_states=[1,2],
                                  qualified_states=[2]), [0,1,0], [0,.1])
    assert result['formation_cumulative'][0] == 1.
    assert result['qualified_cumulative'][0] == 0.
    assert result['present_crack_probability'][1] < 1.
    assert result['formation_cumulative'][1] == pytest.approx(1.)


def test_qualification_already_present_is_recorded_at_initial_time():
    model = joint_crack_records(np.zeros((3,3)), formation_states=[1,2], qualified_states=[2])
    result = evolve_joint_records(model, [.2,.3,.5], [0,2])
    np.testing.assert_allclose(result['formation_cumulative'], [.8,.8])
    np.testing.assert_allclose(result['qualified_cumulative'], [.5,.5])
    assert result['physical_clock'] is None
    assert model['material_labels_verified'] is False


def test_direct_qualification_counts_formation_once():
    g = generator(3, [(0,2,.7), (2,0,2)])
    result = evolve_joint_records(joint_crack_records(g, formation_states=[1,2],
                                  qualified_states=[2]), [1,0,0], [0,.5,2])
    expected = 1-np.exp(-.7*np.array([0,.5,2]))
    np.testing.assert_allclose(result['formation_cumulative'], expected, atol=2e-15)
    np.testing.assert_allclose(result['qualified_cumulative'], expected, atol=2e-15)


def test_first_event_flux_is_the_cdf_derivative():
    g = generator(3, [(0,1,2), (1,0,3), (1,2,.4), (2,1,9)])
    model = joint_crack_records(g, formation_states=[1,2], qualified_states=[2])
    state = np.zeros(9); state[0] = 1
    derivative = (model['generator']@expm(.3*model['generator'].toarray())@state).reshape(3,3)
    expected = [derivative[1:].sum(), derivative[2].sum()]
    got = evolve_joint_records(model, [1,0,0], [.3])
    np.testing.assert_allclose(got['first_event_flux'][0], expected, atol=1e-14)


@pytest.mark.parametrize('formation,qualified', [([1],[2]), ([1,1],[1]), ([1],[1.]), ([],[1])])
def test_invalid_or_non_nested_labels_are_rejected(formation, qualified):
    with pytest.raises(ValueError):
        joint_crack_records(np.zeros((3,3)), formation_states=formation, qualified_states=qualified)
