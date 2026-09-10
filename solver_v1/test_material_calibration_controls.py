import json

import numpy as np
import pytest

from .material_calibration_controls import append_fixed_pair, select_profile_result
from .tail_constrained_material import convex_profile
from .vector_material_calibration import MaterialObservation


def test_pair_control_is_exact_and_never_a_new_fitted_target():
    matrix = np.eye(4)
    rows = [MaterialObservation(str(i), float(i+1), 1., 'test', 'fit') for i in range(4)]
    control = dict(coefficients=[.4, .8], calibration_sha256='synthetic_test_only')
    M, obs = append_fixed_pair(matrix, rows, control)
    assert len(rows) == 4 and len(obs) == 6
    assert [o.role for o in obs[-2:]] == ['exact', 'exact']
    first = convex_profile(M, obs, nonnegative=(0, 1, 2, 3))
    np.testing.assert_allclose(first['coefficients'], [.4, .8, 3., 4.], atol=1e-14)
    assert first['selected_rows'] == [0, 1, 2, 3]
    for invalid in ([0., 1.], [1., -1.], [1., np.nan]):
        with pytest.raises(ValueError, match='positive'):
            append_fixed_pair(matrix, rows, dict(control, coefficients=invalid))
    with pytest.raises(ValueError, match='binding'):
        append_fixed_pair(matrix, rows, dict(coefficients=[1., 1.]))


def test_best_closure_does_not_supplant_positive_LJ_candidate():
    physical = dict(squared_loss=20., strictly_positive_LJ=True)
    closure = dict(squared_loss=.01, strictly_positive_LJ=False)
    assert select_profile_result([closure, physical]) == (physical, True)
    assert select_profile_result([closure]) == (closure, False)
    with pytest.raises(ValueError):
        select_profile_result([])


def test_no_admissible_result_is_preserved_but_cannot_load_as_material(tmp_path):
    from .run_tail_constrained_calibration import finalize_profile_run
    from .run_vector_registry_audit import save_json
    from .validate_tail_calibration import load_material
    row = dict(squared_loss=0., strictly_positive_LJ=False, predictions=[1.],
               residuals=[0.], decays=[3., 5., 7.], coefficients=[0., 1.])
    obs = [MaterialObservation('test', 1., 1., 'test', 'fit')]
    finalize_profile_run(tmp_path, [row], [], [], obs, elapsed_seconds=0., source_sha256='test')
    saved = json.loads((tmp_path/'calibration.json').read_text())
    assert saved['completed'] and not saved['admissible_solution_found']
    assert saved['best_is_diagnostic_closure_only'] and not saved['material_accepted']
    save_json(tmp_path/'definition.json', {})
    with pytest.raises(ValueError, match='closure'):
        load_material(tmp_path)


def test_deterministic_search_rejects_infeasible_trial_without_finite_penalty():
    from .material_calibration_controls import feasible_simplex_search

    def residual(x):
        if x[0] > .031:
            raise ArithmeticError('synthetic infeasible exact-constraint domain')
        return np.array([x[0]-.03, 2*(x[1]-.07)])

    options = dict(initial=[0., 0.], bounds=[[-1., -1.], [1., 1.]], max_evaluations=180)
    first = feasible_simplex_search(residual, **options)
    second = feasible_simplex_search(residual, **options)
    assert first.success and first.rejected_profiles
    assert first.fun < 1e-12 and first.x[0] <= .031
    np.testing.assert_array_equal(first.x, second.x)
    assert first.fun == second.fun


def test_actual_negative_quartic_cancellation_is_resolved_on_exact_face():
    from pathlib import Path
    path = (Path(__file__).resolve().parents[1]/'results/fcc111_active_interface/'
            'even_environment_v16/boundary_reproduction/actual_qp.npz')
    with np.load(path, allow_pickle=False) as data:
        M, target, scale = data['matrix'], data['target'], data['scale']
        obs = [MaterialObservation(str(i), t, s, 'original_scale', str(r))
               for i, (t, s, r) in enumerate(zip(target, scale, data['role']))]
        original = data['returned_coefficients']
        assert original[9] < 0  # exact untouched pre-fix reproduction
        result = convex_profile(M, obs, data['inequalities'],
                                nonnegative=tuple(data['nonnegative']))
    assert result['coefficients'][9] == 0.
    assert 9 in result['exact_nonnegative_boundary_resolve']
    assert result['kkt_residual'] < 1e-6
    assert result['exact_residual'] < 1e-8
    np.testing.assert_allclose(result['predictions'], M@original, rtol=2e-10, atol=2e-9)
