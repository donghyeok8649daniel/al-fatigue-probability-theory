'''Protocol tests only; actual trajectory validation is separate.'''
import pytest
from .run_kinetic_record_refinement import study_plan


def test_record_length_and_cutoff_are_independent_axes():
    plan = study_plan(8192)
    fixed = [p for p in plan if p[0].startswith('prefix') and p[3] == 126]
    assert [p[2] for p in fixed] == [1024, 2048, 4096, 8192]
    assert [p[3] for p in plan if p[0].startswith('prefix_8192_')] == [126, 254, 510, 1022]
    assert all((stop-start)//blocks >= 2*(lag+1) for _, start, stop, lag, blocks in plan)


def test_quarter_windows_are_disjoint():
    windows = [p for p in study_plan(8192) if p[0].startswith('quarter')]
    assert [(p[1], p[2]) for p in windows] == [(0, 2048), (2048, 4096), (4096, 6144), (6144, 8192)]


def test_old_record_keeps_original_full_record_comparison():
    assert study_plan(1024) == [('prefix_1024_lag_126', 0, 1024, 126, 4)]


@pytest.mark.parametrize('frames', [0, 1023, 1024.5])
def test_protocol_refuses_short_or_fractional_records(frames):
    with pytest.raises(ValueError):
        study_plan(frames)
