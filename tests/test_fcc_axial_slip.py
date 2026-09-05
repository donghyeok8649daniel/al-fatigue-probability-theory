import numpy as np

from theory.fcc_axial_slip import (
    fcc_axial_slip_systems,
    fcc_slip_system_indices,
    maximum_schmid_factor,
)


def test_fcc_system_generator_has_twelve_valid_unique_systems():
    systems = fcc_slip_system_indices()
    assert len(systems) == 12
    assert len(set(systems)) == 12
    for plane, direction in systems:
        assert sorted(abs(value) for value in plane) == [1, 1, 1]
        assert sorted(abs(value) for value in direction) == [0, 1, 1]
        assert np.dot(plane, direction) == 0


def test_001_axial_loading_has_known_fcc_maximum_schmid_factor():
    systems = fcc_axial_slip_systems(0, 0, 1)
    assert len(systems) == 12
    assert np.isclose(systems[0].schmid_factor, 1.0 / np.sqrt(6.0))
    assert np.isclose(maximum_schmid_factor(0, 0, 1), 1.0 / np.sqrt(6.0))


def test_resolved_shear_retains_sign_and_scales_with_axial_stress():
    system = fcc_axial_slip_systems(1, 2, 3)[0]
    assert np.isclose(system.resolved_shear_mpa(12.0), 12.0 * system.signed_schmid_factor)
    assert np.isclose(system.resolved_shear_mpa(-12.0), -system.resolved_shear_mpa(12.0))
