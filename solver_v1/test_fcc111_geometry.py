import math

import numpy as np
import pytest

from solver_v1.fcc111_geometry import (
    DIRECT_110,
    SHOCKLEY_112,
    fcc111_geometry,
    fcc111_geometry_from_b,
    full_atomic_distance_squared,
    homogeneous_layer_registry,
    in_plane_lattice_norm_squared,
    registry_path,
)


def test_fcc111_basis_and_surface_cell_geometry_are_exact():
    geometry = fcc111_geometry(4.05)
    basis = np.stack((geometry.e1, geometry.e2, geometry.e3))
    np.testing.assert_allclose(basis @ basis.T, np.eye(3), atol=2.0e-16)
    np.testing.assert_allclose(np.cross(geometry.e1, geometry.e2), geometry.e3)
    assert geometry.b == pytest.approx(4.05 / math.sqrt(2.0))
    assert geometry.h111 == pytest.approx(4.05 / math.sqrt(3.0))
    assert np.linalg.norm(geometry.a1) == pytest.approx(geometry.b)
    assert np.linalg.norm(geometry.a2) == pytest.approx(geometry.b)
    assert geometry.a1 @ geometry.a2 == pytest.approx(0.5 * geometry.b**2)
    assert geometry.atomic_cell_area == pytest.approx(
        math.sqrt(3.0) * 4.05**2 / 4.0
    )


def test_reciprocal_basis_satisfies_duality_and_first_shell_has_six_vectors():
    geometry = fcc111_geometry_from_b(1.0)
    direct = (geometry.a1, geometry.a2)
    reciprocal = (geometry.reciprocal_b1, geometry.reciprocal_b2)
    for i, ai in enumerate(direct):
        for j, bj in enumerate(reciprocal):
            assert ai @ bj == pytest.approx(2.0 * math.pi * (i == j), abs=2e-15)
    first_shell = {
        tuple(np.round(geometry.reciprocal_vector(h, k), 13))
        for h in range(-2, 3)
        for k in range(-2, 3)
        if h * h - h * k + k * k == 1
    }
    assert len(first_shell) == 6


def test_abc_translation_repeats_after_three_planes_modulo_lattice():
    geometry = fcc111_geometry_from_b(1.0)
    np.testing.assert_allclose(3.0 * geometry.tau, geometry.a1 + geometry.a2)
    np.testing.assert_allclose(geometry.abc_shift(0), np.zeros(2))
    np.testing.assert_allclose(geometry.abc_shift(1), geometry.tau)
    np.testing.assert_allclose(geometry.abc_shift(2), 2.0 * geometry.tau)
    np.testing.assert_allclose(geometry.abc_shift(3), np.zeros(2))
    np.testing.assert_allclose(
        homogeneous_layer_registry(3, 0.0, geometry=geometry), np.zeros(2), atol=1e-15
    )


def test_positive_abc_actual_cubic_orientation_not_silently_negative_abc():
    geometry = fcc111_geometry(4.05)
    basis = geometry.plane_basis_in_stacked_cubic_axes()
    assert np.linalg.det(basis) == pytest.approx(1.)
    for layer in range(-3, 4):
        for m, n in ((0, 0), (-1, 1), (2, -3)):
            xy = geometry.lattice_vector(m,n)+geometry.abc_shift(layer)
            local = np.r_[xy,layer*geometry.h111]
            cubic = local@basis/(geometry.lattice_constant/2)
            np.testing.assert_allclose(cubic,np.rint(cubic),atol=3e-15)
            assert int(np.rint(cubic).sum()) % 2 == 0
    # Demonstrate the previously tempting mapping is a different orientation.
    local = np.r_[geometry.tau,geometry.h111]
    wrong = local@np.stack((geometry.e1,geometry.e2,geometry.e3))/(geometry.lattice_constant/2)
    assert np.max(abs(wrong-np.rint(wrong))) > .3


def test_named_registry_paths_are_explicit_and_have_distinct_periods():
    direct = registry_path(DIRECT_110)
    shockley = registry_path(SHOCKLEY_112)
    assert np.linalg.norm(direct.direction()) == pytest.approx(1.0)
    assert np.linalg.norm(shockley.direction()) == pytest.approx(1.0)
    assert direct.period_over_b == pytest.approx(1.0)
    assert shockley.period_over_b == pytest.approx(math.sqrt(3.0))
    assert shockley.partial_increment_over_b == pytest.approx(1.0 / math.sqrt(3.0))
    assert "held-out" in shockley.calibration_status


def test_direct_lattice_norm_distance_and_self_exclusion():
    geometry = fcc111_geometry_from_b(1.0)
    for m, n in ((1, 0), (0, 1), (1, 1), (-2, 1)):
        vector = geometry.lattice_vector(m, n)
        assert vector @ vector == pytest.approx(
            in_plane_lattice_norm_squared(m, n, geometry.b)
        )
    a = math.sqrt(2.0 / 3.0)
    distance = full_atomic_distance_squared(
        0, 0, 1, a, 0.0, geometry=geometry
    )
    assert distance == pytest.approx(1.0)  # ideal FCC nearest neighbour
    with pytest.raises(ValueError, match="self term"):
        full_atomic_distance_squared(0, 0, 0, a, 0.0, geometry=geometry)


@pytest.mark.parametrize("path_id", (DIRECT_110, SHOCKLEY_112))
def test_homogeneous_registry_is_periodic_along_each_named_path(path_id):
    geometry = fcc111_geometry_from_b(1.0)
    period = registry_path(path_id).period_over_b * geometry.b
    for layer in range(1, 7):
        left = homogeneous_layer_registry(
            layer, 0.137, geometry=geometry, path_id=path_id
        )
        right = homogeneous_layer_registry(
            layer, 0.137 + period, geometry=geometry, path_id=path_id
        )
        direct = np.column_stack((geometry.a1, geometry.a2))
        lattice_coordinates = np.linalg.solve(direct, right - left)
        np.testing.assert_allclose(
            lattice_coordinates, np.rint(lattice_coordinates), atol=3.0e-15
        )
