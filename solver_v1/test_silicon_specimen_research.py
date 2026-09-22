"""Physical invariants for the optional Si specimen-to-atom boundary bridge."""
import unittest

import numpy as np
from numpy.testing import assert_allclose

from .silicon_specimen_research import (
    AnisotropicModeI, SI_111_FRAME, atomic_crack_boundary, cubic_tensor,
    specimen_K, specimen_stress,
)
from .silicon_crystal_boundary import with_internal_relaxation


class SpecimenBridgeTests(unittest.TestCase):
    def test_specimen_units_and_inverse(self):
        # 100 MPa at a 10 um half-crack, infinite-plate Y=1.
        k = specimen_K(100., 10e-6, geometry_factor=1.)
        self.assertAlmostEqual(k, .5604991216397929)
        self.assertAlmostEqual(specimen_stress(k, 10e-6, geometry_factor=1.), 100.)
        assert_allclose(specimen_K(100., [1e-6, 4e-6], geometry_factor=1.), [k/np.sqrt(10), 2*k/np.sqrt(10)])

    def test_invalid_geometry_is_not_silently_defaulted(self):
        for value in [0., -1., np.nan, np.inf]:
            with self.assertRaises(ValueError):
                specimen_K(100., value, geometry_factor=1.)
            with self.assertRaises(ValueError):
                specimen_stress(.7, 1e-6, geometry_factor=value)
        with self.assertRaises(TypeError):
            specimen_K(100., 1e-6)

    def test_elastic_energy_is_rotation_invariant(self):
        c = cubic_tensor(165.7, 63.9, 79.6, np.eye(3))
        rotated = cubic_tensor(165.7, 63.9, 79.6)
        strain = np.array([[.012, .003, -.007], [.003, -.005, .002], [-.007, .002, .009]])
        local = SI_111_FRAME@strain@SI_111_FRAME.T
        self.assertAlmostEqual(np.einsum('ij,ijkl,kl', strain, c, strain),
                               np.einsum('ij,ijkl,kl', local, rotated, local), places=14)

    def test_bad_stiffness_and_coupled_orientation_are_rejected(self):
        with self.assertRaises(ValueError):
            AnisotropicModeI(100., 110., 50.)
        angle = .31
        rot = np.array([[1, 0, 0], [0, np.cos(angle), -np.sin(angle)], [0, np.sin(angle), np.cos(angle)]])
        with self.assertRaises(ValueError):
            AnisotropicModeI(165.7, 63.9, 79.6, frame=rot@SI_111_FRAME)

    def test_ahead_stress_normalization_and_face_tractions(self):
        for stiffness in [(151.42449, 76.42212, 56.44934), (142.51243, 75.36411, 69.01294), (180., 60., 60.)]:
            with self.subTest(stiffness=stiffness):
                field = AnisotropicModeI(*stiffness)
                _, _, stress = field.field([[12., 0.], [-12., 0.], [-12., -0.]], .7)
                self.assertAlmostEqual(stress[0, 1], 70/np.sqrt(2*np.pi*12), places=10)
                self.assertAlmostEqual(stress[0, 2], 0., places=11)
                assert_allclose(stress[1:, 1:], 0, atol=1e-11)

    def test_displacement_derivative_and_equilibrium(self):
        field = AnisotropicModeI(165.7, 63.9, 79.6)
        xy = np.array([[4., 3.], [-8., -4.], [30., -9.]])
        _, grad, _ = field.field(xy, .65)
        h = 1e-4
        derivative = []
        dstress = []
        for direction in np.eye(2):
            up, _, sp = field.field(xy+h*direction, .65)
            um, _, sm = field.field(xy-h*direction, .65)
            derivative.append((up-um)/(2*h))
            dstress.append((sp-sm)/(2*h))
        assert_allclose(np.stack(derivative, axis=-1), grad, atol=2e-10, rtol=2e-9)
        assert_allclose(dstress[0][:, 0]+dstress[1][:, 2], 0., atol=2e-9)
        assert_allclose(dstress[0][:, 2]+dstress[1][:, 1], 0., atol=2e-9)

    def test_isotropic_griffith_limit(self):
        young, nu = 170., .23
        mu = young/(2*(1+nu))
        lam = young*nu/((1+nu)*(1-2*nu))
        field = AnisotropicModeI(lam+2*mu, lam, mu)
        self.assertTrue(field.isotropic)
        self.assertAlmostEqual(field.H_GPa_inv, (1-nu**2)/young, places=14)
        k = float(field.griffith_K(2.8))
        self.assertAlmostEqual(k, np.sqrt(2.8*young/(1-nu**2)/1000), places=14)
        self.assertAlmostEqual(float(field.energy_release_J_m2(k)), 2.8, places=13)

    def test_j_integral_and_crack_opening(self):
        for constants in [(151.42449, 76.42212, 56.44934), (142.51243, 75.36411, 69.01294), (180., 60., 60.)]:
            with self.subTest(constants=constants):
                field = AnisotropicModeI(*constants)
                predicted = float(field.energy_release_J_m2(.65))
                for radius in [8., 32., 128.]:
                    self.assertAlmostEqual(field.contour_J(.65, radius_A=radius, points=1024), predicted, places=10)
                u = field.field([[-20., 0.], [-20., -0.]], .65)[0]
                self.assertAlmostEqual(u[0, 1]-u[1, 1], 8*field.H_GPa_inv*65*np.sqrt(20/(2*np.pi)), places=11)

    def test_load_and_length_scaling(self):
        field = AnisotropicModeI(165.7, 63.9, 79.6)
        points = np.array([[2., 3.], [-4., -5.]])
        u, grad, stress = field.field(points, .4)
        u2, grad2, stress2 = field.field(4*points, .8)
        assert_allclose(u2, 4*u, atol=1e-12)
        assert_allclose(grad2, grad, atol=1e-12)
        assert_allclose(stress2, stress, atol=1e-12)
        self.assertAlmostEqual(field.energy_release_J_m2(.8)/field.energy_release_J_m2(.4), 4.)

    def test_boundary_is_actual_lattice_with_independent_free_atoms(self):
        boundary = atomic_crack_boundary(5.43095, radius_A=28., grip_width_A=8., front_repeats=4)
        r = boundary.reference
        self.assertTrue(np.all(np.linalg.norm(r[:, :2], axis=1) < 28.))
        self.assertTrue(np.any(boundary.fixed) and np.any(~boundary.fixed))
        assert_allclose(r[boundary.bonds[:, 1], 1]-r[boundary.bonds[:, 0], 1],
                        5.43095*np.sqrt(3)/4, atol=1e-8)
        field = AnisotropicModeI(151.42449, 76.42212, 56.44934)
        assert_allclose(boundary.displaced(field, 0.), r, atol=1e-12)
        loaded = boundary.displaced(field, .65)
        assert_allclose(loaded[:, 2], r[:, 2], atol=0)
        self.assertTrue(np.all(loaded[r[:, 1] > 0, 1] > r[r[:, 1] > 0, 1]))

    def test_optical_boundary_preserves_zero_and_linearity(self):
        b = atomic_crack_boundary(5.43095, radius_A=28., front_repeats=4)
        response = np.zeros((3, 6))
        response[np.arange(3), [3, 4, 5]] = -.6
        corrected = with_internal_relaxation(b, 5.43095, response)
        # Every crossing shuffle pair connects distinct diamond sublattices.
        assert_allclose(corrected.sublattice_sign[b.bonds[:, 0]], -1.)
        assert_allclose(corrected.sublattice_sign[b.bonds[:, 1]], 1.)
        field = AnisotropicModeI(151.42449, 76.42212, 56.44934)
        assert_allclose(corrected.displaced(field, 0.), b.reference, atol=1e-12)
        d1 = corrected.displaced(field, .4)-b.reference
        d2 = corrected.displaced(field, .8)-b.reference
        assert_allclose(d2, 2*d1, atol=1e-12)
        self.assertGreater(np.max(abs(corrected.displaced(field, .4)-b.displaced(field, .4))), 1e-3)
        with self.assertRaises(ValueError):
            with_internal_relaxation(b, 5.1, response)

    def test_optical_response_uses_cubic_not_crack_axes(self):
        b = atomic_crack_boundary(5.43095, radius_A=28., front_repeats=4)
        response = np.zeros((3, 6))
        response[0, 0] = 2.
        corrected = with_internal_relaxation(b, 5.43095, response)
        # A constant local xx strain transforms to eps_cubic_00=1/6.
        class Uniform:
            def field(self, xy, k):
                grad = np.zeros((len(xy), 2, 2))
                grad[:, 0, 0] = k
                return np.zeros_like(xy), grad, np.zeros((len(xy), 3))
        shift = corrected.displaced(Uniform(), 1.)-b.reference
        expected = corrected.sublattice_sign[:, None]/6*SI_111_FRAME[:, 0]
        assert_allclose(shift, expected, atol=3e-15)


if __name__ == '__main__':
    unittest.main()
