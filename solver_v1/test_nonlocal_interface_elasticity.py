"""Independent half-space elimination; no empirical plasticity or kinetics."""
import numpy as np
import pytest
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

from solver_v1.fcc111_geometry import fcc111_geometry_from_b
from solver_v1.nonlocal_interface_elasticity import (
    cubic_elastic_tensor, rotate_elastic_tensor, halfspace_impedance,
    interface_kernel, screw_elastic_factor,
)


def plane_tensor(c11=110., c12=60., c44=35.):
    g = fcc111_geometry_from_b(1.)
    return rotate_elastic_tensor(cubic_elastic_tensor(c11, c12, c44), [g.e1, g.e2, g.e3])


@pytest.mark.parametrize("angle", [0., .37, 1.2])
def test_isotropic_degenerate_modes_have_exact_jump_kernel(angle):
    lam, mu = 40., 30.
    nu = lam/(2*(lam+mu))
    t = np.array([np.cos(angle), np.sin(angle)])
    c = cubic_elastic_tensor(lam+2*mu, lam, mu)
    result = halfspace_impedance(c, t)
    wave = np.r_[t, 0.]; screw = np.array([-t[1], t[0], 0.]); normal = np.array([0., 0., 1.])
    exact = mu/2*np.outer(screw, screw) + mu/(2*(1-nu))*(np.outer(wave, wave)+np.outer(normal, normal))
    np.testing.assert_allclose(result.jump_per_wave_number, exact, atol=4e-12)
    assert result.hermitian_residual < 2e-13
    assert result.riccati_residual < 2e-13


def test_anisotropic_screw_and_wave_scaling_no_qzero_double_count():
    c = plane_tensor()
    z = halfspace_impedance(c, [0, 1])
    np.testing.assert_allclose(z.jump_per_wave_number[0, 0], screw_elastic_factor(c)/2, rtol=2e-13)
    np.testing.assert_allclose(interface_kernel(c, [0, -.6]), interface_kernel(c, [0, .3]).conj()*2, atol=1e-12)
    np.testing.assert_array_equal(interface_kernel(c, [0, 0]), np.zeros((3, 3)))
    assert np.linalg.eigvalsh(z.jump_per_wave_number)[0] > 0


def finite_depth_schur(c, direction, elements):
    """Independent linear-element energy quadrature, z in [0,16], |q|=1."""
    t = np.r_[np.asarray(direction)/np.linalg.norm(direction), 0.]
    a = np.einsum("ijkl,j,l->ik", c, t, t)
    b = np.einsum("ikl,l->ik", c[:, 2, :, :], t)
    d = c[:, 2, :, 2]
    dz = 16/elements
    matrix = lil_matrix((3*(elements+1),)*2, dtype=complex)
    for cell in range(elements):
        block = np.zeros((6, 6), dtype=complex)
        for gauss in (-1/np.sqrt(3), 1/np.sqrt(3)):
            n = np.array([(1-gauss)/2, (1+gauss)/2]); dn = np.array([-1, 1])/dz
            for i in range(2):
                for j in range(2):
                    block[3*i:3*i+3, 3*j:3*j+3] += dz/2*(a*n[i]*n[j]+d*dn[i]*dn[j]
                        +1j*(b*dn[i]*n[j]-b.T*n[i]*dn[j]))
        start = 3*cell
        matrix[start:start+6, start:start+6] += block
    matrix = matrix[:-3, :-3].tocsc()  # independent fixed far boundary
    return matrix[:3, :3].toarray()-matrix[:3, 3:]@spsolve(matrix[3:, 3:], matrix[3:, :3].toarray())


def test_halfspace_matches_independent_bulk_energy_minimization():
    c = plane_tensor()/100
    t = [.8, .6]
    exact = halfspace_impedance(c, t).upper
    errors = [np.max(abs(finite_depth_schur(c, t, n)-exact)) for n in (64, 128, 256)]
    assert errors[1] < .28*errors[0]
    assert errors[2] < .28*errors[1]
    assert errors[-1] < .0006


def test_common_translation_is_eliminated_by_compliance_sum():
    z = halfspace_impedance(plane_tensor(), [.7, .2])
    jump = np.array([.2+.1j, -.3, .4])
    force = z.jump_per_wave_number@jump
    upper = np.linalg.solve(z.upper, force)
    lower = -np.linalg.solve(z.lower, force)
    np.testing.assert_allclose(upper-lower, jump, atol=1e-14)
    before = (np.vdot(upper, z.upper@upper)+np.vdot(lower, z.lower@lower)).real/2
    after = np.vdot(jump, z.jump_per_wave_number@jump).real/2
    assert before == pytest.approx(after, rel=1e-14)


def test_invalid_geometry_or_unstable_moduli_are_not_repaired():
    with pytest.raises(ValueError):
        cubic_elastic_tensor(40, 60, 30)
    with pytest.raises(ValueError):
        rotate_elastic_tensor(cubic_elastic_tensor(100, 40, 30), np.ones((3, 3)))
    with pytest.raises(ValueError):
        halfspace_impedance(plane_tensor(), [0, 0])
    unstable=cubic_elastic_tensor(3,1,1)-1.8*np.einsum("ij,kl->ijkl",np.eye(3),np.eye(3))
    with pytest.raises(ValueError,match="bulk elastic energy"):
        halfspace_impedance(unstable,[1,0])
