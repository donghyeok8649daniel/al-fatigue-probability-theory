"""Independent real-row / real-3D / finite-layer checks, no parameter fitting."""
from itertools import product

import numpy as np
import pytest

from solver_v1.angular_environment_reference import traceless_second, traceless_third
from solver_v1.discrete_fcc_screw import (
    exponential_row_transform_derivatives, row_exponential_slip_derivatives,
    row_lj_slip_curvature, FCCScrewRowHessian,
)
from solver_v1.run_low_stress_cyclic_diagnostic import build_surface
from solver_v1.run_nonlocal_interface_reference import load_prepared
from solver_v1.static_bulk_stability import StaticBulkHessian
from solver_v1.nonlocal_interface_elasticity import cubic_elastic_tensor, rotate_elastic_tensor


def direct_row(y, z, delta, *, b, kappa, amplitude=1., n=192):
    x = b*np.arange(-n, n+1)+delta
    rvec = np.column_stack((x, np.full(x.shape, y), np.full(x.shape, z)))
    r = np.linalg.norm(rvec, axis=1)
    direction = rvec/r[:, None]
    weight = amplitude*np.exp(-kappa*r)
    xx = np.sum(weight*(kappa**2*direction[:, 0]**2-kappa/r*(1-direction[:, 0]**2)))
    tensors = {}
    for rank in (1, 2, 3):
        raw = np.empty((3,)*rank)
        for indices in product(range(3), repeat=rank):
            factors = np.prod(rvec[:, indices], axis=1)
            derivative = -kappa*direction[:, 0]*factors
            for slot, index in enumerate(indices):
                if index == 0:
                    remaining = indices[:slot]+indices[slot+1:]
                    derivative += np.prod(rvec[:, remaining], axis=1)
            raw[indices] = np.sum(weight*derivative)
        tensors[rank] = raw if rank == 1 else traceless_second(raw) if rank == 2 else traceless_third(raw)
    return xx, tensors


@pytest.mark.parametrize("y,z,delta,b,kappa", [(.7,.8,.173,1.,1.3),
    (.2,.47,.41,.83,4.9), (1.7,-.4,-.217,1.13,.8), (.866,0.,.5,1.,3.)])
def test_infinite_exponential_row_tensor_derivatives_against_direct(y,z,delta,b,kappa):
    exact = row_exponential_slip_derivatives(y,z,delta,b=b,kappa=kappa,amplitude=1.7)
    direct, tensors = direct_row(y,z,delta,b=b,kappa=kappa,amplitude=1.7)
    assert exact["density_xx"] == pytest.approx(direct, abs=4e-13, rel=2e-11)
    for rank in (1,2,3):
        np.testing.assert_allclose(exact["moments"][rank],tensors[rank],atol=2e-13,rtol=2e-10)
    assert exact["last_mode_absolute_envelope"] < 1e-13


def test_transform_high_derivatives_and_phase_periodicity():
    radius,k,kappa = .813,1.3,2.4
    h = 2e-5
    center = exponential_row_transform_derivatives(radius,k,kappa)
    plus = exponential_row_transform_derivatives(radius,k+h,kappa)
    minus = exponential_row_transform_derivatives(radius,k-h,kappa)
    for order in range(3):
        assert (plus[order]-minus[order])/(2*h) == pytest.approx(center[order+1],rel=2e-8,abs=1e-11)
    x = row_exponential_slip_derivatives(.3,.8,.271,b=1.,kappa=2.)
    y = row_exponential_slip_derivatives(.3,.8,1.271,b=1.,kappa=2.)
    for rank in (1,2,3):
        np.testing.assert_allclose(x["moments"][rank],y["moments"][rank],atol=2e-15)


def test_lj_row_curvature_and_self_exclusion():
    for radius,delta in ((.51,.23),(.82,.5),(1.7,0.)):
        x = np.arange(-256,257)+delta
        r = np.sqrt(x*x+radius*radius)
        ep,sigma = .137,.82
        prime = 4*ep*(-12*sigma**12/r**13+6*sigma**6/r**7)
        second = 4*ep*(156*sigma**12/r**14-42*sigma**6/r**8)
        direct = sum(second*x*x/r**2+prime/r*(1-x*x/r**2))
        actual,_,_ = row_lj_slip_curvature(radius,delta,b=1.,epsilon=ep,sigma=sigma)
        assert actual == pytest.approx(direct,rel=3e-12,abs=3e-11)
    with pytest.raises(ValueError):
        row_lj_slip_curvature(0.,0.,b=1.,epsilon=1.,sigma=1.)


@pytest.fixture(scope="module")
def actual():
    surface,units,meta = build_surface(tolerance=2e-11)
    return surface,units,meta,FCCScrewRowHessian.from_surface(surface)


def test_row_geometry_is_the_existing_fcc_not_a_new_reduced_lattice(actual):
    surface,_,_,rows = actual
    g = surface.interface.bulk.geometry
    for j,l in ((-3,-2),(-1,1),(0,3),(2,4)):
        # A row's origin can differ by an integer a1 translation.
        true = j*g.a2+l*g.tau
        assert true[1] == pytest.approx(rows.d*(j+l/3))
        assert true[0] == pytest.approx((j+l)*rows.b/2)
    assert not np.any((rows.j==0)&(rows.l==0))


def test_rigid_step_reconstructs_existing_full_interface_hessian(actual):
    surface,_,_,rows = actual
    reconstructed = rows.rigid_step_stiffness()
    assert reconstructed["total"] == pytest.approx(surface.packed(surface.h,0.)[5],abs=8e-11)
    assert reconstructed["angular"] != 0
    # Per-site sum before squaring: not sum of separate plane embedding terms.
    assert rows.diagnostics.omitted_validation_ring_envelope < 1e-11


def test_symbol_matches_independent_three_dimensional_hessian(actual):
    surface,_,_,rows = actual
    errors = []
    for cutoff in (8.,32.):
        direct = StaticBulkHessian(surface.interface.bulk,cutoff=cutoff,
            D1=surface.vector_amplitude_ev,D2=surface.quadrupole_amplitude_ev,
            D3=surface.amplitude_ev,angular_decay=surface.angular.kappa)
        discrepancy = []
        for q in ([0,.4,.6],[0,1.2,.3],[0,0,1.5],[0,.01,.005]):
            result = direct.evaluate(q)
            exact = rows.symbol(q[1],q[2]*rows.h+q[1]*rows.d/3)
            discrepancy.append(abs(exact-result["matrix"][0,0]))
            assert result["scalar_density"][0,0] < 1e-24
            np.testing.assert_allclose(result["matrix"][0,1:],0.,atol=1e-12)
        errors.append(max(discrepancy))
    assert errors[1] < errors[0]/5
    assert errors[1] < 1e-6


def test_transverse_refinement_translation_and_bloch_symmetry(actual):
    surface,_,_,rows = actual
    q = np.linspace(-3.,3.,17)
    theta = np.linspace(.021,2.7,17)
    fine = FCCScrewRowHessian.from_surface(surface,fixed_ring=10)
    coarse = FCCScrewRowHessian.from_surface(surface,fixed_ring=3)
    intermediate = FCCScrewRowHessian.from_surface(surface,fixed_ring=5)
    exact = fine.symbol(q,theta)
    assert max(abs(coarse.symbol(q,theta)-exact)) <= coarse.diagnostics.omitted_validation_ring_envelope
    assert max(abs(intermediate.symbol(q,theta)-exact)) <= intermediate.diagnostics.omitted_validation_ring_envelope
    assert max(abs(intermediate.symbol(q,theta)-exact)) < max(abs(coarse.symbol(q,theta)-exact))/100
    np.testing.assert_allclose(rows.symbol(q,theta),exact,atol=2e-11,rtol=1e-12)
    np.testing.assert_allclose(rows.symbol(q,theta),rows.symbol(-q,-theta),atol=2e-13)
    np.testing.assert_allclose(rows.symbol(q,theta),rows.symbol(q,theta+2*np.pi),atol=4e-13)
    assert rows.symbol(0.,0.) == 0.


def test_discrete_schur_matches_independent_finite_layer_matrix(actual):
    _,_,_,rows = actual
    n,q = 64,.8
    eye = np.eye(n,dtype=complex)
    matrix = np.zeros_like(eye)
    operators = {rank:np.zeros((3**rank,n,n),complex) for rank in (1,2,3)}
    for index,(j,l) in enumerate(zip(rows.j,rows.l)):
        shift = np.exp(1j*q*rows.d*j)*np.roll(eye,-int(l),axis=0)
        matrix += rows.effective_pair[index]*(eye-(shift+shift.conj().T)/2)
        for rank in (1,2,3):
            operators[rank] += rows.moments[rank][index,:,None,None]*(shift-eye)
    for rank in (1,2,3):
        matrix += 2*rows.amplitudes[rank]*sum(op.conj().T@op for op in operators[rank])
    np.testing.assert_allclose(matrix,matrix.conj().T,atol=2e-13)
    jump = np.zeros(n); jump[1]=1; jump[0]=-1
    finite = 1/(jump@np.linalg.solve(matrix,jump)).real
    assert finite == pytest.approx(rows.jump_stiffness(q)["stiffness"],rel=2e-11)


def test_discrete_jump_resolution_and_constraint_are_explicit(actual):
    _,_,_,rows = actual
    values = [rows.jump_stiffness(.005,layer_phase_points=n)["stiffness"] for n in (1024,4096,8192)]
    assert abs(values[2]-values[1]) < abs(values[1]-values[0])/100
    rigid = rows.rigid_step_stiffness()["total"]
    zero = rows.jump_stiffness(0)["stiffness"]
    assert 0 < zero < rigid
    assert rows.jump_stiffness(0,constraint="spectral_same_y")["stiffness"] == pytest.approx(zero)
    assert abs(rows.jump_stiffness(.3)["stiffness"]-
               rows.jump_stiffness(.3,constraint="spectral_same_y")["stiffness"]) > .1
    with pytest.raises(ValueError):
        rows.jump_stiffness(.1,constraint="guessed")


def test_acoustic_pole_derives_the_long_wave_discrete_slope(actual):
    surface,_,_,rows = actual
    _,_,meta = load_prepared()
    constants = meta["elastic_constants_GPa"]
    tensor = cubic_elastic_tensor(*(constants[f"C{k}_GPa"]*1e9 for k in (11,12,44)))
    tensor = rotate_elastic_tensor(tensor,surface.interface.bulk.geometry.plane_basis_in_stacked_cubic_axes())
    tensor *= meta["elastic_conversion_Pa_to_eV_L0cubed"]
    zero = rows.jump_stiffness(0)["stiffness"]
    for constraint in ("logical_rows","spectral_same_y"):
        predicted = rows.long_wave_jump_slope(tensor,zero,constraint=constraint)["slope"]
        q=.002
        actual_slope = (rows.jump_stiffness(q,layer_phase_points=8192,constraint=constraint)["stiffness"]-zero)/q
        assert actual_slope == pytest.approx(predicted,rel=.006)


def test_exhausted_reciprocal_controls_do_not_silently_accept():
    with pytest.raises(ArithmeticError):
        row_exponential_slip_derivatives(.1,.2,.19,b=1.,kappa=2.,max_modes=2)
    with pytest.raises(ArithmeticError):
        row_lj_slip_curvature(.2,.19,b=1.,epsilon=1.,sigma=1.,max_modes=2)


def test_automatic_quadrature_certification_requires_actual_refinement(actual):
    _,_,_,rows = actual
    result = rows.converged_jump_stiffness(.3)
    assert result["quadrature_converged"]
    assert len(result["refinement_history"]) >= 3
    with pytest.raises(ArithmeticError):
        rows.converged_jump_stiffness(.001,initial_points=32,max_points=128)
