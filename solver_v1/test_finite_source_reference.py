import numpy as np
import pytest

from .fcc111_geometry import fcc111_geometry_from_b
from .finite_source_reference import line_energy_coefficient, pinned_source_branch, solve_pinned_graph, swept_area_to_shear_strain
from .nonlocal_interface_elasticity import cubic_elastic_tensor, rotate_elastic_tensor


@pytest.fixture(scope='module')
def source_tensor():
    # Hypothetical/static-material tests, never actual specimen/kinetic calibration.
    return rotate_elastic_tensor(cubic_elastic_tensor(114e9,62e9,32e9),
        fcc111_geometry_from_b(1).plane_basis_in_stacked_cubic_axes())


def test_isotropic_exact_screw_edge_and_full_angle():
    mu=30e9; nu=.3; lam=2*mu*nu/(1-2*nu); b=2.86e-10
    line=line_energy_coefficient(cubic_elastic_tensor(lam+2*mu,lam,mu),[b,0,0])
    t=np.linspace(-np.pi,np.pi,91)
    exact=mu*b*b/(4*np.pi)*(np.cos(t)**2+np.sin(t)**2/(1-nu))
    np.testing.assert_allclose(line.evaluate(t),exact,rtol=7e-14,atol=1e-24)
    expected=mu*b*b/(4*np.pi)*(np.cos(t)**2+np.sin(t)**2/(1-nu)+
        2*(1/(1-nu)-1)*np.cos(2*t))
    # Second spectral derivatives amplify roundoff in the Schur samples by N^2.
    floor=4*np.finfo(float).eps*line.samples**2*max(abs(exact))
    np.testing.assert_allclose(line.stiffness(t),expected,rtol=0,atol=floor)


def test_anisotropic_angular_refinement_and_burgers_scaling(source_tensor):
    b=np.array([2.86e-10,0.,0.]); line=line_energy_coefficient(source_tensor,b,samples=32)
    fine=line_energy_coefficient(source_tensor,b,samples=64)
    doubled=line_energy_coefficient(source_tensor,2*b,samples=64)
    assert b.flags.writeable  # no mutation of caller ownership
    t=np.linspace(-1.49,1.31,57)
    np.testing.assert_allclose(line.evaluate(t),fine.evaluate(t),rtol=0,atol=1e-21)
    np.testing.assert_allclose(line.stiffness(t),fine.stiffness(t),rtol=0,atol=3e-19)
    np.testing.assert_allclose(doubled.evaluate(t),4*fine.evaluate(t),rtol=3e-14)


def test_constant_tension_circle_independent_formula():
    # Isotropic nu=0 makes the line energy orientation-independent.
    mu=25e9; b=3e-10
    line=line_energy_coefficient(cubic_elastic_tensor(2*mu,0,mu),[b,0,0])
    theta=1.1; L=2e-6; logratio=7.
    result=pinned_source_branch(line,theta,span_m=L,outer_log_ratio=logratio)
    radius=L/(2*np.sin(theta)); gamma=mu*b*b/(4*np.pi)*logratio
    np.testing.assert_allclose(result['x_m'],radius*np.sin(result['theta']),atol=2e-19)
    np.testing.assert_allclose(result['y_m'],radius*(np.cos(result['theta'])-np.cos(theta)),atol=2e-19)
    assert result['applied_shear_Pa']==pytest.approx(gamma/(b*radius),rel=1e-13)
    assert result['experimental_yield_prediction'] is None


def test_independent_discrete_force_and_shape_converge(source_tensor):
    line=line_energy_coefficient(source_tensor,[2.86e-10,0,0])
    analytic=pinned_source_branch(line,1.,span_m=1e-6,outer_log_ratio=8.)
    errors=[]
    for n in (32,64,128):
        solve=solve_pinned_graph(line,span_m=1e-6,outer_log_ratio=8.,
            shear_Pa=analytic['applied_shear_Pa'],segments=n)
        assert solve['converged'] and solve['dimensionless_force_residual']<2e-10
        errors.append(abs(solve['y_m'].max()-analytic['y_m'].max()))
        assert not solve['physical_time_available']
    assert errors[1]<.3*errors[0] and errors[2]<.3*errors[1]


def test_scale_law_is_geometry_not_fitted_strength_factor(source_tensor):
    line=line_energy_coefficient(source_tensor,[3e-10,0,0])
    first=pinned_source_branch(line,1.,span_m=1e-6,outer_log_ratio=8.)
    longer=pinned_source_branch(line,1.,span_m=10e-6,outer_log_ratio=8.)
    assert longer['applied_shear_Pa']==pytest.approx(first['applied_shear_Pa']/10)
    np.testing.assert_allclose(longer['y_m'],10*first['y_m'],atol=1e-20)
    assert not first['core_energy_included'] and not first['finite_source_atomistically_validated']


def test_swept_area_is_not_zero_point_two_percent_yield_by_itself():
    b=3e-10; area=4e-13; volume=8e-12
    strain=swept_area_to_shear_strain(area,b,volume)
    assert strain==pytest.approx(1.5e-11) and strain<.002
    assert swept_area_to_shear_strain(2*area,b,volume)==2*strain


def test_at_or_above_fold_graph_not_certified(source_tensor):
    line=line_energy_coefficient(source_tensor,[3e-10,0,0])
    critical=pinned_source_branch(line,np.pi/2,span_m=1e-6,outer_log_ratio=8.)['critical_shear_outer_only_Pa']
    with pytest.raises(ValueError):
        solve_pinned_graph(line,span_m=1e-6,outer_log_ratio=8.,shear_Pa=critical)


@pytest.mark.parametrize('b',[[0,0,0],[1e-10,0,1e-10],[np.nan,0,0]])
def test_invalid_burgers_refused(source_tensor,b):
    with pytest.raises(ValueError): line_energy_coefficient(source_tensor,b)
