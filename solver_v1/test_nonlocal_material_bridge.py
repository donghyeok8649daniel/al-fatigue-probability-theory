"""Bridge to the actual saved infinite-LJ/Bessel research energy, not a fit."""
import numpy as np
import pytest

from solver_v1.nonlocal_interface_elasticity import cubic_elastic_tensor, rotate_elastic_tensor, halfspace_impedance
from solver_v1.run_low_stress_cyclic_diagnostic import build_surface
from solver_v1.run_nonlocal_interface_reference import load_prepared
from solver_v1.static_bulk_stability import StaticBulkHessian


@pytest.fixture(scope="module")
def actual():
    return build_surface(tolerance=2e-11)


def test_periodic_series_replays_actual_analytic_forces_not_only_energy(actual):
    surface,units,meta=actual
    series,_,saved=load_prepared()
    assert meta["parameter_sha256"]==saved["parameter_sha256"]
    area=units.atomic_cell_area_m2/units.length_scale_m**2
    for s in (.0361,.1742,.2947,.4961):
        exact=surface.packed(surface.h,s)
        for order,index,tol in ((0,0,1e-11),(1,2,2e-10),(2,5,2e-8)):
            assert series.evaluate(s,order)==pytest.approx(exact[index]/area,abs=tol)


def test_actual_bessel_moduli_match_independent_finite_q_with_correct_cubic_frame(actual):
    surface,units,_=actual
    _,_,meta=load_prepared()
    g=surface.interface.bulk.geometry
    c=cubic_elastic_tensor(*(meta["elastic_constants_GPa"][f"C{k}_GPa"]*1e9 for k in (11,12,44)))
    correct=rotate_elastic_tensor(c,g.plane_basis_in_stacked_cubic_axes())
    wrong=rotate_elastic_tensor(c,[g.e1,g.e2,g.e3])
    direct=StaticBulkHessian(surface.interface.bulk,cutoff=20.,D3=surface.amplitude_ev,
        D1=surface.vector_amplitude_ev,D2=surface.quadrupole_amplitude_ev,angular_decay=surface.angular.kappa)
    q=np.array([.005,0.,0.])
    factor=g.atomic_cell_area*surface.h*meta["elastic_conversion_Pa_to_eV_L0cubed"]
    actual_q=direct.evaluate(q)["matrix"]
    expected=factor*np.einsum("ijkl,j,l->ik",correct,q,q)
    erroneous=factor*np.einsum("ijkl,j,l->ik",wrong,q,q)
    assert np.linalg.norm(actual_q-expected)/np.linalg.norm(expected)<8e-5
    assert np.linalg.norm(actual_q-erroneous)/np.linalg.norm(expected)>.3


def test_elastic_unit_change_is_exact_and_no_new_energy_scale_is_fitted():
    c=cubic_elastic_tensor(110.,60.,35.)
    base=halfspace_impedance(c,[.8,.6]).jump_per_wave_number
    np.testing.assert_allclose(halfspace_impedance(c*1e9,[.8,.6]).jump_per_wave_number,
                               base*1e9,rtol=3e-13,atol=1e-3)
