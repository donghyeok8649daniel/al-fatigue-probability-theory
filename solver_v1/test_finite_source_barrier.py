import numpy as np
import pytest

from .finite_source_barrier import pinned_line_barrier,normal_mode_second_variation
from .finite_source_reference import line_energy_coefficient,pinned_source_branch
from .fcc111_geometry import fcc111_geometry_from_b
from .nonlocal_interface_elasticity import cubic_elastic_tensor,rotate_elastic_tensor


@pytest.fixture(scope='module')
def constant_tension():
    # Deterministic analytic geometry test, NOT an Al parameter source.
    mu=25e9;b=3e-10
    return line_energy_coefficient(cubic_elastic_tensor(2*mu,0,mu),[b,0,0]),mu,b


def test_constant_tension_major_and_minor_circular_arc_barrier(constant_tension):
    line,mu,b=constant_tension;L=2e-6;log=7.
    Gamma=mu*b*b/(4*np.pi)*log;critical=2*Gamma/(b*L)
    for fraction in (.2,.6,.9):
        theta=np.arcsin(fraction);p=fraction*critical*b;R=Gamma/p
        result=pinned_line_barrier(line,span_m=L,outer_log_ratio=log,shear_Pa=fraction*critical)
        expected=Gamma**2/p*(np.pi-2*theta-np.sin(2*theta))
        area=R**2*(np.pi-2*theta+np.sin(2*theta))
        assert result['outer_only_barrier_J']==pytest.approx(expected,rel=2e-10)
        assert result['additional_swept_area_m2']==pytest.approx(area,rel=2e-10)
        assert result['outer_only_stress_derivative_m3']==pytest.approx(b*area,rel=2e-10)
        assert not result['finite_source_atomistically_validated']
        assert result['kinetic_rate'] is None and result['production_probability'] is None


def test_anisotropic_stress_derivative_is_envelope_not_chosen_volume():
    rotation=fcc111_geometry_from_b(1).plane_basis_in_stacked_cubic_axes()
    tensor=rotate_elastic_tensor(cubic_elastic_tensor(114e9,62e9,32e9),rotation)
    b=2.86e-10;line=line_energy_coefficient(tensor,[b,0,0])
    options=dict(span_m=1e-6,outer_log_ratio=8.)
    critical=pinned_source_branch(line,np.pi/2,**options)['critical_shear_outer_only_Pa']
    tau=.63*critical
    result=pinned_line_barrier(line,shear_Pa=tau,**options)
    errors=[]
    for step in (critical*2e-4,critical*1e-4):
        plus=pinned_line_barrier(line,shear_Pa=tau+step,**options)
        minus=pinned_line_barrier(line,shear_Pa=tau-step,**options)
        derivative=-(plus['outer_only_barrier_J']-minus['outer_only_barrier_J'])/(2*step)
        errors.append(abs(derivative/result['outer_only_stress_derivative_m3']-1))
    assert errors[-1]<1e-7 and errors[-1]<.3*errors[0]
    assert abs(result['total_energy_identity_residual_J'])<1e-9*result['outer_only_barrier_J']


def test_finite_geometric_scaling_and_fold_exponents(constant_tension):
    line,mu,b=constant_tension;log=7.;L=1e-6
    Gamma=mu*b*b/(4*np.pi)*log;critical=2*Gamma/(b*L)
    a=pinned_line_barrier(line,span_m=L,outer_log_ratio=log,shear_Pa=.7*critical)
    larger=pinned_line_barrier(line,span_m=2*L,outer_log_ratio=log,shear_Pa=.35*critical)
    assert larger['outer_only_barrier_J']==pytest.approx(2*a['outer_only_barrier_J'],rel=2e-10)
    assert larger['outer_only_stress_derivative_m3']==pytest.approx(4*a['outer_only_stress_derivative_m3'],rel=2e-10)
    values=[pinned_line_barrier(line,span_m=L,outer_log_ratio=log,shear_Pa=(1-e)*critical)
            for e in (1e-3,5e-4)]
    assert values[1]['outer_only_barrier_J']/values[0]['outer_only_barrier_J']==pytest.approx(2**-1.5,rel=.003)
    assert values[1]['outer_only_stress_derivative_m3']/values[0]['outer_only_stress_derivative_m3']==pytest.approx(2**-.5,rel=.003)
    # The floating-point Schur/Fourier coefficient need not reproduce the
    # symbolic circle's fold bit-for-bit. Test the implementation's OWN fold.
    computed_fold=pinned_source_branch(line,np.pi/2,span_m=L,outer_log_ratio=log)['critical_shear_outer_only_Pa']
    for invalid in (0.,-1.,computed_fold,2*critical):
        with pytest.raises(ValueError):pinned_line_barrier(line,span_m=L,outer_log_ratio=log,shear_Pa=invalid)


def test_saddle_index_from_independent_perturbed_curve_energy():
    rotation=fcc111_geometry_from_b(1).plane_basis_in_stacked_cubic_axes()
    tensor=rotate_elastic_tensor(cubic_elastic_tensor(114e9,62e9,32e9),rotation)
    line=line_energy_coefficient(tensor,[2.86e-10,0,0]);log=8.;L=1e-6
    for endpoint in (1.,np.pi-1.):
        g=lambda t,d=0:line.evaluate(t,d)*log
        Q=lambda t:g(t)*np.sin(t)+g(t,1)*np.cos(t)
        r=lambda t:-g(t)*np.cos(t)+g(t,1)*np.sin(t)
        p=2*Q(endpoint)/L
        for mode in (1,2):
            exact=normal_mode_second_variation(endpoint_angle=endpoint,pressure_J_m2=p,mode=mode)
            errors=[]
            for segments in (256,512):
                theta=np.linspace(-endpoint,endpoint,segments+1)
                coordinates=np.column_stack([Q(theta)/p,(r(endpoint)-r(theta))/p])
                normal=np.column_stack([np.sin(theta),np.cos(theta)])
                shape=np.sin(mode*np.pi*(theta+endpoint)/(2*endpoint))
                amplitude=1e-4*L
                def energy(amp):
                    v=coordinates+amp*shape[:,None]*normal
                    dv=np.diff(v,axis=0);angle=np.arctan2(dv[:,1],dv[:,0])
                    area=np.sum((v[1:,1]+v[:-1,1])*dv[:,0]/2)
                    return np.sum(g(angle)*np.linalg.norm(dv,axis=1))-p*area
                numerical=(energy(amplitude)+energy(-amplitude)-2*energy(0.))/amplitude**2
                errors.append(abs(numerical/exact-1))
            assert errors[-1]<2e-4 and errors[-1]<.35*errors[0]
            if endpoint>np.pi/2 and mode==1:assert exact<0
            else:assert exact>0
