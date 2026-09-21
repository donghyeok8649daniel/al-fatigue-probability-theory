import numpy as np
from solver_v1.silicon_memory_reduction import (
    positive_memory_quadrature,quadrature_memory,quadrature_release,
)
from solver_v1.silicon_conditional_research import (
    AMU_EV_PS2_A2,conditional_harmonic,harmonic_release,
)


def test_positive_quadrature_preserves_moments_and_full_rank_kernel():
    values=np.array([.2,.8,1.7,3.2,5.])
    weights=np.array([.6,.1,.3,.7,.2])
    reduced=positive_memory_quadrature(values,weights,3)
    for power in range(6):
        np.testing.assert_allclose(reduced['weights_eV_A2']@reduced['nodes_eV_A2']**power,
            weights@values**power,rtol=3e-14)
    full=positive_memory_quadrature(values,weights,5)
    times=np.linspace(0,10,301)
    np.testing.assert_allclose(quadrature_memory(full['nodes_eV_A2'],full['weights_eV_A2'],times),
        quadrature_memory(values,weights,times),atol=1e-12)
    assert reduced['orthogonality_error']<5e-15


def test_oscillator_embedding_has_exact_static_curvature_and_cartesian_response():
    h=np.array([[4.,1.,.4],[1.,2.,0.],[.4,0.,3.]])
    b=np.eye(3)[:,1:];d=np.eye(3)[:,:1]
    reduced=conditional_harmonic(h,b,d)
    nodes=np.array([2.,3.]);weights=np.array([1/2,.16/3])
    mass_amu=1/AMU_EV_PS2_A2
    times=np.linspace(0,3,101)
    expected=harmonic_release(h,reduced.relaxed_lift[:,0],[1.,0.,0.],times,mass_amu=mass_amu)
    actual=quadrature_release(nodes,weights,reduced.relaxed_curvature[0,0],1.,times,mass_amu=mass_amu)
    np.testing.assert_allclose(actual,expected,atol=4e-15)
