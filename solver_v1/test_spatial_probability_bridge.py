import numpy as np
import pytest
from .spatial_probability_bridge import aggregate_surface_history, patch_tractions_mpa


def aggregate(p, area, **kwargs):
    return aggregate_surface_history(p, area, correlation_area_mm2=.1,
        numerical_floor=0., convergence_certified=True, **kwargs)


def test_uniform_mesh_subdivision_invariance():
    p = np.array([[0.], [1.e-16], [.2], [1.]])
    a = aggregate(p, [2.])
    for count in (2, 7, 100):
        b = aggregate(np.repeat(p, count, axis=1), np.full(count, 2/count))
        np.testing.assert_allclose(a.mathematical_extrapolation, b.mathematical_extrapolation, rtol=2e-14)
    with np.errstate(divide='ignore'):
        expected = -np.expm1(20*np.log1p(-p[:,0]))
    np.testing.assert_allclose(a.mathematical_extrapolation, expected, rtol=2e-14)


def test_heterogeneous_patch_split_and_zero_area():
    a = aggregate([[.1,.3]], [1.,2.])
    b = aggregate([[.1,.1,.3,1.]], [.4,.6,2.,0.])
    np.testing.assert_allclose(a.mathematical_extrapolation,b.mathematical_extrapolation)


def test_noise_amplification_never_certifies():
    a = aggregate_surface_history([[1e-16]], [1.], correlation_area_mm2=1e-14,
        numerical_floor=1e-15, convergence_certified=True)
    assert .009 < a.mathematical_extrapolation[0] < .011
    assert np.isnan(a.numerically_resolved_extrapolation[0])
    assert not a.all_active_patches_resolved[0]


def test_missing_patch_certificate_withholds_whole_result():
    a = aggregate_surface_history([[.1,.2]], [1.,1.], correlation_area_mm2=1.,
        numerical_floor=0., convergence_certified=[True,False])
    assert np.isnan(a.numerically_resolved_extrapolation[0])


@pytest.mark.parametrize('p', [[[.2],[.1]], [[-1e-20]], [[1.1]], [[float('nan')]]])
def test_invalid_probabilities_not_clipped(p):
    with pytest.raises(ValueError): aggregate(p,[1.])


def test_tensor_traction_and_rotation_covariance():
    stress=np.array([[100.,20.,30.],[20.,50.,40.],[30.,40.,150.]])
    n=np.array([0.,0.,1.]); m=np.array([1.,0.,0.])
    a=patch_tractions_mpa([stress],[n],[m])
    np.testing.assert_allclose(a,[[150.,30.,40.]])
    angle=.73
    rotation=np.array([[np.cos(angle),0,np.sin(angle)],[0,1,0],[-np.sin(angle),0,np.cos(angle)]])
    b=patch_tractions_mpa([rotation@stress@rotation.T],[rotation@n],[rotation@m])
    np.testing.assert_allclose(a,b,atol=1e-13)


def test_uniform_tension_and_pure_shear():
    n=[[0,0,1]]; m=[[1,0,0]]
    np.testing.assert_allclose(patch_tractions_mpa([np.diag([0,0,100])],n,m),[[100,0,0]])
    np.testing.assert_allclose(patch_tractions_mpa([[[0,0,25],[0,0,0],[25,0,0]]],n,m),[[0,25,0]])
