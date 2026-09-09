from itertools import product

import numpy as np
import pytest

from .fcc111_geometry import fcc111_geometry_from_b
from .isolated_screw_core import IsolatedScrewCore,ScrewFarField,row_environment
from .nonlinear_fcc_screw import NonlinearScrewRows,stf_basis
from .nonlocal_interface_elasticity import cubic_elastic_tensor,rotate_elastic_tensor
from .range_resolved_material import build_range_surface,RangeBulkValidation
from .run_low_stress_cyclic_diagnostic import build_surface
from .static_bulk_stability import StaticBulkHessian
from .vector_fcc_rows import VectorRowKernel


@pytest.fixture(scope='module')
def setup():
    surface,_,_=build_surface(tolerance=2e-11)
    # Unitless values suffice for the geometry/force tests, not a physical fit.
    basis=fcc111_geometry_from_b(1.).plane_basis_in_stacked_cubic_axes()
    tensor=rotate_elastic_tensor(cubic_elastic_tensor(114.,62.,32.),basis)
    center=(np.sqrt(3)/12,np.sqrt(2/3)/2)
    return surface,ScrewFarField(tensor,1.,center)


def test_far_field_sign_affine_traction_and_anisotropic_equation(setup):
    _,far=setup
    theta=np.linspace(0,2*np.pi,257)
    q=np.column_stack([np.zeros_like(theta),far.center[0]+3*np.cos(theta),far.center[1]+3*np.sin(theta)])
    u=far.displacement(q)[:,0]
    assert np.unwrap(2*np.pi*u)[-1]-np.unwrap(2*np.pi*u)[0]==pytest.approx(2*np.pi)
    np.testing.assert_allclose(far.displacement(q,burgers_sign=-1),-far.displacement(q))
    affine=far.displacement(q,shear_traction=.1,burgers_sign=0)[:,0]
    gamma=np.linalg.lstsq(q[:,1:]-far.center,affine,rcond=None)[0]
    np.testing.assert_allclose(far.matrix@gamma,[0,.1],atol=1e-15)
    point=np.array([[0.,2.43,1.91]]); step=2e-4
    f=lambda p:float(far.displacement(p)[0,0])
    dy=np.array([[0,step,0]]); dz=np.array([[0,0,step]])
    yy=(f(point+dy)+f(point-dy)-2*f(point))/step**2
    zz=(f(point+dz)+f(point-dz)-2*f(point))/step**2
    yz=(f(point+dy+dz)-f(point+dy-dz)-f(point-dy+dz)+f(point-dy-dz))/(4*step**2)
    A,B,D=far.matrix[0,0],far.matrix[0,1],far.matrix[1,1]
    assert abs(A*yy+2*B*yz+D*zz)<2e-7
    assert far.log_energy_coefficient>0
    # Independent angular quadrature of the continuum strain energy in a
    # circular annulus verifies the logarithmic coefficient and row repeat.
    angle=2*np.pi*np.arange(4096)/4096
    p=(-B+1j*np.sqrt(A*D-B*B))/D
    zeta=np.cos(angle)+p*np.sin(angle)
    gradient=far.b/(2*np.pi)*np.column_stack([np.imag(1/zeta),np.imag(p/zeta)])
    coefficient=far.b*np.pi*np.mean(np.einsum('ni,ij,nj->n',gradient,far.matrix,gradient))
    assert coefficient==pytest.approx(far.log_energy_coefficient,rel=2e-14)


def test_perfect_crystal_zero_force_and_environment_closure(setup):
    surface,far=setup
    core=IsolatedScrewCore(surface,far,free_radius=1.1,ring=1,burgers_sign=0)
    out=core.evaluate(core.initial)
    assert out['energy']==pytest.approx(0,abs=2e-14)
    assert np.max(abs(out['gradient']))<2e-13
    energy=set(map(tuple,core.indices[core.energy_ids]))
    for j,l in core.indices[core.free_ids]:
        for dj,dl in core.offsets:
            assert (j+dj,l+dl) in energy
    assert len(core.energy_ids)>len(core.free_ids)
    assert not out['all_channel_infinite_tail_certified']


def test_actual_per_site_energy_gradient_hessian_and_translation_balance(setup):
    surface,far=setup
    core=IsolatedScrewCore(surface,far,free_radius=1.1,ring=1)
    q=core.initial+.003*np.sin(np.arange(core.initial.size)).reshape(core.initial.shape)
    out,hessian=core.linearize(q)
    v=np.cos(np.arange(q.size)*.7).reshape(q.shape); v/=np.linalg.norm(v)
    w=np.sin(np.arange(q.size)*.3+.7).reshape(q.shape)
    assert abs(np.sum(w*hessian(v))-np.sum(v*hessian(w)))<2e-12
    errors=[]
    for step in (2e-5,1e-5):
        plus,minus=core.evaluate(q+step*v),core.evaluate(q-step*v)
        errors.append(np.max(abs((plus['gradient']-minus['gradient'])/(2*step)-hessian(v))))
        assert abs((plus['energy']-minus['energy'])/(2*step)-np.sum(out['gradient']*v))<3e-8
    assert errors[1]<errors[0]*.4+1e-9
    assert errors[1]<2e-6
    assert np.max(abs(out['all_site_gradient'].sum(axis=0)))<2e-13
    # Summed linear work is independent of every free displacement. The new
    # annular diagnostic changes no core force or total fixed-boundary energy.
    baseline=core.evaluate(core.initial)
    assert out['linear_reference_site_work'].sum()==pytest.approx(
        baseline['linear_reference_site_work'].sum(),abs=2e-14)
    assert out['quadratic_remainder_site_energy'].sum()==pytest.approx(out['energy'],abs=2e-13)
    assert np.max(abs(out['linear_reference_site_work']))>1e-3
    assert core.boundary_winding(q)==pytest.approx(1.)
    # Whole-row x translations by b are registry-equivalent, not clipping.
    translated=q.copy(); translated[0,0]+=core.rows.b
    shifted=core.evaluate(translated)
    assert shifted['energy']==pytest.approx(out['energy'],abs=2e-13)
    np.testing.assert_allclose(shifted['gradient'],out['gradient'],atol=3e-12)


def test_no_dislocation_relaxes_to_perfect_static_state(setup):
    surface,far=setup
    core=IsolatedScrewCore(surface,far,free_radius=.75,ring=1,burgers_sign=0)
    q=core.initial+.002*np.sin(np.arange(core.initial.size)).reshape(core.initial.shape)
    out=core.relax(q,max_iterations=65,force_tolerance=3e-7)
    assert out['force_converged']
    assert abs(out['boundary_winding'])<1e-12
    assert np.max(abs(out['field']))<2e-7


def test_independently_ranged_row_channels_against_direct_atoms():
    coefficients=[.1,.3,4.,3.,1.2,15.,2.,9.]
    model=build_range_surface(2.1,5.5,3.2,coefficients)
    rows=row_environment(model.surface); kernel=VectorRowKernel(rows)
    point=np.array([[.173,.67,.31]])
    got=kernel.evaluate(point,order=2); x,y,z=point[0]
    atoms=np.column_stack([np.arange(-180,181)+x,np.full(361,y),np.full(361,z)])
    radius=np.linalg.norm(atoms,axis=1); expected=[]
    for rank,invariant in [(1,model.surface.vector),(2,model.surface.quadrupole),(3,model.surface.angular)]:
        weight=invariant.amplitude*np.exp(-invariant.kappa*radius)
        raw=np.array([np.sum(weight*np.prod(atoms[:,indices],axis=1))
                      for indices in product(range(3),repeat=rank)])
        expected.extend(raw@stf_basis(rank))
    np.testing.assert_allclose(got['value'][0,2:],expected,atol=2e-13,rtol=2e-11)
    step=2e-6
    for axis in range(3):
        delta=np.eye(3)[axis]*step
        plus,minus=kernel.evaluate(point+delta),kernel.evaluate(point-delta)
        np.testing.assert_allclose((plus['gradient']-minus['gradient'])/(2*step),
                                   got['hessian'][...,axis],rtol=2e-6,atol=2e-7)


def test_split_bulk_validator_equal_range_no_pair_double_count():
    model=build_range_surface(2.1,5.5,5.5,[.1,.3,4.,3.,1.2,15.,2.,9.])
    new=RangeBulkValidation(model,cutoff=4.)
    old=StaticBulkHessian(model.face.bulk,cutoff=4.,D3=15.,D1=2.,D2=9.,angular_decay=5.5)
    for q in ([.12,.29,.43],[1.2,-.7,.35]):
        np.testing.assert_allclose(new.evaluate(q)['matrix'],old.evaluate(q)['matrix'],atol=2e-14)


def test_common_range_fft_cannot_silently_accept_new_independent_ranges():
    model=build_range_surface(2.1,5.5,3.2,[.1,.3,4.,3.,1.2,15.,2.,9.])
    with pytest.raises(ValueError,match='range-aware vector'):
        NonlinearScrewRows(model.surface)


def test_row_affine_curvature_matches_same_full_bulk_elastic_tensor():
    from .run_isolated_screw_core import material
    surface,tensor,_=material('historical')
    far=ScrewFarField(tensor,1.,(np.sqrt(3)/12,surface.h/2))
    errors=[]
    for ring in (2,3,4):
        core=IsolatedScrewCore(surface,far,free_radius=.75,ring=ring)
        errors.append(np.max(abs(core.affine_antiplane_hessian()-far.matrix)))
    assert errors[2]<errors[1]<errors[0]
    assert errors[-1]<2e-8


@pytest.mark.parametrize('kwargs',[dict(free_radius=0),dict(ring=1.5),dict(burgers_sign=2)])
def test_invalid_core_controls(setup,kwargs):
    surface,far=setup
    with pytest.raises(ValueError):
        IsolatedScrewCore(surface,far,**kwargs)
