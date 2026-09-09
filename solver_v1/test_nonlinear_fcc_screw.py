"""Nonlinear energy checks do not certify an Al core or physical kinetics."""
from itertools import product

import numpy as np
import pytest

from .discrete_fcc_screw import FCCScrewRowHessian
from .nonlinear_fcc_screw import NonlinearScrewRows, stf_basis
from .nonlinear_screw_transverse_audit import omitted_transverse_forces, transverse_row_coefficients
from .run_low_stress_cyclic_diagnostic import build_surface


@pytest.fixture(scope="module")
def model():
    surface,units,meta=build_surface(tolerance=2e-11)
    return surface,units,meta,NonlinearScrewRows(surface,tolerance=2e-12)


@pytest.mark.parametrize("rank",[1,2,3])
def test_stf_compression_is_orthonormal(rank):
    basis=stf_basis(rank)
    np.testing.assert_allclose(basis.T@basis,np.eye(2*rank+1),atol=2e-15)
    tensor=(basis@np.arange(2*rank+1)).reshape((3,)*rank)
    if rank>1:
        # Eigenspace projection has norm-scaled floating arithmetic error.
        roundoff=8*np.finfo(float).eps*np.linalg.norm(tensor)
        np.testing.assert_allclose(np.trace(tensor,axis1=-2,axis2=-1),0.,atol=roundoff)


def test_nonlinear_row_values_against_independent_direct_atoms(model):
    surface,_,_,rows=model
    for j,l,shift in ((1,0,.173),(0,1,-.217),(-1,2,.431)):
        index=np.flatnonzero((rows.j==j)&(rows.l==l))[0]
        y,z=rows.d*(j+l/3),rows.h*l
        delta=rows.b*(j+l)/2
        def real_sum(phase):
            x=rows.b*np.arange(-384,385)+phase
            r=np.sqrt(x*x+y*y+z*z)
            p=rows.bulk.p
            pair=np.sum(4*p.epsilon*((p.sigma_lj/r)**12-(p.sigma_lj/r)**6))
            density=np.sum(rows.bulk.density_params.C_rho*np.exp(-rows.bulk.density_params.kappa*r))
            moments=[]
            weight=surface.angular.amplitude*np.exp(-surface.angular.kappa*r)
            for rank in (1,2,3):
                raw=np.array([np.sum(weight*x**a.count(0)*y**a.count(1)*z**a.count(2))
                              for a in product(range(3),repeat=rank)])
                moments.extend(raw@stf_basis(rank))
            return np.r_[pair,density,moments]
        np.testing.assert_allclose(rows.row_changes(index,shift),real_sum(delta+shift)-real_sum(delta),
                                   atol=5e-13,rtol=2e-10)


def test_fft_is_exact_real_index_convolution(model):
    rows=model[-1]; cell=rows.periodic_cell((7,8))
    j,l=np.indices(cell.shape); u=.137*np.cos(.7*j+.9*l)
    exact=cell.direct_channels(u,.023)
    value=cell.evaluate(u,.023,fields=True)
    np.testing.assert_allclose(value['channels'],exact,atol=2e-14,rtol=3e-11)
    assert value['minimum_site_density']>0


def test_nonlinear_gradient_hessian_and_affine_work(model):
    cell=model[-1].periodic_cell((7,8)); j,l=np.indices(cell.shape)
    u=.077*np.cos(.7*j+.9*l); v=np.sin(.3*j+.7*l); w=np.cos(.8*j+.4*l)
    gamma=.011; dg=.017; dw=-.023; step=2e-6
    center=cell.evaluate(u,gamma,direction=v,gamma_direction=dg)
    plus=cell.evaluate(u+step*v,gamma+step*dg)
    minus=cell.evaluate(u-step*v,gamma-step*dg)
    derivative=(plus['energy']-minus['energy'])/(2*step)
    assert derivative==pytest.approx(np.sum(center['gradient']*v)+dg*center['affine_derivative'],rel=2e-7,abs=2e-8)
    np.testing.assert_allclose((plus['gradient']-minus['gradient'])/(2*step),
                              center['hessian_vector'],atol=2e-8,rtol=3e-7)
    assert (plus['affine_derivative']-minus['affine_derivative'])/(2*step)==pytest.approx(
        center['affine_hessian_vector'],rel=2e-7,abs=3e-8)
    other=cell.evaluate(u,gamma,direction=w,gamma_direction=dw)
    assert (np.sum(w*center['hessian_vector'])+dw*center['affine_hessian_vector'])==pytest.approx(
        np.sum(v*other['hessian_vector'])+dg*other['affine_hessian_vector'],abs=2e-12)
    assert abs(center['gradient'].sum())<2e-12


def test_exact_harmonic_limit_matches_verified_row_symbol(model):
    surface,_,_,rows=model; cell=rows.periodic_cell((9,10))
    harmonic=FCCScrewRowHessian.from_surface(surface,tolerance=2e-12)
    j,l=np.indices(cell.shape); qj,ql=2*np.pi*2/9,2*np.pi*3/10
    mode=np.cos(qj*j+ql*l)
    result=cell.evaluate(np.zeros(cell.shape),direction=mode)
    np.testing.assert_allclose(result['hessian_vector'],harmonic.symbol(qj/rows.d,ql)*mode,
                              atol=3e-12,rtol=2e-12)
    # Affine gamma changes x by gamma*z: independent harmonic row moment check.
    expected=.5*np.sum(harmonic.effective_pair*(harmonic.l*rows.h)**2)
    for rank in (1,2,3):
        collective=(harmonic.moments[rank]*(harmonic.l*rows.h)[:,None]).sum(axis=0)
        expected+=2*harmonic.amplitudes[rank]*(collective@collective)
    affine=cell.evaluate(np.zeros(cell.shape),direction=np.zeros(cell.shape),gamma_direction=1.)
    assert affine['affine_hessian_vector']/cell.size==pytest.approx(expected,abs=3e-12)


def test_full_nonlinear_rigid_cut_matches_independent_plane_energy(model):
    surface,_,_,rows=model; cell=rows.periodic_cell((4,40))
    for slip in (.073,.23,.41):
        u=np.zeros(cell.shape);u[:,10:30]=slip
        actual=cell.evaluate(u)
        independent=surface.packed(surface.h,slip)
        assert actual['energy']/(2*cell.shape[0])==pytest.approx(independent[0],abs=3e-10)
        assert actual['gradient'][:,10:30].sum()/(2*cell.shape[0])==pytest.approx(independent[2],abs=2e-9)


def test_translation_integer_relabeling_and_topology(model):
    rows=model[-1];cell=rows.periodic_cell((24,24))
    u=cell.dipole_seed(7.3); winding=cell.winding(u)
    assert (winding['charge']==1).sum()==1
    assert (winding['charge']==-1).sum()==1
    assert winding['integer_residual']<1e-14
    j,l=np.indices(cell.shape);shift=rows.b*((j+2*l)%3-1)
    old=cell.evaluate(u); new=cell.evaluate(u+shift+.213)
    assert new['energy']==pytest.approx(old['energy'],abs=3e-12)
    np.testing.assert_allclose(new['gradient'],old['gradient'],atol=3e-12)
    np.testing.assert_array_equal(cell.winding(u+shift)['charge'],winding['charge'])
    # A b-valued array is identical atomic registry, not a vortex/dislocation.
    np.testing.assert_array_equal(cell.winding(shift)['charge'],0)
    assert abs(cell.evaluate(shift)['energy'])<1e-12


def test_nonlinear_tail_refinement(model):
    surface,_,_,fine=model
    j,l=np.indices((7,8));u=.21*np.cos(.71*j+.43*l)
    target=fine.periodic_cell(u.shape).evaluate(u,.003)
    errors=[]
    for ring in (2,4,6):
        rows=NonlinearScrewRows(surface,fixed_ring=ring)
        result=rows.periodic_cell(u.shape).evaluate(u,.003)
        errors.append(max(abs(result['energy']-target['energy']),np.max(abs(result['gradient']-target['gradient']))))
    assert errors[2]<errors[1]<errors[0]
    assert errors[2]<2e-8


def test_no_hidden_cutoff_or_acceptance_on_exhaustion(model):
    surface=model[0]
    with pytest.raises(ArithmeticError):
        NonlinearScrewRows(surface,max_modes=3)
    with pytest.raises(ArithmeticError):
        NonlinearScrewRows(surface,max_ring=5,tolerance=1e-20)


def test_registry_shear_identity_and_integer_gauge(model):
    rows=model[-1];cell=rows.periodic_cell((24,24));u=cell.dipole_seed(7.3)
    gamma=.017;j,l=np.indices(cell.shape)
    before=cell.registry_shear(u,gamma)
    after=cell.registry_shear(u+rows.b*((j+2*l)%3-1),gamma)
    for result in (before,after):
        assert abs(result['decomposition_residual'])<2e-16
    assert after['registry_shear']==pytest.approx(before['registry_shear'],abs=2e-16)
    assert after['intrabond_shear']==pytest.approx(before['intrabond_shear'],abs=2e-16)
    # A b-only relabeling does not create macroscopic new slip.
    assert cell.registry_shear(rows.b*((j+2*l)%3-1))['registry_shear']==0.


def test_stress_control_perfect_loading_and_unloading(model):
    cell=model[-1].periodic_cell((6,6));u=np.zeros(cell.shape);gamma=0.
    conversion=1e6*model[1].length_scale_m**3/1.602176634e-19
    for tau in (0.,4.,50.,-50.,0.):
        out=cell.relax(u,gamma=gamma,shear_stress=tau*conversion)
        assert out['force_converged']
        assert out['maximum_force_residual']<2e-8
        assert out['internal_shear_stress']/conversion==pytest.approx(tau,abs=2e-5)
        u=out['displacement'];gamma=out['affine_shear']
        assert cell.registry_shear(u,gamma)['registry_shear']==0.
        assert not out['physical_time_available']
    assert abs(gamma)<1e-9
    assert np.max(abs(u))<1e-12


def test_defect_zero_stress_is_not_zero_strain_and_new_slip_is_separate(model):
    cell=model[-1].periodic_cell((24,24))
    initial=cell.dipole_seed(7.3)
    fixed=cell.evaluate(initial,0.)
    assert abs(fixed['affine_derivative'])>1.
    out=cell.relax(initial,shear_stress=0.)
    assert out['force_converged']
    assert out['maximum_force_residual']<2e-8
    assert abs(out['affine_shear'])>1e-3  # pre-existing defect content, not generated plasticity
    assert abs(out['internal_shear_stress'])<1e-9
    shear=cell.registry_shear(out['displacement'],out['affine_shear'])
    assert abs(shear['decomposition_residual'])<2e-16
    assert abs(shear['registry_shear'])>1e-3
    assert (out['winding']['charge']==1).sum()==1
    assert (out['winding']['charge']==-1).sum()==1


def test_minimum_curvature_excludes_artificial_translation_shift(model):
    cell=model[-1].periodic_cell((4,4));u=np.zeros(cell.shape)
    # Independent dense Hessian including the exact translation null mode.
    columns=[]
    for i in range(cell.size):
        direction=np.eye(cell.size)[i].reshape(cell.shape)
        columns.append(cell.evaluate(u,direction=direction)['hessian_vector'].ravel())
    eigenvalues=np.linalg.eigvalsh(np.array(columns).T)
    out=cell.minimum_curvature(u)
    assert eigenvalues[1]>1.  # makes a blindly reported gauge eigenvalue wrong
    assert out['minimum_eigenvalue']==pytest.approx(eigenvalues[1],rel=2e-8)
    assert out['eigen_residual']<1e-6
    assert out['translation_overlap']<1e-6


def test_transverse_channel_derivatives_include_zero_mode(model):
    rows=model[-1];j=np.array([1,0,-1]);l=np.array([0,1,2])
    y=rows.d*(j+l/3);z=rows.h*l;step=1e-6
    for mode in (0,1,2,3):
        value,derivatives=transverse_row_coefficients(rows,j,l,mode)
        for axis in (0,1):
            plus=transverse_row_coefficients(rows,j,l,mode,y=y+(step if axis==0 else 0),
                                               z=z+(step if axis==1 else 0))[0]
            minus=transverse_row_coefficients(rows,j,l,mode,y=y-(step if axis==0 else 0),
                                                z=z-(step if axis==1 else 0))[0]
            np.testing.assert_allclose((plus-minus)/(2*step),derivatives[...,axis],atol=2e-8,rtol=2e-7)
        assert np.all(np.isfinite(value))


def test_omitted_normal_force_matches_independent_rigid_plane(model):
    surface,_,_,rows=model;cell=rows.periodic_cell((4,48));u=np.zeros(cell.shape);u[:,12:36]=.21
    result=omitted_transverse_forces(cell,u,ring=16)
    # Separate the two interfaces. Translate the lower/upper quarter relative
    # to each other at the first cut; use forces localized around that cut.
    independent=surface.packed(surface.h,.21)
    assert result['gradient'][:,12:24,1].sum()/cell.shape[0]==pytest.approx(independent[1],abs=2e-6)
    np.testing.assert_allclose(result['net_gradient'],0.,atol=2e-12)
    perfect=omitted_transverse_forces(cell,np.zeros(cell.shape),ring=10)
    assert np.max(abs(perfect['gradient']))<2e-12
