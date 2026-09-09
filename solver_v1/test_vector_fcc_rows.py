"""Vector infinite-row identities; static research, not Al strength certification."""
from itertools import product

import numpy as np
import pytest

from .discrete_fcc_screw import exponential_row_transform_derivatives
from .lattice_bessel import _mode_coefficient_and_deda, _mode_coefficient_second_deda
from .nonlinear_fcc_screw import NonlinearScrewRows, stf_basis
from .nonlinear_screw_transverse_audit import exponential_transform_radial_derivatives, omitted_transverse_forces
from .run_low_stress_cyclic_diagnostic import build_surface
from .vector_fcc_rows import VectorPeriodicCell, VectorRowKernel, exponential_radial, power_mode_radial
from .vector_fcc_validation import direct_row_channels


@pytest.fixture(scope="module")
def rows():
    return NonlinearScrewRows(build_surface(tolerance=2e-11)[0],tolerance=2e-12)


def direct_channels(rows,vector,images=384):
    """Independent real atom sum, validation only (no production cutoff)."""
    x,y,z=vector
    atoms=np.column_stack([rows.b*np.arange(-images,images+1)+x,
                          np.full(2*images+1,y),np.full(2*images+1,z)])
    radius=np.linalg.norm(atoms,axis=1);p=rows.bulk.p
    pair=np.sum(4*p.epsilon*((p.sigma_lj/radius)**12-(p.sigma_lj/radius)**6))
    den=rows.bulk.density_params
    density=np.sum(den.C_rho*np.exp(-den.kappa*radius))
    weight=rows.surface.angular.amplitude*np.exp(-rows.surface.angular.kappa*radius)
    moments=[]
    for rank in (1,2,3):
        raw=np.array([np.sum(weight*np.prod(atoms[:,indices],axis=1))
                      for indices in product(range(3),repeat=rank)])
        moments.extend(raw@stf_basis(rank))
    return np.r_[pair,density,moments]


@pytest.mark.parametrize("p",[3,6])
def test_exact_half_integer_polynomial(p):
    radius=np.array([.56,.91,1.71,4.83])
    for mode in (1,2,5):
        actual=power_mode_radial(radius,p,mode,1.17)
        reference=(*_mode_coefficient_and_deda(radius,p,mode,1.17),
                   _mode_coefficient_second_deda(radius,p,mode,1.17))
        np.testing.assert_allclose(actual,reference,rtol=7e-15,atol=1e-28)


def test_integer_exponential_radial_derivatives():
    radius=np.array([.63,.92,1.35,4.12]);step=2e-6
    for g in (0.,2.1,8.3):
        value,first,second=exponential_radial(radius,g,1.287)
        reference=np.stack(exponential_row_transform_derivatives(radius,g,1.287),axis=-1)
        np.testing.assert_allclose(value,reference,rtol=2e-14,atol=1e-15)
        np.testing.assert_allclose(first,exponential_transform_radial_derivatives(radius,g,1.287),rtol=2e-14,atol=2e-15)
        fd=(exponential_radial(radius+step,g,1.287)[1]-exponential_radial(radius-step,g,1.287)[1])/(2*step)
        np.testing.assert_allclose(second,fd,rtol=3e-8,atol=1e-8)


def test_vector_row_direct_atoms_and_truncation(rows):
    vectors=np.array([[.173,.73,.19],[-.217,.29,.79],[.431,-.63,1.37],[-.121,2.41,-.83]])
    kernel=VectorRowKernel(rows)
    actual=kernel.evaluate(vectors,order=2)
    direct=np.array([direct_channels(rows,v) for v in vectors])
    np.testing.assert_allclose(actual['value'],direct,rtol=2e-10,atol=3e-13)
    tighter=VectorRowKernel(rows,tolerance=2e-14).evaluate(vectors,order=2)
    for field in ('value','gradient','hessian'):
        np.testing.assert_allclose(actual[field],tighter[field],atol=5e-13,rtol=2e-13)
    assert actual['modes_used']<40
    assert actual['maximum_last_mode_envelope']<kernel.tolerance


def test_vector_row_analytic_gradient_hessian_and_parity(rows):
    vectors=np.array([[.173,.73,.19],[-.217,.29,.79],[.431,-.63,1.37]])
    kernel=VectorRowKernel(rows);center=kernel.evaluate(vectors,order=2);step=2e-6
    for axis in range(3):
        delta=np.eye(3)[axis]*step
        plus=kernel.evaluate(vectors+delta);minus=kernel.evaluate(vectors-delta)
        np.testing.assert_allclose((plus['value']-minus['value'])/(2*step),
                                  center['gradient'][...,axis],rtol=8e-7,atol=1e-8)
        np.testing.assert_allclose((plus['gradient']-minus['gradient'])/(2*step),
                                  center['hessian'][...,axis],rtol=2e-6,atol=2e-7)
        direct=(np.array([direct_channels(rows,v+delta) for v in vectors])-
                np.array([direct_channels(rows,v-delta) for v in vectors]))/(2*step)
        np.testing.assert_allclose(direct,center['gradient'][...,axis],rtol=8e-7,atol=1e-8)
    np.testing.assert_array_equal(center['hessian'],center['hessian'].swapaxes(-1,-2))
    reverse=kernel.evaluate(-vectors,order=2)
    for field,order in (('value',0),('gradient',1),('hessian',2)):
        parity=kernel.parity.reshape((1,17)+(1,)*order)*(-1)**order
        np.testing.assert_allclose(reverse[field],parity*center[field],atol=2e-13,rtol=1e-13)
    shifted=kernel.evaluate(vectors+np.array([rows.b,0.,0.]),order=2)
    np.testing.assert_allclose(center['value'],shifted['value'],atol=5e-14)


def test_independent_direct_atom_hessian(rows):
    vectors=np.array([[.173,.73,.19],[-.217,.29,.79],[.431,-.63,1.37],[.117,3.13,-.27]])
    reciprocal=VectorRowKernel(rows).evaluate(vectors,order=2)
    references=[direct_row_channels(rows,v) for v in vectors]
    for field in ('value','gradient','hessian'):
        np.testing.assert_allclose(reciprocal[field],np.array([r[field] for r in references]),
                                  atol=2e-12,rtol=2e-10)


def test_vector_kernel_rejects_unresolved_series(rows):
    with pytest.raises(ArithmeticError):
        VectorRowKernel(rows,max_modes=3).evaluate([[.17,.73,.19]])
    with pytest.raises(ValueError):
        VectorRowKernel(rows).evaluate([[.17,0.,0.]])


def vector_field(shape,amplitude=.013):
    j,l=np.indices(shape)
    return amplitude*np.stack([np.sin(.7*j+.8*l),np.cos(.6*j+.2*l),np.sin(.3*j+.7*l)],axis=-1)


def test_per_site_environment_against_direct_atom_assembly(rows):
    cell=VectorPeriodicCell(rows,(4,4),ring=2)
    field=vector_field(cell.shape);gamma=.006
    direct=np.zeros((cell.size,17));flat=field.reshape((cell.size,3))
    for i in range(cell.size):
        for neighbor,reference in zip(cell.destination[i],cell.reference):
            vector=reference+flat[neighbor]-flat[i]
            vector[0]+=gamma*reference[2]
            difference=direct_channels(rows,vector)-direct_channels(rows,reference)
            direct[i]+=difference
            direct[neighbor]+=cell.kernel.parity*difference
    out=cell.evaluate(field,gamma,fields=True)
    np.testing.assert_allclose(out['channels'].reshape((cell.size,17)),direct,atol=5e-13,rtol=2e-10)
    rho=rows.rho_bulk+direct[:,1]
    energy=.5*direct[:,0].sum()+np.sum(rows.embedding.value(rho)-rows.embedding.value(rows.rho_bulk))
    energy+=np.sum(rows.angular_weights*direct[:,2:]**2)
    assert out['energy']==pytest.approx(energy,abs=2e-12)


def test_scalar_subspace_and_omitted_transverse_force_reproduced(rows):
    # Same retained transverse neighborhood in the two independent assemblers.
    ring=6;fixed=NonlinearScrewRows(rows.surface,fixed_ring=ring)
    cell=VectorPeriodicCell(fixed,(7,8),ring=ring)
    scalar=fixed.periodic_cell(cell.shape)
    u=vector_field(cell.shape,.07)[...,0];field=np.zeros(cell.shape+(3,));field[...,0]=u
    gamma=.019
    exact=scalar.evaluate(u,gamma,fields=True);vector=cell.evaluate(field,gamma,fields=True)
    np.testing.assert_allclose(vector['channels'],exact['channels'],atol=3e-13,rtol=3e-10)
    assert vector['energy']==pytest.approx(exact['energy'],abs=8e-12)
    np.testing.assert_allclose(vector['gradient'][...,0],exact['gradient'],atol=8e-12,rtol=3e-10)
    omitted=omitted_transverse_forces(scalar,u,gamma,ring=ring)
    np.testing.assert_allclose(vector['gradient'][...,1:],omitted['gradient'],atol=8e-12,rtol=3e-10)
    assert vector['affine_derivative']==pytest.approx(exact['affine_derivative'],abs=1e-11)


def test_full_vector_gradient_hessian_and_stress_work(rows):
    cell=VectorPeriodicCell(rows,(5,6),ring=4)
    u=vector_field(cell.shape);v=vector_field(cell.shape,.7)
    w=np.roll(v,1,axis=0)*np.array([.9,-1.2,.7])
    gamma=.012;dg=.019;dw=-.031;step=2e-6
    center,apply=cell.linearize(u,gamma);hv,hg=apply(v,dg)
    plus=cell.evaluate(u+step*v,gamma+step*dg);minus=cell.evaluate(u-step*v,gamma-step*dg)
    assert (plus['energy']-minus['energy'])/(2*step)==pytest.approx(
        np.sum(center['gradient']*v)+dg*center['affine_derivative'],abs=2e-8,rel=3e-7)
    np.testing.assert_allclose((plus['gradient']-minus['gradient'])/(2*step),hv,atol=8e-8,rtol=2e-6)
    assert (plus['affine_derivative']-minus['affine_derivative'])/(2*step)==pytest.approx(hg,abs=8e-8,rel=2e-6)
    hw,hwg=apply(w,dw)
    assert (np.sum(w*hv)+dw*hg)==pytest.approx(np.sum(v*hw)+dg*hwg,abs=3e-12)
    np.testing.assert_allclose(center['gradient'].sum(axis=(0,1)),0.,atol=2e-12)
    for axis in range(3):
        shift=np.zeros_like(u);shift[...,axis]=1.
        a,b=apply(shift)
        np.testing.assert_allclose(a,0.,atol=1e-13);assert abs(b)<1e-13


def test_integer_line_relabeling_and_vector_translation(rows):
    cell=VectorPeriodicCell(rows,(5,6),ring=4);u=vector_field(cell.shape)
    j,l=np.indices(cell.shape);other=u+np.array([.123,.059,-.043])
    other[...,0]+=rows.b*((j+2*l)%3-1)
    a=cell.evaluate(u,.012);b=cell.evaluate(other,.012)
    assert a['energy']==pytest.approx(b['energy'],abs=4e-12)
    np.testing.assert_allclose(a['gradient'],b['gradient'],atol=3e-12)


def test_perfect_vector_equilibrium_and_hessian_symmetry(rows):
    cell=VectorPeriodicCell(rows,(4,4),ring=4);u=np.zeros(cell.shape+(3,))
    out,apply=cell.linearize(u)
    assert abs(out['energy'])<2e-13
    np.testing.assert_allclose(out['gradient'],0.,atol=2e-13)
    columns=[]
    for direction in np.eye(3*cell.size):
        columns.append(apply(direction.reshape(u.shape))[0].ravel())
    H=np.array(columns).T
    np.testing.assert_allclose(H,H.T,atol=2e-13)
    eigen=np.linalg.eigvalsh(H)
    np.testing.assert_allclose(eigen[:3],0.,atol=8e-14)
    assert eigen[3]>0
    result=cell.minimum_curvature(u)
    assert result['minimum_eigenvalue']==pytest.approx(eigen[3],abs=3e-8)
    assert result['translation_overlap']<1e-6


def test_vector_force_relaxation_and_static_unload(rows):
    cell=VectorPeriodicCell(rows,(4,4),ring=4);u=vector_field(cell.shape,.015);gamma=0.
    for tau in (0.,.007329590081665426,0.):  # hypothetical imposed eV/L0^3, no Hz
        result=cell.relax(u,gamma=gamma,shear_stress=tau,force_tolerance=2e-7)
        assert result['force_converged']
        assert result['maximum_force_residual']<2e-7
        assert result['transverse_cell_strains_fixed']
        assert not result['physical_time_available']
        u=result['displacement'];gamma=result['affine_shear']
    assert abs(gamma)<2e-7
    assert np.max(abs(u))<2e-7


def test_pair_zero_mode_tail_bound_is_not_all_channel_certificate(rows):
    field=vector_field((5,6),.07);previous=None
    for ring in (4,6,8,12):
        cell=VectorPeriodicCell(rows,field.shape[:2],ring=ring)
        bound=cell.pair_zero_mode_tail_bounds(field)
        assert not bound['all_channel_infinite_tail_certified']
        assert bound['pair_zero_mode_force_norm_bound']>0
        if previous is not None:
            assert bound['pair_zero_mode_force_norm_bound']<previous
        previous=bound['pair_zero_mode_force_norm_bound']
        zero=cell.pair_zero_mode_tail_bounds(np.zeros_like(field))
        assert zero['pair_zero_mode_energy_bound_ev']==0.
        assert zero['pair_zero_mode_force_norm_bound']==0.


def test_released_small_pair_does_not_infer_full_topology_from_x_winding(rows):
    # A real nonlinear solve, not an assertion against a saved endpoint.
    # This supplied close pair is not a predicted defect population or yield.
    cell=VectorPeriodicCell(rows,(8,8),ring=4)
    scalar=cell.scalar.relax(cell.scalar.dipole_seed(3.),shear_stress=0.)
    assert scalar['force_converged']
    field=np.zeros(cell.shape+(3,));field[...,0]=scalar['displacement']
    before=cell.evaluate(field,scalar['affine_shear'])
    assert np.max(abs(before['gradient'][...,1:]))>.1
    out=cell.relax(field,gamma=scalar['affine_shear'],shear_stress=0.,max_iterations=130)
    assert out['force_converged']
    assert out['energy']<before['energy']
    assert not np.any(out['x_winding']['charge'])
    # Actual counterexample: x winding vanished, but a cell-spanning vector
    # registry fault remains. Do NOT force the test to a desired perfect state.
    assert np.max(abs(out['displacement'][...,1:]))>.1
    layers=cell.layer_registry(out['displacement'],out['affine_shear'])
    fault=min(layers,key=lambda r:r['distance_to_extra_tau'])
    assert fault['distance_to_extra_tau']<.04
    assert fault['minimum_scalar_x_partition_margin']<1e-6
    assert fault['scalar_x_partition_is_not_vector_plastic_strain']
    assert max(r['x_phase_dispersion'] for r in layers)<1e-6
    assert cell.minimum_curvature(out['displacement'],out['affine_shear'],stress_control=True)['minimum_eigenvalue']>0


def test_vector_registry_labels_are_invariant_to_row_integer_relabeling(rows):
    cell=VectorPeriodicCell(rows,(4,5),ring=4);field=vector_field(cell.shape)
    other=field.copy();j,l=np.indices(cell.shape);other[...,0]+=rows.b*((j+2*l)%3-1)
    a=cell.layer_registry(field,.017);b=cell.layer_registry(other,.017)
    for left,right in zip(a,b):
        for key in ('distance_to_perfect_extra_registry','distance_to_extra_tau','distance_to_extra_2tau',
                    'minimum_scalar_x_partition_margin','x_phase_dispersion'):
            assert left[key]==pytest.approx(right[key],abs=3e-15)


def test_zero_mode_tail_bound_covers_independent_shell_difference(rows):
    field=vector_field((5,6),.05);p=rows.bulk.p;results=[]
    for ring in (4,14):
        cell=VectorPeriodicCell(rows,field.shape[:2],ring=ring)
        u=field.reshape((cell.size,3))
        r=cell.reference[None,:,:]+u[cell.destination]-u[:,None,:]
        actual=np.linalg.norm(r[...,1:],axis=-1)
        reference=np.linalg.norm(cell.reference[:,1:],axis=-1)
        energy=np.zeros_like(actual);derivative=np.zeros_like(actual)
        for power,scale in ((6,4*p.epsilon*p.sigma_lj**12),(3,-4*p.epsilon*p.sigma_lj**6)):
            value,grad,_=power_mode_radial(actual,power,0,rows.b)
            base=power_mode_radial(reference,power,0,rows.b)[0]
            energy+=scale*(value-base);derivative+=scale*grad
        bond_force=derivative[...,None]*r[...,1:]/actual[...,None]
        gradient=-bond_force.sum(axis=1)+cell._scatter(bond_force)
        results.append((energy.sum(),gradient,cell.pair_zero_mode_tail_bounds(field)))
    a,b=results
    assert abs(a[0]-b[0])<=a[2]['pair_zero_mode_energy_bound_ev']+b[2]['pair_zero_mode_energy_bound_ev']
    assert np.max(np.linalg.norm(a[1]-b[1],axis=-1))<=a[2]['pair_zero_mode_force_norm_bound']+b[2]['pair_zero_mode_force_norm_bound']
