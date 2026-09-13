import numpy as np
import pytest
from .polynomial_exponential_density import QuadraticEnvelopeDensity
from .quadratic_density_interface import QuadraticScalarEnvironment, QuadraticDensityInterface, bulk_scalar_columns
from .coordination_screening import CoordinationScreenedBulk

SHAPE=(4.9,10.,11.,1000.,8.4,1.)
KERNEL=QuadraticEnvelopeDensity(1.,8.94,.6843,.0894)


@pytest.fixture(scope='module')
def env():
    return QuadraticScalarEnvironment(KERNEL)


def test_infinite_bulk_normalization_direct_and_elastic_derivatives(env):
    bulk=CoordinationScreenedBulk(SHAPE,radius=10.,law='power')
    k=env.kernel; r=bulk.r; e=bulk.e
    f,fp,fpp=[k.radial(r,n) for n in (0,1,2)]
    assert f.sum() == pytest.approx(env.reference_density,abs=2e-12)
    rr=np.einsum('ni,nj->nij',e,e)
    H=(fpp-fp/r)[:,None,None]*rr+(fp/r)[:,None,None]*np.eye(3)
    normal=np.column_stack([np.zeros(len(r)),np.zeros(len(r)),bulk.R[:,2]])
    shear=np.column_stack([bulk.R[:,2],np.zeros(len(r)),np.zeros(len(r))])
    result=env.bulk()
    assert result[1] == pytest.approx(np.einsum('ni,ni->',fp[:,None]*e,normal),abs=2e-11)
    assert result[2] == pytest.approx(np.einsum('ni,nij,nj->',normal,H,normal),abs=2e-10)
    assert result[4] == pytest.approx(np.einsum('ni,nij,nj->',shear,H,shear),abs=2e-10)
    one,two=env.bulk_columns(8e-4),env.bulk_columns(4e-4)
    np.testing.assert_allclose(one,two,atol=2e-7,rtol=2e-8)
    hydro=(fpp*r*r).sum()/result[0]
    assert two[2,1] == pytest.approx(hydro,abs=2e-7)


def test_interface_force_hessian_and_reference():
    model=QuadraticDensityInterface(SHAPE,np.ones(10),KERNEL)
    assert abs(model.evaluate((model.h,0.,0.)).energy)<1e-12
    q=np.array([.879,.267,.113]);value=model.evaluate(q);step=1e-5
    hi=[model.evaluate(q+step*e) for e in np.eye(3)]
    lo=[model.evaluate(q-step*e) for e in np.eye(3)]
    g=[(a.energy-b.energy)/(2*step) for a,b in zip(hi,lo)]
    H=np.column_stack([(a.gradient-b.gradient)/(2*step) for a,b in zip(hi,lo)])
    np.testing.assert_allclose(g,value.gradient,atol=1e-6,rtol=2e-6)
    np.testing.assert_allclose(H,value.hessian,atol=3e-5,rtol=3e-6)
    np.testing.assert_array_equal(value.hessian,value.hessian.T)
    repeat=model.evaluate(q+[0.,1.,0.])
    np.testing.assert_allclose(repeat.hessian,value.hessian,atol=2e-10)


def test_bulk_bloch_direct_actual_embedding_energy(env):
    bulk=CoordinationScreenedBulk(SHAPE,radius=10.,law='power');q=(1/3,)*3
    columns,tails=bulk_scalar_columns(bulk,env,q)
    R=bulk.R;layers=np.rint(R[:,2]/bulk.geometry.h111)
    def energy(amplitude,v):
        values=[]
        for site in range(6):
            shift=amplitude*(np.cos(2*np.pi*(site+layers)/3)-np.cos(2*np.pi*site/3))
            rho=np.sum(env.kernel.radial(np.linalg.norm(R+shift[:,None]*v,axis=1)))/env.reference_density
            values.append([-np.sqrt(rho),rho-1,(rho-1)**2])
        return np.mean(values,axis=0)
    for v in np.eye(3):
        slopes=[2*(energy(h,v)+energy(-h,v)-2*energy(0.,v))/h**2 for h in (2e-4,1e-4)]
        exact=np.einsum('i,cij,j->c',v,columns,v)
        np.testing.assert_allclose((4*slopes[1]-slopes[0])/3,exact,atol=2e-6)
        assert np.all(tails>=0)


def test_calibration_matrix_replays_actual_interface_jets():
    from .run_quadratic_density_material_v32 import interface_matrix
    from .interface_tangent_calibration import TangentCalibrationProblem
    from .run_current_material_core import ROOT
    problem=TangentCalibrationProblem(ROOT/'results/fcc111_active_interface/coordination_screening_v19/joint_refinement')
    coefficients=np.array([.06,.5,4.,-2.,.8,3.,.4,-.2,.7,9.])
    model=QuadraticDensityInterface(SHAPE,coefficients,KERNEL)
    matrix=interface_matrix(problem,SHAPE,model.scalar_environment)
    # Every nonbulk observation is independently evaluated by the full model,
    # not by reconstructing the matrix under test.
    states={o.state for o in problem.raw_observations if o.state is not None}
    values={state:model.evaluate(state) for state in states}
    for i,o in enumerate(problem.raw_observations):
        if o.state is None:
            continue
        v=values[o.state]
        jet=np.r_[v.energy,v.gradient,v.hessian[0],v.hessian[1,1:],v.hessian[2,2]]
        assert matrix[i]@coefficients == pytest.approx(np.asarray(o.jet_weights)@jet,abs=3e-9,rel=3e-11)


def test_density_search_coordinates_fix_only_amplitude_gauge():
    from .run_quadratic_density_joint_v32 import encode,decode
    original=QuadraticEnvelopeDensity(7.,8.94,.6843,.0894)
    replay=decode(encode(original))
    assert replay.C==1.
    np.testing.assert_allclose([replay.k,replay.center,replay.width],
                               [original.k,original.center,original.width])
    with pytest.raises(ValueError):
        decode([1.,np.nan,2.])


def test_joint_shape_decay_uses_actual_new_scalar_kernel():
    from .run_quadratic_density_joint_v32 import encode,joint_decode
    shape=(4.,6.,8.,100.,3.,-1.)
    # Old placeholder 4 would pass 2*3-4>0, but actual density decay9 fails.
    with pytest.raises(ValueError,match='decay'):
        joint_decode(encode(KERNEL),shape)
    with pytest.raises(ValueError,match='decay'):
        QuadraticDensityInterface(shape,np.ones(10),KERNEL)
    x=np.r_[encode(KERNEL),np.log(SHAPE[1:5]),SHAPE[5]]
    kernel,replay=joint_decode(x,SHAPE,True)
    np.testing.assert_allclose(replay,SHAPE)
    assert kernel.k==pytest.approx(KERNEL.k)
