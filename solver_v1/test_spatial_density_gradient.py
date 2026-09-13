import numpy as np
import pytest
from .spatial_density_gradient import SpatialDensityGradient,gradient_plane,gradient_bulk_column
from .coordination_screening import CoordinationScreenedInterface,CoordinationScreenedBulk,screened_site_jet
from .vector_interface_reference import HESSIAN_INDICES

SHAPE=(4.9,10.,11.,1000.,8.4,1.)


def test_plane_gradient_third_tensor_against_direct_lattice():
    environment=SpatialDensityGradient(SHAPE[0]);kernel=environment.kernel;g=environment.geometry
    d=.81;delta=np.array([.17,-.23]);jet,_=gradient_plane(kernel,d,delta,g)
    m,n=np.meshgrid(np.arange(-16,17),np.arange(-16,17),indexing='ij')
    xy=m.ravel()[:,None]*g.a1+n.ravel()[:,None]*g.a2+delta
    R=np.column_stack([np.full(len(xy),d),xy]);r=np.linalg.norm(R,axis=1);e=R/r[:,None]
    f,fp,fpp,fppp=[kernel.radial(r,i) for i in range(4)]
    H=(fpp-fp/r)[:,None,None]*e[:,:,None]*e[:,None,:]+(fp/r)[:,None,None]*np.eye(3)
    T=np.einsum('n,ni,nj,nk->nijk',fppp-3*fpp/r+3*fp/r**2,e,e,e)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                T[:,i,j,k]+=(fpp/r-fp/r**2)*((i==j)*e[:,k]+(i==k)*e[:,j]+(j==k)*e[:,i])
    np.testing.assert_allclose(jet[0],-np.sum(fp[:,None]*e,axis=0),atol=2e-12)
    np.testing.assert_allclose(jet[1:4],-H.sum(axis=0),atol=3e-12)
    np.testing.assert_allclose(jet[HESSIAN_INDICES],-T.sum(axis=0),atol=3e-11)
    Q=np.sum(f[:,None]*R,axis=0)
    assert np.linalg.norm(jet[0]-Q*(Q@jet[0])/(Q@Q))>1e-4


def test_actual_density_gradient_interface_energy_derivatives_and_reference():
    environment=SpatialDensityGradient(SHAPE[0]);base=CoordinationScreenedInterface(SHAPE,np.ones(10),law='power')
    def jet(q):
        density,_,_=base.site_inputs(tuple(q));gradients,_=environment.site_jets(tuple(q))
        return 2*screened_site_jet(density,gradients,-1.,law='power')
    q=np.array([.879,.267,.113]);v=jet(q);h=1e-5
    hi=[jet(q+h*e) for e in np.eye(3)];lo=[jet(q-h*e) for e in np.eye(3)]
    np.testing.assert_allclose([(a[0]-b[0])/(2*h) for a,b in zip(hi,lo)],v[1:4],atol=2e-7,rtol=2e-6)
    np.testing.assert_allclose(np.column_stack([(a[1:4]-b[1:4])/(2*h) for a,b in zip(hi,lo)]),
                               v[HESSIAN_INDICES],atol=3e-6,rtol=3e-6)
    perfect=jet((base.h,0.,0.));assert abs(perfect[0])<1e-24
    np.testing.assert_allclose(perfect[1:4],0.,atol=1e-12)
    np.testing.assert_allclose(jet(q+[0.,1.,0.]),v,atol=2e-11)


def test_spatial_gradient_bloch_matches_actual_site_energy_and_acoustic_limit():
    bulk=CoordinationScreenedBulk(SHAPE,radius=12.,law='power')
    env=SpatialDensityGradient(SHAPE[0]);kernel=env.kernel;R=bulk.R
    layers=np.rint(R[:,2]/bulk.geometry.h111);q=(1/3,)*3
    H,tail=gradient_bulk_column(bulk,q)
    def energy(amplitude,v):
        values=[]
        for site in range(6):
            shift=amplitude*(np.cos(2*np.pi*(site+layers)/3)-np.cos(2*np.pi*site/3))
            pos=R+shift[:,None]*v;r=np.linalg.norm(pos,axis=1)
            density=kernel.radial(r).sum();g=-(kernel.radial(r,1)/r)@pos
            values.append((g@g)/density)
        return np.mean(values)
    for v in np.eye(3):
        slopes=[2*(energy(h,v)+energy(-h,v)-2*energy(0.,v))/h**2 for h in (2e-4,1e-4)]
        assert abs((4*slopes[1]-slopes[0])/3-v@H@v)<3e-6+tail
    one,_=gradient_bulk_column(bulk,(1e-4,1e-4,0.))
    half,_=gradient_bulk_column(bulk,(5e-5,5e-5,0.))
    assert np.linalg.norm(one)/np.linalg.norm(half)==pytest.approx(16.,rel=1e-6)
