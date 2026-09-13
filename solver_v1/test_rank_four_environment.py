import numpy as np
import pytest
from .rank_four_environment import RankFourEnvironment,stf_fourth,fourth_basis
from .coordination_screening import CoordinationScreenedBulk
from .vector_interface_reference import HESSIAN_INDICES


def direct_moment(R,k,C):
    w=C*np.exp(-k*np.linalg.norm(R,axis=1))
    raw=np.einsum('n,ni,nj,nk,nl->ijkl',w,R,R,R,R)
    return raw.ravel()@fourth_basis()


def test_rank_four_projection_and_rotation():
    B=fourth_basis();np.testing.assert_allclose(B.T@B,np.eye(9),atol=2e-15)
    np.testing.assert_allclose(np.einsum('c i i k l->c k l',B.T.reshape(9,3,3,3,3)),0,atol=1e-15)
    R=np.array([[1.,2.,.3],[-.4,1.,.8],[.7,-.6,1.2]])
    angle=.37;rotation=np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1.]])
    left=direct_moment(R,4.9,1.);right=direct_moment(R@rotation.T,4.9,1.)
    assert left@left==pytest.approx(right@right,rel=2e-14)
    single=direct_moment(R[:1],4.9,1.);r=np.linalg.norm(R[0])
    assert single@single==pytest.approx((8/35)*r**8*np.exp(-2*4.9*r),rel=2e-14)
    raw=np.arange(81.).reshape(3,3,3,3)
    projected=stf_fourth(raw)
    np.testing.assert_allclose(stf_fourth(projected),projected,atol=3e-14)


def test_rank_four_plane_poisson_and_analytic_derivatives():
    model=RankFourEnvironment(4.9);g=model.geometry
    m,n=np.meshgrid(np.arange(-18,19),np.arange(-18,19),indexing='ij')
    xy=m.ravel()[:,None]*g.a1+n.ravel()[:,None]*g.a2
    q=np.array([.87,.19,-.21])
    def evaluate(q):return model.plane(q[0],tuple(q[1:]))[0]
    jet=evaluate(q)
    R=np.column_stack([xy+q[1:],np.full(len(xy),q[0])])
    np.testing.assert_allclose(jet[0],direct_moment(R,model.k,model.C),atol=2e-12)
    h=1e-5
    for i,e in enumerate(np.eye(3)):
        upper,lower=evaluate(q+h*e),evaluate(q-h*e)
        np.testing.assert_allclose((upper[0]-lower[0])/(2*h),jet[1+i],atol=2e-8,rtol=3e-6)
        np.testing.assert_allclose((upper[1:4]-lower[1:4])/(2*h),jet[HESSIAN_INDICES[:,i]],atol=2e-7,rtol=3e-6)


def test_rank_four_bulk_background_and_elastic_mapping():
    model=RankFourEnvironment(4.9)
    R=CoordinationScreenedBulk((4.9,10.,11.,1000.,8.4,1.),radius=12.,law='power').R
    Q=model.bulk()[0];assert np.linalg.norm(Q[0])>1e-4
    np.testing.assert_allclose(Q[0],direct_moment(R,model.k,model.C),atol=2e-11)
    column=model.bulk_column();refined=model.bulk_column(2e-4)
    np.testing.assert_allclose(column,refined,atol=2e-7,rtol=1e-5)
    def energy(strain,mode):
        transformed=R.copy()
        if mode=='hydro':transformed*=1+strain
        if mode=='normal':transformed[:,2]*=1+strain
        if mode=='shear':transformed[:,0]+=strain*R[:,2]
        value=direct_moment(transformed,model.k,model.C)
        return value@value
    step=2e-4
    for i,mode in enumerate(('hydro','normal','shear')):
        x={j:energy(j*step,mode) for j in (-2,-1,0,1,2)}
        second=(-x[2]+16*x[1]-30*x[0]+16*x[-1]-x[-2])/(12*step**2)
        assert column[2+i]==pytest.approx(second,abs=3e-7,rel=2e-5)
    assert column[1]==pytest.approx(-energy(0.,'normal'),abs=1e-13)


def test_rank_four_interface_keeps_bulk_and_pressure_work():
    model=RankFourEnvironment(4.9)
    perfect,_=model.interface_jet((model.h,0.,0.))
    assert abs(perfect[0])<1e-25
    assert perfect[1]==pytest.approx(model.bulk_column()[0]/model.h,abs=3e-11)
    q=np.array([.879,.267,.113]);v,_=model.interface_jet(tuple(q));h=1e-5
    for i,e in enumerate(np.eye(3)):
        upper=model.interface_jet(tuple(q+h*e))[0];lower=model.interface_jet(tuple(q-h*e))[0]
        assert (upper[0]-lower[0])/(2*h)==pytest.approx(v[1+i],abs=2e-8,rel=2e-5)
        np.testing.assert_allclose((upper[1:4]-lower[1:4])/(2*h),v[HESSIAN_INDICES[:,i]],atol=2e-7,rtol=3e-6)
    repeat=model.interface_jet(tuple(q+[0.,1.,0.]))[0]
    np.testing.assert_allclose(repeat,v,atol=3e-11)


def test_rank_four_neighbor_analytic_jets():
    from .rank_four_environment import neighbor_fourth_jets
    R=np.array([[.2,-.8,1.1],[1.,0.,0.],[0.,.7,.3]])
    k=4.9;C=2.;value,J,H=neighbor_fourth_jets(R,k,C)
    for j,row in enumerate(R):
        np.testing.assert_allclose(value[j],direct_moment(row[None,:],k,C),atol=1e-16)
    step=1e-5
    for i,e in enumerate(np.eye(3)):
        hi=neighbor_fourth_jets(R+step*e,k,C);lo=neighbor_fourth_jets(R-step*e,k,C)
        np.testing.assert_allclose((hi[0]-lo[0])/(2*step),J[:,i],atol=2e-10,rtol=2e-6)
        np.testing.assert_allclose((hi[1]-lo[1])/(2*step),H[:,:,i],atol=2e-9,rtol=2e-6)


def test_rank_four_bloch_matches_actual_energy_with_background():
    from .rank_four_environment import rank_four_bulk_column
    model=RankFourEnvironment(4.9)
    bulk=CoordinationScreenedBulk((4.9,10.,11.,1000.,8.4,1.),radius=10.,law='power')
    H,tail=rank_four_bulk_column(model,bulk,(1/3,)*3)
    R=bulk.R;layer=np.rint(R[:,2]/bulk.geometry.h111)
    def energy(amplitude,direction):
        energies=[]
        for site in range(3):
            displacement=amplitude*(np.cos(2*np.pi*(site+layer)/3)-np.cos(2*np.pi*site/3))
            Q=direct_moment(R+displacement[:,None]*direction,model.k,model.C)
            energies.append(Q@Q)
        return np.mean(energies)
    for direction in np.eye(3):
        reference=energy(0.,direction)
        slopes=[2*(energy(h,direction)+energy(-h,direction)-2*reference)/h**2 for h in (4e-4,2e-4)]
        assert abs((4*slopes[1]-slopes[0])/3-direction@H@direction)<2e-6+tail
    fine=CoordinationScreenedBulk((4.9,10.,11.,1000.,8.4,1.),radius=12.,law='power')
    H2,tail2=rank_four_bulk_column(model,fine,(1/3,)*3)
    assert np.linalg.norm(H-H2)<=tail+tail2
