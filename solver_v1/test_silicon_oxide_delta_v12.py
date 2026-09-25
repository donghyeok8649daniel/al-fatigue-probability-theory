import numpy as np
from numpy.testing import assert_allclose
from solver_v1.silicon_oxide_delta_v12 import correction_features,radial_kernels


def test_pair_force_is_energy_gradient_and_conserves_force_torque():
    r=np.array([[0.,0.,0.],[1.6,.2,.1],[.1,1.8,-.1]])
    z=np.array([14,8,8]);cell=np.eye(3)*12;c=np.arange(1,9)/13
    energy,force=correction_features(z,r,cell,False)
    errors=[]
    for step in (1e-3,5e-4):
        fd=np.zeros_like(r)
        for i in range(3):
            for j in range(3):
                delta=np.zeros_like(r);delta[i,j]=step
                ep,_=correction_features(z,r+delta,cell,False);em,_=correction_features(z,r-delta,cell,False)
                fd[i,j]=-(ep-em)@c/(2*step)
        errors.append(np.max(abs(fd-force@c)))
    assert errors[1]<errors[0]/3.8
    assert_allclose((force@c).sum(axis=0),0,atol=2e-15)
    assert_allclose(np.cross(r,force@c).sum(axis=0),0,atol=2e-15)


def test_periodic_replication_and_rigid_rotation():
    from ase import Atoms
    a=Atoms(numbers=[14,8,8],positions=[[0,0,0],[1.5,.1,0],[.1,1.6,0]],cell=np.eye(3)*4,pbc=True)
    e,f=correction_features(a.numbers,a.positions,a.cell,a.pbc)
    repeated=a.repeat((2,1,1));er,fr=correction_features(repeated.numbers,repeated.positions,repeated.cell,repeated.pbc)
    assert_allclose(er,2*e,atol=1e-13)
    assert_allclose(fr,np.tile(f,(2,1,1)),atol=1e-13)
    theta=.71;c=np.cos(theta);s=np.sin(theta);rot=np.array([[c,-s,0],[s,c,0],[0,0,1.]])
    ez,fz=correction_features(a.numbers,a.positions@rot,a.cell.array@rot,a.pbc)
    assert_allclose(ez,e,atol=1e-13)
    assert_allclose(fz,np.einsum('nik,ij->njk',f,rot),atol=1e-13)


def test_kernel_support_and_pure_si_control():
    value,derivative=radial_kernels(np.array([.5,.6,3.5,3.6]))
    assert np.max(abs(value))==0 and np.max(abs(derivative))==0
    e,f=correction_features([14,14],[[0,0,0],[2.35,0,0]],np.eye(3)*8,False)
    assert np.max(abs(e))==0 and np.max(abs(f))==0
