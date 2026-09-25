import numpy as np
import pytest
from numpy.testing import assert_allclose
from solver_v1.silicon_oxide_angular_delta_v12 import angular_features,combined_features


@pytest.mark.parametrize('numbers',[[14,8,8],[8,14,14],[14,14,8]])
def test_angular_force_energy_consistency(numbers):
    r=np.array([[0.,0.,0.],[1.6,.2,.1],[.1,1.8,-.2]])
    c=np.linspace(-.3,.7,9);cell=np.eye(3)*10
    e,f=angular_features(numbers,r,cell,False);errors=[]
    for h in (1e-3,5e-4):
        numerical=np.zeros_like(r)
        for i in range(3):
            for k in range(3):
                d=np.zeros_like(r);d[i,k]=h
                ep,_=angular_features(numbers,r+d,cell,False);em,_=angular_features(numbers,r-d,cell,False)
                numerical[i,k]=-(ep-em)@c/(2*h)
        errors.append(np.max(abs(numerical-f@c)))
    assert errors[1]<errors[0]/3.8
    assert_allclose((f@c).sum(axis=0),0,atol=1e-14)
    assert_allclose(np.cross(r,f@c).sum(axis=0),0,atol=1e-14)


def test_angular_periodic_rotation_replication_and_atom_order():
    from ase import Atoms
    a=Atoms(numbers=[14,8,14,8],positions=[[0,0,0],[1.4,.2,0],[.1,2.2,.1],[2.,2.,.2]],cell=np.eye(3)*4,pbc=True)
    e,f=combined_features(a.numbers,a.positions,a.cell,a.pbc)
    b=a.repeat((2,1,1));er,fr=combined_features(b.numbers,b.positions,b.cell,b.pbc)
    assert_allclose(er,2*e,atol=1e-12);assert_allclose(fr,np.tile(f,(2,1,1)),atol=1e-12)
    perm=np.array([2,0,3,1]);ep,fp=combined_features(a.numbers[perm],a.positions[perm],a.cell,a.pbc)
    assert_allclose(ep,e,atol=1e-13);assert_allclose(fp,f[perm],atol=1e-13)
    theta=.53;c=np.cos(theta);s=np.sin(theta);rot=np.array([[1.,0,0],[0,c,-s],[0,s,c]])
    ez,fz=combined_features(a.numbers,a.positions@rot,a.cell.array@rot,a.pbc)
    assert_allclose(ez,e,atol=1e-12);assert_allclose(fz,np.einsum('nik,ij->njk',f,rot),atol=1e-12)


@pytest.mark.parametrize('change',[np.diag([-1,1,1]),np.array([[1,1,0],[0,1,0],[0,0,1]])])
def test_equivalent_left_handed_and_sheared_periodic_bases(change):
    z=np.array([14,8,14,8]);r=np.array([[0,0,0],[1.4,.2,0],[.1,2.2,.1],[2.,2.,.2]])
    cell=np.array([[4.,0,0],[.5,4.1,0],[.2,.3,4.2]])
    e,f=combined_features(z,r,cell,True)
    changed_e,changed_f=combined_features(z,r,change@cell,True)
    assert_allclose(changed_e,e,atol=1e-12)
    assert_allclose(changed_f,f,atol=1e-12)
