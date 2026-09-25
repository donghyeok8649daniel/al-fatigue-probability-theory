import numpy as np
import pytest
from solver_v1.silicon_structure_diagnostics_v12 import local_affine,edge_bottleneck,neighbor_lists,pair_set,reference_sections
from results.silicon_wafer_feasibility.check_prism_dense_hessian_v12 import constraint_basis


def test_affine_and_frame_covariance():
    r=np.array([[0,0,0],[1,1,1],[-1,-1,1],[-1,1,-1],[1,-1,-1]],float)
    nb=[[1,2,3,4],[],[],[],[]]
    a=np.array([[1.12,.07,0],[0,.93,.02],[.01,0,1.04]])
    q=r@a+np.array([4.,-3.,7.])
    result=local_affine(r,q,nb)[0]
    assert result['D2_A2']<1e-28
    assert np.allclose(result['green_principal_strains'],np.linalg.eigvalsh((a@a.T-np.eye(3))/2))
    q[1]+=[.2,-.1,.05]
    reference=local_affine(r,q,nb)[0]
    theta=.713;c=np.cos(theta);s=np.sin(theta);rot=np.array([[c,-s,0],[s,c,0],[0,0,1.]])
    changed=local_affine(r@rot+3,q@rot-2,nb)[0]
    assert changed['D2_A2']==pytest.approx(reference['D2_A2'],rel=1e-12)
    assert np.allclose(changed['green_principal_strains'],reference['green_principal_strains'])


def test_rank_deficient_neighbors_have_no_strain():
    r=np.array([[0,0,0],[1,0,0],[0,1,0]],float)
    result=local_affine(r,r,[[1,2],[0],[0]])
    assert result[0]['rank']==2
    assert result[0]['green_principal_strains'] is None


def test_bottleneck_and_atom_permutation():
    pairs={(0,1),(1,3),(0,2),(2,3)};low=np.array([1,0,0,0],bool);high=low[::-1]
    assert edge_bottleneck(4,pairs,low,high)['edge_disjoint_paths']==2
    assert edge_bottleneck(4,pairs-{(1,3)},low,high)['edge_disjoint_paths']==1
    assert edge_bottleneck(4,pairs-{(1,3),(2,3)},low,high)['edge_disjoint_paths']==0
    permutation=np.array([2,0,3,1]);remapped={tuple(sorted((permutation[i],permutation[j]))) for i,j in pairs}
    l=np.zeros(4,bool);u=l.copy();l[permutation[0]]=True;u[permutation[3]]=True
    assert edge_bottleneck(4,remapped,l,u)['edge_disjoint_paths']==2


def test_pair_count_and_sections():
    r=np.array([[0,0,0],[0,0,1],[0,0,2],[0,0,3]],float);pairs=pair_set(r,1.1)
    assert pairs=={(0,1),(1,2),(2,3)}
    assert [len(x) for x in neighbor_lists(4,pairs)]==[1,2,2,1]
    q=r.copy();q[:,2]*=1.5
    rows=reference_sections(r,q,pairs,pairs,np.ones(4,bool))
    assert all(x['current_mean_layer_gap_A']==1.5 for x in rows)
    assert all(x['lost_initial_pairs']==0 for x in rows)


def test_grip_basis_energy_and_force_chain_rule():
    lower=np.array([1,1,0,0,0,0],bool);upper=np.array([0,0,0,0,1,1],bool);free=~(lower|upper)
    b,norm=constraint_basis(free,lower,upper,'force')
    assert np.allclose(b.T@b,np.eye(7))
    rng=np.random.default_rng(54);x=rng.normal(size=18);a=rng.normal(size=(18,18));h=a.T@a
    direction=rng.normal(size=7);shift=b@direction;step=1e-4
    reduced=b.T@h@b
    fd=b.T@(h@(x+step*shift)-h@(x-step*shift))/(2*step)
    assert np.allclose(fd,reduced@direction,rtol=1e-10,atol=1e-10)
    grip=(b[:,-1]*norm).reshape(-1,3)
    assert np.all(grip[lower,2]==-.5) and np.all(grip[upper,2]==.5)
    assert np.all(grip[free]==0)
    fixed,_=constraint_basis(free,lower,upper,'displacement')
    assert np.array_equal(fixed,b[:,:-1])
