"""Bloch Hessian versus independent conservative force response."""
import numpy as np
from ase.build import bulk
from results.silicon_wafer_feasibility.run_static_probe import source_parameters
from .silicon_sw_phonon import SWBlochMatrix
from .silicon_sw_material_fit import force_features


def test_translation_hermiticity_and_full_wavevector_angular_invariance():
    p,_=source_parameters();plain=SWBlochMatrix(p,5.431);modified=SWBlochMatrix(p,5.431,angle_beta=.375)
    gamma=plain.evaluate([0.,0.,0.])
    for axis in range(3):
        translation=np.zeros(6);translation[axis]=translation[axis+3]=1.
        np.testing.assert_allclose(gamma@translation,0.,atol=3e-12)
    for q in [[0.,0.,0.],[0.,.5,.5],[.5,.5,.5],[.11,.19,.27]]:
        matrix=plain.evaluate(q)
        np.testing.assert_allclose(matrix,matrix.conj().T,atol=3e-12)
        np.testing.assert_allclose(modified.evaluate(q),matrix,atol=3e-9)
        values,_=plain.spectrum(q)
        assert np.min(values)>-1e-11


def test_commensurate_bloch_matrix_matches_supercell_force_derivative():
    p,_=source_parameters();lattice=5.431;model=SWBlochMatrix(p,lattice,angle_beta=.375)
    primitive=bulk('Si','diamond',a=lattice);atoms=primitive.repeat((4,4,4));x=atoms.positions
    np.testing.assert_allclose(primitive.cell.array,model.bulk.cell,atol=1e-14)
    basis=np.tile([0,1],64);u=np.array([1.,.3,-.2,-.7,.11,.4])+1j*np.array([.1,-.2,.3,.4,-.17,-.1])
    for q in [[0.,0.,0.],[0.,.5,.5],[.5,.5,.5]]:
        k=np.array(q)@(2*np.pi*np.linalg.inv(primitive.cell.array).T);phase=np.exp(1j*(x@k));field=phase[:,None]*u.reshape(2,3)[basis]
        h=2e-5;response=[]
        for component in [field.real,field.imag]:
            fp,fm=[force_features(x+sign*h*component,atoms.cell.array,atoms.pbc,p,angle_beta=.375)[1].sum(axis=2) for sign in [1,-1]]
            response.append(-(fp-fm)/(2*h))
        complex_force=response[0]+1j*response[1]
        projected=np.array([(complex_force[basis==b]*phase[basis==b,None].conj()).mean(axis=0) for b in range(2)]).ravel()
        np.testing.assert_allclose(projected,model.evaluate(q)@u,atol=2e-7,rtol=2e-7)


def test_longer_support_bloch_response_keeps_angular_bulk_background():
    p,_=source_parameters();p=dict(p,a=4.0/p['sigma']);lattice=5.461
    plain=SWBlochMatrix(p,lattice);modified=SWBlochMatrix(p,lattice,angle_beta=.375)
    q=[0.,.5,.5]
    assert np.max(abs(plain.evaluate(q)-modified.evaluate(q)))>1e-6
    primitive=bulk('Si','diamond',a=lattice);atoms=primitive.repeat((4,4,4));x=atoms.positions
    basis=np.tile([0,1],64);u=np.array([1.,.3,-.2,-.7,.11,.4]);k=np.array(q)@(2*np.pi*np.linalg.inv(primitive.cell.array).T)
    phase=np.exp(1j*(x@k));field=phase[:,None]*u.reshape(2,3)[basis];h=2e-5;responses=[]
    for component in [field.real,field.imag]:
        fp,fm=[force_features(x+sign*h*component,atoms.cell.array,atoms.pbc,p,angle_beta=.375)[1].sum(axis=2) for sign in [1,-1]]
        responses.append(-(fp-fm)/(2*h))
    complex_response=responses[0]+1j*responses[1]
    projected=np.array([(complex_response[basis==b]*phase[basis==b,None].conj()).mean(axis=0) for b in range(2)]).ravel()
    np.testing.assert_allclose(projected,modified.evaluate(q)@u,atol=3e-7,rtol=3e-7)
