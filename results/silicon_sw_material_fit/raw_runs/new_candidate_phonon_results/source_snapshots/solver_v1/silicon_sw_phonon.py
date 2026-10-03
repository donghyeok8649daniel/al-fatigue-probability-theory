"""Exact local harmonic Bloch matrix from centered-site energy Hessians.

No MD, friction or overdamped time conversion. Optical/acoustic frequencies
are Newtonian harmonic diagnostics with the explicitly stated atomic mass.
Neighbor displacements keep the original SW support and all basis phases.
"""
import numpy as np
from .silicon_environment_research import DiamondCell,_state,_finish
from .silicon_sw_anharmonic import angular_site

SI_MASS_AMU=28.0855
EV_A2_AMU_TO_THZ=np.sqrt(1.602176634e-19/(1e-20*1.66053906660e-27))/(2*np.pi*1e12)


class SWBlochMatrix:
    def __init__(self,parameters,lattice,*,angle_beta=0.):
        self.bulk=DiamondCell(parameters,lattice);self.beta=angle_beta;self.sites=[]
        for center,environment in enumerate(self.bulk.environments):
            entries=[(v,center+sign) for v,sign in environment if np.linalg.norm(v)<parameters['sigma']*parameters['a']]
            vectors=np.array([v for v,_ in entries]);q=_state(vectors.ravel(),vectors.size,True)
            jet=_finish(angular_site([q[i:i+3] for i in range(0,len(q),3)],parameters,angle_beta),vectors.size,True)
            self.sites.append((center,entries,jet.hessian))

    def evaluate(self,fractional_q):
        q=np.asarray(fractional_q,float)
        if q.shape!=(3,) or not np.all(np.isfinite(q)):raise ValueError('finite primitive reciprocal q required')
        k=q@(2*np.pi*np.linalg.inv(self.bulk.cell).T);result=np.zeros((6,6),complex)
        for center,entries,hessian in self.sites:
            b=np.zeros((3*len(entries),6),complex)
            for j,(v,basis) in enumerate(entries):
                b[3*j:3*j+3,3*basis:3*basis+3]+=np.exp(1j*k@v)*np.eye(3)
                b[3*j:3*j+3,3*center:3*center+3]-=np.eye(3)
            result+=b.conj().T@hessian@b
        if np.max(abs(result-result.conj().T))>1e-10:
            raise FloatingPointError('Bloch harmonic matrix is not Hermitian')
        return result

    def spectrum(self,q,*,mass_amu=SI_MASS_AMU):
        if not np.isfinite(mass_amu) or mass_amu<=0:raise ValueError('positive atomic mass required')
        eigenvalues=np.linalg.eigvalsh(self.evaluate(q))
        # Signed square roots expose negative curvature; no stability clipping.
        frequencies=np.sign(eigenvalues)*np.sqrt(abs(eigenvalues)/mass_amu)*EV_A2_AMU_TO_THZ
        return eigenvalues,frequencies
