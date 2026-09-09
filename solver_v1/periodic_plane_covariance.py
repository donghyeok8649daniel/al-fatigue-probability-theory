"""Equal-time normalization of measured periodic (111) plane averages.

This is STATIC harmonic equipartition, not a clock or new dynamics. Source
tabulated EAM is a validator only; the analytic LJ/Bessel model is untouched.
"""
from types import SimpleNamespace

import numpy as np

from .aluminum_calibration import EV_J
from .static_bulk_stability import neighbors_in_plane_frame


KB_J_K=1.380649e-23


class SourceEAMBlochHessian:
    """Published finite-cutoff EAM harmonic matrix [eV/Angstrom^2].

    All neighbors and the per-atom F'' collective variation are included.
    No atomic mass appears and eigenvalues are NOT called frequencies.
    """
    def __init__(self,reference):
        self.reference=reference
        bulk=SimpleNamespace(geometry=reference.geometry,a0=reference.h)
        self.R=neighbors_in_plane_frame(bulk,reference.r[-1]*(1-1e-12))
        r=np.linalg.norm(self.R,axis=1); e=self.R/r[:,None]
        rr=e[:,:,None]*e[:,None,:]
        z=reference._rphi(r);zp=reference._rphi(r,1);zpp=reference._rphi(r,2)
        phi1=zp/r-z/r**2;phi2=zpp/r-2*zp/r**2+2*z/r**3
        rho=reference._rho(r);rho1=reference._rho(r,1);rho2=reference._rho(r,2)
        self.rho_bulk=float(np.sum(rho))
        if abs(self.rho_bulk-reference._rho_bulk)>2e-12:
            raise ArithmeticError('source density counting differs from plane sum')
        first=float(reference._F(self.rho_bulk,1))
        self.second=float(reference._F(self.rho_bulk,2))
        effective1=phi1+2*first*rho1;effective2=phi2+2*first*rho2
        self.pair_hessian=effective2[:,None,None]*rr+(effective1/r)[:,None,None]*(np.eye(3)-rr)
        self.density_gradient=rho1[:,None]*e

    def evaluate(self,q_cubic_per_angstrom):
        q=self.reference.geometry.plane_basis_in_stacked_cubic_axes()@np.asarray(q_cubic_per_angstrom,float)
        phase=self.R@q
        g=np.sin(phase)@self.density_gradient
        H=np.einsum('n,nij->ij',2*np.sin(phase/2)**2,self.pair_hessian)+self.second*np.outer(g,g)
        return (H+H.T)/2

    def direct_sinusoidal_energy(self, *, planes, mode, polarization_plane, amplitude_angstrom):
        """Independent actual source EAM energy of a periodic plane wave.

        Average per-atom site energies over all plane classes, evaluating F
        AFTER each full neighbor density. The source's published cutoff is
        used; it is not a new approximation of canonical infinite LJ energy.
        A neighbor buffer encloses every possible shifted source interaction.
        """
        polarization=np.asarray(polarization_plane,float)
        if (int(planes)!=planes or planes<3 or int(mode)!=mode or not 0<mode<planes
                or polarization.shape!=(3,) or np.any(~np.isfinite(polarization))
                or not np.isclose(polarization@polarization,1.,rtol=1e-12,atol=1e-12)
                or not np.isfinite(amplitude_angstrom)):
            raise ValueError('periodic integer mode, unit polarization and finite amplitude required')
        source=self.reference
        bulk=SimpleNamespace(geometry=source.geometry,a0=source.h)
        R=neighbors_in_plane_frame(bulk,source.r[-1]+2*abs(amplitude_angstrom)+1e-5)
        layer=np.rint(R[:,2]/source.h).astype(int)
        site=np.arange(planes)[:,None]
        phase=2*np.pi*mode/planes
        difference=amplitude_angstrom*(np.cos(phase*(site+layer))-np.cos(phase*site))
        displaced=R[None,:,:]+difference[:,:,None]*polarization
        r=np.linalg.norm(displaced,axis=-1)
        inside=r<source.r[-1]
        density=np.zeros_like(r);pair=np.zeros_like(r)
        density[inside]=source._rho(r[inside])
        pair[inside]=source._rphi(r[inside])/r[inside]
        return float(np.mean(.5*pair.sum(axis=1)+source.F(density.sum(axis=1))))


def plane_mode_equipartition(hessians_J_m2, *, planes, atoms_per_plane, temperature_K):
    r"""C_{difference,k}=4sin^2(pi k/N) kBT/(Np) H(k)^-1 [m^2].

    N is the number of periodic plane classes, Np their atom count. With a
    unitary N-point DFT the harmonic energy is Np/2 sum u_k^* H(k)u_k.
    Thus real-space equal-time covariance is sum C_difference,k / N.
    The translation k=0 drops out; it is never inverted or regularized.
    """
    H=np.asarray(hessians_J_m2,float)
    if (int(planes)!=planes or planes<2 or int(atoms_per_plane)!=atoms_per_plane or atoms_per_plane<1
            or not np.isfinite(temperature_K) or temperature_K<=0
            or H.shape!=(planes-1,3,3) or np.any(~np.isfinite(H))):
        raise ValueError('explicit periodic counts, positive T and nonzero-mode 3x3 Hessians required')
    if not np.allclose(H,H.swapaxes(1,2),rtol=1e-12,atol=1e-12):
        raise ValueError('symmetric physical Hessians required')
    if np.min(np.linalg.eigvalsh(H))<=0:
        raise ValueError('stable nonzero harmonic modes required')
    weight=4*np.sin(np.pi*np.arange(1,planes)/planes)**2
    modes=weight[:,None,None]*KB_J_K*temperature_K/atoms_per_plane*np.linalg.inv(H)
    return modes,np.sum(modes,axis=0)/planes


def empirical_plane_modes(coordinates_m):
    q=np.asarray(coordinates_m,float)
    if q.ndim!=3 or q.shape[1]<2 or q.shape[2]!=3 or len(q)<4 or np.any(~np.isfinite(q)):
        raise ValueError('finite (time,periodic plane,normal/slip/transverse) displacements required')
    q=q-q.mean(axis=0,keepdims=True)
    modes=np.fft.fft(q,axis=1,norm='ortho')
    C=np.einsum('tki,tkj->kij',modes,modes.conj())/len(q)
    return C,np.einsum('tki,tkj->ij',q,q)/(len(q)*q.shape[1])


def source_predicted_covariance(reference, *, repeats, temperature_K):
    operator=SourceEAMBlochHessian(reference)
    box_angstrom=repeats*reference.geometry.lattice_constant
    H=np.array([operator.evaluate(np.ones(3)*2*np.pi*k/box_angstrom) for k in range(1,repeats)])
    # Actual +ABC stacked cubic axes give (-e1,-e2,e3) as the plane basis.
    # Experimental projection order is (e3,+e1,+e2), hence these two signs.
    transform=np.array([[0.,0.,1.],[-1.,0.,0.],[0.,-1.,0.]])
    H=np.einsum('ai,kij,bj->kab',transform,H,transform)*EV_J/1e-20
    modes,local=plane_mode_equipartition(H,planes=repeats,
        atoms_per_plane=4*repeats**2,temperature_K=temperature_K)
    return dict(mode_covariances_m2=modes,local_covariance_m2=local,
        hessians_J_m2=H,atoms_per_plane=4*repeats**2,
        source_harmonic_temperature_effect='geometry at source finite-T box; Hessian is static, no anharmonic renormalization')
