"""External Mishin reference: generate LIKE-FOR-LIKE static targets only.

This tabulated EAM is NOT an energy model selectable by the probability PDE.
The canonical analytic LJ pair is not replaced. Its published finite cutoff
defines the *reference potential*, not a cutoff approximation of our LJ sum.
Coordinates here are dimensional angstrom; energies are eV. No kinetics.
"""
from __future__ import annotations

import hashlib
from pathlib import Path
import argparse

import numpy as np
from scipy.interpolate import CubicSpline

from .fcc111_geometry import fcc111_geometry, registry_path, DIRECT_110


SOURCE_URL = ("https://www.ctcms.nist.gov/potentials/Download/"
    "1999--Mishin-Y-Farkas-D-Mehl-M-J-Papaconstantopoulos-D-A--Al/2/Al99.eam.alloy")
SOURCE_SHA256 = "60c8a085be79d273324ab421f5b1447578fef55c1acfc6492c0999f15ee8a284"
DEFAULT_CACHE = Path(__file__).resolve().parents[1]/".cache/al-reference/Al99.eam.alloy"


class MishinRigidFCCReference:
    """Published source data on the SAME rigid FCC half-crystal geometry.

    Interpolated source functions are used only in this isolated reference
    evaluator. They are not a newly fitted model or a tabulated production path.
    """

    def __init__(self, path=DEFAULT_CACHE, *, lattice_angstrom=4.05,
                 interpolation="cubic", verify_hash=True):
        path=Path(path)
        data=path.read_bytes()
        self.sha256=hashlib.sha256(data).hexdigest()
        if verify_hash and self.sha256!=SOURCE_SHA256:
            raise ValueError("reference file checksum does not match documented NIST release")
        lines=data.decode("ascii").splitlines()
        nrho,drho,nr,dr,cutoff=map(float,lines[4].split())
        nrho,nr=int(nrho),int(nr)
        values=np.fromstring(" ".join(lines[6:]),sep=" ")
        if len(values)!=nrho+2*nr or lines[3].split()!=["1","Al"]:
            raise ValueError("expected single-Al NIST setfl layout")
        self.cutoff=float(cutoff)
        self.r=np.arange(nr)*dr
        self.rho=np.arange(nrho)*drho
        self.embedding=values[:nrho]
        self.density=values[nrho:nrho+nr]
        self.r_pair=values[nrho+nr:]
        self.interpolation=interpolation
        if interpolation not in ("cubic","linear"):
            raise ValueError("unknown reference interpolation")
        self._F=CubicSpline(self.rho,self.embedding,extrapolate=False)
        self._rho=CubicSpline(self.r,self.density,extrapolate=False)
        self._rphi=CubicSpline(self.r,self.r_pair,extrapolate=False)
        self.geometry=fcc111_geometry(lattice_angstrom)
        self.h=self.geometry.h111
        # Guaranteed to enclose the published cutoff even for a cell shift.
        radius=int(np.ceil(2*self.cutoff/self.geometry.b))+4
        m,n=np.meshgrid(np.arange(-radius,radius+1),
                         np.arange(-radius,radius+1),indexing="ij")
        self.R=m.ravel()[:,None]*self.geometry.a1+n.ravel()[:,None]*self.geometry.a2
        self._rho_bulk,self._pair_bulk=self._bulk_components()

    def F(self,rho):
        x=np.asarray(rho,dtype=float)
        if np.any((x<0)|(x>self.rho[-1])):
            raise ValueError("reference embedding density outside tabulated range")
        if self.interpolation=="linear":
            return np.interp(x,self.rho,self.embedding)
        return self._F(x)

    def plane(self,d,delta,*,exclude_self=False):
        radial=np.sqrt(float(d)**2+np.sum((self.R+delta)**2,axis=1))
        select=radial<self.r[-1]
        if exclude_self:
            select &= radial>1e-12
        r=radial[select]
        if np.any(r<=0):
            raise ValueError("self atom included in reference pair sum")
        if self.interpolation=="linear":
            den=np.interp(r,self.r,self.density)
            pair=np.interp(r,self.r,self.r_pair)/r
        else:
            den=self._rho(r)
            pair=self._rphi(r)/r
        return float(np.sum(pair)),float(np.sum(den))

    def _bulk_components(self):
        pair,density=self.plane(0,np.zeros(2),exclude_self=True)
        pair*=.5
        for k in range(1,int(np.ceil(self.cutoff/self.h))+1):
            wp,rho=self.plane(k*self.h,self.geometry.abc_shift(k))
            pair+=wp; density+=2*rho
        return density,pair

    def bulk_energy(self):
        return self._pair_bulk+float(self.F(self._rho_bulk))

    def cohesion(self):
        return float(self.F(0))-self.bulk_energy()

    def interface_energy(self,a_angstrom,s_angstrom,*,path_id=DIRECT_110):
        """One cut; explicit per-depth atom density, per interface-cell energy."""
        if a_angstrom<=0:
            raise ValueError("positive physical interface spacing required")
        slip=float(s_angstrom)*registry_path(path_id).direction()
        count=int(np.ceil((self.cutoff+abs(float(a_angstrom)-self.h))/self.h))+1
        differences=[]; pair=0.
        for k in range(1,count+1):
            baseline=self.plane(k*self.h,self.geometry.abc_shift(k))
            active=self.plane(a_angstrom+(k-1)*self.h,
                              self.geometry.abc_shift(k)+slip)
            pair+=k*(active[0]-baseline[0])
            differences.append(active[1]-baseline[1])
        delta=np.cumsum(differences[::-1])[::-1]
        embedding=2*np.sum(self.F(self._rho_bulk+delta)-self.F(self._rho_bulk))
        return float(pair+embedding)

    def work_separation(self):
        return self.interface_energy(self.cutoff+self.h,0)

    def _plane_derivatives(self,d,delta,direction):
        """Analytic derivatives of the interpolated SOURCE, validation only.

        Return value,a,s,aa,as,ss for the pair and density plane sums. This is
        not a finite-difference drift or a replacement for analytic LJ/Bessel.
        """
        if self.interpolation!="cubic":
            raise ValueError("source derivative checks require cubic interpolation")
        xy=self.R+delta
        radii=np.sqrt(float(d)**2+np.sum(xy*xy,axis=1))
        select=(radii>0)&(radii<self.r[-1])
        r=radii[select]; along=xy[select]@direction
        ra=d/r; rs=along/r
        raa=1/r-d*d/r**3; ras=-d*along/r**3; rss=1/r-along**2/r**3
        z=self._rphi(r); zp=self._rphi(r,1); zpp=self._rphi(r,2)
        def pack(value,first,second):
            return np.array([np.sum(value),np.sum(first*ra),np.sum(first*rs),
                np.sum(second*ra*ra+first*raa),np.sum(second*ra*rs+first*ras),
                np.sum(second*rs*rs+first*rss)])
        return (pack(z/r,zp/r-z/r**2,zpp/r-2*zp/r**2+2*z/r**3),
                pack(self._rho(r),self._rho(r,1),self._rho(r,2)))

    def interface_derivatives(self,a_angstrom,s_angstrom,*,path_id=DIRECT_110):
        """Source-only [E,E_a,E_s,E_aa,E_as,E_ss], eV and angstrom.

        This intentionally does not implement the probability energy-surface
        interface. Source interpolation has a published cutoff; tests compare
        derivatives away from neighbor-cutoff crossings and refine stencils.
        """
        if not np.isfinite(a_angstrom) or a_angstrom<=0:
            raise ValueError("finite positive interface spacing required")
        direction=registry_path(path_id).direction()
        slip=float(s_angstrom)*direction
        count=int(np.ceil((self.cutoff+abs(float(a_angstrom)-self.h))/self.h))+1
        pair=np.zeros(6); changes=[]
        for k in range(1,count+1):
            baseline=self.plane(k*self.h,self.geometry.abc_shift(k))
            wp,rho=self._plane_derivatives(a_angstrom+(k-1)*self.h,
                self.geometry.abc_shift(k)+slip,direction)
            wp[0]-=baseline[0]; rho[0]-=baseline[1]
            pair+=k*wp; changes.append(rho)
        depths=np.cumsum(np.asarray(changes)[::-1],axis=0)[::-1]
        emb=np.zeros(6)
        for drho in depths:
            density=self._rho_bulk+drho[0]
            first=float(self._F(density,1)); second=float(self._F(density,2))
            emb+=2*np.array([float(self.F(density)-self.F(self._rho_bulk)),
                first*drho[1],first*drho[2],second*drho[1]**2+first*drho[3],
                second*drho[1]*drho[2]+first*drho[4],
                second*drho[2]**2+first*drho[5]])
        return pair+emb


def download_reference(path=DEFAULT_CACHE):
    """Explicit opt-in source download with checksum and no overwrite."""
    from urllib.request import urlopen
    destination=Path(path)
    if destination.exists():
        if hashlib.sha256(destination.read_bytes()).hexdigest()!=SOURCE_SHA256:
            raise ValueError("existing reference differs; preserved, not overwritten")
        return destination
    with urlopen(SOURCE_URL,timeout=45) as response:
        data=response.read()
    if hashlib.sha256(data).hexdigest()!=SOURCE_SHA256:
        raise ValueError("download checksum mismatch; no file written")
    destination.parent.mkdir(parents=True,exist_ok=True)
    with destination.open("xb") as stream:
        stream.write(data)
    return destination


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download",action="store_true",help="fetch and verify NIST source; never replace local data")
    args=parser.parse_args()
    if args.download:
        download_reference()
    reference=MishinRigidFCCReference()
    print("Source-only reference SHA256:",reference.sha256)
    print("Cohesion [eV/atom]:",reference.cohesion())
    print("Rigid work of separation [J/m2]:",reference.work_separation()*16.02176634/reference.geometry.atomic_cell_area)
