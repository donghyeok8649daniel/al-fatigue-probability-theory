"""Auxiliary pair/EAM gauge audit, NOT fitting a replacement pair potential.

The source density column is an idealized target-only comparison; it is never
used as a canonical kernel. An auxiliary-gauge fit is not material validation.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.optimize import lsq_linear
from .reference_eam_targets import MishinRigidFCCReference,SOURCE_URL
from .polynomial_exponential_density import QuadraticEnvelopeDensity
from .run_vector_registry_audit import save_json
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_low_frequency_forcing_v29 import sha


def run(density_fit,out):
    density_fit,out=Path(density_fit),Path(out)
    if out.exists():
        raise FileExistsError('fresh gauge diagnostic required')
    source=MishinRigidFCCReference();L=source.geometry.b
    fits=json.loads((density_fit/'fits.json').read_bytes())
    best=min((p for p in fits if p['name'].startswith('quadratic_')),key=lambda p:p['squared_loss'])
    kernel=QuadraticEnvelopeDensity(**best['parameters'])
    train=np.arange(2.1,6.101,.25);validation=np.arange(2.225,6.1,.25)
    def phi(r):
        return source._rphi(r)/r
    def matrix(r,kind,order=0):
        x=r/L
        powers=(x**-12,-x**-6) if order==0 else (-12*x**-13/L,6*x**-7/L) if order==1 else (156*x**-14/L**2,-42*x**-8/L**2)
        columns=list(powers)
        if kind=='source_density_gauge':
            columns.append(2*source._rho(r,order))
        if kind=='quadratic_density_gauge':
            columns.append(2*kernel.radial(r,order))
        return np.column_stack(columns)
    target=phi(train);scale=np.maximum(.05*abs(target),.005)
    out.mkdir(parents=True)
    save_json(out/'definition.json',dict(source_sha256=source.sha256,source_url=SOURCE_URL,
        auxiliary_density_fit_sha256=sha(density_fit/'fits.json'),fit_radii_A=train,validation_radii_A=validation,
        objective='auxiliary pair values only, scale=max(.05*abs(phi_source),.005 eV); NOT experimental uncertainty',
        gauge_identity='phi_LJ+2 B rho with F_new=F_source+B rho represents source energy only if phi_source=phi_LJ+2 B rho',
        derivative_units=['eV','eV/Angstrom','eV/Angstrom^2'],source_density_is_target_only=True,
        material_calibration=False,production_changed=False))
    profiles=[];rows=[]
    for name in ('LJ_only','source_density_gauge','quadratic_density_gauge'):
        M=matrix(train,name);J=M/scale[:,None];D=1/np.linalg.norm(J,axis=0)
        lower=np.r_[0.,0.,[-np.inf] if M.shape[1]==3 else []]
        fit=lsq_linear(J*D,target/scale,bounds=(lower,np.full(M.shape[1],np.inf)),
            tol=1e-12,max_iter=1000)
        c=D*fit.x
        profiles.append(dict(name=name,success=bool(fit.success),coefficients=c,
            squared_loss=np.sum(((M@c-target)/scale)**2),singular_values=np.linalg.svd(J*D,compute_uv=False),
            strictly_positive_float_LJ=bool(np.all(c[:2]>0)),optimizer_active_mask=fit.active_mask,
            LJ_away_from_optimizer_bound=bool(np.all(fit.active_mask[:2]==0)),
            coefficient_units='c12,c6 in eV at length L0; optional B in eV/source auxiliary density'))
        for role,radii in [('fit',train),('excluded',validation)]:
            z=source._rphi(radii);zp=source._rphi(radii,1);zpp=source._rphi(radii,2)
            references=(z/radii,zp/radii-z/radii**2,zpp/radii-2*zp/radii**2+2*z/radii**3)
            for order,ref in enumerate(references):
                for r,t,p in zip(radii,ref,matrix(radii,name,order)@c):
                    rows.append(dict(model=name,role=role,order=order,radius_A=r,reference=t,prediction=p,error=p-t))
    save_json(out/'profiles.json',profiles);write_csv(out/'radial_comparison.csv',rows)
    save_json(out/'summary.json',dict(completed=True,profiles=profiles,whole_energy_family_impossibility_proved=False,
        material_accepted=False,production_changed=False))
    print(profiles)


if __name__=='__main__':
    p=argparse.ArgumentParser(__doc__)
    p.add_argument('--density-fit',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();run(a.density_fit,a.out)
