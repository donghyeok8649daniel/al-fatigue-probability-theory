"""Actual one-site atomistic tractions for the declared affine screw far field.

F=I+ex tensor (gamma_y ey+gamma_z ez), det F=1. P_iJ=dE_atom/dF_iJ/Omega.
For J=y,z the corresponding Cauchy traction columns equal these Piola columns.
Reference pressure and nonlinear normal reactions are reported, not clipped.
The input MPa is the nominal LINEAR-elastic far-field traction, not a hidden
exact nonlinear stress-control algorithm. No A_c or kinetic scale appears.
"""
import numpy as np


def affine_row_traction(core,gamma):
    gamma=np.asarray(gamma,float)
    if gamma.shape!=(2,) or np.any(~np.isfinite(gamma)):
        raise ValueError('finite two-component anti-plane engineering strain required')
    if getattr(core,'reference_only',False):
        transform=np.array([[core.rows.d,core.rows.d/3],[0.,core.rows.h]])
        if np.linalg.svd(transform,compute_uv=False).min()*(core.ring+1)<=core.source.r[-1]/core.length_angstrom:
            raise ValueError('source affine neighborhood misses published cutoff neighbors')
    R=core.reference; vectors=R.copy();vectors[:,0]+=R[:,1:]@gamma
    row=core.kernel.evaluate(vectors)
    z=(row['value']-core.reference_channels).sum(axis=0,keepdims=True)
    energy,weights,_=core._site_response(z,second=False)
    volume=core.rows.b*core.rows.d*core.rows.h
    traction=np.einsum('c,rci,rj->ij',weights[0],row['gradient'],R[:,1:])/volume
    return dict(energy_change_per_atom_eV=float(energy[0]),
                Piola_yz_columns_eV_L0cubed=traction,atomic_volume_over_L0cubed=volume,
                determinant=1.,gamma=gamma,transverse_ring=core.ring)


def nominal_shear_audit(core,shear_MPa,*,length_scale_m,energy_scale_J=1.602176634e-19):
    if not np.all(np.isfinite([shear_MPa,length_scale_m,energy_scale_J])) or min(length_scale_m,energy_scale_J)<=0:
        raise ValueError('finite nominal shear and positive dimensional scales required')
    scale=energy_scale_J/length_scale_m**3/1e6
    gamma=np.linalg.solve(core.far_field.matrix,[0.,shear_MPa/scale])
    actual=affine_row_traction(core,gamma);zero=affine_row_traction(core,[0.,0.])
    stress=actual['Piola_yz_columns_eV_L0cubed']*scale
    before=zero['Piola_yz_columns_eV_L0cubed']*scale
    return dict(nominal_sigma_xz_MPa=float(shear_MPa),gamma=gamma,
        actual_traction_columns_MPa=stress,zero_load_columns_MPa=before,
        traction_increment_columns_MPa=stress-before,
        sigma_xz_error_MPa=float(stress[0,1]-shear_MPa),
        actual_energy_change_eV_atom=actual['energy_change_per_atom_eV'],
        nonlinear_stress_control=False,physical_Hz=False)
