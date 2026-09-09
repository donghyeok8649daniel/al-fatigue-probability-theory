"""Provenance-checked initialization of a larger STATIC atomistic core disk.

Only the initial displacement is transferred. All target forces, environments,
fixed boundary sites and final Morse checks are recomputed at target resolution.
No interpolation, strain scaling, force repair or fitted boundary is used.
"""
import numpy as np


def continue_core_field(core, prior, saved, *, parameter_sha256, length_scale_m,
                        shear_mpa, allow_expansion=False, require_stable=False,
                        prior_summary=None):
    """Preserve common logical rows; initialize newly free rows from far field.

    For inherited rows u_new^0=u_old+u_affine(tau_new-tau_old). Newly free rows
    start at u_Volterra+u_affine(tau_new). The outside field remains the target
    Dirichlet data. This is continuation, not proof of infinite-domain accuracy.
    """
    if not np.isfinite(length_scale_m) or length_scale_m<=0:
        raise ValueError('positive finite physical length scale required')
    if prior['parameter_sha256'] != parameter_sha256:
        raise ValueError('continuation cannot switch the material potential')
    if not np.isclose(prior['length_scale_m'],length_scale_m,rtol=1e-13,atol=0):
        raise ValueError('continuation cannot change the physical length scale')
    if not np.allclose(prior['center_over_L0'],core.far_field.center,rtol=0,atol=1e-13):
        raise ValueError('continuation must preserve the prescribed far-field center')
    radius=float(prior['radius_over_L0'])
    if not np.isfinite(radius) or radius<=0:
        raise ValueError('invalid source free radius')
    if radius != core.free_radius and not allow_expansion:
        raise ValueError('domain change requires explicit expansion option')
    if core.free_radius < radius:
        raise ValueError('this continuation only permits equal or larger free disks')
    if require_stable:
        probe=(prior_summary or {}).get('final_stability_probe',{})
        if not probe.get('stable_on_tested_fixed_boundary',False):
            raise ValueError('source is not a verified fixed-boundary stable core')
    indices=core.indices[core.free_ids]
    positions=core.xyz[core.free_ids]
    inner=np.linalg.norm(positions[:,1:]-core.far_field.center,axis=1)<radius
    expected={tuple(i) for i in indices[inner]}
    if set(saved)!=expected or len(saved)!=prior['free_sites']:
        raise ValueError('saved rows do not exactly cover the source free disk')
    if any(np.shape(v)!=(3,) or not np.all(np.isfinite(v)) for v in saved.values()):
        raise ValueError('finite unscaled three-component source displacements required')
    prior_shear=float(prior['shear_traction_MPa'])
    if not np.isfinite(shear_mpa) or not np.isfinite(prior_shear):
        raise ValueError('finite static shear tractions required')
    conversion=1e6*length_scale_m**3/1.602176634e-19
    affine=core.far_field.displacement(positions,
        shear_traction=(shear_mpa-prior_shear)*conversion,burgers_sign=0)
    seed=core.initial.copy()
    for k,index in enumerate(indices):
        key=tuple(index)
        if key in saved:
            seed[k]=np.asarray(saved[key])+affine[k]
    return seed,dict(source_radius_over_L0=radius,target_radius_over_L0=core.free_radius,
        inherited_free_rows=len(saved),new_free_rows=len(indices)-len(saved),
        domain_expansion=core.free_radius>radius,
        source_stability_required=bool(require_stable),
        initialization_only=True,interpolation_used=False,
        exterior_boundary_recomputed=True,domain_convergence_certified=False)
