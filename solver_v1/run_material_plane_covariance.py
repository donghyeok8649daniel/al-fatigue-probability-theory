"""Condition-matched static material/MD covariance validation, not a refit."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from .aluminum_calibration import EV_J
from .isotropic_bulk_validation import IsotropicBulkBasis
from .periodic_plane_covariance import plane_mode_equipartition,empirical_plane_modes
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json
from .tail_constrained_material import declared_wavepoints
from .vector_material_calibration import LENGTH_M


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source',type=Path)
    parser.add_argument('--models',type=Path,nargs='+',required=True)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--additional-stretches',type=float,nargs='+',default=[],
        help='independent interpolation points; do not change the material density gauge')
    parser.add_argument('--wavepoint-step',type=float,default=.25,
        help='additional off-axis reciprocal stability grid; not a phonon frequency')
    args=parser.parse_args()
    if args.out.exists():raise FileExistsError('preserve existing validation')
    if np.any(~np.isfinite(args.additional_stretches)) or any(x<=0 for x in args.additional_stretches):
        raise ValueError('positive finite additional geometry stretches required')
    extra_points=declared_wavepoints(args.wavepoint_step)
    began=time.perf_counter()
    meta=json.loads((args.source/'source_metadata.json').read_bytes())
    source_hash=hashlib.sha256((args.source/'source_Al99.eam.alloy').read_bytes()).hexdigest()
    raw=(args.source/'plane_coordinates.npz').read_bytes()
    with np.load(args.source/'plane_coordinates.npz',allow_pickle=False) as data:q=data['coordinates_m']
    if q.shape!=(meta['frames'],meta['plane_classes'],3):raise ValueError('source projection shape mismatch')
    N=q.shape[1];Np=meta['atoms_per_plane'];T=meta['temperature_K']
    measured_modes,measured=empirical_plane_modes(q)
    block_values=[empirical_plane_modes(q[i*len(q)//4:(i+1)*len(q)//4]) for i in range(4)]
    blocks=[value[1] for value in block_values]
    empirical_error=np.max(abs(np.array(blocks)-measured),axis=0)
    mode_error=np.max(abs(np.array([value[0] for value in block_values])-measured_modes),axis=0)
    stretch=meta['lattice_angstrom']/(np.sqrt(2)*LENGTH_M/1e-10)
    transform=np.array([[0.,0.,1.],[-1.,0.,0.],[0.,-1.,0.]])
    results=[];wave_rows=[];derivatives=[];model_records=[];mode_rows=[];stability=[];bindings=[]
    for directory in args.models:
        data=json.loads((directory/'calibration.json').read_bytes())
        definition=json.loads((directory/'definition.json').read_bytes())
        if definition['source_sha256']!=source_hash:
            raise ValueError('static-target source potential differs from the actual MD potential')
        if not data['completed'] or not definition.get('quartic_angular_extension'):
            raise ValueError('completed quartic-family fit required')
        if definition.get('density_mixture_extension') or definition.get('cubic_density_extension'):
            raise ValueError('this validator does not implement these different site laws')
        name=directory.name;c=np.array(data['best']['coefficients']);decays=data['best']['decays'][:3]
        bindings.append(dict(model=name,calibration_sha256=hashlib.sha256((directory/'calibration.json').read_bytes()).hexdigest()))
        cross=definition.get('density_angular_cross_extension',False)
        saturation=data['best']['decays'][3] if definition.get('rational_angular_extension') else 0.
        for dilation in sorted(set([.99,1.,stretch,1.01]+args.additional_stretches)):
            previous=None
            for radius in (12.,16.):
                basis=IsotropicBulkBasis(decays,stretch=dilation,radius=radius,include_cross=cross)
                # Eleven actual MD plane modes, plus a fixed off-axis static
                # wedge test. No phonon frequency/mass/kinetic fit is performed.
                plane_points=[tuple(np.ones(3)*k/N) for k in range(1,N)]
                points=list(dict.fromkeys(plane_points+[tuple(v) for v in extra_points]))
                matrices=[];tails=[];worst=float('inf')
                for point in points:
                    columns,error=basis.evaluate(point)
                    H=np.einsum('c,cij->ij',c,columns);bound=error@abs(c)
                    margin=np.linalg.eigvalsh(H)[0]-bound;worst=min(worst,float(margin))
                    matrices.append(H);tails.append(bound)
                    wave_rows.append(dict(model=name,stretch=dilation,radius=radius,
                        qx=point[0],qy=point[1],qz=point[2],minimum_H=float(np.linalg.eigvalsh(H)[0]),
                        tail_bound=bound,robust_margin=margin,normalized_density=basis.x,
                        normalized_density_truncation=basis.density_series_tail))
                matrices=np.array(matrices[:N-1]);change=None
                stability.append(dict(model=name,stretch=dilation,radius=radius,
                    wavepoints=len(points),minimum_sampled_robust_margin=worst,
                    all_sampled_positive_above_tail=bool(worst>0),whole_zone_proof=False))
                if previous is not None:change=float(np.max(np.linalg.norm(matrices-previous,ord=2,axis=(-2,-1))))
                previous=matrices
                if dilation==stretch:
                    physical=transform@matrices@transform.T*EV_J/LENGTH_M**2
                    if np.any(np.linalg.eigvalsh(matrices)[:,0]<=np.array(tails[:N-1])):
                        model_records.append(dict(model=name,radius=radius,stable_plane_modes=False,
                                                  covariance_not_computed='unstable or unresolved mode'))
                        continue
                    modes,local=plane_mode_equipartition(physical,planes=N,atoms_per_plane=Np,temperature_K=T)
                    model_records.append(dict(model=name,radius=radius,stable_plane_modes=True,
                        local_covariance_m2=local,normalized_density=basis.x,
                        minimum_sampled_robust_margin=worst,radius_matrix_change=change,
                        fixed_reference_density=basis.reference_density,
                        physical_length_scale_m=LENGTH_M,physical_energy_scale_J=EV_J,
                        density_gauge_reset=False,atom_count_fitted=False))
                    for j,label in enumerate(('normal','slip_110','transverse')):
                        results.append(dict(model=name,radius=radius,coordinate=label,
                            fixed_model_static_RMS_m=np.sqrt(local[j,j]),actual_MD_RMS_m=np.sqrt(measured[j,j]),
                            variance_relative_difference=local[j,j]/measured[j,j]-1,
                            MD_empirical_block_variance_range_m2=empirical_error[j,j],
                            source_temperature_K=T,source_lattice_angstrom=meta['lattice_angstrom'],
                            kinetic_calibration=False))
                        for k,model_mode in enumerate(modes,start=1):
                            observed=float(measured_modes[k,j,j].real)
                            mode_rows.append(dict(model=name,radius=radius,coordinate=label,
                                mode_index=k,planes=N,q_cubic_fraction=k/N,
                                predicted_mode_variance_m2=model_mode[j,j],
                                actual_MD_mode_variance_m2=observed,
                                mode_variance_relative_difference=model_mode[j,j]/observed-1,
                                MD_empirical_block_mode_range_m2=mode_error[k,j,j],
                                conjugate_mode_index=N-k,conjugate_modes_independent=False,
                                source_temperature_K=T,kinetic_calibration=False))
                if dilation==stretch and radius==12.:
                    for j,v in enumerate(np.eye(3)):
                        H=matrices[3]  # k=4/N, an explicitly fixed nonzero mode
                        keyword=dict(planes=N,mode=4,polarization_plane=v,
                                     validation_radius=radius,saturation=saturation)
                        zero=basis.direct_sinusoidal_energy(c,amplitude=0.,**keyword)
                        fd=[]
                        for amplitude in (4e-4,2e-4):
                            plus=basis.direct_sinusoidal_energy(c,amplitude=amplitude,**keyword)
                            minus=basis.direct_sinusoidal_energy(c,amplitude=-amplitude,**keyword)
                            fd.append(2*(plus+minus-2*zero)/amplitude**2)
                        derivatives.append(dict(model=name,polarization_plane_axis=j,
                            analytic_H=float(v@H@v),central_4e4=fd[0],central_2e4=fd[1],
                            fourth_order_richardson=(4*fd[1]-fd[0])/3,
                            absolute_richardson_error=abs((4*fd[1]-fd[0])/3-v@H@v),
                            density_truncation=basis.density_series_tail))
        print('actual fixed-parameter MD-box validation',name,flush=True)
    write_csv(args.out/'material_MD_variance.csv',results)
    write_csv(args.out/'strained_static_modes.csv',wave_rows)
    write_csv(args.out/'direct_energy_derivative_check.csv',derivatives)
    write_csv(args.out/'material_MD_mode_variance.csv',mode_rows)
    write_csv(args.out/'strain_stability_summary.csv',stability)
    save_json(args.out/'scope.json',dict(completed=True,source_doi=meta['doi'],
        source_potential_sha256=source_hash,
        projection_sha256=hashlib.sha256(raw).hexdigest(),reference_length_scale_m=LENGTH_M,
        actual_box_stretch=stretch,temperature_K=T,atoms_per_plane=Np,plane_classes=N,
        model_records=model_records,parameter_refit=False,
        model_bindings=bindings,additional_stretches=args.additional_stretches,
        wavepoint_grid_step=args.wavepoint_step,
        anharmonic_thermal_renormalization=False,physical_time_calibrated=False,
        production_model_changed=False,whole_zone_stability_proof=False,
        elapsed_seconds=time.perf_counter()-began))
    print('static plane-covariance validation complete',flush=True)


if __name__=='__main__':main()
