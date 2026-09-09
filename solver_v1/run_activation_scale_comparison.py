"""Actually compute MPa mixed-traction stationary barrier sensitivities.

0K static per-cell barrier derivatives and 293K experimental apparent rate
slopes are different observables. This is a scale/normalization diagnosis,
NOT a joint calibration, guessed activation volume or specimen A_c estimate.
"""
import argparse
import hashlib
from pathlib import Path
import time

import numpy as np

from .activation_stress_sensitivity import barrier_and_derivative, apparent_rate_slope_from_barrier
from .run_vector_material_calibration import source_and_targets
from .validate_tail_calibration import load_material
from .vector_material_calibration import UNITS,LENGTH_M
from .run_vector_registry_audit import save_json
from .run_low_stress_cyclic_diagnostic import write_csv


ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--additional-model',type=Path,action='append',default=[],
        help='completed research calibration to compare; not production promotion')
    args=parser.parse_args()
    if args.out.exists():raise FileExistsError('preserve completed experiments')
    began=time.perf_counter();source,_,states=source_and_targets()
    old=ROOT/'results/fcc111_active_interface/tail_calibration_v14/quartic_exact_bulk'
    new=ROOT/'results/fcc111_active_interface/interface_development_v15/joint_quartic'
    models=[('Mishin_Al99_source_only',source),('v14_bulk_only_interface_unvalidated',load_material(old)[2]),
            ('v15_joint_interface_rejected',load_material(new)[2])]
    bindings=[]
    for path in args.additional_model:
        _,definition,model,_,_=load_material(path)
        if definition['source_sha256']!=source.reference.sha256:raise ValueError('comparison source mismatch')
        name=path.name+'_research_only'
        if any(label==name for label,_ in models):raise ValueError('duplicate comparison model name')
        models.append((name,model))
        bindings.append(dict(model=name,calibration_sha256=hashlib.sha256((path/'calibration.json').read_bytes()).hexdigest()))
    rows=[];checks=[]
    for name,model in models:
        for normal in (-25.,0.,25.):
            minimum=states['perfect'];saddle=states['saddle']
            for shear in (0.,.5,1.,2.,4.,6.,10.):
                tractions=np.array([normal,shear,0.])
                result=barrier_and_derivative(model,minimum,saddle,traction_MPa=tractions,units=UNITS)
                minimum=result['minimum']['q'];saddle=result['saddle']['q']
                volume=result['minus_barrier_traction_derivative_m3']
                rows.append(dict(model=name,normal_MPa=normal,shear1_MPa=shear,shear2_MPa=0.,
                    barrier_eV_per_cell=result['barrier_eV_cell'],
                    normal_derivative_m3=volume[0],shear1_derivative_m3=volume[1],shear2_derivative_m3=volume[2],
                    shear1_derivative_over_model_b3=volume[1]/LENGTH_M**3,
                    barrier_only_rate_slope_per_MPa_at_293K=apparent_rate_slope_from_barrier(volume[1],293.),
                    min_a=minimum[0],min_s1=minimum[1],min_s2=minimum[2],
                    saddle_a=saddle[0],saddle_s1=saddle[1],saddle_s2=saddle[2],
                    min_force_residual=result['minimum']['force_residual'],
                    saddle_force_residual=result['saddle']['force_residual'],
                    minimum_hessian_eigenvalue=result['minimum']['eigenvalues'][0],
                    saddle_negative_eigenvalue=result['saddle']['eigenvalues'][0],
                    static_reference_temperature_K=0.,physical_rate_available=False))
                if normal==0. and shear in (0.,4.,10.):
                    for step in (.02,.01):
                        energies=[]
                        for sign in (-1,1):
                            load=tractions.copy();load[1]+=sign*step
                            energies.append(barrier_and_derivative(model,minimum,saddle,
                                traction_MPa=load,units=UNITS)['barrier_eV_cell'])
                        fd=-(energies[1]-energies[0])/(2*step)
                        exact=result['minus_barrier_traction_derivative_eV_per_MPa'][1]
                        checks.append(dict(model=name,shear_MPa=shear,fd_step_MPa=step,
                            envelope_eV_per_MPa=exact,finite_difference_eV_per_MPa=fd,
                            absolute_error=abs(fd-exact),relative_error=abs(fd/exact-1)))
            print(name,'normal',normal,'completed static MPa sweep',flush=True)
    write_csv(args.out/'cell_barrier_stress_derivatives.csv',rows)
    write_csv(args.out/'envelope_validation.csv',checks)
    save_json(args.out/'scope.json',dict(completed=True,elapsed_seconds=time.perf_counter()-began,
        actual_static_cases=len(rows),source_sha256=source.reference.sha256,
        additional_model_bindings=bindings,
        static_temperature_K=0.,comparison_experimental_analysis_temperature_K=293.,
        geometry='uniform infinite half-crystal displacement; energy/atomic interface cell',
        shear='local resolved direct_110 shear1; shear2=0; three-state stationary solve',
        source_and_model_saddle='continued local index-one saddle, not global minimum-path proof',
        finite_source_or_loop_activation_solved=False,experimental_activation_fitted=False,
        apparent_activation_equal_to_specimen_correlation_area=False,
        production_pde_changed=False,physical_seconds_available=False))


if __name__=='__main__':main()
