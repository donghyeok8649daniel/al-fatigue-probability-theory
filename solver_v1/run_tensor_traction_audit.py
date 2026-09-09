"""Actual static multiaxial scenarios; NOT a tensor-input production PDE.

Stress is in corrected +ABC cubic axes, with explicit interface normal and
two slip directions. No phenomenological chi, specimen area or physical time
is inferred. Quasistatic return to zero is not a dynamic residual-plasticity test.
"""
import argparse
import hashlib
from pathlib import Path
import time

import numpy as np

from .fcc111_geometry import fcc111_geometry
from .interface_static_scenarios import resolved_tensor_tractions
from .run_vector_material_calibration import source_and_targets
from .validate_tail_calibration import load_material
from .vector_material_calibration import UNITS,LENGTH_M
from .vector_registry_audit import stationary_state
from .run_vector_registry_audit import save_json
from .run_low_stress_cyclic_diagnostic import write_csv


def declared_tensor_scenarios():
    """Named physical MPa examples, not a selected specimen orientation fit."""
    plane=fcc111_geometry(np.sqrt(2)).plane_basis_in_stacked_cubic_axes()
    m,m2,n=plane
    cases={}
    for j,axis in enumerate('xyz'):
        tensor=np.zeros((3,3));tensor[j,j]=50.;cases['normal_'+axis+'x50MPa']=tensor
    for i,j,label in ((0,1,'xy'),(0,2,'xz'),(1,2,'yz')):
        tensor=np.zeros((3,3));tensor[i,j]=tensor[j,i]=10.;cases['shear_'+label+'10MPa']=tensor
    cases['biaxial_opposite25MPa']=np.diag([25.,-25.,0.])
    cases['mixed_MPa']=np.array([[40.,4.,-3.],[4.,-15.,2.],[-3.,2.,5.]])
    cases['pure_interface_normal50MPa']=50.*np.outer(n,n)
    cases['pure_interface_shear1_4MPa']=4.*(np.outer(n,m)+np.outer(m,n))
    cases['pure_interface_shear2_4MPa']=4.*(np.outer(n,m2)+np.outer(m2,n))
    return cases,n,m,m2


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--models',type=Path,nargs='+',required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    if args.out.exists():raise FileExistsError('preserve prior static scenarios')
    began=time.perf_counter();source,_,states=source_and_targets()
    models=[('source_Al99_only',source)];bindings=[]
    for path in args.models:
        _,definition,model,_,_=load_material(path)
        if definition['source_sha256']!=source.reference.sha256:raise ValueError('static source mismatch')
        if any(name==path.name for name,_ in models):raise ValueError('unique model names required')
        models.append((path.name,model));bindings.append(dict(model=path.name,
            calibration_sha256=hashlib.sha256((path/'calibration.json').read_bytes()).hexdigest()))
    cases,n,m,m2=declared_tensor_scenarios();rows=[];returns=[]
    for name,model in models:
        zero=stationary_state(model,states['perfect'],expected_index=0)
        if not zero['valid']:raise ArithmeticError('stable unloaded reference required')
        for case,tensor in cases.items():
            q=zero['q'];valid=True
            for fraction in (0.,.5,1.,.5,0.,-.5,-1.,-.5,0.):
                sigma=fraction*tensor;traction=resolved_tensor_tractions(sigma,n,m)
                force=UNITS.traction_mpa_to_force(traction)
                result=stationary_state(model,q,force=force,expected_index=0)
                valid=bool(valid and result['valid']);q=result['q']
                delta=q-zero['q']
                rows.append(dict(model=name,case=case,load_fraction=fraction,
                    sigma_xx_MPa=sigma[0,0],sigma_yy_MPa=sigma[1,1],sigma_zz_MPa=sigma[2,2],
                    sigma_xy_MPa=sigma[0,1],sigma_xz_MPa=sigma[0,2],sigma_yz_MPa=sigma[1,2],
                    local_normal_MPa=traction[0],local_shear1_MPa=traction[1],local_shear2_MPa=traction[2],
                    a_over_L0=q[0],s1_over_L0=q[1],s2_over_L0=q[2],
                    delta_a_angstrom=delta[0]*LENGTH_M/1e-10,
                    delta_s1_angstrom=delta[1]*LENGTH_M/1e-10,
                    delta_s2_angstrom=delta[2]*LENGTH_M/1e-10,
                    local_normal_opening_over_h=delta[0]/model.h,
                    conjugate_work_eV_cell=float(force@delta),
                    force_residual=result['force_residual'],minimum_H=result['eigenvalues'][0],
                    local_equilibrium_verified=result['valid'],dynamic_plasticity_inferred=False))
                if not valid:break  # do not jump to a new basin and call it continuation
            returns.append(dict(model=name,case=case,path_verified=valid,
                static_return_state_error_over_L0=float(np.linalg.norm(q-zero['q'])) if valid else None,
                dynamic_hold_performed=False,residual_plasticity_validated=False))
        print(name,'actual cubic tensor/load-unload states completed',flush=True)
    write_csv(args.out/'static_tensor_states.csv',rows)
    write_csv(args.out/'static_return_summary.csv',returns)
    save_json(args.out/'scope.json',dict(completed=True,actual_states=len(rows),
        source_sha256=source.reference.sha256,model_bindings=bindings,
        cubic_interface_normal=n,cubic_slip1=m,cubic_slip2=m2,
        physical_stress_unit='MPa',physical_displacement_unit='Angstrom',
        orientation='explicit +ABC cubic basis, not an inferred calibrated specimen orientation',
        shear_work='A_atomic_cell*L0*[Tn,tau1,tau2] dot deltaq; no fitted chi',
        temperature_scope='static0K energy surfaces',model_parameters_changed=False,
        production_tensor_input_connected=False,probability_PDE_run=False,
        dynamic_residual_or_yield_inferred=False,physical_PDE_Hz=False,
        elapsed_seconds=time.perf_counter()-began))


if __name__=='__main__':main()
