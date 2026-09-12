"""Optional-engine synthetic oscillator check; NO Al mass or mobility fit."""
import argparse
from pathlib import Path
import numpy as np
from .collective_forcing import harmonic_response
from .run_vector_registry_audit import save_json


def run(output):
    from lammps import lammps
    output=Path(output)
    if output.exists():raise FileExistsError('fresh result directory required')
    omega=2*np.pi*.05;force=.01;drag=.002
    exact=1/(1-omega**2+1j*omega*drag)
    rows=[]
    for dt in (.005,.0025,.00125):
        engine=lammps(cmdargs=['-log','none','-screen','none'])
        try:
            for command in ('units lj','atom_style atomic','atom_modify map array',
                'boundary f f f','region box block -5 5 -5 5 -5 5',
                'create_box 1 box',f'create_atoms 1 single {force*exact.real:.17g} 0 0',
                'mass 1 1','pair_style zero 1.0','pair_coeff * *',
                f'velocity all set {-omega*force*exact.imag:.17g} 0 0',
                f'variable drive atom -x+{force}*cos({omega:.17g}*time)',
                'fix forcing all addforce v_drive 0 0',f'fix damping all viscous {drag}',
                'fix dynamics all nve',f'timestep {dt}','thermo 1000000','run 0'):
                engine.command(command)
            positions=[];times=[]
            for i in range(16001):
                if i:engine.command(f'run {int(round(.025/dt))} pre no post no')
                positions.append(engine.gather_atoms('x',1,3)[0]);times.append(i*.025)
            t=np.asarray(times);mask=(t>=40)&(t<400)
            fit=harmonic_response(t[mask],np.asarray(positions)[mask],.05,force)
            measured=fit['real_A2_eV']+1j*fit['imag_A2_eV']
            rows.append(dict(dt=dt,measured_real=float(measured.real),measured_imag=float(measured.imag),
                complex_error=float(abs(measured-exact)),phase_error_rad=float(np.angle(measured/exact))))
        finally:engine.close()
    output.mkdir(parents=True)
    save_json(output/'summary.json',dict(completed=True,synthetic_dimensionless_oscillator=True,
        exact_real=float(exact.real),exact_imag=float(exact.imag),rows=rows,
        numerical_parameters_are_not_Al_calibration=True))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',required=True)
    run(p.parse_args().out)
