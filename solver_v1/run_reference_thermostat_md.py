"""Independent Al99 reference MD for kinetic inference; never a PDE backend.

LAMMPS is optional and imported only by this explicitly invoked research tool.
The published potential is hash-bound. All ensembles start from the same saved
equilibrated state; NVE does not have a thermostat. Output is collective plane
coordinates, not Monte Carlo trajectory counts or a fatigue-life fit.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from .public_aluminum_kinetics import FCCPlaneProjection
from .run_vector_registry_audit import save_json


SOURCE_MD5 = 'eb0f0b204ea40787274efcf8e44d6de2'


def plane_velocity_statistics(velocity_angstrom_ps, plane_labels, planes, *, mass_amu):
    """Al99 reference velocities: per-atom KE, streaming KE removed separately.

    1 Angstrom/ps =100 m/s. This is source MD kinetic energy, not a mass-derived
    production overdamped time. Membership is the fixed FCC site projection.
    """
    v=np.asarray(velocity_angstrom_ps,float);labels=np.asarray(plane_labels)
    if (v.ndim!=2 or v.shape[1]!=3 or labels.shape!=(len(v),)
            or not np.all(np.isfinite(v)) or np.any(labels!=np.rint(labels))
            or np.any(labels<0) or np.any(labels>=planes)
            or not np.isfinite(mass_amu) or mass_amu<=0):
        raise ValueError('finite source velocities and valid plane labels required')
    labels=labels.astype(int);counts=np.bincount(labels,minlength=planes)
    if np.any(counts==0):raise ValueError('all represented planes must contain atoms')
    mean=np.stack([np.bincount(labels,weights=v[:,j],minlength=planes)/counts for j in range(3)],axis=1)
    factor=.5*mass_amu*1.66053906660e-27*100**2/1.602176634e-19
    raw=np.bincount(labels,weights=np.sum(v*v,axis=1),minlength=planes)/counts*factor
    internal=np.bincount(labels,weights=np.sum((v-mean[labels])**2,axis=1),minlength=planes)/counts*factor
    return mean,raw,internal


def sample_protocol(dt_ps, frame_ps, duration_ps):
    values = np.array([dt_ps, frame_ps, duration_ps], float)
    if np.any(~np.isfinite(values)) or np.min(values) <= 0:
        raise ValueError('positive finite MD timing required')
    stride = int(round(frame_ps / dt_ps))
    intervals = int(round(duration_ps / frame_ps))
    if (stride < 1 or intervals < 1
            or not np.isclose(stride * dt_ps, frame_ps, rtol=1e-12, atol=0)
            or not np.isclose(intervals * frame_ps, duration_ps, rtol=1e-12, atol=0)):
        raise ValueError('duration/frame and frame/dt must be integral')
    return stride, intervals


def run(potential, output, *, ensemble, dt_ps=.005, duration_ps=25.,
        frame_ps=.025, repeats=12, restart=None, damping_ps=1., threads=2,
        lattice_angstrom=4.065, thermal_observables=False, drive=None):
    from lammps import lammps

    stride, intervals = sample_protocol(dt_ps, frame_ps, duration_ps)
    potential, output = Path(potential).resolve(), Path(output)
    if hashlib.md5(potential.read_bytes()).hexdigest() != SOURCE_MD5:
        raise ValueError('source Al99 potential checksum mismatch')
    if output.exists():
        raise FileExistsError('fresh output directory required')
    if (ensemble not in ('nve', 'nvt') or repeats < 3 or damping_ps <= 0
            or not np.isfinite(lattice_angstrom) or lattice_angstrom<=0):
        raise ValueError('declared ensemble/geometry/damping required')
    output.mkdir(parents=True)
    started = time.perf_counter()
    args=['-log', str((output/'lammps.log').resolve()), '-screen', 'none']
    if threads>1:
        args+=['-pk','omp',str(threads),'-sf','omp']
    lmp = lammps(cmdargs=args)
    commands = []

    def command(line):
        commands.append(line)
        lmp.command(line)

    def gather_positions():
        return np.ctypeslib.as_array(lmp.gather_atoms('x', 1, 3),
                                    shape=(4*repeats**3*3,)).reshape(-1, 3).copy()

    try:
        command('units metal')
        command('atom_style atomic')
        command('atom_modify map array')
        if restart is None:
            command('boundary p p p')
            command(f'lattice fcc {lattice_angstrom:.17g}')
            command(f'region box block 0 {repeats} 0 {repeats} 0 {repeats} units lattice')
            command('create_box 1 box')
            command('create_atoms 1 box')
            command('mass 1 26.98')
        else:
            command(f'read_restart "{Path(restart).resolve().as_posix()}"')
        command('pair_style eam/alloy')
        command(f'pair_coeff * * "{potential.as_posix()}" Al')
        # EAM pair_coeff sets the element mass from its header, overriding
        # an earlier mass command. Query the engine, never assume 26.98.
        atomic_mass_amu=float(lmp.extract_atom('mass')[1])
        command('neighbor 2.0 bin')
        command('neigh_modify delay 0 every 1 check yes')
        command('thermo_style custom step time temp press pe ke etotal')
        command('thermo 1000')
        command(f'timestep {dt_ps:.17g}')
        ids = np.arange(1, 4*repeats**3+1)
        lo, hi, xy, yz, xz, periodic, flag = lmp.extract_box()
        if any(abs(v) > 1e-12 for v in (xy, yz, xz)) or not all(periodic):
            raise ValueError('only fixed periodic cubic source geometry supported')
        bounds = np.array([lo, hi]).T
        lengths = np.diff(bounds).ravel()
        projection = FCCPlaneProjection(ids, (gather_positions()-lo)/lengths,
            bounds, repeats=repeats, lattice_angstrom=lattice_angstrom)
        if restart is None:
            command('velocity all create 600 28459 rot yes dist gaussian mom yes')
            command('fix equil all nvt temp 300 300 1.0')
            command(f'run {int(round(25./dt_ps))}')
            command('unfix equil')
            command(f'write_restart "{(output/"equilibrated.restart").resolve().as_posix()}"')
        command('reset_timestep 0')
        if drive is not None:
            from .collective_forcing import conjugate_atom_weights
            axis=int(drive['axis']);amplitude=float(drive['force_eV_A'])
            frequency=float(drive['frequency_per_ps'])
            if axis not in (0,1,2) or not np.isfinite(amplitude+frequency) or amplitude==0 or frequency<=0:
                raise ValueError('explicit finite conjugate drive required')
            weights=conjugate_atom_weights(projection.plane)
            for plane,sign in ((0,-1),(1,1)):
                members=ids[projection.plane==plane]
                command(f'group drive{plane} id '+ ' '.join(map(str,members)))
                variables=[]
                for j in range(3):
                    name=f'drive{plane}_{j}'
                    value=sign*amplitude*projection.basis[axis,j]/len(members)
                    command(f'variable {name} equal {value:.17g}*cos(2*PI*{frequency:.17g}*time)')
                    variables.append('v_'+name)
                command(f'fix forcing{plane} drive{plane} addforce '+ ' '.join(variables))
        if ensemble == 'nve':
            command('fix dynamics all nve')
        else:
            command(f'fix dynamics all nvt temp 300 300 {damping_ps:.17g}')
        command('run 0')
        q = np.empty((intervals+1, repeats, 3))
        thermo = np.empty((intervals+1, 5))
        extra={}
        if drive is not None:
            extra.update(conjugate_force_eV_A=np.empty(intervals+1),
                external_power_eV_ps=np.empty(intervals+1),
                center_velocity_angstrom_ps=np.empty((intervals+1,3)))
        if thermal_observables:
            extra.update(plane_velocity_angstrom_ps=np.empty_like(q),
                plane_kinetic_eV_per_atom=np.empty((intervals+1,repeats)),
                plane_internal_kinetic_eV_per_atom=np.empty((intervals+1,repeats)))
        for index in range(intervals+1):
            if index:
                lmp.command(f'run {stride} pre no post no')
            q[index] = projection.evaluate(ids, (gather_positions()-lo)/lengths, bounds)
            thermo[index] = [lmp.get_thermo(k) for k in ('temp','pe','ke','etotal','press')]
            if thermal_observables or drive is not None:
                velocity=np.ctypeslib.as_array(lmp.gather_atoms('v',1,3),shape=(len(ids)*3,)).reshape(-1,3).copy()
            if drive is not None:
                force=amplitude*np.cos(2*np.pi*frequency*index*frame_ps)
                extra['conjugate_force_eV_A'][index]=force
                extra['external_power_eV_ps'][index]=force*(weights@velocity@projection.basis[axis])
                extra['center_velocity_angstrom_ps'][index]=velocity.mean(axis=0)
            if thermal_observables:
                mean,raw,internal=plane_velocity_statistics(velocity,projection.plane,repeats,mass_amu=atomic_mass_amu)
                extra['plane_velocity_angstrom_ps'][index]=mean@projection.basis.T
                extra['plane_kinetic_eV_per_atom'][index]=raw
                extra['plane_internal_kinetic_eV_per_atom'][index]=internal
            if index % 1000 == 0 or index == intervals:
                np.savez_compressed(output/'checkpoint.npz', coordinates_m=q[:index+1],
                    time_seconds=np.arange(index+1)*frame_ps*1e-12,
                    thermo=thermo[:index+1],**{key:value[:index+1] for key,value in extra.items()})
                print(json.dumps(dict(frame=index, frames=intervals+1,
                    elapsed_seconds=time.perf_counter()-started,
                    temperature_K=float(thermo[index,0]))), flush=True)
        np.savez_compressed(output/'plane_coordinates.npz', coordinates_m=q,
            time_seconds=np.arange(intervals+1)*frame_ps*1e-12, thermo=thermo,**extra)
        save_json(output/'summary.json', dict(completed=True, ensemble=ensemble,
            dt_ps=dt_ps, frame_ps=frame_ps, duration_ps=duration_ps,
            repeats=repeats, atom_count=len(ids), atoms_per_plane=4*repeats**2,
            lattice_angstrom=lattice_angstrom, potential_md5=SOURCE_MD5,
            source='Independent reference MD using Zenodo 10014454 Al99 potential',
            published_trajectory_reproduction=False,
            restart_sha256=None if restart is None else hashlib.sha256(Path(restart).read_bytes()).hexdigest(),
            thermostat_damping_ps=damping_ps if ensemble=='nvt' else None,
            mean_temperature_K=float(thermo[:,0].mean()),
            total_energy_range_eV=float(np.ptp(thermo[:,3])),
            total_energy_change_eV=float(thermo[-1,3]-thermo[0,3]),
            thermo_columns=['temperature_K','potential_eV','kinetic_eV','total_eV','pressure_bar'],
            lammps_version=lmp.version(), elapsed_seconds=time.perf_counter()-started,
            initialization_commands=commands,
            thermal_observables=thermal_observables,
            conjugate_drive=drive,
            drive_energy_convention='thermo energy excludes external potential; compare internal energy change with integral F*qdot' if drive else None,
            atomic_mass_amu=atomic_mass_amu,
            plane_kinetic_postprocess_mass_amu=atomic_mass_amu if thermal_observables else None,
            production_clock_calibrated=False))
    finally:
        lmp.close()


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--potential', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    p.add_argument('--ensemble', choices=['nve','nvt'], required=True)
    p.add_argument('--dt-ps', type=float, default=.005)
    p.add_argument('--duration-ps', type=float, default=25.)
    p.add_argument('--repeats', type=int, default=12)
    p.add_argument('--restart', type=Path)
    p.add_argument('--damping-ps', type=float, default=1.)
    p.add_argument('--threads', type=int, default=2)
    p.add_argument('--lattice-angstrom', type=float, default=4.065)
    p.add_argument('--thermal-observables', action='store_true')
    p.add_argument('--drive-axis',type=int,choices=(0,1,2))
    p.add_argument('--drive-force-eV-A',type=float)
    p.add_argument('--drive-frequency-per-ps',type=float)
    a=p.parse_args()
    values=(a.drive_axis,a.drive_force_eV_A,a.drive_frequency_per_ps)
    if any(v is not None for v in values) and not all(v is not None for v in values):
        p.error('all three drive arguments must be supplied together')
    drive=None if a.drive_axis is None else dict(axis=a.drive_axis,
        force_eV_A=a.drive_force_eV_A,frequency_per_ps=a.drive_frequency_per_ps)
    run(a.potential,a.out,ensemble=a.ensemble,dt_ps=a.dt_ps,
        duration_ps=a.duration_ps,repeats=a.repeats,restart=a.restart,
        damping_ps=a.damping_ps,threads=a.threads,lattice_angstrom=a.lattice_angstrom,
        thermal_observables=a.thermal_observables,drive=drive)
