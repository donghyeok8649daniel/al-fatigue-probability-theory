"""Independent raw-force, work and constraint replay of force-controlled states.

No potential calls. A solved generalized axial load is distinguished from the
finite residual of the full grip support force/torque. No material modulus,
strength, nucleation rate or physical clock is certified by this calculation.
"""
from __future__ import annotations
import argparse, csv, hashlib, json
from pathlib import Path
import numpy as np


def main(args):
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError('fresh output required')
    args.output.mkdir(parents=True, exist_ok=True)
    protocol = json.loads((args.run/'protocol.json').read_text(encoding='utf-8'))
    if hashlib.sha256(args.geometry.read_bytes()).hexdigest() != protocol['geometry_sha256']:
        raise ValueError('geometry hash differs from executed protocol')
    with np.load(args.geometry) as data:
        reference = data['reference']; lower = data['lower']; upper = data['upper']; free = data['free']
        area = float(data['area_A2']); gauge = float(data['gauge_length_A'])
    fixed = lower | upper; records = []; sources = []
    for path in sorted(args.run.glob('state_*/result.json')):
        row = json.loads(path.read_text(encoding='utf-8'))
        with np.load(path.parent/'raw.npz') as data:
            positions = data['positions']; forces = data['forces']; energy = float(data['energy'])
            extension = float(data['x'][-1]); saved_gradient = data['gradient']; enthalpy = float(data['enthalpy'])
            free_coordinate_error = float(np.max(abs(data['x'][:-1].reshape(-1,3)-positions[free])))
        exact = reference.copy(); exact[lower,2] -= extension/2; exact[upper,2] += extension/2
        grip_error = float(np.max(abs(exact[fixed]-positions[fixed])))
        target = row['target_nominal_GPa']*area/160.2176634
        conjugate = .5*(forces[lower,2].sum()-forces[upper,2].sum())
        gradient = np.r_[-forces[free].ravel(), conjugate-target]
        gradient_error = float(np.max(abs(gradient-saved_gradient)))
        reaction = -forces[fixed].sum(axis=0); free_resultant = forces[free].sum(axis=0)
        torque = -np.cross(positions[fixed],forces[fixed]).sum(axis=0)
        free_torque = np.cross(positions[free],forces[free]).sum(axis=0)
        calculated = dict(target_force_eV_A=target, conjugate_force_eV_A=float(conjugate),
            axial_force_residual_eV_A=float(conjugate-target), actual_nominal_GPa=float(conjugate/area*160.2176634),
            extension_A=extension, grip_centre_distance_A=gauge+extension, energy_eV=energy,
            enthalpy_eV=energy-target*extension, free_force_max_eV_A=float(np.linalg.norm(forces[free],axis=1).max()),
            reaction_force_residual_eV_A=float(np.linalg.norm(reaction)), reaction_torque_eV=float(np.linalg.norm(torque)),
            total_internal_force_eV_A=float(np.linalg.norm(forces.sum(axis=0))),
            total_internal_torque_eV=float(np.linalg.norm(np.cross(positions,forces).sum(axis=0))))
        observable_error = max(abs(calculated[k]-row[k]) for k in calculated)
        enthalpy_error = abs(enthalpy-calculated['enthalpy_eV'])
        if max(observable_error, enthalpy_error, gradient_error, grip_error, free_coordinate_error) > 1e-11:
            raise ValueError('raw state/gradient/constraint replay failed')
        records.append(dict(state=row['state'], target_nominal_GPa=row['target_nominal_GPa'],
            converged=row['converged'], **calculated,
            axial_target_residual_MPa=float((conjugate-target)/area*160217.6634),
            support_resultant_equivalent_MPa=float(np.linalg.norm(reaction)/area*160217.6634),
            support_axial_sum_equivalent_MPa=float(reaction[2]/area*160217.6634),
            lower_axial_reaction_eV_A=float(-forces[lower,2].sum()), upper_axial_reaction_eV_A=float(-forces[upper,2].sum()),
            grip_coordinate_error_A=grip_error, raw_gradient_error_eV_A=gradient_error,
            maximum_observable_replay_difference=observable_error,
            reaction_free_force_identity_error_eV_A=float(np.linalg.norm(reaction-free_resultant)),
            reaction_free_torque_identity_error_eV=float(np.linalg.norm(torque-free_torque))))
        sources.append(dict(file=path.parent.name+'/raw.npz', sha256=hashlib.sha256((path.parent/'raw.npz').read_bytes()).hexdigest()))
    if not records:
        raise ValueError('no completed raw state to audit')
    with (args.output/'force_control_replay.csv').open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(records[0])); writer.writeheader(); writer.writerows(records)
    comparison = None
    if len(records) == 2 and all(r['converged'] for r in records):
        a,b = records; delta = b['extension_A']-a['extension_A']; length = a['grip_centre_distance_A']
        if delta == 0:
            raise ValueError('zero response under two distinct requested loads')
        comparison = dict(extension_difference_A=delta, strain_relative_to_zero_axial_load=delta/length,
            secant_stress_strain_ratio_GPa=(b['actual_nominal_GPa']-a['actual_nominal_GPa'])/(delta/length),
            zero_load_contraction_relative_to_bulk_grip_length=a['extension_A']/gauge,
            energy_difference_eV=b['energy_eV']-a['energy_eV'],
            coarse_trapezoid_work_eV=.5*(a['conjugate_force_eV_A']+b['conjugate_force_eV_A'])*delta,
            scope='two-state secant of a finite bare prism with rigid grips; not bulk Young modulus or a measured wafer response')
    summary = dict(states_verified=len(records), records=records, comparison=comparison,
        source_raw_files=sources, geometry_sha256=protocol['geometry_sha256'],
        full_support_equilibrium_exact=False,
        equilibrium_status='axial generalized force target and all free-atom residuals reported; finite support resultants are not hidden',
        new_potential_calls=0, new_DFT=0, new_MD=0, initiation_probability=None, physical_clock=None)
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    if comparison:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,axes=plt.subplots(1,2,figsize=(10,4.5),layout='constrained')
        strain=[100*(r['extension_A']-records[0]['extension_A'])/records[0]['grip_centre_distance_A'] for r in records]
        axes[0].plot(strain,[1000*r['actual_nominal_GPa'] for r in records],'o-',color='#3f799c')
        axes[0].set_xlabel('Extension relative to zero axial load (%)'); axes[0].set_ylabel('Nominal axial stress (MPa)')
        axes[0].set_title('0 and 100 MPa prescribed load states')
        x=np.arange(len(records))
        axes[1].bar(x-.17,[abs(r['axial_target_residual_MPa']) for r in records],.34,label='Axial target error')
        axes[1].bar(x+.17,[r['support_resultant_equivalent_MPa'] for r in records],.34,label='Full support resultant / area')
        axes[1].set_xticks(x,['0 MPa','100 MPa']); axes[1].set_ylabel('Residual expressed as stress (MPa)'); axes[1].legend(fontsize=8)
        axes[1].set_title('Finite convergence errors remain visible')
        fig.suptitle('Bare 360-atom prism: force control replay, no crack-initiation claim')
        fig.savefig(args.output/'force_control_response.png',dpi=180); plt.close(fig)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--geometry',type=Path,required=True); parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args())
