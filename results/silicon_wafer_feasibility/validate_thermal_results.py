"""Reevaluate saved SW states, constraints, hashes and timestep refinement.

Checks operate on raw outputs; passing them does not certify equilibrium,
physical mobility, an exhaustive transition network or wafer calibration.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from solver_v1.silicon_crack_research import RelaxedCoordinates,SpatialSW
from solver_v1.silicon_conditional_research import AMU_EV_PS2_A2
from solver_v1.silicon_thermal_research import bounded_force_control
from .run_local_crack_audit import save_json
from .run_static_probe import source_parameters


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def reference_path(recorded,root):
    path=Path(recorded)
    if path.exists():return path
    archived=root/'references'/path.name
    if not archived.exists():raise FileNotFoundError('recorded and archived references are missing')
    return archived


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path('results/silicon_thermal_v4'))
    args=parser.parse_args();root=args.root
    parameters,_=source_parameters();ensemble_rows=[];md_rows=[];direct_rows=[];hash_rows=[]
    for directory in ['hmc_pilot','mean_force_profile','front8_ensembles','multibasin_box15']:
        folder=root/directory
        metadata_file=folder/'summary.json'
        if not metadata_file.exists():metadata_file=folder/'interrupted_summary.json'
        if not metadata_file.exists():continue
        summary=read(metadata_file)
        for filename,expected in summary['code_sha256'].items():
            hash_rows.append(dict(calculation=directory,file=filename,recorded=expected,
                current=digest(filename),equal=expected==digest(filename)))
        references=reference_path(summary['source_references'],root)
        metadata=summary['reference_summary']
        model=SpatialSW(parameters,front_period=metadata['geometry']['front_period_A'])
        direct_done=set()
        for run in summary['records']:
            source=np.load(references/(run['reference']+'.npz'),allow_pickle=False)
            data=np.load(folder/(run['label']+'.npz'),allow_pickle=False)
            coords=RelaxedCoordinates(source['reference_positions'],source['fixed'],bond=source['bond'])
            r=data['final_positions_A'];result=model.evaluate(r)
            gradient,reaction=coords.pullback(result.gradient);x=coords.encode(r)
            control=bounded_force_control(x,gradient,source['response'],run['temperature_K'],run['bath_box_halfwidth_A'])
            expected=np.array([reaction[0],reaction[0]-control,np.sum(result.site_energy-source['sites'])])
            error=float(np.max(abs(expected-data['observations'][-1,:3])))
            fixed_error=float(np.max(abs(r[source['fixed']]-source['reference_positions'][source['fixed']]),initial=0))
            gap_error=float(abs(r[source['bond'][1],1]-r[source['bond'][0],1]-run['gap_A']))
            if error>2e-9 or fixed_error>1e-13 or gap_error>1e-13:
                raise AssertionError('saved ensemble observations/constraints do not reproduce')
            ensemble_rows.append(dict(calculation=directory,label=run['label'],observation_max_error=error,
                fixed_error_A=fixed_error,gap_error_A=gap_error,
                reference_hash_match=digest(references/(run['reference']+'.npz'))==run['reference_sha256']))
            group=(run['reference'],run['temperature_K'])
            if group not in direct_done:
                independent=model.evaluate(r,representation='direct')
                direct_rows.append(dict(calculation=directory,label=run['label'],
                    energy_error_eV=float(abs(np.sum(independent.site_energy-result.site_energy))),
                    force_max_error_eV_A=float(np.max(abs(independent.gradient-result.gradient)))))
                direct_done.add(group)
    for directory in ['thermal_dynamics','front8_dynamics','timestep_half','reconstruction_dynamics']:
        folder=root/directory
        if not (folder/'summary.json').exists():continue
        summary=read(folder/'summary.json')
        references=reference_path(summary.get('reference_source','.cache/si-thermal-v4/front4'),root)
        metadata=read(references/'summary.json')
        model=SpatialSW(parameters,front_period=metadata['geometry']['front_period_A'])
        for filename,expected in summary['code_sha256'].items():
            hash_rows.append(dict(calculation=directory,file=filename,recorded=expected,
                current=digest(filename),equal=expected==digest(filename)))
        for run in summary['records']:
            source=np.load(references/(run.get('reference','initial')+'.npz'),allow_pickle=False)
            data=np.load(folder/(run['label']+'.npz'),allow_pickle=False)
            coords=RelaxedCoordinates(source['reference_positions'],source['fixed'],bond=source['bond'])
            r=data['final_positions_A'];result=model.evaluate(r)
            reaction=coords.pullback(result.gradient)[1][0]
            initial=model.evaluate(coords.positions(data['initial_bath_coordinates_A'],[float(source['gap_A'])]))
            kinetic=.5*run['mass_amu']*AMU_EV_PS2_A2*np.sum(data['final_bath_velocities_A_ps']**2)
            errors=[abs(reaction-data['observations'][-1,0]),
                abs(np.sum(result.site_energy-initial.site_energy)-data['potential_difference_eV'][-1]),
                abs(kinetic-data['kinetic_energy_eV'][-1])]
            if max(errors)>2e-9:raise AssertionError('saved dynamics does not reproduce')
            np.testing.assert_allclose(coords.positions(data['final_bath_coordinates_A'],[float(source['gap_A'])]),r,rtol=0,atol=1e-13)
            md_rows.append(dict(calculation=directory,label=run['label'],force_energy_max_error=max(errors),
                gap_error_A=float(abs(r[source['bond'][1],1]-r[source['bond'][0],1]-float(source['gap_A']))),
                fixed_error_A=float(np.max(abs(r[source['fixed']]-source['reference_positions'][source['fixed']]),initial=0))))
            if directory=='reconstruction_dynamics' or run['reflection_count']>0:
                independent=model.evaluate(r,representation='direct')
                direct_rows.append(dict(calculation=directory,label=run['label'],
                    energy_error_eV=float(abs(np.sum(independent.site_energy-result.site_energy))),
                    force_max_error_eV_A=float(np.max(abs(independent.gradient-result.gradient)))))
    timestep=[]
    if (root/'timestep_half/summary.json').exists():
        for run in read(root/'timestep_half/summary.json')['records']:
            half=np.load(root/'timestep_half'/(run['label']+'.npz'),allow_pickle=False)
            coarse=np.load(root/'thermal_dynamics'/(run['preparation_label']+'_dt0.00100.npz'),allow_pickle=False)
            n=len(half['time_ps'])
            np.testing.assert_array_equal(half['initial_bath_coordinates_A'],coarse['initial_bath_coordinates_A'])
            np.testing.assert_array_equal(half['initial_bath_velocities_A_ps'],coarse['initial_bath_velocities_A_ps'])
            np.testing.assert_allclose(half['time_ps'],coarse['time_ps'][:n],atol=1e-14)
            coarse_error=float(np.max(abs(coarse['energy_residual_eV'][:n])))
            fine_error=float(np.max(abs(half['energy_residual_eV'])))
            short=half['time_ps']<=.1
            timestep.append(dict(label=run['preparation_label'],comparison_duration_ps=run['duration_ps'],
                same_preparation=True,coarse_max_energy_residual_eV=coarse_error,
                fine_max_energy_residual_eV=fine_error,observed_error_ratio=coarse_error/fine_error,
                short_time_force_max_difference_eV_A=float(np.max(abs(
                    coarse['observations'][:n,0][short]-half['observations'][short,0])))))
    if any(not r['reference_hash_match'] for r in ensemble_rows):raise AssertionError('reference checksum mismatch')
    if any(r['energy_error_eV']>2e-9 or r['force_max_error_eV_A']>2e-10 for r in direct_rows):
        raise AssertionError('independent triple evaluation disagrees')
    save_json(root/'raw_validation.json',dict(ensemble_states=ensemble_rows,md_states=md_rows,
        direct_triple_checks=direct_rows,source_hash_checks=hash_rows,timestep_refinement=timestep,
        source_hash_mismatches=[r for r in hash_rows if not r['equal']],
        incomplete_directories=[name for name in ['mean_force_profile','front8_dynamics','timestep_half','multibasin_box15']
                                if not (root/name/'summary.json').exists()],
        implementation_checks_passed=True,material_calibration_passed=False,physical_clock_calibrated=False,
        code_sha256=digest(__file__)))
    print('states',len(ensemble_rows),len(md_rows),'independent triples',len(direct_rows),
          'dt ratios',[r['observed_error_ratio'] for r in timestep],flush=True)


if __name__=='__main__':main()
