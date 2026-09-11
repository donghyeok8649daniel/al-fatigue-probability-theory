"""Compare completed MD controls, dimensional mobility, and external linewidth.

All estimates remain scoped to the measured periodic-plane PMF and the DHO
extrapolation. No automatic PDE calibration or selected best cutoff is emitted.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from .mode_kinetic_calibration import gradient_drag_hypothesis
from .run_modal_mobility_projection import scalar_mixture_mobility
from .run_vector_registry_audit import save_json
from .run_low_stress_cyclic_diagnostic import write_csv


def summarize(analyses, output, *, sources_root=None,phonon_reference=None):
    output=Path(output)
    if output.exists():raise FileExistsError('fresh comparison output required')
    rows=[];modal=[];gradient=[];controls=[]
    for directory in map(Path,analyses):
        report=json.loads((directory/'summary.json').read_bytes())
        if not report['completed']:raise ValueError('completed analysis required')
        meta=report['source_metadata'];N=meta['repeats'];Np=meta['atoms_per_plane']
        tag=directory.name
        controls.append(dict(study=tag,ensemble=meta['ensemble'],dt_ps=meta['dt_ps'],
            planes=N,lattice_angstrom=meta['lattice_angstrom'],
            temperature_K=report['mean_temperature_K'],
            retained_ps=report['retained_duration_ps'],
            energy_slope_eV_per_atom_ps=report['energy_slope_eV_per_atom_ps'],
            energy_range_eV_per_atom=report['energy_range_eV_per_atom'],
            max_modal_cross_covariance=report['max_normalized_modal_offdiagonal']))
        if sources_root is not None:
            source_file=Path(sources_root)/tag/'plane_coordinates.npz'
            import hashlib
            if hashlib.sha256(source_file.read_bytes()).hexdigest()!=report['projection_sha256']:
                raise ValueError('thermodynamic summary must match analyzed trajectory')
            with np.load(source_file,allow_pickle=False) as data:
                thermo=data['thermo'][int(round(25/meta['frame_ps'])):]
            controls[-1].update(mean_pressure_bar=float(thermo[:,4].mean()),
                temperature_standard_deviation_K=float(thermo[:,0].std()))
        records=report['records']
        for cutoff in sorted({r['cutoff_ps'] for r in records}):
            for half in (0,1):
                for axis in ('normal','direct110','transverse','transverse_symmetry_pooled'):
                    subset=sorted([r for r in records if r['cutoff_ps']==cutoff and r['coordinate']==axis],key=lambda r:r['mode'])
                    if len(subset)!=N//2:raise ValueError('complete independent mode set required')
                    variance=np.array([r['training_variance_m2'] for r in subset])
                    g=np.array([r['damping_per_ps'] for r in subset])
                    w=np.array([r['frequency_rad_ps'] for r in subset])
                    temperature=subset[0]['temperature_used_K']
                    if half:
                        variance*=np.array([r['variance_ratio_second_to_first'] for r in subset])
                        g=np.array([r['second_half_damping_per_ps'] for r in subset])
                        w=np.array([r['second_half_frequency_rad_ps'] for r in subset])
                        temperature=subset[0]['second_half_temperature_K']
                    weight=np.full(len(subset),2.);weight[-1]=1. if N%2==0 else 2.
                    estimate=scalar_mixture_mobility(variance,2*g/(g*g+w*w)*1e-12,
                        multiplicities=weight,planes=N,temperature_K=temperature)
                    rows.append(dict(study=tag,cutoff_ps=cutoff,half=half,axis=axis,
                        temperature_K=temperature,**estimate,
                        cell_extensivity_hypothesis_mobility_m2_per_J_s=Np*estimate['mobility_m2_per_J_s'],
                        production_certified=False))
                    for index,r in enumerate(subset):
                        modal.append(dict(study=tag,cutoff_ps=cutoff,half=half,axis=axis,
                            mode=r['mode'],wavevector_xi=r['mode']/N,
                            damping_per_ps=g[index],envelope_time_ps=1/g[index],
                            linewidth_cycle_THz=g[index]/np.pi,
                            heldout_correlation_rmse=r['heldout_rmse'],
                            overdamped_control_rmse=r['heldout_overdamped_rmse']))
                    if len(subset)>3:
                        check=gradient_drag_hypothesis([r['mode'] for r in subset],g,
                            planes=N,fit_modes=[1,2,3])
                        gradient.append(dict(study=tag,cutoff_ps=cutoff,half=half,axis=axis,
                            fit_modes='1,2,3',heldout_modes='4,5,6',
                            coefficient_per_ps=check['coefficient_per_ps'],
                            heldout_rmse_per_ps=check['heldout_rmse_per_ps'],
                            minimum_per_mode_coefficient=float(np.min(check['per_mode_coefficient_per_ps'])),
                            maximum_per_mode_coefficient=float(np.max(check['per_mode_coefficient_per_ps']))))
    output.mkdir(parents=True)
    write_csv(output/'mobility_projection.csv',rows)
    write_csv(output/'modal_damping.csv',modal)
    write_csv(output/'controls.csv',controls)
    if gradient:write_csv(output/'gradient_drag_hypothesis.csv',gradient)
    if phonon_reference is not None:
        reference=json.loads(Path(phonon_reference).read_bytes())
        experimental=reference['records'][0]
        dft=reference['dft_comparison_records'][0]
        external=[]
        for control in controls:
            values=[r['linewidth_cycle_THz'] for r in modal if r['study']==control['study']
                    and r['axis']=='transverse_symmetry_pooled' and r['wavevector_xi']==.5]
            external.append(dict(study=control['study'],md_temperature_K=control['temperature_K'],
                md_linewidth_THz_min=min(values),md_linewidth_THz_max=max(values),
                experimental_temperature_K=experimental['temperature_K_from_figure'],
                experimental_linewidth_THz=experimental['linewidth_THz_from_figure'],
                dft_temperature_K=dft['temperature_K_from_figure'],
                dft_linewidth_THz=dft['linewidth_THz_from_figure'],
                md_to_experimental_min=min(values)/experimental['linewidth_THz_from_figure'],
                md_to_experimental_max=max(values)/experimental['linewidth_THz_from_figure'],
                source_doi=reference['doi'],experimental_confidence_interval_available=False))
        write_csv(output/'external_linewidth_validation.csv',external)
    save_json(output/'summary.json',dict(completed=True,controls=controls,
        mobility_records=rows,mode_records=modal,
        uncertainties_are_sensitivity_not_confidence_intervals=True,
        frequency_dependent_memory_and_cell_projection_not_certified=True,
        production_M_a_phys=None,production_M_s_phys=None,production_t0_seconds=None))
    print('Completed MD control and dimensional mobility comparison',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('analyses',nargs='+',type=Path);p.add_argument('--out',required=True,type=Path)
    p.add_argument('--sources-root',type=Path)
    p.add_argument('--phonon-reference',type=Path)
    a=p.parse_args();summarize(a.analyses,a.out,sources_root=a.sources_root,phonon_reference=a.phonon_reference)
