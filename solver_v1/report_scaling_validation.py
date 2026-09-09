"""Plot completed static/covariance/source-scale checks, with scope labels."""
import argparse
from pathlib import Path

import numpy as np

from .report_public_calibration_evidence import read_csv,save_figure
from .run_vector_registry_audit import save_json
from .run_low_stress_cyclic_diagnostic import write_csv
from .aluminum_calibration import EV_J
from .vector_material_calibration import LENGTH_M,UNITS,IDEAL_H


def atomic_unit_ledger():
    """Declared full-FCC static units, NOT a yield or clock calibration.

    The crystallographic per-atom volume is exactly A_atomic_cell*h; it is
    only a bulk energy-density unit, not an inserted activation volume.
    """
    area=UNITS.atomic_cell_area_m2;h=IDEAL_H*LENGTH_M
    entries=[
        ('reference_b_L0',LENGTH_M,'m','declared FCC4.05Angstrom reference'),
        ('reference_h111',h,'m','crystallographic plane spacing'),
        ('atomic_interface_cell_area',area,'m^2','microscopic crystallography; NOT A_c'),
        ('crystallographic_per_atom_volume',area*h,'m^3','FCC atom volume; NOT activation volume'),
        ('energy_unit',EV_J,'J','one electronvolt'),
        ('one_eV_per_interface_cell',EV_J/area,'J/m^2','surface/GSF conversion'),
        ('one_eV_per_reduced_coordinate_traction',EV_J/(area*LENGTH_M)/1e9,'GPa','unit conversion; NOT ideal or actual strength'),
        ('one_eV_per_reduced_coordinate_force',EV_J/LENGTH_M,'N','conjugate coordinate force'),
        ('one_eV_per_reduced_coordinate_squared',EV_J/LENGTH_M**2,'J/m^2','physical coordinate Hessian'),
        ('one_eV_per_atom_energy_density',EV_J/(area*h)/1e9,'GPa','bulk energy-density unit; NOT yield'),
        ('one_eV_per_L0_line_energy',EV_J/LENGTH_M,'J/m','line energy; NOT finite activation energy'),
        ('mobility_time_scale_numerator',LENGTH_M**2/EV_J,'m^2/J','M_phys=M_star*this/t0; t0 unavailable'),
        ('one_MPa_interface_force',float(UNITS.traction_mpa_to_force(1.)),'eV/reduced_coordinate','explicit microscopic work conversion; not production kappa replacement'),
        ('source_300K_thermal_energy',1.380649e-23*300/EV_J,'eV','actual MD temperature; not default production kT'),
    ]
    return [dict(quantity=name,value=value,units=units,scope=scope,
                 actual_strength_calibrated=False,physical_PDE_time_calibrated=False)
            for name,value,units,scope in entries]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results',type=Path)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    if args.out.exists():raise FileExistsError('preserve existing scale report')
    args.out.mkdir(parents=True)
    write_csv(args.out/'atomic_unit_ledger.csv',atomic_unit_ledger())
    import matplotlib.pyplot as plt
    model=read_csv(args.results/'fixed_material_modes/material_MD_mode_variance.csv')
    source=read_csv(args.results/'plane_normalization/mode_variance.csv')
    fig,axes=plt.subplots(1,3,figsize=(12,4),layout='constrained')
    for ax,coordinate in zip(axes,('normal','slip_110','transverse')):
        data=[r for r in source if r['coordinate']==coordinate and int(r['mode'])<=6]
        k=np.array([int(r['mode']) for r in data])
        ax.errorbar(k,[float(r['actual_MD_variance_m2'])*1e24 for r in data],
            yerr=[float(r['empirical_four_block_range_m2'])*1e24 for r in data],
            fmt='ko',ms=3,capsize=2,label='MD, empirical block range (not CI)')
        ax.plot(k,[float(r['harmonic_source_variance_m2'])*1e24 for r in data],
                'k--',label='Source static Al99')
        for name in ('quartic_exact_bulk','joint_quartic','spatial_stiffness_fit','cross_verified_polish'):
            rows=[r for r in model if r['model']==name and float(r['radius'])==16
                  and r['coordinate']==coordinate and int(r['mode_index'])<=6]
            ax.plot([int(r['mode_index']) for r in rows],
                    [float(r['predicted_mode_variance_m2'])*1e24 for r in rows],label=name)
        ax.set(title=coordinate,xlabel='Periodic plane mode k (N=12)',ylabel='Variance [10^-24 m²]')
        ax.grid(alpha=.2)
    axes[1].legend(fontsize=6)
    fig.suptitle('Actual 300 K MD box: static normalization test, not a mobility fit')
    save_figure(fig,args.out/'mode_resolved_variance.svg')

    cases=read_csv(args.results/'finite_source_refined/hypothetical_finite_source_barriers.csv')
    fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    for span in (.2,1.,5.):
        rows=[r for r in cases if np.isclose(float(r['span_um']),span) and float(r['core_radius_over_b'])==2]
        x=[float(r['load_fraction']) for r in rows]
        axes[0].semilogy(x,[float(r['outer_barrier_eV']) for r in rows],'o-',label=f'L={span:g} micrometre')
        axes[1].semilogy(x,[float(r['derived_stress_derivative_over_b3']) for r in rows],'o-',label=f'L={span:g} micrometre')
    axes[0].set(ylabel='Outer-only barrier [eV]')
    axes[1].set(ylabel='Derived -d(barrier)/d(shear) [b³]')
    for ax in axes:
        ax.set(xlabel='Applied shear / outer critical shear')
        ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Hypothetical pinned-source geometry; NOT calibrated yield or an activated rate')
    save_figure(fig,args.out/'finite_source_barrier_scales.svg')
    save_json(args.out/'scope.json',dict(completed=True,
        actual_MD_observed=True,material_parameters_refitted=False,
        source_geometry_hypothetical=True,physical_time_of_PDE_calibrated=False,
        activated_rate_or_probability_computed=False,whole_Brillouin_zone_certified=False))
    print('completed static scale plots; physical interpretation flags retained')


if __name__=='__main__':main()
