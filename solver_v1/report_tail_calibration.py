"""Report completed v14 optimizations and run low-MPa static checks.

No optimization, kinetics, fatigue-life inference or production registration.
The selected candidate has a scoped bulk calibration, NOT full Al validation.
Independent curve/held-out failures remain visible in tables and plots.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import time
from io import StringIO

import numpy as np
from scipy.optimize import brentq

from .run_vector_material_calibration import source_and_targets
from .run_vector_registry_audit import save_json,stress_probes,point_row,derivative_rows
from .run_low_stress_cyclic_diagnostic import write_csv
from .validate_tail_calibration import load_material
from .vector_registry_audit import stationary_state
from .vector_material_calibration import UNITS,LENGTH_M
from .symmetry_resolved_material import quadrupole_channel_gauge


def read_csv(path):
    with path.open(encoding='utf-8',newline='') as stream:
        return list(csv.DictReader(stream))


def opening_extrema(model, samples):
    """Bracket ALL observed curvature sign changes; not a unimodal peak fit.

    Independent 129/257-point brackets check the extrema on [h,5h]. No claim
    is made about unobserved double roots or the entire unbounded domain.
    """
    def at(a): return model.evaluate([float(a),0.,0.])
    grid=np.linspace(model.h,5*model.h,samples)
    values=[at(a) for a in grid];rows=[]
    for a,b,va,vb in zip(grid[:-1],grid[1:],values[:-1],values[1:]):
        if va.hessian[0,0]*vb.hessian[0,0]<0:
            root=brentq(lambda z:at(z).hessian[0,0],a,b,xtol=2e-12,rtol=2e-12)
            v=at(root)
            rows.append(dict(bracket_samples=samples,a_over_h=root/model.h,
                energy_J_m2=float(UNITS.energy_to_surface(v.energy)),
                traction_MPa=float(UNITS.force_to_traction_mpa(v.gradient[0])),
                curvature_eV_L0sq=v.hessian[0,0]))
    return rows


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',type=Path,help='completed selected calibration')
    parser.add_argument('--validation-name',default='validation_refined')
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();out=args.out;directory=args.directory
    if out.exists(): raise FileExistsError('preserve existing report')
    start=time.perf_counter()
    data,definition,model,_,_=load_material(directory)
    if (not definition.get('quartic_angular_extension') or definition.get('rational_angular_extension')
            or definition.get('density_mixture_extension') or definition.get('quadrupole_saturation_extension')):
        raise ValueError('this scoped report requires the unsaturated quartic candidate, not another ablation')
    validation=directory/args.validation_name
    decision=json.loads((validation/'decision.json').read_bytes())
    if not decision['completed']: raise ValueError('independent validation required first')
    source,_,states=source_and_targets();best=data['best']
    if source.reference.sha256!=definition['source_sha256']: raise ValueError('source binding mismatch')
    summary=[]
    for sub in sorted(directory.parent.iterdir()):
        if not sub.is_dir() or not (sub/'definition.json').exists(): continue
        completed=(sub/'calibration.json').exists()
        value=json.loads((sub/('calibration.json' if completed else 'checkpoint.json')).read_bytes())
        config=json.loads((sub/'definition.json').read_bytes())
        row=value.get('best')
        if row is None: row=min(value['profiles'],key=lambda p:p['squared_loss'])
        summary.append(dict(run=sub.name,completed=completed,profiles=len(value['profiles']),
            successful_starts=sum(bool(p['success']) for p in value['optimizers']),
            numerical_stops=sum(bool(p.get('numerical_stop')) for p in value['optimizers']),
            exact_bulk=config.get('exact_bulk_stage',False),parameters=len(row['coefficients'])+len(row['decays'])-2,
            fit_loss=row['squared_loss'],heldout_normalized_rms=row['heldout_normalized_rms'],
            C11_GPa=row['cubic_GPa'][0],C12_GPa=row['cubic_GPa'][1],C44_GPa=row['cubic_GPa'][2],
            minimum_training_q_margin=row['minimum_robust_margin'],
            full_material_accepted=False,
            interpretation='completed optimization' if completed else 'incomplete numerical attempt; not a fitted result'))
    write_csv(out/'family_comparison.csv',summary)
    # Re-evaluate source and candidate roots rather than labeling source-point
    # energies as the candidate's own relaxed fault/saddle.
    points=[]
    for label,index in [('perfect',0),('fault',0),('saddle',1)]:
        for name,surface in [('source_Al99',source),('scoped_static_candidate',model)]:
            root=stationary_state(surface,states[label],expected_index=index)
            points.append(dict(valid=root['valid'],**point_row(name,label,root['q'],root['evaluation'],UNITS)))
    write_csv(out/'relaxed_state_comparison.csv',points)
    stresses=[]
    for name,surface in [('source_Al99',source),('scoped_static_candidate',model)]:
        stresses.extend(stress_probes(surface,name,UNITS))
    write_csv(out/'static_stress_scenarios.csv',stresses)
    deriv=derivative_rows(model,'scoped_static_candidate')
    write_csv(out/'derivative_refinement.csv',deriv)
    extrema=[]
    for name,surface in [('source_Al99',source),('scoped_static_candidate',model)]:
        for n in (129,257):
            extrema.extend(dict(model=name,**r) for r in opening_extrema(surface,n))
    write_csv(out/'opening_traction_extrema.csv',extrema)
    # Store physical unit conventions separately from dimensionless density gauges.
    c=np.asarray(best['coefficients']);decays=best['decays']
    length_angstrom=LENGTH_M/1e-10
    epsilon=c[1]**2/(4*c[0]);sigma=(c[0]/c[1])**(1/6)
    save_json(out/'scoped_parameter_set.json',dict(
        name='Al99_matched_bulk_quartic_v14',model_class=type(model).__name__,
        status='bulk-calibrated static research candidate; interface validation failed',
        production_registered=False,full_material_accepted=False,
        source_title='Mishin et al. (1999), Al99 EAM matched 0 K reference',
        source_doi='10.1103/PhysRevB.59.3393',source_sha256=source.reference.sha256,
        calibration_result_sha256=hashlib.sha256((directory/'calibration.json').read_bytes()).hexdigest(),
        target_temperature_K=0.,energy_unit='eV',length_scale_m=LENGTH_M,
        atomic_interface_area_m2=float(np.sqrt(3)/2*LENGTH_M**2),
        coefficient_names=['u','v','A','B','C','D3','D1','D2','D_E','K3'],
        coefficient_unit='eV; coordinates and per-site environment moments reduced by L0',
        coefficients=c,decays=decays,density_scale_convention='rho_ref=perfect bulk density for each fixed radial channel',
        quadrupole_gauge=quadrupole_channel_gauge(decays[1],decays[2]),
        epsilon_LJ_eV=float(epsilon),sigma_LJ_L0=float(sigma),sigma_LJ_angstrom=float(sigma*length_angstrom),
        lattice_constant_angstrom=length_angstrom*np.sqrt(2),
        bulk_cohesive_energy_eV_atom=decision['cohesive_eV_atom'],
        bulk_cubic_GPa=decision['cubic_GPa'],whole_Brillouin_zone_proved=False,
        physical_mobility_calibrated=False,physical_seconds=False,physical_Hz=False,
        experimental_yield_validated=False,
        permitted_use='scoped static comparisons only; not a production PDE or specimen-strength model'))
    curves=read_csv(validation/'interface_curves.csv')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    figure,axes=plt.subplots(1,3,figsize=(13,3.7))
    for ax,path,title in zip(axes,['shockley','direct110','opening'],
        ['Rigid Shockley path (partial translation)','Rigid direct <110> path','Rigid opening']):
        rows=[r for r in curves if r['path']==path and float(r['a_over_h'])<=4.]
        xkey='a_over_h' if path=='opening' else 'fraction'
        for label,key,style in [('Al99 reference','source_J_m2','k--'),('Scoped candidate','candidate_J_m2','C0-')]:
            ax.plot([float(r[xkey]) for r in rows],[float(r[key]) for r in rows],style,label=label)
        ax.set(xlabel='a / h' if path=='opening' else 'Path fraction',ylabel='Energy [J / m²]',title=title)
        ax.grid(alpha=.2)
    axes[0].legend(fontsize=8);figure.suptitle('Bulk calibrated; full interface NOT accepted')
    figure.tight_layout()
    svg=StringIO();figure.savefig(svg,format='svg')
    (out/'interface_comparison.svg').write_text(
        '\n'.join(line.rstrip() for line in svg.getvalue().splitlines())+'\n',encoding='utf-8')
    figure.savefig(out/'interface_comparison.png',dpi=140);plt.close(figure)
    selected=[r for r in stresses if r['model']=='scoped_static_candidate']
    save_json(out/'report_status.json',dict(completed=True,elapsed_seconds=time.perf_counter()-start,
        stress_states=len(stresses),maximum_static_force_residual_eV_L0=max(r['equilibrium_residual_eV_L0'] for r in selected),
        maximum_static_unload_registry_L0=max(r['du_norm_L0'] for r in selected if r['static_unload']),
        maximum_gradient_FD_error=max(r['gradient_max_abs_error_eV_L0'] for r in deriv),
        maximum_hessian_FD_error=max(r['hessian_max_abs_error_eV_L0sq'] for r in deriv),
        dynamic_hold_performed=False,PDE_changed=False,material_accepted=False,
        physical_time_calibrated=False))
    print('completed static calibration report and actual MPa scenarios',flush=True)


if __name__=='__main__': main()
