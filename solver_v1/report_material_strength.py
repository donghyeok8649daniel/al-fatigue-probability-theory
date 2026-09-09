"""Summarize executed v10 studies; no optimization or fabricated strength.

Run after both validation workflows have actually completed. Figure readers
must not interpret differently defined stresses as a matched prediction error.
"""
import csv
import hashlib
from io import StringIO
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from .run_low_stress_cyclic_diagnostic import FIT, write_csv
from .run_vector_material_calibration import source_and_targets
from .run_vector_registry_audit import save_json
from .strength_validation import StrengthConditions, compare_strength
from .vector_material_calibration import (
    LENGTH_M, UNITS, VectorCoefficientBasis, observation_matrix, fit_coefficients,
)
from .full_fcc_calibration_audit import cubic_constants_gpa


ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/fcc111_active_interface/material_strength_v10'


def read_csv(path):
    with path.open(newline='',encoding='utf-8') as stream:
        return list(csv.DictReader(stream))


def main():
    original=json.loads(FIT.read_bytes())['angular_monotone_opening']
    fits=[('old_candidate',original),
          ('joint_fit',json.loads((OUT/'joint_fit/calibration.json').read_text())['best']),
          ('opening_exchange',json.loads((OUT/'opening_exchange/calibration.json').read_text())['best'])]
    for name in ('joint_fit','opening_exchange'):
        if not json.loads((OUT/name/'validation_status.json').read_text())['completed']:
            raise RuntimeError('validation not completed')
    if not json.loads((OUT/'opening_exchange/independent_verification.json').read_text())['completed']:
        raise RuntimeError('independent calculations not completed')
    parameters=[]
    for name,fit in fits:
        c=fit['coefficients']; sigma=(c[0]/c[1])**(1/6); eps=c[1]**2/(4*c[0])
        parameters.append(dict(model=name,scalar_decay=fit['scalar_decay'],angular_decay=fit['angular_decay'],
            **dict(zip(('u','v','A','B','C','D3','D1','D2'),c)),
            epsilon_LJ_eV=eps,sigma_LJ_L0=sigma,sigma_LJ_angstrom=sigma*LENGTH_M/1e-10,
            isolated_atom_embedding_eV=-c[3]+c[4],production_promoted=False))
    write_csv(OUT/'parameter_sets.csv',parameters)
    residuals=read_csv(OUT/'joint_fit/fit_and_heldout_residuals.csv')
    summaries=[]
    for name,selection in [('old_candidate',[r for r in residuals if r['model']=='old_candidate']),
                           ('joint_fit',[r for r in residuals if r['model']=='best_feasible_same_family']),
                           ('opening_exchange',read_csv(OUT/'opening_exchange/fit_and_heldout_residuals.csv'))]:
        for role in ('fit','heldout'):
            v=np.array([float(r['normalized_residual']) for r in selection if r['role']==role])
            summaries.append(dict(model=name,role=role,count=len(v),sum_squared_normalized_error=float(v@v),
                                  normalized_RMS=float(np.sqrt(np.mean(v*v))),max_abs_normalized_error=float(max(abs(v)))))
    write_csv(OUT/'fit_quality_summary.csv',summaries)
    # Recompute source targets after the empty-neighbor identity fix, not replay.
    source,observations,_=source_and_targets()
    saved=json.loads((OUT/'joint_fit/target_definition.json').read_text())['observations']
    changes=[abs(o.target-r['target']) for o,r in zip(observations,saved)]
    if len(observations)!=len(saved) or any(o.name!=r['name'] for o,r in zip(observations,saved)):
        raise RuntimeError('source target identity changed')
    if max(changes)>1e-11:
        raise RuntimeError('source target values changed; refit required')
    # The SOURCE itself is not monotone in opening. Audit this rather than
    # presenting positivity as a universal law or silently editing the target.
    source_extrema=[]
    for count in (181,361):
        grid=np.linspace(1.001,2.8,count)*source.h
        curv=[source.evaluate((a,0.,0.)).hessian[0,0] for a in grid]
        for l,r,cl,cr in zip(grid[:-1],grid[1:],curv[:-1],curv[1:]):
            if cl*cr>=0:
                continue
            a=brentq(lambda a:source.evaluate((a,0.,0.)).hessian[0,0],l,r,xtol=2e-12)
            v=source.evaluate((a,0.,0.))
            independent=source.reference.interface_energy(a*source.length,0.)
            for step in (2.5e-6,1.25e-6):
                fd=(source.evaluate((a+step,0.,0.)).energy-source.evaluate((a-step,0.,0.)).energy)/(2*step)
                source_extrema.append(dict(samples=count,a_over_h=a/source.h,
                    normal_traction_MPa=float(UNITS.force_to_traction_mpa(v.gradient[0])),
                    energy_eV_cell=v.energy,independent_scalar_energy_error_eV=abs(v.energy-independent),
                    step_L0=step,FD_traction_error_MPa=float(UNITS.force_to_traction_mpa(fd-v.gradient[0]))))
    write_csv(OUT/'source_opening_extrema_audit.csv',source_extrema)
    # Quantify sensitivity to the monotonic shape prior at the SAME decays.
    fit=fits[-1][1]; basis=VectorCoefficientBasis(fit['scalar_decay'],fit['angular_decay'])
    matrix=observation_matrix(basis,observations)
    unconstrained=fit_coefficients(matrix,observations)
    save_json(OUT/'opening_shape_prior_ablation.json',dict(
        radial_decays=[fit['scalar_decay'],fit['angular_decay']],
        result_without_monotonicity=unconstrained,
        bulk_constants=cubic_constants_gpa(matrix[:5]@unconstrained['coefficients']),
        observations=[dict(name=o.name,role=o.role,reference=o.target,prediction=p,units=o.units,
                           normalized_residual=(p-o.target)/o.scale)
                      for o,p in zip(observations,matrix@unconstrained['coefficients'])],
        not_adopted=True,full_stability_not_validated=True,
        interpretation='monotonicity is a declared shape prior; SOURCE curve has a resolved negative-traction lobe'))
    ident=json.loads((OUT/'joint_fit/identifiability.json').read_text())
    names=['A','B','C','D3','D1','D2','log_scalar_decay','log_angular_decay']
    cos=np.array(ident['column_cosines'])
    correlations=[dict(left=names[i],right=names[j],cosine=float(cos[i,j]))
                  for i in range(8) for j in range(i+1,8) if abs(cos[i,j])>.90]
    write_csv(OUT/'strong_sensitivity_correlations.csv',correlations)
    experimental=read_csv(OUT/'experimental_benchmark/digitized_flow_points.csv')
    checks=[]
    fold=read_csv(OUT/'opening_exchange/fold_independent_refinement.csv')[-1]
    for row in experimental:
        if not row['stress_mpa']:
            continue
        reference=StrengthConditions('stress_at_plastic_strain','resolved_shear',float(row['plastic_shear_strain']),
            '99.99% as-cast Al single-crystal wire103um',None,'room temperature','300nm/s tensile; source figure2b')
        prediction=StrengthConditions('uniform_ideal_fold','resolved_shear',None,
            'rigid perfect FCC half-crystals',None,'0 K','static pure-x shear traction')
        checks.append(dict(experimental_flow_MPa=float(row['stress_mpa']),
            experimental_plastic_shear_strain=float(row['plastic_shear_strain']),
            uniform_ideal_MPa=float(fold['traction_x_MPa']),
            **compare_strength(reference,prediction,float(row['stress_mpa']),float(fold['traction_x_MPa']))))
    save_json(OUT/'strength_comparison_status.json',dict(comparisons=checks,
        experimental_strength_prediction_available=False,finite_source_validated=False,
        source_length_inferred_from_wire_diameter=False,physical_mobility_calibrated=False))
    save_json(OUT/'material_decision.json',dict(
        status='not accepted: bulk elastic and held-out interface discrepancies remain',
        entire_family_impossibility_proved=False,global_optimum_claimed=False,
        new_terms_added=False,production_promoted=False,source_targets_recomputed_max_change=max(changes),
        old_parameter_sha256=hashlib.sha256(FIT.read_bytes()).hexdigest(),
        source_kernel_sha256=hashlib.sha256((ROOT/'solver_v1/vector_interface_reference.py').read_bytes()).hexdigest(),
        source_minimum_opening_traction_MPa=min(r['normal_traction_MPa'] for r in source_extrema),
        opening_monotonicity_is_shape_prior_not_universal_physical_law=True,
        numerical_tests_are_not_material_acceptance=True,finite_source_strength_available=False,
        physical_seconds=False,physical_Hz=False))
    plots(source,experimental)
    print(json.dumps(dict(summary=summaries,target_change=max(changes),correlations=correlations)),flush=True)


def plots(source,experiment):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['svg.fonttype']='none'
    plt.rcParams['svg.hashsalt']='material_strength_v10'
    def save_svg(fig,path):
        buffer=StringIO()
        fig.savefig(buffer,format='svg',metadata={'Date':None})
        # Mechanical whitespace normalization of generated SVG, not image editing.
        path.write_text('\n'.join(line.rstrip() for line in buffer.getvalue().splitlines())+'\n',encoding='utf-8')
    bulk=read_csv(OUT/'opening_exchange/bulk_elastic_comparison.csv')
    barrier=read_csv(OUT/'opening_exchange/barrier_comparison.csv')
    opening=read_csv(OUT/'opening_exchange/opening_curve.csv')
    fig,ax=plt.subplots(1,3,figsize=(13,3.8),layout='constrained')
    for i,(label,vals) in enumerate([
        ('0 K source',[114,62,32]),('previous',[float(bulk[0][k]) for k in ('C11_GPa','C12_GPa','C44_GPa')]),
        ('v10 research',[float(bulk[1][k]) for k in ('C11_GPa','C12_GPa','C44_GPa')])]):
        ax[0].bar(np.arange(3)+(i-1)*.24,vals,.24,label=label)
    ax[0].set_xticks(range(3),['C11','C12','C44']); ax[0].set_ylabel('GPa'); ax[0].set_title('Bulk: not adequately matched')
    # Recompute source barrier values; saved states are initial guesses only.
    from .vector_registry_audit import stationary_state
    target=json.loads((OUT/'joint_fit/target_definition.json').read_text())['source_states']
    roots={name:stationary_state(source,target[name],expected_index=1 if name=='saddle' else 0)
           for name in ('saddle','fault')}
    if not all(r['valid'] for r in roots.values()):
        raise RuntimeError('source stationary state not validated for report')
    values={name:r['evaluation'].energy for name,r in roots.items()}
    source_vals=UNITS.energy_to_surface([values['saddle'],values['fault'],values['saddle']-values['fault']])
    for i,(label,vals) in enumerate([('source',source_vals),('previous',[float(barrier[0][k]) for k in ('forward_J_m2','fault_J_m2','reverse_J_m2')]),
                                    ('v10 research',[float(barrier[1][k]) for k in ('forward_J_m2','fault_J_m2','reverse_J_m2')])]):
        ax[1].bar(np.arange(3)+(i-1)*.24,vals,.24,label=label)
    ax[1].set_xticks(range(3),['forward','fault','reverse']); ax[1].set_ylabel('J/m$^2$'); ax[1].set_title('Connected stationary states')
    for name,label in [('old_candidate','previous'),('best_feasible_same_family','v10 research')]:
        rows=[r for r in opening if r['model']==name]
        ax[2].plot([float(r['a_over_h']) for r in rows],[float(r['normal_traction_MPa'])/1000 for r in rows],label=label)
    grid=np.linspace(1.,6.,100)
    vals=[float(UNITS.force_to_traction_mpa(source.evaluate((r*source.h,0.,0.)).gradient[0]))/1000 for r in grid]
    ax[2].plot(grid,vals,'--',label='0 K source'); ax[2].set_xlim(1,6)
    ax[2].set_xlabel('a/h'); ax[2].set_ylabel('normal traction [GPa]'); ax[2].set_title('Cleavage curve, not specimen strength')
    for a in ax: a.legend(fontsize=7); a.grid(alpha=.2)
    save_svg(fig,OUT/'material_comparison.svg')
    fig.savefig(ROOT/'.cache/material_strength_v10_comparison.png',dpi=130)
    plt.close(fig)
    known=[r for r in experiment if r['stress_mpa']]
    fig,ax=plt.subplots(figsize=(6,3.8),layout='constrained')
    ax.errorbar([float(r['plastic_shear_strain']) for r in known],[float(r['stress_mpa']) for r in known],
        xerr=[float(r['strain_digitization_halfwidth']) for r in known],
        yerr=[float(r['stress_digitization_halfwidth_mpa']) for r in known],fmt='o',capsize=3)
    ax.set(xlabel='plastic shear strain',ylabel='resolved shear stress [MPa]',
        title='Krebs2017 Fig2b: 103 micrometre Al wire\nDigitized flow points; NOT yield, no model prediction')
    ax.grid(alpha=.2); save_svg(fig,OUT/'experimental_flow_points.svg')
    fig.savefig(ROOT/'.cache/material_strength_v10_experiment.png',dpi=130)
    plt.close(fig)


if __name__=='__main__':
    main()
