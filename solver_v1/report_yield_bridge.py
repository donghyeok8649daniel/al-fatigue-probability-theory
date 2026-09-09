"""Re-evaluate material states, summarize completed source calculations/plots.

Static stationary points and FD validation are executed here; calibration and
pinned-source branches are read from their separately completed runners.
"""
import csv
import io
import json
from pathlib import Path
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from .run_low_stress_cyclic_diagnostic import FIT, write_csv
from .run_vector_material_calibration import source_and_targets
from .run_vector_registry_audit import save_json
from .vector_material_calibration import build_coefficient_surface, VectorCoefficientBasis, observation_matrix, UNITS, IDEAL_H
from .vector_registry_audit import stationary_state
from .full_fcc_calibration_audit import cubic_constants_gpa

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/fcc111_active_interface/yield_bridge_v11'


def rows(path):
    with path.open(newline='') as f: return list(csv.DictReader(f))


def save_plot(fig,name):
    plt.rcParams['svg.fonttype']='none'; plt.rcParams['svg.hashsalt']='yield-bridge-v11'
    stream=io.StringIO(); fig.savefig(stream,format='svg',metadata={'Date':None},bbox_inches='tight')
    (OUT/(name+'.svg')).write_text('\n'.join(line.rstrip() for line in stream.getvalue().splitlines())+'\n',encoding='utf-8')
    cache=ROOT/'.cache/yield_bridge_v11'; cache.mkdir(parents=True,exist_ok=True)
    fig.savefig(cache/(name+'.png'),dpi=120,bbox_inches='tight'); plt.close(fig)


def main():
    started=time.perf_counter(); source,obs,states=source_and_targets()
    study=json.loads((OUT/'material_metric_refined/calibration.json').read_text())
    fits={'historical_candidate':json.loads(FIT.read_bytes())['angular_monotone_opening']}
    fits.update(study['best'])
    table=[]; stationary=[]; derivative=[]
    for name,fit in fits.items():
        c=np.asarray(fit['coefficients']); model=build_coefficient_surface(fit['scalar_decay'],fit['angular_decay'],c)
        matrix=observation_matrix(VectorCoefficientBasis(fit['scalar_decay'],fit['angular_decay']),obs)
        y=matrix@c; moduli=cubic_constants_gpa(y[:5])
        ref=model.evaluate((IDEAL_H,0.,0.))
        roots={}
        for label,guess,index in [('perfect',[IDEAL_H,0,0],0),('fault',states['fault'],0),('saddle',states['saddle'],1)]:
            result=stationary_state(model,guess,expected_index=index)
            roots[label]=result
            state=result['q']; value=result['evaluation']
            stationary.append(dict(model=name,state=label,valid=result['valid'],a=state[0],x=state[1],y=state[2],
                energy_J_m2=float(UNITS.energy_to_surface(value.energy)),
                gradient_max_eV_L0=float(max(abs(value.gradient))),
                eigenvalues=str(np.linalg.eigvalsh(value.hessian).tolist())))
        table.append(dict(model=name,**moduli,cohesion_eV_atom=y[1],
            zero_state_force_max=float(max(abs(ref.gradient))),minimum_interface_eigenvalue=np.linalg.eigvalsh(ref.hessian)[0],
            opening_40h_J_m2=float(UNITS.energy_to_surface(model.evaluate((40*IDEAL_H,0.,0.)).energy)),
            source_opening_40h_J_m2=float(UNITS.energy_to_surface(source.evaluate((40*IDEAL_H,0.,0.)).energy)),
            fault_energy_J_m2=float(UNITS.energy_to_surface(roots['fault']['evaluation'].energy)) if roots['fault']['valid'] else None,
            saddle_energy_J_m2=float(UNITS.energy_to_surface(roots['saddle']['evaluation'].energy)) if roots['saddle']['valid'] else None,
            heldout_normalized_rms=float(np.sqrt(np.mean([((pred-o.target)/o.scale)**2 for pred,o in zip(y,obs) if o.role=='heldout']))),
            adopted=False))
        q=np.array([1.09*IDEAL_H,.23,.11]); exact=model.evaluate(q)
        for step in (4e-5,2e-5):
            numerical_gradient=[]; numerical_hessian=[]
            for axis in range(3):
                d=np.eye(3)[axis]*step; plus=model.evaluate(q+d); minus=model.evaluate(q-d)
                numerical_gradient.append((plus.energy-minus.energy)/(2*step))
                numerical_hessian.append((plus.gradient-minus.gradient)/(2*step))
            derivative.append(dict(model=name,step_L0=step,
                gradient_error_eV_L0=float(np.max(abs(numerical_gradient-exact.gradient))),
                hessian_error_eV_L0sq=float(np.max(abs(np.array(numerical_hessian).T-exact.hessian)))))
    write_csv(OUT/'material_summary.csv',table)
    write_csv(OUT/'material_stationary_states.csv',stationary)
    write_csv(OUT/'material_derivative_check.csv',derivative)
    # Actual executed leading-log source results, not fitted yield curves.
    threshold=rows(OUT/'finite_source/hypothetical_source_thresholds.csv')
    shape=rows(OUT/'finite_source/source_shapes.csv')
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    for key,label in [('historical_candidate','Historical LJ/Bessel candidate'),
                      ('v11_cubic_5pct_not_adopted','New C-metric candidate (rejected)'),
                      ('source_elastic_comparator_0K','Source elasticity comparator')]:
        chosen=[r for r in threshold if r['model']==key and float(r['core_radius_over_b'])==1.]
        axes[0].loglog([float(r['span_um']) for r in chosen],[float(r['critical_outer_only_MPa']) for r in chosen],'o-',label=label)
    axes[0].set(xlabel='Hypothetical pin spacing L [micrometre]',ylabel='Outer-elastic bow-out threshold [MPa]',
                title='R=L, r_core=b; core/finite part omitted')
    axes[0].legend(fontsize=7); axes[0].grid(True,alpha=.25)
    for angle in (.35,.8,1.2,1.4):
        chosen=[r for r in shape if r['model']=='historical_candidate' and float(r['endpoint_angle'])==angle]
        x=np.array([float(r['x_over_span']) for r in chosen]); y=np.array([float(r['y_over_span']) for r in chosen])
        axes[1].plot(np.r_[-x[::-1],x],np.r_[y[::-1],y],label=f"{float(chosen[0]['applied_outer_only_MPa']):.2f} MPa")
    axes[1].set(xlabel='x/L',ylabel='y/L',title='Finite pinned-line equilibrium; L=1 micrometre')
    axes[1].set_aspect('equal'); axes[1].legend(fontsize=7); axes[1].grid(True,alpha=.25)
    fig.suptitle('Long-wave mechanism diagnostic — NOT measured Al yield',fontsize=11)
    fig.tight_layout(); save_plot(fig,'finite_source_diagnostic')
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    x=np.arange(3); width=.2
    for index,row in enumerate(table):
        axes[0].bar(x+index*width,[row[f'C{k}_GPa'] for k in (11,12,44)],width,label=row['model'])
    axes[0].plot(x+width,[114,62,32],'k_',markersize=16,label='0 K target')
    axes[0].set(xticks=x+width,xticklabels=['C11','C12','C44'],ylabel='GPa',title='Independent elasticity; fits not adopted')
    axes[0].legend(fontsize=7)
    data=rows(OUT/'experimental_benchmark/yield_reference_points.csv')
    axes[1].errorbar([float(r['inverse_diameter_b_over_D']) for r in data],[float(r['CRSS_over_G']) for r in data],
        xerr=[float(r['inverse_diameter_reading_halfwidth']) for r in data],
        yerr=[float(r['normalized_stress_reading_halfwidth']) for r in data],fmt='s',markersize=4,color='red')
    axes[1].set(xlabel='b/D (NOT b/source length)',ylabel='CRSS/G at 0.2% plastic shear',
        title='7 isolated measured points: Krebs et al.2017')
    axes[1].ticklabel_format(style='sci',axis='both',scilimits=(0,0))
    axes[1].grid(True,alpha=.25)
    fig.tight_layout(); save_plot(fig,'yield_target_and_material')
    save_json(OUT/'decision.json',dict(completed=True,material_accepted=False,
        experimental_yield_reproduced=False,physical_yield_MPa=None,
        derived_finite_source_outer_elastic_reference=True,finite_source_atomistic_core_validated=False,
        exact_bulk_does_not_certify_interface=True,physical_seconds=False,physical_Hz=False,
        source_lengths_are_hypotheses_not_measurements=True,
        empirical_yield_or_hardening_added=False,production_PDE_or_UI_changed=False,
        elapsed_seconds=time.perf_counter()-started,
        next='derive/calibrate independent core and finite-source geometry; joint material compatibility still fails'))


if __name__=='__main__': main()
