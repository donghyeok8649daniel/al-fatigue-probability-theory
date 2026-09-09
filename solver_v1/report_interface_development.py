"""Summarize actual v15 optimizations and independent validation, not refit.

Numerical termination, exact-bulk constraints, held-out tests and material
acceptance are distinct. Incomplete/failed starts stay visible in the report.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from .report_public_calibration_evidence import read_csv,save_figure
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json


def summarize_replay(rows):
    """Keep interface and optional Bloch heldouts separate; no heldout loss."""
    groups={name:[] for name in ('fit','exact','interface_heldout','bloch_heldout')}
    for row in rows:
        role=row['role']
        group=('bloch_heldout' if row['observable'].startswith('bloch_') else 'interface_heldout') if role=='heldout' else role
        if group not in groups:raise ValueError('unknown observation role')
        target=float(row['target']);pred=float(row['prediction']);scale=float(row['scale'])
        if not np.isfinite([target,pred,scale]).all() or scale<=0:raise ValueError('finite values and positive scale required')
        residual=(pred-target)/scale
        if not np.isclose(residual,float(row['normalized_residual']),rtol=2e-10,atol=2e-10):
            raise ValueError('stored residual differs from stated units/scale')
        groups[group].append(residual)
    result={}
    for name,values in groups.items():
        values=np.asarray(values)
        result[name+'_count']=len(values)
        result[name+'_normalized_rms']=float(np.sqrt(np.mean(values**2))) if len(values) else None
        result[name+'_normalized_max']=float(np.max(abs(values))) if len(values) else None
    result['replayed_fit_loss']=float(np.dot(groups['fit'],groups['fit']))
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',type=Path)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    if args.out.exists():raise FileExistsError('preserve existing report')
    comparisons=[];optimizers=[];residuals=[];parameters=[];svd=[];stationary=[];curves={};bindings=[]
    for directory in sorted(args.directory.iterdir()):
        if not (directory/'calibration.json').exists():continue
        raw=(directory/'calibration.json').read_bytes();data=json.loads(raw)
        definition=json.loads((directory/'definition.json').read_bytes());best=data['best']
        row=dict(run=directory.name,runner_completed=bool(data['completed']),
            profiles=len(data['profiles']),optimizer_starts=len(data['optimizers']),
            successful_stopping_criteria=sum(bool(o['success']) for o in data['optimizers']),
            numerical_stops=sum(bool(o.get('numerical_stop')) for o in data['optimizers']),
            failed_profile_evaluations=len(data.get('failed_profiles',[])),
            stability_wavepoint_count=len(definition['wavepoints_cubic']),
            derivative_scheme=definition.get('numerical_derivative_scheme','2-point (historical runner)'),
            derivative_step=definition.get('numerical_relative_step',2e-4),
            exact_bulk=definition.get('exact_bulk_stage',False),
            fit_loss=best['squared_loss'],profile_heldout_normalized_rms=best['heldout_normalized_rms'],
            C11_GPa=best['cubic_GPa'][0],C12_GPa=best['cubic_GPa'][1],C44_GPa=best['cubic_GPa'][2],
            independent_validation=False,full_material_accepted=False,
            physical_time_calibrated=False)
        for i,opt in enumerate(data['optimizers']):
            optimizers.append(dict(run=directory.name,start_index=i,
                **{k:json.dumps(v) if isinstance(v,(list,dict)) else v for k,v in opt.items()}))
        if not definition.get('quartic_angular_extension'):
            raise ValueError('this v15 report requires an explicitly defined quartic-family basis')
        names=['u','v','A','B','C','D3','D1','D2','D_E','K3']
        if definition.get('density_angular_cross_extension'):names.append('E_xI')
        if len(names)!=len(best['coefficients']):raise ValueError('unknown coefficient basis')
        for name,value in zip(names,best['coefficients']):
            parameters.append(dict(run=directory.name,parameter=name,value=value,units='eV per normalized site basis'))
        for i,value in enumerate(best['decays']):
            parameters.append(dict(run=directory.name,parameter='shape_'+str(i),value=value,
                units='see definition.json: reduced exponential decay or declared extension shape'))
        validation=directory/'validation_exact_tangent'
        if (validation/'decision.json').exists():
            decision=json.loads((validation/'decision.json').read_bytes())
            if not decision['completed']:raise ValueError('incomplete independent validation')
            replay=read_csv(validation/'fit_replay.csv');row.update(summarize_replay(replay))
            row.update(independent_validation=True,
                finite_q_sample_positive=decision['finite_q_grid_positive_above_tail'],
                continuous_q_search_positive=decision['continuous_searches_positive_above_tail'])
            for r in replay:residuals.append(dict(run=directory.name,**r))
            ident=json.loads((validation/'exact_stage_identifiability.json').read_bytes())
            row.update(exact_tangent_rank=ident['rank'],exact_matrix_rank=ident['exact_matrix_rank'],
                exact_tangent_condition=ident['condition_number'],
                exact_derivative_residual=ident['exact_derivative_residual'])
            for i,s in enumerate(ident['singular_values']):
                svd.append(dict(run=directory.name,index=i,singular_value=s,
                    exact_rows=len(ident['exact_rows']),condition_number=ident['condition_number'],
                    active_inequality_cones_imposed=False,statistical_confidence_claimed=False))
            stationary.extend(dict(run=directory.name,**r) for r in read_csv(validation/'stationary_states.csv'))
            curves[directory.name]=read_csv(validation/'interface_curves.csv')
        comparisons.append(row)
        bindings.append(dict(run=directory.name,source_sha256=definition['source_sha256'],
            calibration_sha256=hashlib.sha256(raw).hexdigest()))
    if not comparisons:raise ValueError('no actual completed calibration runs found')
    # Some runs have no independent replay. A CSV still has the same columns
    # for every run, with missing diagnostics explicitly blank, never zero.
    def rectangular(rows):
        keys=list(dict.fromkeys(k for r in rows for k in r))
        return [{k:r.get(k) for k in keys} for r in rows]
    for name,rows in [('model_comparison',comparisons),('optimizer_terminations',optimizers),
                      ('parameter_sets',parameters),('independent_residuals',residuals),
                      ('exact_stage_svd',svd),('stationary_states',stationary)]:
        write_csv(args.out/(name+'.csv'),rectangular(rows))
    import matplotlib.pyplot as plt
    chosen=[name for name in ('joint_quartic','spectral_refined_fit','spectral_cross_fit','joint_bulk_interface') if name in curves]
    fig,axes=plt.subplots(1,3,figsize=(13,4),layout='constrained')
    for axis,path,title in zip(axes,('shockley','direct110','opening'),
                               ('Rigid Shockley path','Rigid direct <110> path','Rigid opening')):
        for j,name in enumerate(chosen):
            rows=[r for r in curves[name] if r['path']==path and float(r['a_over_h'])<=4]
            x=np.array([float(r['a_over_h'] if path=='opening' else r['fraction']) for r in rows])
            if j==0:axis.plot(x,[float(r['source_J_m2']) for r in rows],'k--',label='Al99 source (0 K)')
            axis.plot(x,[float(r['candidate_J_m2']) for r in rows],label=name)
        axis.set(title=title,xlabel='a / h' if path=='opening' else 'Path fraction',ylabel='Energy [J / m²]')
        axis.grid(alpha=.2)
    axes[0].legend(fontsize=6)
    fig.suptitle('Actual v15 fits: no candidate accepted as a full Al interface')
    save_figure(fig,args.out/'interface_tradeoffs.svg')
    save_json(args.out/'report_scope.json',dict(completed=True,bindings=bindings,
        optimization_rerun_by_report=False,global_minimum_proved=False,full_material_accepted=False,
        production_model_changed=False,kinetic_time_calibrated=False,
        interpretation='ftol/xtol termination is not material validation; failed starts and heldout residuals retained'))
    print('completed actual calibration/validation summary',flush=True)


if __name__=='__main__':main()
