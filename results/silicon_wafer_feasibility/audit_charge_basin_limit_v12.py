"""Hard first-entry basins and fast charge do not reduce as finite sinks.

An exact small synthetic counterexample in the existing deterministic law.
No silicon energy, charge rate, fracture criterion, or physical time is fitted.
"""
from __future__ import annotations
import argparse,csv,json,sys
from pathlib import Path
import numpy as np
from scipy.linalg import expm,solve
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from solver_v1.silicon_initiation_probability import initiation_collectors,first_passage_moments,evolve_first_passage,check_generator


def generator(rate,occupation,spatial):
    g=np.zeros((6,6))
    def edge(source,destination,value):
        g[destination,source]+=value;g[source,source]-=value
    for sector in (0,1):
        for q in (0,1):
            edge(3*sector+q,3*sector+q+1,spatial)
            edge(3*sector+q+1,3*sector+q,spatial)
    for q in range(3):
        edge(q,3+q,rate*occupation);edge(3+q,q,rate*(1-occupation))
    check_generator(g);return g


def main(args):
    if args.output.exists():raise ValueError('fresh output required')
    args.output.mkdir(parents=True)
    rows=[];boundary=[];matrix_errors=[];mass_errors=[];maximum_moment_error=0.
    times=np.array([0.,.2,2.,8.,20.]);spatial=.4
    spatial_generator=np.array([[-spatial,spatial,0],[spatial,-2*spatial,spatial],[0,spatial,-spatial]])
    common=initiation_collectors(spatial_generator,{'target_projection':[2]})
    limit=evolve_first_passage(common,[1.,0.],times)
    exact_limit_MFPT=3/spatial
    if abs(first_passage_moments(common)['mean_first_passage_time'][0]-exact_limit_MFPT)>1e-12:
        raise ValueError('independent reflecting-chain MFPT formula failed')
    for occupation in (.02,.2,.8):
        for rate in (.01,.1,1.,10.,100.,1000.,10000.,100000.):
            g=generator(rate,occupation,spatial)
            model=initiation_collectors(g,{'q2_sector1':[5]})
            # q0 is intact in both charge sectors, initially in charge equilibrium.
            initial=np.array([1-occupation,0.,0.,occupation,0.])
            moments=first_passage_moments(model)
            mfpt=float(initial@moments['mean_first_passage_time'])
            # Independent backward equations written in physical state order.
            # a_q: sector0 MFPT, b_q: sector1 MFPT; b_2=0.
            u=rate*occupation;v=rate*(1-occupation);r=spatial
            coefficients=np.array([[r+u,-r,0,-u,0],[-r,2*r+u,-r,0,-u],
                [0,-r,r+u,0,0],[-v,0,0,r+v,-r],[0,-v,0,-r,2*r+v]])
            direct=solve(coefficients,np.ones(5))
            maximum_moment_error=max(maximum_moment_error,float(np.max(abs(direct-moments['mean_first_passage_time']))))
            augmented=model['generator'].toarray();full0=np.r_[initial,0.]
            for index,time in enumerate(times):
                # Fixed six-dimensional Padé exponential handles stiffness without
                # a large Krylov iteration count; it is a separate evaluator.
                full=expm(time*augmented)@full0;survival=float(full[:-1].sum())
                mass_errors.append(abs(float(full.sum()-1)))
                rows.append(dict(occupation=occupation,charge_rate=rate,time_model=float(time),
                    survival=survival,projected_hard_basin_limit_survival=float(limit['survival'][index]),
                    difference_from_limit=survival-float(limit['survival'][index]),
                    MFPT_model=mfpt,hard_basin_limit_MFPT_model=exact_limit_MFPT,
                    boundary_charge_rate=rate*occupation,spatial_rate=spatial))
            if rate in (.01,1.,100.):
                independent=evolve_first_passage(model,initial,times)
                dense=np.array([expm(t*augmented)@full0 for t in times])
                matrix_errors.append(float(np.max(abs(dense[:,:-1]-independent['transient_mass']))))
        # At q2 with spatial motion disabled, surviving sector0 is killed at
        # nu*w_B. Renormalizing its charge weights would erase that probability.
        for rate in (.1,10.,1000.):
            for scaled_time in (.1,1.,10.):
                time=scaled_time/(rate*occupation)
                original=generator(rate,occupation,0.)
                model=initiation_collectors(original,{'q2_sector1':[5]})
                initial=np.array([0.,0.,1.,0.,0.])
                actual=evolve_first_passage(model,initial,[time])['survival'][0]
                expected=np.exp(-scaled_time)
                if abs(actual-expected)>5e-13:raise ValueError('exact boundary-layer exponential failed')
                boundary.append(dict(occupation=occupation,rate=rate,time_model=time,
                    survival=float(actual),exact_survival=float(expected),error=float(abs(actual-expected))))
    if max(mass_errors)>2e-9 or max(matrix_errors)>1e-11:raise ValueError('synthetic numerical balance unresolved')
    for filename,records in [('fast_charge_hard_basin.csv',rows),('boundary_layer_exact.csv',boundary)]:
        with (args.output/filename).open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
    result=dict(complete=True,synthetic=True,common_law_unchanged=True,states=6,charge_occupations=[.02,.2,.8],
        hard_target='q2 in sector1 only',initial='q0, local charge-equilibrium weights',
        cases=24,stored_time_points=len(rows),boundary_exact_checks=len(boundary),
        maximum_mass_error=max(mass_errors),maximum_Pade_vs_Krylov_error=max(matrix_errors),
        maximum_independent_backward_equation_error=maximum_moment_error,
        exact_boundary_survival='exp(-nu*w_B*t) for initial q2/sector0 and no spatial motion',
        hard_fast_limit='q2 is entirely absorbing when nu*w_B dominates spatial rates; both sectors must be removed at this projected coordinate',
        finite_rate_condition='nu*w_B, not nu alone; a small equilibrium target weight can delay the fast boundary limit',
        caveat='fixed-grid synthetic limit; does not identify actual Si initiating basins or commute automatically with spatial mesh refinement',
        new_DFT=0,new_MD=0,new_MACE=0,physical_clock=None,Si_initiation_calibrated=False)
    (args.output/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(10.5,4.5),layout='constrained')
    for occupation in (.02,.2,.8):
        selected=[row for row in rows if row['occupation']==occupation and row['time_model']==8.]
        axes[0].loglog([x['charge_rate'] for x in selected],[x['difference_from_limit'] for x in selected],'o-',label=f'w_B={occupation}')
        axes[1].loglog([x['boundary_charge_rate']/spatial for x in selected],
            [x['MFPT_model']-exact_limit_MFPT for x in selected],'o-',label=f'w_B={occupation}')
    axes[0].set(xlabel='Charge reset rate nu (model units)',ylabel='Survival minus hard-limit value at t=8')
    axes[1].set(xlabel='Boundary charge rate nu*w_B / spatial rate',ylabel='Mean-time excess above 7.5 (model units)')
    for ax in axes:ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.suptitle('Synthetic hard first-entry basin with fast charge\nNo silicon rate, fracture criterion, seconds or Hz')
    fig.savefig(args.output/'charge_basin_limit.png',dpi=180);plt.close(fig)
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args())
