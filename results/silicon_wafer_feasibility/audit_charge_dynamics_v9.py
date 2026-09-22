"""Deterministic fast-charge, spatial-grid, absorption and quenched controls."""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import sys
import numpy as np
from scipy.linalg import expm
from scipy.special import logsumexp

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from solver_v1.silicon_charge_dynamics import coupled_charge_generator, reversible_propagator, sg_chain


def setup(cells, rate):
    dx = 2./cells
    q = -1+dx*(np.arange(cells)+.5)
    free = np.array([.5*(q+.35)**2, .8*(q-.25)**2+.03])
    r = coupled_charge_generator(free, excess_counts=[0,1], chemical_potential=.02,
        thermal_energy=.12, spacing=dx, mobilities=[1.,.3], charge_rate=rate)
    p = np.exp(-((q+.15)/.2)**2); p /= p.sum()
    return q, p, r


def save_csv(path, rows):
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def audit(out):
    out.mkdir(parents=True, exist_ok=True)
    rates, grids, absorption, quenched, independent = [], [], [], [], []
    mass_errors, minimum_probabilities, stationarity, balance = [], [], [], []
    for cells in [24,48,96]:
        for rate in [0.,3.,30.,300.,3000.,30000.]:
            q, p, r = setup(cells, rate)
            l, pi = r['generator'], r['equilibrium']
            conservation = float(np.max(abs(np.asarray(l.sum(axis=0)))))
            stationary = float(np.max(abs(l@pi)))
            db = l.toarray()*pi
            db = float(np.max(abs(db-db.T)))
            stationarity.append(stationary); balance.append(conservation)
            for initial_kind in ['local_charge_equilibrium','all_sector_0']:
                pjoint = r['lift']@p if initial_kind == 'local_charge_equilibrium' else np.r_[p,np.zeros(cells)]
                for t in [.03,.3,1.]:
                    actual_joint = reversible_propagator(l,pi,t)@pjoint
                    actual = r['projection']@actual_joint
                    limit = reversible_propagator(r['coarse'],r['coarse_equilibrium'],t)@p
                    mass_error = float(abs(actual_joint.sum()-1))
                    mass_errors.append(mass_error); minimum_probabilities.append(float(actual_joint.min()))
                    rates.append(dict(cells=cells, charge_rate_model=rate, initial=initial_kind,
                        time_model=t, L1_marginal_error=float(sum(abs(actual-limit))),
                        L1_charge_disequilibrium=float(sum(abs(actual_joint-r['lift']@actual))),
                        probability_mass_error=mass_error, minimum_probability=float(actual_joint.min()),
                        generator_column_sum_max_abs=conservation, stationary_residual=stationary,
                        detailed_balance_residual=db))
            if cells == 24 and rate in [0.,30.,30000.]:
                spectral = reversible_propagator(l,pi,.3)
                dense = expm(.3*l.toarray())
                independent.append(dict(cells=cells, rate=rate,
                    spectral_vs_pade_max_abs=float(np.max(abs(spectral-dense)))))
    save_csv(out/'charge_rate_convergence.csv', rates)
    for cells in [12,24,48,96,192]:
        q,p,r = setup(cells,30.)
        for t in [.03,.3,1.]:
            exact = reversible_propagator(r['coarse'],r['coarse_equilibrium'],t)@p
            continuous = reversible_propagator(r['continuum'],r['coarse_equilibrium'],t)@p
            grids.append(dict(cells=cells, spacing=2/cells, time_model=t,
                L1_finite_grid_vs_continuum_SG=float(sum(abs(exact-continuous))),
                coarse_mass_error=float(abs(exact.sum()-1)),
                continuum_mass_error=float(abs(continuous.sum()-1))))
    save_csv(out/'spatial_grid_convergence.csv',grids)
    # Synthetic sector-dependent sink. This validates the algebra of killing,
    # not a Si crack dividing surface, barrier, or measured fracture rate.
    for rate in [0.,3.,30.,300.,3000.,30000.]:
        q,p,r = setup(24,rate)
        sink = np.vstack([.5*(q>.25),2.*(q>.25)])
        killing = r['generator'].toarray()-np.diag(sink.ravel())
        average_sink = np.sum(sink*r['weights'],axis=0)
        coarse_killing = r['coarse'].toarray()-np.diag(average_sink)
        # Augmented absorbing counter integrates the actual sink flux exactly.
        augmented = np.zeros((49,49)); augmented[:48,:48] = killing
        augmented[-1,:48] = sink.ravel()
        reduced = np.zeros((25,25)); reduced[:24,:24] = coarse_killing
        reduced[-1,:24] = average_sink
        for t in [.03,.3,1.]:
            actual = expm(t*augmented)@np.r_[r['lift']@p,0.]
            limit = expm(t*reduced)@np.r_[p,0.]
            marginal = r['projection']@actual[:-1]
            absorption.append(dict(rate_model=rate,time_model=t,
                full_absorbed=float(actual[-1]),coarse_absorbed=float(limit[-1]),
                absorption_error=abs(float(actual[-1]-limit[-1])),
                L1_survivor_error=float(sum(abs(marginal-limit[:-1]))),
                mass_plus_absorbed_error=abs(float(sum(actual)-1)),
                minimum_augmented_probability=float(min(actual))))
    save_csv(out/'absorption_convergence.csv',absorption)
    cells,kt = 48,.12
    q = -1+2/cells*(np.arange(cells)+.5)
    f = [.5*(q+.4)**2,.5*(q-.4)**2]
    ls = [sg_chain(a,thermal_energy=kt,spacing=2/cells,face_mobility=1).toarray() for a in f]
    pi = [np.exp(-a/kt-logsumexp(-a/kt)) for a in f]
    p = np.exp(-(q/.2)**2); p /= p.sum()
    def mix(t):
        return .5*reversible_propagator(ls[0],pi[0],t)+.5*reversible_propagator(ls[1],pi[1],t)
    for t in [.03,.1,.3,1.]:
        one = mix(2*t)@p
        reset = mix(t)@mix(t)@p
        average = expm(t*(ls[0]+ls[1]))@p
        quenched.append(dict(time_per_half=t, L1_semigroup_failure=float(sum(abs(one-reset))),
            L1_average_generator_failure=float(sum(abs(one-average))),
            actual_mass_error=abs(float(one.sum()-1))))
    save_csv(out/'quenched_mixture_counterexample.csv',quenched)
    main_rates = [r for r in rates if r['cells']==48 and r['initial']=='local_charge_equilibrium' and r['time_model']==.3]
    main_grids = [r for r in grids if r['time_model']==.3]
    report = dict(status='synthetic deterministic generator validation; actual Si electronic rate unavailable',
        units='arbitrary model coordinates, energies and time; not physical seconds or Hz',
        energy_definition=['F0=0.5(q+0.35)^2','F1=0.8(q-0.25)^2+0.03'],
        kT=.12, mu=.02, sector_counts=[0,1],mobilities=[1.,.3],
        maximum_mass_error=max(mass_errors), minimum_probability=min(minimum_probabilities),
        max_stationary_residual=max(stationarity), max_column_sum_residual=max(balance),
        independent_exponential_comparisons=independent,
        main_rate_slice=main_rates,main_grid_slice=main_grids,
        charge_rate_cases=len(rates), spatial_cases=len(grids), absorption_cases=len(absorption),
        absorption_mass_balance_max=max(r['mass_plus_absorbed_error'] for r in absorption),
        new_DFT_runs=0,new_MD_runs=0)
    (out/'dynamics_summary.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes = plt.subplots(2,2,figsize=(10,7),layout='constrained')
    rows=[r for r in main_rates if r['charge_rate_model']>0]
    axes[0,0].loglog([r['charge_rate_model'] for r in rows],[r['L1_marginal_error'] for r in rows],'o-')
    axes[0,0].set(xlabel='Charge rate (model units)',ylabel='Marginal L1 error',title='Fast charge: fixed 48-cell grid')
    axes[0,1].loglog([r['spacing'] for r in main_grids],[r['L1_finite_grid_vs_continuum_SG'] for r in main_grids],'o-')
    axes[0,1].set(xlabel='Grid spacing',ylabel='L1 error',title='Projected grid vs grand-potential SG')
    rows=[r for r in absorption if r['time_model']==.3 and r['rate_model']>0]
    axes[1,0].loglog([r['rate_model'] for r in rows],[r['absorption_error'] for r in rows],'o-')
    axes[1,0].set(xlabel='Charge rate (model units)',ylabel='Absorbed-mass difference',title='Sector-dependent synthetic sink')
    axes[1,1].plot([r['time_per_half'] for r in quenched],[r['L1_semigroup_failure'] for r in quenched],'o-',label='Reset hidden configuration')
    axes[1,1].plot([r['time_per_half'] for r in quenched],[r['L1_average_generator_failure'] for r in quenched],'s-',label='Average generator')
    axes[1,1].set(xlabel='Half interval (model units)',ylabel='L1 error',title='Quenched mixture retains memory')
    axes[1,1].legend(fontsize=8)
    for ax in axes.flat: ax.grid(alpha=.2)
    fig.suptitle('Synthetic validation only: no silicon rate or fracture calibration')
    fig.savefig(out/'charge_reduction.png',dpi=170)
    plt.close(fig)
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    audit(p.parse_args().output)
