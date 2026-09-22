"""Exact fixed-total charge vs fixed-mean charge vs reservoir ensemble.

Synthetic single-occupation spin-resolved orbitals, no Coulomb term or Si data. Enumerates
small finite occupation spaces deterministically; no trajectory sampling.
"""
from itertools import combinations
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.special import expit,logsumexp


def coefficients(size):
    index=np.arange(size)
    return .018*np.cos(index+.3),.4+.07*index,.1+.015*index


def exact_canonical(q,number,kt):
    q=np.asarray(q,float);offset,g,h=coefficients(len(q))
    energies=offset+g*q+.5*h*q*q;gradient=g+h*q
    configs=[]
    for selected in combinations(range(len(q)),number):
        n=np.zeros(len(q));n[list(selected)]=1.;configs.append(n)
    configs=np.array(configs)
    values=configs@energies
    logz=logsumexp(-values/kt);w=np.exp(-values/kt-logz)
    mean=w@configs
    centered=configs-mean
    covariance=(centered.T*w)@centered
    return dict(free_energy=-kt*logz,gradient=mean*gradient,
        hessian=np.diag(mean*h)-gradient[:,None]*covariance*gradient[None,:]/kt,
        mean_occupations=mean,occupation_covariance=covariance,states=len(configs))


def fixed_mean(q,number,kt,mu_fixed=None):
    q=np.asarray(q,float);offset,g,h=coefficients(len(q))
    energies=offset+g*q+.5*h*q*q;gradient=g+h*q
    if mu_fixed is None:
        mu=brentq(lambda m:sum(expit((m-energies)/kt))-number,
            min(energies)-100*kt,max(energies)+100*kt,xtol=1e-14)
    else:mu=mu_fixed
    n=expit((mu-energies)/kt);variance=n*(1-n)
    omega=-kt*sum(np.logaddexp(0,(mu-energies)/kt))
    hessian=np.diag(n*h-variance*gradient**2/kt)
    if mu_fixed is None:
        v=variance*gradient/kt
        hessian+=np.outer(v,v)/(sum(variance)/kt)
    return dict(free_energy=omega+mu*number if mu_fixed is None else omega,
                gradient=n*gradient,hessian=hessian,mu=mu,mean_occupations=n)


def main(out):
    out.mkdir(parents=True,exist_ok=True)
    records=[];snapshots=[]
    for sites in [2,4,8,16]:
        number=sites//2
        for kt in [.01,.025,.10]:
            q=.015*np.sin(np.arange(sites)+.4)
            exact=exact_canonical(q,number,kt)
            mean=fixed_mean(q,number,kt)
            reservoir=fixed_mean(q,number,kt,mu_fixed=mean['mu'])
            for ensemble,reference,fn in [
                    ('exact_total_N',exact,lambda x:exact_canonical(x,number,kt)),
                    ('fixed_mean_N',mean,lambda x:fixed_mean(x,number,kt)),
                    ('fixed_mu_reservoir',reservoir,lambda x:fixed_mean(x,number,kt,mu_fixed=mean['mu']))]:
                for step in [1e-3,1e-4,1e-5]:
                    fd=np.empty((sites,sites))
                    gradient_fd=np.empty(sites)
                    for j in range(sites):
                        qp,qm=q.copy(),q.copy();qp[j]+=step;qm[j]-=step
                        positive,negative=fn(qp),fn(qm)
                        fd[:,j]=(positive['gradient']-negative['gradient'])/(2*step)
                        gradient_fd[j]=(positive['free_energy']-negative['free_energy'])/(2*step)
                    error=float(np.max(abs(fd-reference['hessian'])))
                    gradient_error=float(np.max(abs(gradient_fd-reference['gradient'])))
                    records.append(dict(sites=sites,total_electrons=number,kT=kt,ensemble=ensemble,
                        step=step,Hessian_max_abs_FD_error=error,gradient_max_abs_FD_error=gradient_error))
                if error>3e-4 or gradient_error>1e-5:raise AssertionError('free energy derivative finite-difference failure')
            if abs(sum(exact['mean_occupations'])-number)>2e-12:raise AssertionError('canonical number not conserved')
            if np.max(abs(exact['occupation_covariance'].sum(axis=0)))>2e-12:raise AssertionError('canonical covariance constraint failed')
            cross=exact['hessian'].copy();np.fill_diagonal(cross,0)
            snapshots.append(dict(sites=sites,total_electrons=number,kT=kt,
                exact_free_energy=exact['free_energy'],fixed_mean_free_energy=mean['free_energy'],
                exact_vs_mean_Hessian_max_difference=float(np.max(abs(exact['hessian']-mean['hessian']))),
                exact_cross_Hessian_max_abs=float(np.max(abs(cross))),
                exact_covariance_row_sum_residual=float(np.max(abs(exact['occupation_covariance'].sum(axis=0)))),
                canonical_configurations=exact['states'],
                exact_mean_occupations=exact['mean_occupations'].tolist(),
                fixed_mean_occupations=mean['mean_occupations'].tolist(),
                exact_Hessian=exact['hessian'].tolist(),fixed_mean_Hessian=mean['hessian'].tolist(),
                reservoir_Hessian=reservoir['hessian'].tolist()))
    with (out/'finite_difference_checks.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
    (out/'ensemble_comparisons.json').write_text(json.dumps(snapshots,indent=2)+'\n',encoding='utf-8')
    report=dict(status='synthetic finite global-charge ensemble validation',cases=len(snapshots),
        site_definition='single-occupation spin-resolved electronic orbital; not an atomic site',
        derivative_checks=len(records),maximum_finest_Hessian_error=max(r['Hessian_max_abs_FD_error'] for r in records if r['step']==1e-5),
        maximum_finest_gradient_error=max(r['gradient_max_abs_FD_error'] for r in records if r['step']==1e-5),
        globally_fixed_total_charge_is_not_fixed_local_mean_charge=True,
        electrostatics_included=False,Si_calibration=False,physical_time_calibrated=False)
    (out/'global_charge_summary.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)
