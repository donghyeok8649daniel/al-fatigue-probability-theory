"""Joint spatial/charge resolution of a specified synthetic first-entry boundary.

Two equal diffusivities on [0,L], charge rates nu*p and nu*(1-p).
Both sectors reflect at 0; sector 0 reflects and sector 1 absorbs at L.
No actual Si rates, basins, dopant fractions, or physical time are supplied.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, sys
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from solver_v1.silicon_initiation_probability import check_generator, initiation_collectors


def generator(intervals, occupation, rate, diffusion=.4, length=1.):
    n=intervals+1;h=length/intervals
    volume=np.full(n,h);volume[[0,-1]]=h/2
    left=np.arange(n-1);right=left+1
    forward=diffusion/(h*volume[left]);reverse=diffusion/(h*volume[right])
    chain=sparse.coo_matrix((np.r_[forward,-forward,reverse,-reverse],
        (np.r_[right,left,left,right],np.r_[left,left,right,right])),shape=(n,n)).tocsc()
    eye=sparse.eye(n,format='csc')
    g=sparse.block_diag([chain,chain],format='csc')+sparse.bmat(
        [[-rate*occupation*eye,rate*(1-occupation)*eye],
         [rate*occupation*eye,-rate*(1-occupation)*eye]],format='csc')
    return check_generator(g),volume


def exact_continuum(p,nu,diffusion=.4,length=1.):
    k=np.sqrt(nu/diffusion)
    return length**2/(2*diffusion)+(1-p)*length/(p*diffusion*k*np.tanh(k*length))


def exact_discrete(intervals,p,nu,diffusion=.4,length=1.):
    h=length/intervals
    theta=2*np.arcsinh(.5*h*np.sqrt(nu/diffusion))
    return length**2/(2*diffusion)+(1-p)*length*h/(p*diffusion*np.sinh(theta)*np.tanh(intervals*theta))


def main(args):
    if args.output.exists():raise ValueError('fresh continuum-boundary audit required')
    args.output.mkdir(parents=True)
    rows=[];details={};max_residual=0.;max_relative=0.;max_symmetry=0.
    for p in (.02,.2,.8):
        for nu in (.01,1.,100.,10000.,1000000.):
            continuum=exact_continuum(p,nu)
            for intervals in (8,16,32,64,128,256,512,1024):
                g,vol=generator(intervals,p,nu);n=intervals+1
                model=initiation_collectors(g,{'specified_sector1_boundary':[2*n-1]})
                q=model['transient_generator'].tocsc()
                mean=spsolve(q.T,-np.ones(q.shape[0]))
                if np.any(mean<=0) or not np.all(np.isfinite(mean)):
                    raise ValueError('unresolved synthetic mean first-entry time')
                value=(1-p)*mean[0]+p*mean[n]
                reference=exact_discrete(intervals,p,nu)
                error=abs(value-reference)/reference
                residual=float(np.max(abs(q.T@mean+1)))
                equilibrium=np.r_[(1-p)*vol,p*vol]/vol.sum()
                conductance=g@sparse.diags(equilibrium)
                asym=conductance-conductance.T
                symmetry=float(np.max(abs(asym.data))) if asym.nnz else 0.
                scale=float(np.max(abs(conductance.data)))
                max_symmetry=max(max_symmetry,symmetry/scale)
                max_relative=max(max_relative,error);max_residual=max(max_residual,residual)
                if error>1e-7 or symmetry>1e-12*scale:
                    raise ValueError('closed-form or detailed-balance replay failed')
                row=dict(occupation=p,charge_rate=nu,intervals=intervals,diffusion=.4,length=1.,
                    mean_time_numerical=float(value),mean_time_discrete_exact=float(reference),
                    mean_time_continuum_exact=float(continuum),fast_absorbing_mean_time=1.25,
                    discrete_formula_relative_error=float(error),linear_system_max_residual=residual,
                    continuum_relative_error=float((value-continuum)/continuum),
                    excess_over_fast_boundary=float(value-1.25),continuum_excess=float(continuum-1.25),
                    charge_to_grid_rate_ratio=float(nu/(.4*intervals**2)),
                    cells_across_charge_boundary_layer=float(intervals*np.sqrt(.4/nu)))
                rows.append(row)
                # Keep actual solved vectors for representative resolutions,
                # without turning every matrix entry into a redundant artifact.
                if p==.2 and nu in (1.,10000.) and intervals in (16,256,1024):
                    key=f'p02_nu{nu:g}_N{intervals}'
                    details[key]=mean
    path=args.output/'joint_resolution.csv'
    with path.open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    np.savez_compressed(args.output/'selected_backward_solutions.npz',**details)
    # Independent analytic backward equations and mixed boundary conditions.
    formula_checks=[]
    for p in (.02,.2,.8):
        for nu in (.01,1.,100.,10000.):
            d=.4;k=np.sqrt(nu/d);x=np.linspace(0.,1.,129)
            # Stable cosh(k*x)/sinh(k) and sinh(k*x)/sinh(k).
            ratio_c=(np.exp(k*(x-1))+np.exp(-k*(x+1)))/(-np.expm1(-2*k))
            ratio_s=(np.exp(k*(x-1))-np.exp(-k*(x+1)))/(-np.expm1(-2*k))
            delta=ratio_c/(p*d*k);delta_prime=ratio_s/(p*d)
            average=(1-x*x)/(2*d)+(1-p)*delta[-1]
            t0=average+p*delta;t1=average-(1-p)*delta
            t0_second=-1/d+p*k*k*delta;t1_second=-1/d-(1-p)*k*k*delta
            residual0=d*t0_second+nu*p*(t1-t0)+1
            residual1=d*t1_second+nu*(1-p)*(t0-t1)+1
            checks=dict(occupation=p,charge_rate=nu,
                backward_equation_error=float(max(np.max(abs(residual0)),np.max(abs(residual1)))),
                reflected_zero_derivative_error=float(max(abs(p*delta_prime[0]),abs((1-p)*delta_prime[0]))),
                reflected_right_derivative_error=float(abs(-1/d+p*delta_prime[-1])),
                absorbed_right_value_error=float(abs(t1[-1])),
                initial_mean_formula_error=float(abs((1-p)*t0[0]+p*t1[0]-exact_continuum(p,nu))))
            if max(v for key,v in checks.items() if key.endswith('error'))>1e-8:
                raise ValueError('continuum backward equations or boundary conditions failed')
            formula_checks.append(checks)
    result=dict(complete=True,synthetic_backward_systems=len(rows),analytic_boundary_checks=len(formula_checks),
        maximum_discrete_formula_relative_error=max_relative,maximum_backward_system_absolute_residual=max_residual,
        maximum_detailed_balance_relative_error=max_symmetry,formula_checks=formula_checks,
        continuum_mean_time='L^2/(2D) + (1-p)*L*coth(L*sqrt(nu/D))/(p*sqrt(D*nu))',
        discrete_mean_time='L^2/(2D) + (1-p)*L*h*coth(N*theta)/(p*D*sinh(theta)); theta=2*asinh(h*sqrt(nu/D)/2)',
        continuum_fast_excess='(1-p)*L/(p*sqrt(D*nu))',
        fixed_grid_fast_excess='2*(1-p)*L/(p*h*nu)',
        distinction='same absorbing limit, different leading finite-rate corrections; h must resolve sqrt(D/nu)',
        initial_state='x=0 with charge occupation (1-p,p)',spatial_boundary='both reflect at 0; sector0 reflects and sector1 absorbs at L',
        occupation_is_not_dopant_concentration=True,physical_clock=None,new_MD=0,new_DFT=0,new_potential_calls=0,
        source_csv_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        material_approved=False,first_crack_calibrated=False)
    (args.output/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4.7),layout='constrained')
    for nu in (1.,100.,10000.,1000000.):
        subset=[r for r in rows if r['occupation']==.2 and r['charge_rate']==nu]
        axes[0].loglog([r['intervals'] for r in subset],
            [abs(r['continuum_relative_error']) for r in subset],'.-',label='nu='+str(nu))
    axes[0].set(xlabel='Spatial intervals N',ylabel='Relative mean-time error vs continuum');axes[0].legend(fontsize=8)
    rates=np.logspace(-2,7,200)
    axes[1].loglog(rates,[exact_continuum(.2,nu)-1.25 for nu in rates],color='black',label='Continuum exact')
    for n in (16,64,256):
        axes[1].loglog(rates,[exact_discrete(n,.2,nu)-1.25 for nu in rates],label='Fixed N='+str(n))
    axes[1].set(xlabel='Synthetic charge rate nu',ylabel='Mean-time excess above absorbing limit');axes[1].legend(fontsize=8)
    fig.suptitle('Charge-dependent initiation boundary: spatial resolution matters\nSynthetic D=0.4, L=1, charge occupation p=0.2; no physical clock',fontsize=11)
    fig.savefig(args.output/'joint_charge_spatial_resolution.png',dpi=180);plt.close(fig)
    print(json.dumps({k:v for k,v in result.items() if k!='formula_checks'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args())
