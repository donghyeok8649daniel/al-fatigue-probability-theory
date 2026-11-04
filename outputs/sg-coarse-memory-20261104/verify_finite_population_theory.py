"""Focused synthetic checks of the new theory, not material calibration."""
from fractions import Fraction as F
from itertools import product
from math import factorial
from pathlib import Path
import json
import time
import numpy as np
from scipy.linalg import expm

def main():
    start=time.perf_counter(); checks=[]
    def check(name,ok,value=None):
        checks.append({'name':name,'pass':bool(ok),'value':None if value is None else float(value)})
        if not ok:raise AssertionError(name)
    N=3; d=3
    u=[F(0),F(1,7),F(2,5)];g=[F(-1,5),F(1,10),F(2,5)]
    v=F(2,3);K=F(7,5);theta=F(2,7);nu=np.array([1/6,1/3,1/2])
    states=list(product(range(d),repeat=N)); index={x:i for i,x in enumerate(states)}
    counts=sorted(set(tuple(x.count(i) for i in range(d)) for x in states))
    cindex={x:i for i,x in enumerate(counts)}
    H=np.zeros((len(counts),len(states)))
    for j,x in enumerate(states):H[cindex[tuple(x.count(i) for i in range(d))],j]=1
    def energy(x,e):
        m=sum((g[i] for i in x),F(0))/N
        return sum((u[i] for i in x),F(0))+N*v*K*(e-m)**2/2
    def B(z):return 1.0 if abs(z)<1e-14 else z/np.expm1(z)
    def generators(e):
        G=np.zeros((len(states),len(states)));wrong=G.copy();hetero=G.copy()
        for col,x in enumerate(states):
            m=sum((g[i] for i in x),F(0))/N;sig=K*(e-m)
            for r,i in enumerate(x):
                for j in range(d):
                    if abs(i-j)!=1:continue
                    y=list(x);y[r]=j;y=tuple(y);row=index[y];dg=g[j]-g[i]
                    de=energy(y,e)-energy(x,e)
                    formula=u[j]-u[i]-v*sig*dg+v*K*dg*dg/(2*N)
                    assert de==formula
                    sig_after=sig-K*dg/N
                    assert de==u[j]-u[i]-v*(sig+sig_after)*dg/2
                    rate=B(float(de/theta))/nu[i]
                    G[row,col]+=rate;G[col,col]-=rate
                    bad=B(float((u[j]-u[i]-v*sig*dg)/theta))/nu[i]
                    wrong[row,col]+=bad;wrong[col,col]-=bad
                    extra=F(1,4)*(j-i) if r==0 else F(0)
                    bad=B(float((de+extra)/theta))/nu[i]
                    hetero[row,col]+=bad;hetero[col,col]-=bad
        Gc=np.zeros((len(counts),len(counts)))
        for col,n in enumerate(counts):
            m=sum((n[i]*g[i] for i in range(d)),F(0))/N;sig=K*(e-m)
            for i in range(d):
                for j in range(d):
                    if n[i]==0 or abs(i-j)!=1:continue
                    target=list(n);target[i]-=1;target[j]+=1
                    de=u[j]-u[i]-v*sig*(g[j]-g[i])+v*K*(g[j]-g[i])**2/(2*N)
                    rate=n[i]*B(float(de/theta))/nu[i];row=cindex[tuple(target)]
                    Gc[row,col]+=rate;Gc[col,col]-=rate
        return G,Gc,wrong,hetero
    e=F(3,10);G,Gc,wrong,hetero=generators(e)
    check('exact_rational_energy_and_midpoint_all_edges',True)
    pi=np.array([np.exp(-float(energy(x,e)/theta))*np.prod([nu[i] for i in x]) for x in states]);pi/=pi.sum()
    flux=G*pi[None,:]
    check('full_detailed_balance',np.max(abs(flux-flux.T))<2e-13,np.max(abs(flux-flux.T)))
    pic=np.array([factorial(N)/np.prod([factorial(ni) for ni in n])*np.prod(nu**np.array(n))*np.exp(-float(energy(tuple(i for i,ni in enumerate(n) for _ in range(ni)),e)/theta)) for n in counts]);pic/=pic.sum()
    check('multinomial_equilibrium_pushforward',np.max(abs(H@pi-pic))<2e-13,np.max(abs(H@pi-pic)))
    check('occupancy_detailed_balance',np.max(abs(Gc*pic[None,:]-(Gc*pic[None,:]).T))<2e-13)
    check('exact_occupancy_generator_intertwining',np.max(abs(H@G-Gc@H))<2e-13,np.max(abs(H@G-Gc@H)))
    check('before_stress_omission_breaks_Gibbs',np.max(abs(wrong@pi))>1e-5,np.max(abs(wrong@pi)))
    HG=H@hetero;diam=0
    for row in range(len(counts)):
        cols=np.flatnonzero(H[row]);diam=max(diam,float(np.max(np.ptp(HG[:,cols],axis=1))))
    check('heterogeneous_sites_break_untyped_histogram_closure',diam>1e-5,diam)
    rng=np.random.default_rng(20261104);p=rng.random(len(states));p/=p.sum()
    check('full_occupancy_semigroup_agreement',np.max(abs(H@expm(.37*G)@p-expm(.37*Gc)@H@p))<2e-13)
    energies=np.array([float(energy(x,e)) for x in states]);weights=np.array([np.prod([nu[i] for i in x]) for x in states])
    mu=energies+float(theta)*(np.log(p/weights)+1)
    diss=0.0
    for i in range(len(states)):
        for j in range(i):
            a=p[i]*G[j,i];b=p[j]*G[i,j]
            if a>0 and b>0:diss+=float(theta)*(a-b)*np.log(a/b)
    derivative=float(mu@(G@p))
    check('closed_free_energy_identity',abs(derivative+diss)<2e-12,abs(derivative+diss))
    check('closed_dissipation_nonnegative',diss>=0,diss)
    alive=np.array([2 not in x for x in states]);ac=np.array([n[2]==0 for n in counts])
    A=G[np.ix_(alive,alive)];Ac=Gc[np.ix_(ac,ac)];Ha=H[np.ix_(ac,alive)]
    rho=p[alive].copy();rhoc=Ha@rho
    check('killed_generator_intertwining',np.max(abs(Ha@A-Ac@Ha))<2e-13)
    check('killed_CDF_agreement',abs(np.sum(expm(.7*A)@rho)-np.sum(expm(.7*Ac)@rhoc))<2e-13)
    exitrate=-A.sum(axis=0);q=rho/rho.sum();m=np.array([sum(float(g[i]) for i in x)/N for x in np.array(states)[alive]])
    Lint=A+np.diag(exitrate);h=float(exitrate@q)
    direct=float(m@(A@rho)/rho.sum()-(m@q)*(np.sum(A@rho)/rho.sum()))
    rhs=float(m@(Lint@q)-(m*exitrate)@q+(m@q)*h)
    check('survivor_selection_covariance',abs(direct-rhs)<2e-13,abs(direct-rhs))
    Fsur=float((energies[alive]+float(theta)*(np.log(rho/weights[alive])+1))@(A@rho))
    boundary=float((rho*exitrate)@(energies[alive]+float(theta)*(np.log(rho/weights[alive])+1)))
    dint=0.0
    for i in range(len(rho)):
        for j in range(i):
            a=rho[i]*A[j,i];b=rho[j]*A[i,j]
            if a>0 and b>0:dint+=float(theta)*(a-b)*np.log(a/b)
    check('raw_survivor_free_energy_boundary_balance',abs(Fsur+dint+boundary)<2e-12,abs(Fsur+dint+boundary))
    _,G2,_,_=generators(F(-1,10));A2=G2[np.ix_(ac,ac)]
    ta,tb=.2,.3;PA=expm(ta*Ac);PB=expm(tb*A2);cycle=PB@PA
    reverse=PA@PB
    r0=np.zeros(sum(ac))
    r0[np.array(counts)[ac].tolist().index([N,0,0])]=1
    sequence_gap=float(abs(np.sum(cycle@r0)-np.sum(reverse@r0)))
    check('load_order_changes_first_exit_probability',sequence_gap>1e-8,sequence_gap)
    def integral(A,t):
        k=len(A);block=np.zeros((2*k,2*k));block[:k,:k]=A;block[:k,k:]=np.eye(k)
        return expm(t*block)[:k,k:]
    ell=np.ones(len(Ac))@(integral(Ac,ta)+integral(A2,tb)@PA)
    mean=float(ell@np.linalg.solve(np.eye(len(cycle))-cycle,r0))
    z=r0.copy();summed=0.0
    for _ in range(200):summed+=float(ell@z);z=cycle@z
    check('periodic_mean_resolvent_vs_time_integral_sum',abs(mean-summed)<2e-11,abs(mean-summed))
    check('finite_mean_positive',mean>0,mean)
    result={'status':'PASS_DECLARED_FINITE_POPULATION_THEORY_CHECKS','new_checks':len(checks),'checks':checks,
            'synthetic_only':True,'physical_material_calibration':False,'independent_review':False,
            'full_states':len(states),'occupancy_states':len(counts),'alive_full_states':int(sum(alive)),
            'alive_occupancy_states':int(sum(ac)),'finite_population':N,
            'exact_fraction_energy_all_edges':True,'cyclic_first_exit_mean_model_time':mean,
            'wall_seconds':time.perf_counter()-start,'prior_science_suites_repeated':False}
    (Path(__file__).parent/'FINITE_POPULATION_THEORY_CHECKS.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result))

if __name__=='__main__':main()
