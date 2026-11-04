"""New matched-marginal/energy comparison on captured original SG operators."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import json
import sys
import time
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from solver_v1.model import TwoRowLJ, ModelParams
from certify_frozen_table import apply, encode

def main():
    sys.set_int_max_str_digits(50000)
    start=time.perf_counter();root=Path(__file__).parent;checks=[]
    def check(name,ok,value=None):
        checks.append({'name':name,'passed':bool(ok),'value':value})
        if not ok:raise AssertionError(name)
    d=np.load(root/'SOURCE_OPERATOR_PAIR.npz');steps=d['steps'];a=d['a'];s=d['s']
    common=np.logical_and.reduce(d['valid'][:,:-1])
    indices=np.flatnonzero(common);nr=len(a);ns=len(s);size=len(indices)
    aa,ss=np.meshgrid(a,s,indexing='ij')
    model=TwoRowLJ(ModelParams(n_cells=1,kT=.02,mobility_a=1.,mobility_s=.05,chi_axial_projection=.2))
    energy=model.energy_gradient_batch(aa.reshape(-1,1),ss.reshape(-1,1),0.)[0].ravel()
    energy=energy[indices]
    sparse=[csr_matrix(t) for t in steps]
    def response(cycles):
        r=np.zeros(steps.shape[1]);r[-1]=1
        for _ in range(cycles):
            for t in sparse[::-1]:r=t.T@r
        return r[indices]
    responses={n:response(n) for n in [1,2,4,8]}
    rowids=indices//ns;colids=indices%ns
    O=np.vstack([*(rowids==i for i in np.unique(rowids)),*(colids==j for j in np.unique(colids)),
                 (energy-energy.mean())/max(1.,np.ptp(energy))]).astype(float)
    uniform=np.ones(size)/size; y=O@uniform
    # Constraint family frozen before evaluation at 2,4,8 cycles.
    a0=aa.ravel()[indices];s0=ss.ravel()[indices]
    cov=(a0-a0.mean())*(s0-s0.mean())
    even=(a0-a0.mean())*((s0-s0.mean())**2-np.mean((s0-s0.mean())**2))
    constraints={'marginals_energy':O,
                 'plus_simple_covariance':np.vstack([O,cov/max(np.max(abs(cov)),1e-15)]),
                 'plus_even_correlation':np.vstack([O,even/max(np.max(abs(even)),1e-15)]),
                 'plus_one_cycle_event_memory':np.vstack([O,responses[1]/max(np.ptp(responses[1]),1e-15)])}
    def envelope(r,A,b):
        scale=max(float(np.ptp(r)),1e-15);objective=(r-r.min())/scale
        opts={'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9}
        low=linprog(objective,A_eq=A,b_eq=b,bounds=(0,None),method='highs',options=opts)
        high=linprog(-objective,A_eq=A,b_eq=b,bounds=(0,None),method='highs',options=opts)
        if not(low.success and high.success):raise ValueError((low.message,high.message))
        residual=max(np.max(abs(A@low.x-b)),np.max(abs(A@high.x-b)))
        check('LP_feasible',residual<1e-8,float(residual))
        return {'lower':float(r@low.x),'upper':float(r@high.x),'width':float(r@(high.x-low.x)),
                'max_constraint_residual':float(residual)},low.x,high.x
    comparisons=[];extremes={}
    for n,r in responses.items():
        row={'cycles':n,'uniform_reference_cdf':float(r@uniform),'summaries':{}}
        for name,A in constraints.items():
            rec,lo,hi=envelope(r,A,A@uniform);row['summaries'][name]=rec
            if name=='marginals_energy':extremes[n]=(lo,hi)
        comparisons.append(row)
    # A response moment is only useful if it is known accurately enough.
    # Test independent information uncertainty rather than silently treating it as exact.
    noisy=[];r1=responses[1];z1=float(r1@uniform)
    for n in [2,4,8]:
        r=responses[n];scale=max(np.ptp(r),1e-15);objective=(r-r.min())/scale
        for noise in [1e-6,1e-5,1e-4]:
            inequalities=np.vstack([r1,-r1]);bounds=np.array([z1+noise,-z1+noise])
            opts={'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9}
            low=linprog(objective,A_ub=inequalities,b_ub=bounds,A_eq=O,b_eq=y,bounds=(0,None),method='highs',options=opts)
            high=linprog(-objective,A_ub=inequalities,b_ub=bounds,A_eq=O,b_eq=y,bounds=(0,None),method='highs',options=opts)
            check('noisy_probe_LP_solved',low.success and high.success)
            check('noisy_probe_LP_feasible',max(np.max(abs(O@low.x-y)),np.max(abs(O@high.x-y)),np.max(inequalities@low.x-bounds),np.max(inequalities@high.x-bounds))<2e-9)
            noisy.append({'cycles':n,'probe_absolute_uncertainty':noise,
                          'lower':float(r@low.x),'upper':float(r@high.x),
                          'width':float(r@(high.x-low.x))})
    # Exact rational correlation witness: sum of two rectangle circulation vectors.
    # All support lies inside every mask; there is no forced initial mask crossing.
    local={int(j):i for i,j in enumerate(indices)};cycles=[]
    for i,k in combinations(np.unique(rowids),2):
        for j,l in combinations(range(ns),2):
            cells=[i*ns+j,k*ns+l,i*ns+l,k*ns+j]
            if not all(c in local for c in cells):continue
            v=np.zeros(size)
            for c,sign in zip(cells,[1,1,-1,-1]):v[local[c]]=sign
            cycles.append(v)
    V=np.asarray(cycles);ev=V@energy;rv=V@responses[1]
    best=(-1,None,None)
    for i,v in enumerate(V):
        z=ev[:,None]*v[None,:]-ev[i]*V
        mass=np.maximum(z,0).sum(axis=1)
        gap=np.divide(abs(ev*rv[i]-ev[i]*rv),mass,out=np.zeros_like(mass),where=mass>1e-12)
        j=int(np.argmax(gap))
        if gap[j]>best[0]:best=(float(gap[j]),i,j)
    ef=[F.from_float(float(x)) for x in energy]
    v1,v2=V[best[1]].astype(int),V[best[2]].astype(int)
    E1=sum((ef[k]*int(v1[k]) for k in range(size)),F(0))
    E2=sum((ef[k]*int(v2[k]) for k in range(size)),F(0))
    z=[E2*int(v1[k])-E1*int(v2[k]) for k in range(size)]
    mass=sum((max(x,F(0)) for x in z),F(0))
    if mass<=0:raise ValueError('NO_EXACT_WITNESS')
    plus=[max(x,F(0))/mass for x in z];minus=[max(-x,F(0))/mass for x in z]
    if sum(float(plus[k]-minus[k])*responses[1][k] for k in range(size))<0:plus,minus=minus,plus
    for rows in (rowids,colids):
        for level in np.unique(rows):
            check('exact_marginal_equality',sum((plus[k]-minus[k] for k in range(size) if rows[k]==level),F(0))==0)
    check('exact_energy_mean_equality',sum(((plus[k]-minus[k])*ef[k] for k in range(size)),F(0))==0)
    check('exact_unit_mass',sum(plus,F(0))==sum(minus,F(0))==1)
    # To avoid point-mass claims only, mix each witness with the SAME common law.
    # The equality of both marginals and mean energy survives this smoothing exactly.
    smoothing=F(1,2)
    plus=[smoothing*x+(1-smoothing)/size for x in plus]
    minus=[smoothing*x+(1-smoothing)/size for x in minus]
    check('smoothed_preparations_strictly_positive',min(plus+minus)>0)
    check('smoothed_preparations_exact_energy',sum(((plus[k]-minus[k])*ef[k] for k in range(size)),F(0))==0)
    for rows in (rowids,colids):
        for level in np.unique(rows):
            check('smoothed_exact_marginal',sum((plus[k]-minus[k] for k in range(size) if rows[k]==level),F(0))==0)
    preparation={'indices':indices.tolist(),'a':a.tolist(),'s':s.tolist(),
                 'energy_binary64_hex':[float(x).hex() for x in energy],
                 'plus':[[str(x.numerator),str(x.denominator)] for x in plus],
                 'minus':[[str(x.numerator),str(x.denominator)] for x in minus],
                 'meaning':'Declared exact rational initial measures, not measured physical data'}
    (root/'CORRELATION_PREPARATION_EXACT.json').write_text(json.dumps(preparation,indent=2),encoding='utf-8')
    fullplus=[F(0)]*steps.shape[1];fullminus=fullplus.copy()
    for k,i in enumerate(indices):fullplus[i]=plus[k];fullminus[i]=minus[k]
    exact_columns=[];max_repair=F(0)
    for t in steps:
        columns=[]
        for j in range(t.shape[1]):
            entries={int(i):F.from_float(float(t[i,j])) for i in np.flatnonzero(t[:,j])}
            adjustment=1-sum(entries.values(),F(0));anchor=max(entries,key=entries.get)
            entries[anchor]+=adjustment
            check('repair_positive',entries[anchor]>=0)
            max_repair=max(max_repair,abs(adjustment));columns.append(list(entries.items()))
        exact_columns.append(columns)
    xp,xm=fullplus,fullminus
    for col in exact_columns:xp=apply(col,xp);xm=apply(col,xm)
    exact_gap=xp[-1]-xm[-1]
    check('positive_exact_cdf_gap',exact_gap>0)
    pf=np.array([float(x) for x in plus]);mf=np.array([float(x) for x in minus])
    witness={'support':'Intersection of ALL masks, 75 interior cells; both preparations strictly positive on this support',
             'common_mixture_weight':str(smoothing),'initial_absorbed_mass':0,
             'same_a_marginal_exact':True,'same_s_marginal_exact':True,'same_mean_energy_exact':True,
             'corrected_exact_cdf_plus':encode(xp[-1]),'corrected_exact_cdf_minus':encode(xm[-1]),
             'exact_gap':encode(exact_gap),'any_same_summary_predictor_max_error_lower_bound':encode(exact_gap/2),
             'maximum_column_repair':encode(max_repair),
             'energy_target':'Exact-real values of captured binary64 energy, not exact continuum W'}
    for n,r in responses.items():
        witness[f'captured_table_cdf_{n}_cycles']=[float(r@pf),float(r@mf)]
    check('float_exact_corrected_difference_small',abs(float(exact_gap)-float(responses[1]@(pf-mf)))<1e-12)
    for row in comparisons:
        baseline=row['summaries']['marginals_energy']['width']
        check('extra_information_cannot_widen_LP',all(c['width']<=baseline+1e-9 for c in row['summaries'].values()))
    np.savez_compressed(root/'CORRELATION_MEMORY_COMPARISON.npz',indices=indices,energy=energy,O=O,
                        plus=pf,minus=mf,uniform=uniform,cycles=np.array(list(responses)),
                        responses=np.array(list(responses.values())))
    (root/'CORRELATION_WITNESS_EXACT.json').write_text(json.dumps(witness,indent=2),encoding='utf-8')
    result={'status':'PASS_MATCHED_MARGINAL_ENERGY_SOURCE_COMPARISON','dateKST':'2026-11-04',
            'common_interior_states':size,'original_fine_states':steps.shape[1],
            'constraints_frozen_before_heldout_horizons':True,'train_cycles':[1],'heldout_cycles':[2,4,8],
            'heldout_meaning':'Other horizons of same fixed operator cycle, not independent load/material data',
            'comparisons':comparisons,'noisy_probe_comparisons':noisy,'witness':witness,'checks_count':len(checks),
            'checks_all_passed':all(x['passed'] for x in checks),
            'check_count_meaning':'Includes repeated column positivity assertions, not independent experiments',
            'new_source_energy_evaluation':True,
            'old_source_solver_suite_or_MD_repeated':False,'physical_material_or_clock_validation':False,
            'formal_novelty_established':False,'wall_seconds':time.perf_counter()-start}
    (root/'CORRELATION_MEMORY_RESULT.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps({'status':result['status'],'seconds':result['wall_seconds'],'witness_gap':float(exact_gap),
                      'comparisons':comparisons},ensure_ascii=False))

if __name__=='__main__':main()
