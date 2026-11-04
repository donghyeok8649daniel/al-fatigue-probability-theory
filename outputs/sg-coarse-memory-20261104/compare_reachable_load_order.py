"""New load-order probes of the original SG generator on the captured grid.

Continuous-time finite SG chain with killed mask; not continuum path certification.
"""
from pathlib import Path
from types import SimpleNamespace
import json
import sys
import time
import numpy as np
from scipy.linalg import expm

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from solver_v1.model import ModelParams,TwoRowLJ
from solver_v1.probability_pde_2d import _sg_generator_2d

def main():
    start=time.perf_counter();root=Path(__file__).parent
    d=np.load(root/'SOURCE_OPERATOR_PAIR.npz');a,s=d['a'],d['s'];size=len(a)*len(s)
    grid=SimpleNamespace(da=float(a[1]-a[0]),ds=float(s[1]-s[0]))
    model=TwoRowLJ(ModelParams(n_cells=1,kT=.02,mobility_a=1.,mobility_s=.05,chi_axial_projection=.2))
    aa,ss=np.meshgrid(a,s,indexing='ij');dt=.0006096053822325939
    masks={1.:d['valid'][0,:-1],1.4:d['valid'][8,:-1],.8:d['valid'][16,:-1]}
    common=np.logical_and.reduce(list(masks.values()))
    G={};checks=[]
    def check(name,ok,value=None):
        checks.append({'name':name,'pass':bool(ok),'value':value})
        if not ok:raise AssertionError(name)
    for f in masks:
        energy=model.energy_gradient_batch(aa.reshape(-1,1),ss.reshape(-1,1),f)[0].reshape(aa.shape)
        G[f]=_sg_generator_2d(energy,model,grid).toarray()
        check('source_generator_mass',np.max(abs(G[f].sum(axis=0)))<2e-12)
    def jump(mask):
        J=np.eye(size+1);J[:-1,:-1]=np.diag(mask.astype(float));J[-1,:-1]=~mask
        return J
    def killed(f,mask):
        A=np.zeros((size+1,size+1));A[:-1,:-1]=mask[:,None]*G[f]*mask[None,:]
        A[-1,:-1]=-A[:-1,:-1].sum(axis=0)
        # The intact columns have all exiting generator mass routed to cemetery.
        check('killed_generator_mass',np.max(abs(A.sum(axis=0)))<2e-12)
        off=A.copy();np.fill_diagonal(off,0)
        check('killed_offdiagonal_positive',off.min()>=-2e-12,float(off.min()))
        return A
    cases=[]
    x0=np.r_[common.astype(float)/common.sum(),0.]
    for mode in ['original_load_masks','common_fixed_mask']:
        mm={f:(masks[f] if mode=='original_load_masks' else common) for f in masks}
        A={f:killed(f,mm[f]) for f in masks};J={f:jump(mm[f]) for f in masks}
        for h in [4,16,64]:
            t=h*dt
            TA=expm(t*A[.8])@J[.8];TB=expm(t*A[1.4])@J[1.4]
            p=J[1.]@TB@TA@x0;q=J[1.]@TA@TB@x0
            rec={'mode':mode,'preparation_steps_equivalent':h,'each_window_model_time':t,
                 'A_force':.8,'B_force':1.4,'common_probe_force':1.,
                 'preparation_absorbed_AB_BA':[float(p[-1]),float(q[-1])],'probes':[]}
            for duration in [0,8,32,128]:
                C=expm(duration*dt*A[1.]);rp=C@p;rq=C@q
                gap=float(rp[-1]-rq[-1]);newgap=float((rp[-1]-p[-1])-(rq[-1]-q[-1]))
                check('direct_composition_matches_projection',abs(gap-float(C[-1]@(p-q)))<1e-13)
                check('probability_bounds',min(rp.min(),rq.min())>=-2e-13 and max(rp.sum(),rq.sum())<1+2e-12)
                rec['probes'].append({'probe_steps_equivalent':duration,'cdf_AB':float(rp[-1]),'cdf_BA':float(rq[-1]),
                                      'total_gap':gap,'new_absorption_gap_after_preparation':newgap})
            cases.append(rec)
        if mode=='common_fixed_mask':
            comm=A[1.4]@A[.8]-A[.8]@A[1.4]
            # A common future event makes the leading coefficient observable.
            probe=expm(32*dt*A[1.]);coefficient=float(probe[-1]@comm@x0)
            scaling=[]
            for h in [1.,.5,.25,.125]:
                t=h*dt
                gap=float(probe[-1]@(expm(t*A[1.4])@expm(t*A[.8])-expm(t*A[.8])@expm(t*A[1.4]))@x0)
                scaling.append({'window_time':t,'gap':gap,'gap_over_t_squared':gap/t**2,
                                'error_to_commutator_coefficient':abs(gap/t**2-coefficient)})
            check('short_window_commutator_convergence',scaling[-1]['error_to_commutator_coefficient']<scaling[0]['error_to_commutator_coefficient'])
    result={'status':'PASS_SOURCE_SG_REACHABLE_LOAD_ORDER_COMPARISON','dateKST':'2026-11-04',
            'source':'Original N1 SG generator, captured grid and original opening masks',
            'event_scope':'First exit of continuous-time finite SG chain; original continuum and explicit mask timestepping are distinct',
            'initial':'Uniform on75 cells inside intersection of every original load mask, declared numerical preparation',
            'cases':cases,'fixed_mask_commutator_coefficient':coefficient,'short_window_scaling':scaling,
            'checks_count':len(checks),'checks_all_passed':all(c['pass'] for c in checks),
            'new_protocols_not_old_suite':True,'material_or_physical_time_validation':False,
            'independent_review':False,'wall_seconds':time.perf_counter()-start}
    (root/'REACHABLE_LOAD_ORDER_RESULT.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result))

if __name__=='__main__':main()
