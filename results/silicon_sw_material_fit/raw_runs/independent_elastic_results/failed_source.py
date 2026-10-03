"""Energy and strain-work tangents from independent conservative E/F/stress code."""
from pathlib import Path
import sys,json,time,hashlib
import numpy as np
from scipy.optimize import minimize
cache=Path(__file__).resolve().parent;root=cache.parents[2]/'aft-silicon-wafer';sys.path.insert(0,str(root))
from solver_v1.silicon_sw_material_fit import force_features
from results.silicon_wafer_feasibility.run_static_probe import source_parameters
out=cache/'independent_elastic_results';assert not out.exists();out.mkdir();start=time.perf_counter()
reference=json.loads((cache/'results/bulk_diagnostic.json').read_text())['original']
p,_=source_parameters();lattice=reference['lattice_A'];cell=lattice/2*np.array([[0,1,1],[1,0,1],[1,1,0.]])
basis=np.array([[0.,0.,0.],[lattice/4]*3]);volume=abs(np.linalg.det(cell));pbc=np.ones(3,bool)
e0=force_features(basis,cell,pbc,p)[0].sum();calls=1
modes=dict(hydro=np.eye(3),tetragonal=np.diag([1.,-1.,0.]),
           simple_xy=np.array([[0.,1.,0.],[0.,0.,0.],[0.,0.,0.]]))
states=[];tangents=[]
for name,direction in modes.items():
    for step in [.01,.003,.001,.0003]:
        for relaxed in [False,True]:
            current=[]
            for sign in [-1,1]:
                deformation=np.eye(3)+sign*step*direction;cc=cell@deformation.T;xx=basis@deformation.T
                def evaluator(u,full=False):
                    global calls
                    calls+=1;x=xx.copy();x[1]+=u
                    result=force_features(x,cc,pbc,p,strain_derivative=full)
                    return result
                def objective(u):
                    e,f=evaluator(u)
                    return e.sum(),-f[1].sum(axis=1)
                u=np.zeros(3)
                if relaxed:
                    result=minimize(objective,u,jac=True,method='BFGS',options=dict(gtol=1e-10,maxiter=100))
                    u=result.x
                e,f,derivative=evaluator(u,True)
                force=float(np.max(abs(f.sum(axis=2))))
                if relaxed and force>2e-8:raise ValueError('independent internal relaxation unresolved')
                piola=derivative.sum(axis=2)@np.linalg.inv(deformation).T/volume
                row=dict(mode=name,step=step,sign=sign,internally_relaxed=relaxed,
                         energy_eV=float(e.sum()),force_max_eV_A=force,shift_A=u.tolist(),
                         conjugate_Piola_eV_A3=float(np.sum(piola*direction)))
                states.append(row);current.append(row)
            minus,plus=current
            ee=(plus['energy_eV']+minus['energy_eV']-2*e0)/(volume*step**2)*160.2176634
            ss=(plus['conjugate_Piola_eV_A3']-minus['conjugate_Piola_eV_A3'])/(2*step)*160.2176634
            tangents.append(dict(mode=name,step=step,internally_relaxed=relaxed,
                                 energy_tangent_GPa=ee,stress_tangent_GPa=ss,difference_GPa=ee-ss))
constants=[]
for step in [.01,.003,.001,.0003]:
    for relaxed in [False,True]:
        for kind in ['energy','stress']:
            t={r['mode']:r[kind+'_tangent_GPa'] for r in tangents if r['step']==step and r['internally_relaxed']==relaxed}
            s,d=t['hydro']/3,t['tetragonal']/2
            constants.append(dict(step=step,internally_relaxed=relaxed,method=kind,
                                  C11_GPa=(s+2*d)/3,C12_GPa=(s-d)/3,C44_GPa=t['simple_xy']))
published=json.loads((cache/'framework_sources/test-results/model-SW-test-bulk_diamond-properties.json').read_text())
summary=dict(complete=True,new_conservative_E_F_calls=calls,new_DFT=0,new_MD=0,new_LAMMPS=0,
    lattice_A=lattice,reference_jet_C44_GPa=reference['C44_GPa'],
    independent_fine_relaxed=[r for r in constants if r['step']==.0003 and r['internally_relaxed']],
    published_framework_SW_C44_GPa=published['diamond_c44'],
    published_SW_epsilon_eV=2.1675,current_epsilon_eV=p['epsilon'],
    published_finite_fit_protocol_fully_reproduced=False,
    elapsed_seconds=time.perf_counter()-start,
    source_sha256=hashlib.sha256((root/'solver_v1/silicon_sw_material_fit.py').read_bytes()).hexdigest())
for name,object_ in [('states.json',states),('tangents.json',tangents),('constants.json',constants),('summary.json',summary)]:
    (out/name).write_text(json.dumps(object_,indent=2,allow_nan=False)+'\n')
print(json.dumps(summary))
