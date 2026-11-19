import time,json
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply
import fvm2d as m
from cram_action import action

m.set_model('coherent_g4')
start=time.perf_counter()
grid=m.Grid2D.original(0,0,1);op=m.assemble(grid,4.5)
p,_=m.initial_mass(grid)
mass=p[None,:]
exact,F,a=m.advance_uniform_batch(op,mass,np.zeros(1),.3)
n=len(op.keep)
Q=sparse.bmat([[op.matrix,sparse.csr_matrix((n,1))],[sparse.csr_matrix(op.event[None,:]),sparse.csr_matrix((1,1))]],format='csc')
z=np.r_[p[op.keep],0.]
results=[]
for substeps in (1,2):
    t=time.perf_counter();actual,c=action(Q,z,.3,substeps)
    results.append(dict(substeps=substeps,seconds=time.perf_counter()-t,
        L1=float(abs(actual-np.r_[exact[0,op.keep],F[0]]).sum()),CDF=float(abs(actual[-1]-F[0])),
        negative=float(-np.minimum(actual,0).sum()),balance=float(abs(actual.sum()-1.)),check=c))
report=dict(results=results,seconds=time.perf_counter()-start)
report['accepted']=all(r['L1']<2e-12 and r['CDF']<2e-15 and r['negative']<2e-12 and r['balance']<1e-9 for r in results)
(m.ROOT/'outputs/fvm-cell2d/CRAM_ACTION_CHECK.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2),flush=True)
assert report['accepted']
