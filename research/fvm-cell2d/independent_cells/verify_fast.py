import json,time
import numpy as np
import fvm2d as m
from fast_propagation import advance

m.set_model('coherent_g4')
g=m.Grid2D.original();p,_=m.initial_mass(g);op=m.assemble(g,4.5)
t=time.perf_counter()
mass=np.stack([p,p]);F=np.zeros(2)
reference,rf,ra=m.advance_batch(op,mass,F,.3)
actual,af,aa=advance(op,mass,F,.3)
report=dict(L1=float(abs(actual-reference).sum(axis=1).max()),CDF=float(abs(af-rf).max()),
    fast_audit=aa,reference_audit=ra,seconds=time.perf_counter()-t)
report['accepted']=report['L1']<2e-12 and report['CDF']<2e-15
(m.ROOT/'outputs/fvm-cell2d/FAST_ACTION_CHECK.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2),flush=True)
assert report['accepted']
