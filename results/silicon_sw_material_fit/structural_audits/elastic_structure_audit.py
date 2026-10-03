"""Derive nearest-neighbor SW elastic constraint from independent bond geometry."""
import json
from pathlib import Path
import numpy as np
n=np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]])/np.sqrt(3)
strain=[]
for axis in range(3):
    e=np.zeros((3,3));e[axis,axis]=1;strain.append(e)
for a,b in [(1,2),(0,2),(0,1)]:
    e=np.zeros((3,3));e[a,b]=e[b,a]=.5;strain.append(e)
radial=np.zeros((9,9));angular=np.zeros((9,9))
for sign in [1.,-1.]:
    ns=sign*n
    for v in ns:
        g=np.r_[[v@e@v for e in strain],sign*v]
        radial+=.5*np.outer(g,g)
    for i in range(4):
        for j in range(i+1,4):
            v,w=ns[i],ns[j];c=v@w
            g=np.r_[[((w-c*v)@e@v+(v-c*w)@e@w) for e in strain],
                    sign*((w-c*v)+(v-c*w))]
            angular+=2*np.outer(g,g)
def elastic(k,w):
    h=k*radial+w*angular
    relaxed=h[:6,:6]-h[:6,6:]@np.linalg.solve(h[6:,6:],h[6:,:6])
    return relaxed[[0,0,3],[0,1,3]]
data=json.loads((Path(__file__).parent/'results/bulk_diagnostic.json').read_text())
records=[]
for key,value in data.items():
    c=np.array([value['C11_GPa'],value['C12_GPa'],value['C44_GPa']])
    matrix=np.array([[radial[0,0],angular[0,0]],[radial[0,1],angular[0,1]]])
    coefficients=np.linalg.solve(matrix,c[:2])
    pred=elastic(*coefficients)
    records.append(dict(label=key,actual=c.tolist(),coefficients=coefficients.tolist(),
                        geometry_predicted=pred.tolist(),C44_residual=float(pred[2]-c[2])))
target=np.array([153.3,56.3,72.2]);coeff=np.linalg.solve(matrix,target[:2]);pred=elastic(*coeff)
result=dict(radial_hessian=radial.tolist(),angular_hessian=angular.tolist(),
            controls=records,published_DFT_target_GPa=target.tolist(),
            target_matching_C11_C12_predicted_GPa=pred.tolist(),
            target_C44_discrepancy_GPa=float(pred[2]-target[2]))
(Path(__file__).parent/'elastic_structure_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if not k.endswith('hessian')}))
print('pair',radial[:6,:6].round(12));print('angle',angular[:6,:6].round(12))
print('pair coupling',radial[3,6:],radial[6:,6:].round(12))
print('angle coupling',angular[3,6:],angular[6:,6:].round(12))
