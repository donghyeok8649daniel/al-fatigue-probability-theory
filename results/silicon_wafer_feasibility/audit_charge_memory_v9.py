"""Independent Schur/Laplace check for exact finite-rate electronic memory."""
import argparse
import csv
import json
from pathlib import Path
import sys
import numpy as np
from scipy.linalg import expm

sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from solver_v1.silicon_charge_dynamics import charge_memory_blocks,coupled_charge_generator


def main(out):
    out.mkdir(parents=True,exist_ok=True)
    rows=[];kernels=[]
    cells=16;q=np.linspace(-1,1,cells)
    for sectors in [2,3]:
        f=np.array([.5*(q+.2*j)**2+.01*j for j in range(sectors)])
        for rate in [0.,3.,30.,300.,3000.]:
            r=coupled_charge_generator(f,excess_counts=np.arange(sectors),chemical_potential=.02,
                thermal_energy=.12,spacing=q[1]-q[0],mobilities=np.linspace(.3,1,sectors),charge_rate=rate)
            m=charge_memory_blocks(r)
            a,b,c,d=[m[k].toarray() for k in ['A','B','C','D']]
            for omega in [0.,.1,1.,10.,100.]:
                z=.4+omega*1j
                memory=b@np.linalg.solve(z*np.eye(len(d))-d,c)
                reduced=np.linalg.solve(z*np.eye(cells)-a-memory,np.eye(cells))
                full=r['projection']@np.linalg.solve(z*np.eye(cells*sectors)-r['generator'].toarray(),r['lift'].toarray())
                markov=np.linalg.solve(z*np.eye(cells)-a,np.eye(cells))
                rows.append(dict(sectors=sectors,charge_rate_model=rate,laplace_real=.4,laplace_imaginary=omega,
                    exact_memory_vs_full_max_abs=float(np.max(abs(reduced-full))),
                    dropped_memory_max_abs=float(np.max(abs(markov-full))),
                    memory_matrix_norm=float(np.linalg.norm(memory,2))))
            for t in [0.,.001,.01,.1,1.]:
                kernels.append(dict(sectors=sectors,charge_rate_model=rate,time_model=t,
                    kernel_norm=float(np.linalg.norm(b@expm(t*d)@c,2))))
    for name,data in [('resolvent_comparisons.csv',rows),('memory_kernel.csv',kernels)]:
        with (out/name).open('w',encoding='utf-8',newline='') as stream:
            w=csv.DictWriter(stream,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
    report=dict(status='exact fixed-parameter algebra on synthetic charge dynamics',
        cases=len(rows),kernel_cases=len(kernels),
        maximum_resolvent_error=max(r['exact_memory_vs_full_max_abs'] for r in rows),
        maximum_error_if_memory_omitted=max(r['dropped_memory_max_abs'] for r in rows),
        physical_Si_rate=None,physical_time_calibrated=False)
    (out/'memory_summary.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)
