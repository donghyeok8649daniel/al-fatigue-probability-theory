"""Synthetic validation of initiation first-passage bookkeeping; no Si rates."""
import argparse,csv,json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from solver_v1.silicon_charge_dynamics import sg_chain
from solver_v1.silicon_initiation_probability import (
    initiation_collectors,evolve_first_passage,basin_committor,first_passage_moments)


def write_csv(path,rows):
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def main(out):
    if out.exists() and any(out.iterdir()):raise ValueError('fresh output required')
    out.mkdir(parents=True,exist_ok=True)
    g=np.zeros((5,5))
    for s,d,k in [(0,1,2),(1,0,3),(0,2,1),(2,0,4),(1,3,.4),(2,4,.5),(3,1,7),(4,2,9)]:g[d,s]+=k;g[s,s]-=k
    model=initiation_collectors(g,{'surface_demo':[3],'interface_demo':[4]})
    early=initiation_collectors(g,{'all_intermediate_entries':[1,2,3,4]})
    t=np.r_[0.,np.geomspace(.001,30,90)]
    actual=evolve_first_passage(model,[1,0,0],t);premature=evolve_first_passage(early,[1],t)
    rows=[dict(time_model_units=float(ti),survival=float(actual['survival'][i]),
        surface_demo_CDF=float(actual['cause_cumulative'][i,0]),interface_demo_CDF=float(actual['cause_cumulative'][i,1]),
        premature_intermediate_CDF=float(1-premature['survival'][i]),mass_residual=float(actual['mass_residual'][i])) for i,ti in enumerate(t)]
    write_csv(out/'synthetic_first_passage.csv',rows)
    moments=first_passage_moments(model)
    h=basin_committor(g,intact_core=[0],initiated_basins=[3,4])
    refinement=[];D=.2;time=2.;k=np.arange(60);odd=2*k+1
    exact=4/np.pi*np.sum((-1.)**k/odd*np.exp(-D*(odd*np.pi/2)**2*time))
    for count in (31,61,121):
        x=np.linspace(-1,1,count)
        sg=sg_chain(np.zeros(count),thermal_energy=1,spacing=x[1]-x[0],face_mobility=D)
        bm=initiation_collectors(sg,{'left':[0],'right':[count-1]});p=np.zeros(count-2);p[(count-2)//2]=1
        evolved=evolve_first_passage(bm,p,[time]);got=float(evolved['survival'][0])
        mt=first_passage_moments(bm)['mean_first_passage_time']
        refinement.append(dict(grid_points=count,spacing=x[1]-x[0],survival=got,continuum_survival=exact,
            error=abs(got-exact),MFPT_max_error=float(np.max(abs(mt-(1-x[1:-1]**2)/(2*D))))))
    write_csv(out/'diffusion_refinement.csv',refinement)
    summary=dict(scope='synthetic generator and model time; no actual Si kinetics or calibrated initiation probability',
        committor=h['committor'].tolist(),mean_first_passage_time_model_units=moments['mean_first_passage_time'].tolist(),
        eventual_cause_probabilities=moments['eventual_cause_probability'].tolist(),
        maximum_mass_residual=float(np.max(abs(actual['mass_residual']))),
        minimum_probability=actual['minimum_probability'],
        diffusion_error_ratios=[refinement[i]['error']/refinement[i+1]['error'] for i in (0,1)],
        Si_initiation_probability=None,Si_physical_clock=None,new_atomic_calculations=0)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8');print(json.dumps(summary,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
