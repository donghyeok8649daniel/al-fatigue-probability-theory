"""New tests of the reduced model, rather than reused 4D action results."""
from pathlib import Path
import json,time,hashlib
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply
import fvm2d as m
import sys
sys.path.insert(0,str(m.ROOT/'independent_cells'))
from run_family import save_json


def main():
    start=time.perf_counter();m.set_model('coherent_g4')
    points=np.array([[.75,-.03],[.8,.02],[.86,.07],[1.1,.12],[2.,.18]])
    a,s=points.T;d=a-m.A0;j=m.original
    independent=2*j.local(m.A0,0.)+2*j.local(a,s)+j.phi(3.)+j.phi(1.)
    for shift in (1+s,2-s,2+s,1-s):independent+=j.phi(np.hypot(d,shift))
    f=4.5;independent-=2*f*(d+.2*s)
    energy_difference=float(abs(m.energy(points,f)-independent).max())
    derivative_checks=[]
    h=2e-6
    for q in points[:3]:
        value,grad,hess=m.energy_jet(q,f)
        fdgrad=np.array([(m.energy((q+np.eye(2)[i]*h)[None,:],f)[0]-m.energy((q-np.eye(2)[i]*h)[None,:],f)[0])/(2*h) for i in range(2)])
        fdhess=np.column_stack([(m.energy_jet(q+np.eye(2)[i]*h,f)[1]-m.energy_jet(q-np.eye(2)[i]*h,f)[1])/(2*h) for i in range(2)])
        derivative_checks.append(dict(gradient_relative=float(abs(fdgrad-grad).max()/max(1.,abs(grad).max())),
            hessian_relative=float(abs(fdhess-hess).max()/max(1.,abs(hess).max()))))
    projected_mobility=np.linalg.inv(m.PROJECTION.T@np.diag(1/m.MOBILITY)@m.PROJECTION)
    mobility_difference=float(abs(projected_mobility-np.diag(m.MOBILITY[:2]/2)).max())
    grid=m.Grid2D(np.r_[np.linspace(.7,1.45,25),1.8,2.5,3.5,5.,6.],np.linspace(-.08,.25,18))
    p,_=m.initial_mass(grid);mass=np.stack([p,np.roll(p,3)])
    op=m.assemble(grid,4.5)
    htime=.03
    result,F,audit=m.advance_batch(op,mass,np.zeros(2),htime)
    K=op.matrix;n=len(op.keep)
    Q=sparse.bmat([[K,sparse.csr_matrix((n,1))],[sparse.csr_matrix(op.event[None,:]),sparse.csr_matrix((1,1))]],format='csr')
    sweep=mass.sum(axis=1)-mass[:,op.keep].sum(axis=1)
    inputs=np.vstack([mass[:,op.keep].T,sweep[None,:]])
    reference=expm_multiply(Q*htime,inputs,traceA=htime*float(K.diagonal().sum()))
    actual=np.vstack([result[:,op.keep].T,F[None,:]])
    action_L1=float(abs(reference-actual).sum(axis=0).max())
    action_CDF=float(abs(reference[-1]-F).max())
    solo,singleF,_=m.advance_batch(op,mass[:1],np.zeros(1),htime)
    isolation=float(abs(solo-result[:1]).sum()+abs(singleF-F[:1]).sum())
    reflected=m.assemble(grid,-4.5)
    assert len(reflected.keep)==grid.size and not np.any(reflected.event)
    e=m.energy(grid.points,-4.5)
    gibbs=np.exp(-(e-e.min())/m.THETA)*grid.volume.ravel();gibbs/=gibbs.sum()
    stationarity=float(abs(reflected.matrix@gibbs).sum()/max(1.,(-reflected.matrix.diagonal()*gibbs).sum()))
    # A joint density is stored and evolved, never replaced by product marginals.
    pdf=p.reshape(grid.shape);aa,ss=np.meshgrid(*grid.x,indexing='ij')
    covariance=float((pdf*aa*ss).sum()-(pdf*aa).sum()*(pdf*ss).sum())
    product_distance=float(abs(pdf-np.outer(pdf.sum(axis=1),pdf.sum(axis=0))).sum())
    report=dict(model=m.model_metadata(),energy_expansion_absolute_difference=energy_difference,
        finite_difference=derivative_checks,projected_mobility_difference=mobility_difference,
        action_against_independent_augmented_expm_L1=action_L1,
        action_against_independent_augmented_expm_CDF=action_CDF,independent_row_isolation_error=isolation,
        reflecting_gibbs_stationarity_scaled=stationarity,initial_a_s_covariance=covariance,
        initial_joint_vs_product_marginals_L1=product_distance,
        a_s_joint_density_preserved=True,per_step_renormalization=False,
        audit=audit,operator_checks=[op.metadata,reflected.metadata],seconds=time.perf_counter()-start)
    report['accepted']=(energy_difference<1e-12 and mobility_difference<1e-15
        and all(c['gradient_relative']<1e-6 and c['hessian_relative']<1e-6 for c in derivative_checks)
        and action_L1<2e-12 and action_CDF<2e-15 and isolation<2e-12
        and stationarity<2e-12 and product_distance>1e-6)
    save_json(m.ROOT/'outputs/fvm-cell2d/COMPONENT_CHECKS.json',report)
    print(json.dumps(report,indent=2),flush=True)
    return 0 if report['accepted'] else 1


if __name__=='__main__':raise SystemExit(main())
