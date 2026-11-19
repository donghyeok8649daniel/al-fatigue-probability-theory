"""Independent-cell, two-coordinate conservative first-passage FVM.

mass[cell, a_bin * n_s + s_bin] is a table of local distributions. A cell
index is never another coordinate of a joint probability grid. The local
Two explicit representative models are available. coherent_g4 retains every
old G4 energy term on the constrained manifold a1=a2=a, s1=s2=s. n1_local
retains only W. Neither is an exact marginal evolution of the unconstrained
G4 process. Physical-cell independence alone does not establish either
microscopic closure. The selected closure must accompany every result.
"""

from dataclasses import dataclass
from pathlib import Path
import sys

import numpy as np
from scipy import ndimage,sparse
from scipy.optimize import root

ROOT=Path(__file__).resolve().parents[1]
CANDIDATE=ROOT/'outputs/fvm-integration/candidate'
sys.path.insert(0,str(CANDIDATE))
import joint_model as original
from joint_model import A0,THETA,MOBILITY,local,local_jet
from mesh_plan_loader import original_axes

MODEL='coherent_g4'
PROJECTION=np.array([[1.,0.],[0.,1.],[1.,0.],[0.,1.]])


def set_model(model):
    global MODEL
    if model not in ('n1_local','coherent_g4'):
        raise ValueError(model)
    MODEL=model


def model_metadata():
    return dict(model_id=MODEL,active_dimension=2,independent_physical_cells=True,
        energy_formula='W(a,s)-f*((a-A0)+0.2*s)' if MODEL=='n1_local' else 'G4(a,s,a,s;f)',
        mobility=(MOBILITY[:2] if MODEL=='n1_local' else MOBILITY[:2]/2).tolist(),
        theta=THETA,old_four_coordinate_dynamics_exactly_preserved=False,
        microscopic_closure=('one local W representative; explicit mobile-site pairs omitted'
            if MODEL=='n1_local' else 'holonomic coherent motion a1=a2=a, s1=s2=s; all original G4 energy terms retained'),
        physical_calibration_status='not calibrated to seconds, MPa, characteristic volume or 30um initiation')


def energy(points,force):
    if MODEL=='coherent_g4':
        return original.energy(np.asarray(points)@PROJECTION.T,force)
    a,s=np.asarray(points).T
    return local(a,s)-force*((a-A0)+.2*s)


def energy_jet(q,force):
    if MODEL=='coherent_g4':
        value,gradient,hessian=original.jet(PROJECTION@q,force)
        return value,PROJECTION.T@gradient,PROJECTION.T@hessian@PROJECTION
    v,a,s,aa,ass,ss=local_jet(*q)
    return float(v-force*((q[0]-A0)+.2*q[1])),np.array([a-force,s-.2*force]),np.array([[aa,ass],[ass,ss]])


def reference(force):
    def gradient(q):
        return energy_jet(q,force)[1]
    def hessian(q):
        return energy_jet(q,force)[2]
    for guess in ((A0,0.),(.85,.05),(.9,.1)):
        try:
            result=root(gradient,guess,jac=hessian,options={'xtol':1e-11})
        except ValueError:
            continue
        q=result.x
        if .7<q[0]<1.45 and -.08<q[1]<.25 and np.linalg.norm(gradient(q),np.inf)<1e-8:
            if np.linalg.eigvalsh(hessian(q)).min()>0:
                return dict(state=q.tolist(),energy=float(energy(q[None,:],force)[0]),
                            gradient_residual=float(np.linalg.norm(gradient(q),np.inf)))
    raise ArithmeticError(f'No resolved intact minimum for {MODEL} at force={force}')


def dual_weights(edge):
    x=(edge[:-1]+edge[1:])/2
    return np.r_[x[0]-edge[0]+(x[1]-x[0])/2,(x[2:]-x[:-2])/2,
                 edge[-1]-x[-1]+(x[-1]-x[-2])/2]


@dataclass
class Grid2D:
    a_edges:np.ndarray
    s_edges:np.ndarray

    def __post_init__(self):
        self.edges=[np.asarray(self.a_edges),np.asarray(self.s_edges)]
        self.x=[(edge[:-1]+edge[1:])/2 for edge in self.edges]
        self.shape=tuple(len(x) for x in self.x)
        self.size=self.shape[0]*self.shape[1]
        self.coords=np.array(np.unravel_index(np.arange(self.size),self.shape))
        self.points=np.column_stack([x[c] for x,c in zip(self.x,self.coords)])
        self.strides=(self.shape[1],1)
        self.width=[np.diff(edge) for edge in self.edges]
        self.volume=dual_weights(self.edges[0])[:,None]*dual_weights(self.edges[1])[None,:]

    @classmethod
    def original(cls,a_level=0,s_level=0,outer=1):
        return cls(*original_axes(a_level,s_level,outer))


def event_region(grid,force,ref):
    labels,_=ndimage.label(energy(grid.points,force).reshape(grid.shape)<=ref['energy'])
    seeds=np.unique(labels[-1,:]);seeds=seeds[seeds!=0]
    return np.isin(labels,seeds).ravel()


def initial_mass(grid,force=1.):
    ref=reference(force)
    values=energy(grid.points,force).reshape(grid.shape)
    event=event_region(grid,force,ref).reshape(grid.shape)
    seed=tuple(int(np.argmin(abs(x-q))) for x,q in zip(grid.x,ref['state']))
    low,high=float(values[seed]),ref['energy']+5.
    for _ in range(37):
        level=(low+high)/2
        labels,_=ndimage.label(values<level)
        core=labels==labels[seed]
        if np.any(core&event):high=level
        else:low=level
    labels,_=ndimage.label(values<low)
    assert labels[seed]!=0
    core=labels==labels[seed]
    weights=np.exp(-(values-ref['energy'])/THETA)*grid.volume*core
    norm=float(weights.sum())
    assert norm>0 and np.isfinite(norm)
    return (weights/norm).ravel(),dict(reference=ref,normalization_count=1,
        intact_basin_nodes=int(core.sum()),normalizer=norm,barrier_interval=[low,high])


def first_crossing(q0,q1,force,level):
    previous=energy(q0,force)-level
    assert np.all(previous>-2e-12)
    first=np.full(len(q0),-1,dtype=int)
    for k in range(1,65):
        value=energy(q0+(q1-q0)*(k/64),force)-level
        take=(first<0)&(previous>0)&(value<=0)
        first[take]=k
        previous=value
    assert np.all(first>0)
    low,high=(first-1)/64,first/64
    for _ in range(48):
        mid=(low+high)/2
        positive=energy(q0+(q1-q0)*mid[:,None],force)>level
        low=np.where(positive,mid,low);high=np.where(positive,high,mid)
    t=(low+high)/2
    point=q0+(q1-q0)*t[:,None]
    return point,t,float(abs(energy(point,force)-level).max(initial=0.))


def assemble(grid,force,order=8):
    ref=reference(force)
    event_mask=event_region(grid,force,ref)
    keep=np.flatnonzero(~event_mask)
    n=len(keep)
    mapping=np.full(grid.size,-1,dtype=np.int32);mapping[keep]=np.arange(n)
    values=energy(grid.points,force)
    outgoing,events=np.zeros(n),np.zeros(n)
    rows,cols,data=[],[],[]
    nodes,weights=np.polynomial.legendre.leggauss(order)
    max_root=0.
    for axis,mobility in enumerate(MOBILITY[:2]/(2 if MODEL=='coherent_g4' else 1)):
        dual=np.tile(grid.width[axis][grid.coords[axis,keep]],(2,1))
        pending=[]
        for side,sign in enumerate((-1,1)):
            coordinate=grid.coords[axis,keep]
            source=np.flatnonzero((coordinate+sign>=0)&(coordinate+sign<grid.shape[axis]))
            gi=keep[source];gj=gi+sign*grid.strides[axis]
            q0,q1=grid.points[gi],grid.points[gj].copy()
            distance=abs(q1[:,axis]-q0[:,axis])
            target=event_mask[gj]
            if np.any(target):
                cut,fraction,residual=first_crossing(q0[target],q1[target],force,ref['energy'])
                q1[target]=cut;distance[target]*=fraction
                max_root=max(max_root,residual)
            dual[side,source]=distance
            resistance=np.full(len(source),-np.inf)
            for t,weight in zip((nodes+1)/2,weights/2):
                resistance=np.logaddexp(resistance,
                    (energy(q0+(q1-q0)*t,force)-values[gi])/THETA+np.log(weight))
            rate=mobility*THETA/distance*np.exp(-resistance)
            pending.append((source,gj,target,rate))
        for source,gj,target,rate in pending:
            rate=rate/dual.mean(axis=0)[source]
            events[source[target]]+=rate[target]
            outgoing[source]+=rate
            inside=~target
            rows.append(mapping[gj[inside]]);cols.append(source[inside]);data.append(rate[inside])
    matrix=sparse.coo_matrix((np.r_[np.concatenate(data),-outgoing],
        (np.r_[np.concatenate(rows),np.arange(n)],np.r_[np.concatenate(cols),np.arange(n)])),shape=(n,n)).tocsr()
    defect=np.asarray(matrix.sum(axis=0)).ravel()+events
    scaled=float(np.max(abs(defect)/np.maximum(1.,outgoing)))
    assert scaled<2e-14 and max_root<1e-10
    return Operator2D(force,grid.size,keep,matrix,events,
        dict(reference=ref,model=model_metadata(),state_dimension=2,shape=list(grid.shape),states=n,nonzeros=matrix.nnz,
             sparse_bytes=sum(x.nbytes for x in (matrix.data,matrix.indices,matrix.indptr)),
             scaled_column_defect=scaled,maximum_event_root_residual=max_root))


@dataclass
class Operator2D:
    force:float
    common_states:int
    keep:np.ndarray
    matrix:object
    event:np.ndarray
    metadata:dict


def advance_uniform_batch(operator,mass,cdf,duration):
    """Positive Poisson action on multiple local PDFs, with no cell-product grid."""
    sys.path.insert(0,str(CANDIDATE/'reference_kernel/outputs/raw-sg-convergence-20261108'))
    from positive_uniformization import poisson_band,stable_poisson_weights
    assert mass.ndim==2 and mass.shape[1]==operator.common_states
    physical_event=np.ones(operator.common_states,dtype=bool);physical_event[operator.keep]=False
    moving=mass[:,physical_event].sum(axis=1)
    initial_f=cdf+moving
    pk=mass[:,operator.keep].T.copy(order='C')
    nu=1.01*float(np.max(-operator.matrix.diagonal()))
    transition=operator.matrix/nu+sparse.eye(len(operator.keep),format='csr')
    assert transition.data.min()>=0
    rates=operator.event/nu
    kappa=max(1.,float((np.asarray(transition.sum(axis=0)).ravel()+rates).max())+32*np.finfo(float).eps)
    mu=nu*duration
    coarse=poisson_band(mu,kappa,1e-16);fine=poisson_band(mu,kappa,1e-20)
    if fine[1]>200_000:
        return advance_stiff_batch(operator,mass,cdf,duration)
    weights=stable_poisson_weights(mu,fine[0],fine[1])
    yf,yc=np.zeros_like(pk),np.zeros_like(pk)
    absorbed=np.zeros(mass.shape[0]);ff=np.zeros_like(absorbed);fc=np.zeros_like(absorbed)
    for power in range(fine[1]+1):
        if power>=fine[0]:
            w=weights[power-fine[0]]
            yf+=w*pk;ff+=w*absorbed
            if coarse[0]<=power<=coarse[1]:yc+=w*pk;fc+=w*absorbed
        if power==fine[1]:break
        absorbed+=rates@pk
        pk=transition@pk
    result=np.zeros_like(mass);result[:,operator.keep]=yf.T
    final_f=initial_f+ff
    l1=np.abs(yf-yc).sum(axis=0)+abs(ff-fc)
    cdf_error=abs(ff-fc)
    balance=abs(result.sum(axis=1)+final_f-1.)
    negative=-np.minimum(result,0.).sum(axis=1)
    assert np.max(l1)<2e-12 and np.max(cdf_error)<2e-15
    assert np.max(balance)<1e-9 and np.max(negative)<2e-12
    assert np.all(final_f>=cdf-2e-15)
    audit=dict(state_dimension=2,independent_cell_count=mass.shape[0],mass_table_shape=list(result.shape),
        common_operator_shape=list(operator.matrix.shape),poisson_powers=fine[1],
        maximum_successive_L1=float(l1.max()),maximum_successive_CDF=float(cdf_error.max()),
        maximum_mass_residual=float(balance.max()),maximum_negative_mass=float(negative.max()),
        explicit_moving_event_increment=moving.tolist(),first_passage_flux=ff.tolist(),
        poisson_mean=mu,poisson_tail_bound=fine[2],probability_table_bytes=result.nbytes,
        numerical_mass_loss_added_to_event=False,probability_renormalized=False)
    return result,final_f,audit


def advance_stiff_batch(operator,mass,cdf,duration):
    """Full sparse linear BDF action, checked at two independent tolerances.

    First-passage probability is an explicit absorbing state driven by event
    flux. Signed roundoff is retained, without clipping or normalization.
    Sparse LU storage grows with the 2D local grid, never a joint grid of cells.
    """
    from scipy.integrate import solve_ivp
    n=len(operator.keep)
    Q=sparse.bmat([[operator.matrix,sparse.csr_matrix((n,1))],
        [sparse.csr_matrix(operator.event[None,:]),sparse.csr_matrix((1,1))]],format='csc')
    moving=mass.sum(axis=1)-mass[:,operator.keep].sum(axis=1)
    # Sum the event mask directly; subtracting two near-unit masses loses tiny events.
    mask=np.ones(operator.common_states,dtype=bool);mask[operator.keep]=False
    moving=mass[:,mask].sum(axis=1)
    result=np.zeros_like(mass);final_f=np.empty_like(cdf);audits=[]
    cached={}
    for row in range(len(mass)):
        initial=np.r_[mass[row,operator.keep],0.]
        # Identical independent replicas share a linear action, not a joint law.
        key=initial.tobytes()
        if key in cached:
            fine,coarse,nfev=cached[key]
        else:
            solutions=[];nfev=[]
            for rtol,atol in ((2e-13,1e-20),(2.5e-14,1e-22)):
                sol=solve_ivp(lambda t,y:Q@y,(0.,duration),initial,method='BDF',jac=Q,
                    rtol=rtol,atol=atol,t_eval=[duration])
                if not sol.success:raise RuntimeError(sol.message)
                solutions.append(sol.y[:,-1]);nfev.append(sol.nfev)
            coarse,fine=solutions;cached[key]=(fine,coarse,nfev)
        result[row,operator.keep]=fine[:-1]
        final_f[row]=cdf[row]+moving[row]+fine[-1]
        audits.append(dict(L1=float(abs(fine-coarse).sum()),CDF=float(abs(fine[-1]-coarse[-1])),nfev=nfev))
    l1=max(x['L1'] for x in audits);df=max(x['CDF'] for x in audits)
    balance=float(abs(result.sum(axis=1)+final_f-1.).max())
    negative=float((-np.minimum(result,0.)).sum(axis=1).max())
    if not (l1<2e-12 and df<2e-15 and balance<1e-9 and negative<2e-12):
        raise ArithmeticError(dict(L1=l1,CDF=df,balance=balance,negative=negative,solves=audits))
    if not np.all(final_f>=cdf-2e-15):raise ArithmeticError('Nonmonotone first passage')
    return result,final_f,dict(method='full_sparse_BDF_successive_tolerances',state_dimension=2,
        independent_cell_count=len(mass),mass_table_shape=list(mass.shape),
        common_operator_shape=list(operator.matrix.shape),maximum_successive_L1=l1,
        maximum_successive_CDF=df,maximum_mass_residual=balance,maximum_negative_mass=negative,
        explicit_moving_event_increment=moving.tolist(),first_passage_flux=(final_f-cdf-moving).tolist(),
        numerical_mass_loss_added_to_event=False,probability_renormalized=False,
        probability_table_bytes=result.nbytes,linear_solve_checks=audits)


def advance_batch(operator,mass,cdf,duration):
    """Full-domain CRAM48 action: one step versus two half-steps, explicit F flux."""
    from cram_action import action
    n=len(operator.keep)
    Q=sparse.bmat([[operator.matrix,sparse.csr_matrix((n,1))],
        [sparse.csr_matrix(operator.event[None,:]),sparse.csr_matrix((1,1))]],format='csc')
    mask=np.ones(operator.common_states,dtype=bool);mask[operator.keep]=False
    moving=mass[:,mask].sum(axis=1)
    initial=np.vstack([mass[:,operator.keep].T,np.zeros((1,len(mass)))])
    coarse,cc=action(Q,initial,duration,1)
    fine,fc=action(Q,initial,duration,2)
    result=np.zeros_like(mass);result[:,operator.keep]=fine[:-1].T
    final_f=cdf+moving+fine[-1]
    l1=float(abs(fine-coarse).sum(axis=0).max());df=float(abs(fine[-1]-coarse[-1]).max())
    balance=float(abs(result.sum(axis=1)+final_f-1.).max())
    negative=float((-np.minimum(result,0.)).sum(axis=1).max())
    if not (l1<2e-12 and df<2e-15 and balance<1e-9 and negative<2e-12):
        raise ArithmeticError(dict(method='CRAM48',L1=l1,CDF=df,balance=balance,negative=negative))
    if not np.all(final_f>=cdf-2e-15):raise ArithmeticError('Nonmonotone first passage')
    return result,final_f,dict(method='full_sparse_CRAM48_successive_substeps',state_dimension=2,
        independent_cell_count=len(mass),mass_table_shape=list(mass.shape),
        common_operator_shape=list(operator.matrix.shape),maximum_successive_L1=l1,
        maximum_successive_CDF=df,maximum_mass_residual=balance,maximum_negative_mass=negative,
        explicit_moving_event_increment=moving.tolist(),first_passage_flux=fine[-1].tolist(),
        numerical_mass_loss_added_to_event=False,probability_renormalized=False,
        probability_table_bytes=result.nbytes,rational_action_checks=[cc,fc])
