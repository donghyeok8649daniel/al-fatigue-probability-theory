"""Positive Poisson-series action of the identical finite SG generator.

For nu>=max(-diag K), exp(hQ)=exp(-nu*h) sum_k (nu*h)^k/k! (I+Q/nu)^k.
Centered summation avoids underflow. Chernoff bounds control omitted tails.
Signed roundoff in the initial vector is retained, never clipped into F/U.
"""
import math,time,numpy as np
from scipy.special import gammaln

def stable_poisson_weights(mu,low,high):
    # Avoid cancellation of O(mu*log(mu)) terms in the central log PMF.
    mode=int(math.floor(mu));weights=np.empty(high-low+1)
    if mode<16:
        value=math.exp(-mu)
        for k in range(high+1):
            if k>=low:weights[k-low]=value
            value*=mu/(k+1)
        return weights
    m=float(mode);v=1/m;v2=v*v
    stirling=v*(1/12+v2*(-1/360+v2*(1/1260+v2*(-1/1680+v2*(1/1188+v2*(-691/360360+v2*(7/1092-v2*3617/122400)))))))
    delta=(mu-m)/m
    logpm=m*(math.log1p(delta)-delta)-.5*math.log(2*math.pi*m)-stirling
    weights[mode-low]=math.exp(logpm)
    for k in range(mode,low,-1):weights[k-1-low]=weights[k-low]*k/mu
    for k in range(mode,high):weights[k+1-low]=weights[k-low]*mu/(k+1)
    return weights

def poisson_band(mu,kappa,tol):
    m=mu*kappa;s=math.sqrt(2*math.log(16/tol));width=s*math.sqrt(m)+s*s
    for _ in range(20):
        low=max(0,int(math.floor(m-width)));high=int(math.ceil(m+width))
        def logbound(k):
            if k<=0:return -m
            return -m+k-k*math.log(k/m)
        lower=0. if low==0 else math.exp(min(0.,logbound(low-1)+1e-8))
        upper=math.exp(min(0.,logbound(high+1)+1e-8))
        bound=math.exp(mu*(kappa-1))*(lower+upper)
        if bound<=tol:return low,high,bound
        width=width*1.2+1
    raise RuntimeError('Poisson tail bracket failed')

def propagate(K,p,F,U,r,ru,h,threads=4,progress=None,fine_tol=1e-20,peak_out=None,physical_offset=None):
    import torch
    torch.set_num_threads(threads);start=time.perf_counter();n=len(p)
    nu=1.01*float(np.max(-K.diagonal()));assert nu>0 and np.isfinite(nu)
    # The on-disk K is unchanged. Its loaded buffer becomes I+K/nu.
    K.data/=nu;K.setdiag(K.diagonal()+1.);assert K.data.min()>=0
    colsum=np.asarray(K.sum(axis=0)).ravel()+(r+ru)/nu
    kappa=max(1.,float(colsum.max())+32*np.finfo(float).eps)
    mu=nu*h;coarse=poisson_band(mu,kappa,1e-16);fine=poisson_band(mu,kappa,fine_tol)
    if fine[1]>200000:raise RuntimeError(f'Uniformization cap: {fine[1]} powers')
    weights=stable_poisson_weights(mu,fine[0],fine[1])
    P=torch.sparse_csr_tensor(torch.from_numpy(K.indptr),torch.from_numpy(K.indices),torch.from_numpy(K.data),size=K.shape,dtype=torch.float64,check_invariants=True)
    pk=torch.from_numpy(np.ascontiguousarray(p).copy());work=torch.empty_like(pk);yf=torch.zeros_like(pk);yc=torch.zeros_like(pk)
    peak=pk.clone() if peak_out is not None else None
    rf=torch.from_numpy(np.ascontiguousarray(r/nu)) if np.any(r) else None
    uf=torch.from_numpy(np.ascontiguousarray(ru/nu)) if np.any(ru) else None
    Fk=Uk=0.;Ffine=Ufine=Fcoarse=Ucoarse=0.
    with torch.inference_mode():
        for k in range(fine[1]+1):
            if k>=fine[0]:
                w=float(weights[k-fine[0]]);yf.add_(pk,alpha=w);Ffine+=w*Fk;Ufine+=w*Uk
                if coarse[0]<=k<=coarse[1]:yc.add_(pk,alpha=w);Fcoarse+=w*Fk;Ucoarse+=w*Uk
            if k==fine[1]:break
            if rf is not None:Fk+=torch.dot(rf,pk).item()
            if uf is not None:Uk+=torch.dot(uf,pk).item()
            torch.mv(P,pk,out=work);pk,work=work,pk
            if peak is not None:torch.maximum(peak,pk,out=peak)
            if progress and (k+1)%250==0:progress(dict(power=k+1,total_powers=fine[1],nu=nu,mu=mu,seconds=time.perf_counter()-start))
    a=yf.numpy().copy();b=yc.numpy();y=np.r_[a,F+Ffine,U+Ufine]
    d=float(abs(a-b).sum()+abs(Ffine-Fcoarse)+abs(Ufine-Ucoarse));df=abs(Ffine-Fcoarse);du=abs(Ufine-Ucoarse)
    physical=y if physical_offset is None else physical_offset+y
    neg=float(-physical[physical<0].sum());mass=abs(float(y.sum())-(float(p.sum())+F+U)) if physical_offset is None else abs(float(physical.sum())-1.)
    check=dict(method='positive_uniformization',order=fine[1],successive_L1=d,successive_CDF=df,successive_unresolved=du,negative_mass=neg,mass_residual=mass,nu=nu,poisson_mean=mu,coarse_band=list(coarse),fine_band=list(fine),poisson_tail_L1_bound=fine[2]*float(abs(p).sum()),column_norm_upper=kappa,negative_input_mass=float(-p[p<0].sum()),seconds=time.perf_counter()-start,physical_offset_used=physical_offset is not None)
    if not (d<2e-12 and df<2e-15 and du<2e-16 and neg<2e-12 and mass<1e-9):raise RuntimeError(str(check))
    if peak_out is not None:peak_out.append(peak.numpy().copy())
    return y,[check]

if __name__=='__main__':
    from pathlib import Path
    import sys,json
    from scipy import sparse
    from scipy.sparse.linalg import expm_multiply
    here=Path(__file__).parent;case=here/'stencil_small_new_pad1';f=-4.5;h=.3
    with np.load(case/'initial.npz') as z:p0=z['mass']
    with np.load(case/f'phase_{f:+.1f}.npz') as z:keep=z['keep'];r=z['event'];ru=z['unresolved']
    K=sparse.load_npz(case/f'K_{f:+.1f}.npz');p=np.ones(K.shape[0]);p/=p.sum();z=np.r_[p,0.,0.]
    # Explicit zero scalar blocks fix the augmented matrix shape.
    Q=sparse.bmat([[K,sparse.csr_matrix((len(p),2))],[sparse.csr_matrix(np.vstack((r,ru))),sparse.csr_matrix((2,2))]],format='csr')
    ref=expm_multiply(h*Q,z,traceA=h*float(K.diagonal().sum()))
    y,checks=propagate(K,p,0.,0.,r,ru,h,threads=1)
    report=dict(L1_difference=float(abs(y-ref).sum()),CDF_difference=float(abs(y[-2]-ref[-2])),unresolved_difference=float(abs(y[-1]-ref[-1])),checks=checks)
    report['accepted']=report['L1_difference']<2e-12 and report['CDF_difference']<2e-15 and report['unresolved_difference']<2e-15
    (here/'UNIFORMIZATION_CHECK.json').write_text(json.dumps(report,indent=2));print(report);assert report['accepted']
