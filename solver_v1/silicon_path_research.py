"""Nudged elastic-band research search; optimization springs are not physics."""
from __future__ import annotations

import numpy as np


def reparameterize_path(points,images):
    path=np.asarray(points,float)
    if path.ndim!=2 or len(path)<2 or images<3 or not np.isfinite(path).all():
        raise ValueError('finite path and at least three images required')
    distance=np.r_[0.,np.cumsum(np.linalg.norm(np.diff(path,axis=0),axis=1))]
    keep=np.r_[True,np.diff(distance)>1e-14]
    if distance[-1]<=0:
        raise ValueError('distinct endpoints/path points required')
    target=np.linspace(0,distance[-1],images)
    return np.column_stack([np.interp(target,distance[keep],path[keep,j]) for j in range(path.shape[1])])


def neb_forces(images,energies,gradients,*,spring,climb):
    """Improved energy-weighted tangent and optional highest climbing image."""
    x,e,g=np.asarray(images,float),np.asarray(energies,float),np.asarray(gradients,float)
    if (x.ndim!=2 or x.shape!=g.shape or e.shape!=(len(x),) or len(x)<3
            or not all(np.isfinite(a).all() for a in (x,e,g)) or spring<=0):
        raise ValueError('finite band energies/gradients and positive numerical spring required')
    forces=np.zeros_like(x)
    highest=1+int(np.argmax(e[1:-1]))
    for i in range(1,len(x)-1):
        forward,backward=x[i+1]-x[i],x[i]-x[i-1]
        up,down=abs(e[i+1]-e[i]),abs(e[i-1]-e[i])
        if e[i+1]>e[i]>e[i-1]:
            tangent=forward
        elif e[i+1]<e[i]<e[i-1]:
            tangent=backward
        elif e[i+1]>e[i-1]:
            tangent=forward*max(up,down)+backward*min(up,down)
        else:
            tangent=forward*min(up,down)+backward*max(up,down)
        if np.linalg.norm(tangent)<1e-14:
            tangent=forward+backward
        norm=np.linalg.norm(tangent)
        if norm<1e-14:
            raise ValueError('degenerate band tangent')
        tangent=tangent/norm
        parallel=float(g[i]@tangent)
        if climb and i==highest:
            forces[i]=-g[i]+2*parallel*tangent
        else:
            forces[i]=-g[i]+parallel*tangent+spring*(np.linalg.norm(forward)-np.linalg.norm(backward))*tangent
    return forces,highest


def relax_neb(evaluate,initial,*,spring=1.,max_steps=2000,climb_after=200,
              tolerance=2e-4,step=.02,max_step=.15,max_image_move=.08,
              climb_force_threshold=.03,progress=None):
    """FIRE relaxation of the band, preserving both endpoints exactly.

    evaluate(index,x) returns the ACTUAL potential energy and gradient. Numerical
    spring, fictitious FIRE time and displacement limiter are optimizer settings
    and are never included in reported material energies or physical dynamics.
    """
    x=np.asarray(initial,float).copy();velocity=np.zeros_like(x)
    endpoints=x[[0,-1]].copy()
    dt=float(step);alpha=.1;positive=0;history=[];climbing=False
    if max_steps<=climb_after or tolerance<=0 or step<=0 or max_step<step or max_image_move<=0:
        raise ValueError('valid band optimization controls required')
    for iteration in range(max_steps+1):
        evaluated=[evaluate(i,point) for i,point in enumerate(x)]
        energies=np.array([a[0] for a in evaluated]);gradients=np.array([a[1] for a in evaluated])
        started_climbing=False
        if not climbing and iteration>=climb_after:
            ordinary,_=neb_forces(x,energies,gradients,spring=spring,climb=False)
            if np.max(abs(ordinary[1:-1]))<=climb_force_threshold:
                climbing=True;started_climbing=True
        force,highest=neb_forces(x,energies,gradients,spring=spring,climb=climbing)
        residual=float(np.max(abs(force[1:-1])))
        record=dict(iteration=iteration,highest_image=highest,highest_energy=float(energies[highest]),
            neb_force_max=residual,climbing=climbing,fire_step=dt)
        if iteration%25==0 or residual<tolerance or iteration==max_steps:
            history.append(record)
            if progress is not None: progress(record,x,energies)
        if climbing and residual<tolerance:
            break
        if iteration==max_steps:
            break
        if started_climbing:
            velocity[:]=0;dt=step;alpha=.1;positive=0
        velocity+=dt*force
        power=float(np.sum(velocity*force))
        if power>0:
            positive+=1
            fnorm=np.linalg.norm(force);vnorm=np.linalg.norm(velocity)
            if fnorm:
                velocity=(1-alpha)*velocity+alpha*force*(vnorm/fnorm)
            if positive>5:
                dt=min(dt*1.1,max_step);alpha*=.99
        else:
            velocity[:]=0;dt*=.5;alpha=.1;positive=0
        displacement=dt*velocity
        largest=float(np.max(np.linalg.norm(displacement[1:-1],axis=1)))
        if largest>max_image_move:
            displacement*=max_image_move/largest
        x[1:-1]+=displacement[1:-1]
        np.testing.assert_array_equal(x[[0,-1]],endpoints)
    return dict(images=x,energies=energies,gradients=gradients,history=history,
        converged=bool(climbing and residual<tolerance),iterations=iteration,
        highest_image=highest,neb_force_max=residual)
