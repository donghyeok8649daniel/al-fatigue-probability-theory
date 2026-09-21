"""One explicit collective displacement retained beside a physical crack gap."""
from __future__ import annotations

import numpy as np


class RetainedBathDirection:
    """Householder null basis with O(n) application, no dense coordinate map.

    base is an orthonormal bath at fixed gap q. In its bath x, retain u=w.x,
    w.w=1. x=H[:,keep] y+w u. H is an orthogonal Householder reflection.
    Both the retained coordinate and y have Cartesian length units. A plotted
    fraction u/path_length is a label, not a change of energy or thermal scale.
    """
    def __init__(self,base,direction):
        vector=np.asarray(direction,float)
        if vector.shape!=(base.dimension,) or not np.isfinite(vector).all() or np.linalg.norm(vector)<=0:
            raise ValueError('nonzero finite direction in the base bath required')
        self.base=base
        self.path_length=float(np.linalg.norm(vector))
        self.direction=vector/self.path_length
        self.pivot=int(np.argmax(abs(self.direction)))
        reflection=self.direction.copy()
        reflection[self.pivot]+=np.copysign(1.,reflection[self.pivot])
        self.reflection=reflection/np.linalg.norm(reflection)
        self.keep=np.arange(base.dimension)!=self.pivot
        self.dimension=base.dimension-1

    def reflect(self,vector):
        x=np.asarray(vector,float)
        return x-2*self.reflection*(self.reflection@x)

    def lift(self,variables,retained):
        y=np.asarray(variables,float)
        if y.shape!=(self.dimension,) or not np.isfinite(y).all() or not np.isfinite(retained):
            raise ValueError('finite bath variables and retained length required')
        x=np.zeros(self.base.dimension);x[self.keep]=y
        return self.reflect(x)+self.direction*retained

    def positions(self,variables,values):
        q,retained=values
        return self.base.positions(self.lift(variables,retained),[q])

    def encode(self,positions):
        return self.reflect(self.base.encode(positions))[self.keep]

    def pullback(self,gradient):
        bath,reactions=self.base.pullback(gradient)
        return self.reflect(bath)[self.keep],np.r_[reactions,self.direction@bath]

    def restricted_hessian(self,base_bath_hessian):
        a=np.asarray(base_bath_hessian,float)
        v=self.reflection;av=a@v
        rotated=a-2*np.outer(v,av)-2*np.outer(av,v)+4*(v@av)*np.outer(v,v)
        return rotated[np.ix_(self.keep,self.keep)]
