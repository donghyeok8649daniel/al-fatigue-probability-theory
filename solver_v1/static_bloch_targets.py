"""Independent spatial force-constant targets, no mass, frequency or mobility.

The existing fit knew cubic elastic constants and sampled finite-q positivity,
but not finite-q stiffness magnitudes. Seven predeclared source projections
now constrain those magnitudes; seven different q projections remain held out.
This adds observations, NOT a new potential term or fitted fatigue quantity.
"""
from dataclasses import dataclass

import numpy as np

from .periodic_plane_covariance import SourceEAMBlochHessian
from .vector_material_calibration import LENGTH_M,MaterialObservation


@dataclass(frozen=True)
class BlochTarget:
    observation: MaterialObservation
    fractional_cubic: tuple
    polarization_cubic: tuple


def source_bloch_targets(reference):
    operator=SourceEAMBlochHessian(reference)
    if not np.isclose(reference.geometry.lattice_constant,4.05,rtol=0,atol=1e-12):
        raise ValueError('0K static 4.05 Angstrom source state required; do not mix MD thermal box')
    plans=[('fit',(.25,.25,.25)),('fit',(.5,0.,0.)),('fit',(.4,.4,0.)),
           ('heldout',(.5,.5,.5)),('heldout',(.9,0.,0.)),('heldout',(.6,.6,0.))]
    targets=[];basis=reference.geometry.plane_basis_in_stacked_cubic_axes()
    for role,point in plans:
        q=np.asarray(point);long=q/np.linalg.norm(q)
        trans=np.array([1.,-1.,0.]) if q[1] else np.array([0.,1.,0.])
        trans/=np.linalg.norm(trans)
        if abs(trans@long)>1e-12:raise ArithmeticError('nonorthogonal high-symmetry polarization')
        vectors=[('L',long),('T1',trans)]
        if q[1] and not q[2]:vectors.append(('T2',np.cross(long,trans)))
        H=operator.evaluate(q*2*np.pi/reference.geometry.lattice_constant)*(LENGTH_M/1e-10)**2
        for label,direction in vectors:
            v=basis@direction;value=float(v@H@v)
            name='bloch_'+'_'.join(f'{x:g}' for x in point)+'_'+label
            if value<=0:raise ArithmeticError('positive source restoring stiffness required')
            targets.append(BlochTarget(MaterialObservation(name,value,.10*value,
                'eV/L0^2',role),tuple(point),tuple(direction)))
    return tuple(targets)


def append_bloch_problem(matrix,observations,basis,targets):
    columns=[];errors=[]
    rotation=basis.geometry.plane_basis_in_stacked_cubic_axes()
    cache={}
    for target in targets:
        q=target.fractional_cubic
        if q not in cache:cache[q]=basis.evaluate(q)
        H,bound=cache[q];v=rotation@np.asarray(target.polarization_cubic)
        columns.append(np.einsum('i,cij,j->c',v,H,v));errors.append(bound)
    return np.vstack([matrix,columns]),list(observations)+[t.observation for t in targets],np.asarray(errors)
