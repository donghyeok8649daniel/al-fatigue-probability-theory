"""Static Morse audit; small forces alone do not certify a stable core."""
import numpy as np
from scipy.sparse.linalg import LinearOperator,eigsh


def lowest_core_mode(core,field,*,check_energy=True,force_tolerance=2e-6,method='operator'):
    field=np.asarray(field,float); out,apply=core.linearize(field); n=field.size
    operator=LinearOperator((n,n),matvec=lambda v:apply(v.reshape(field.shape)).ravel(),dtype=float)
    if method=='operator':
        values,vectors=eigsh(operator,k=1,which='SA',tol=1e-8,maxiter=4000,
                             v0=np.cos(np.arange(n)*.7))
    elif method=='assembled':
        from scipy.linalg import eigh
        matrix=apply.explicit_matrix()
        asymmetry=float(np.max(abs(matrix-matrix.T)))
        if asymmetry>2e-11*max(1.,float(np.max(abs(matrix)))):
            raise ArithmeticError('assembled analytic Hessian is not symmetric')
        values,vectors=eigh(matrix,subset_by_index=[0,0],check_finite=True)
    else:
        raise ValueError('explicit operator or assembled Morse method required')
    vector=vectors[:,0]; index=int(np.argmax(abs(vector)))
    vector*=1 if vector[index]>=0 else -1
    residual=float(np.linalg.norm(operator@vector-values[0]*vector))
    mode=vector.reshape(field.shape); checks=[]
    if check_energy:
        for step in (2e-3,1e-3):
            plus=core.evaluate(field+step*mode)['energy']
            minus=core.evaluate(field-step*mode)['energy']
            curvature=(plus+minus-2*out['energy'])/step**2
            checks.append(dict(step_over_L0=step,curvature=curvature,
                eigenvalue_error=abs(curvature-values[0]),
                plus_energy_change=plus-out['energy'],minus_energy_change=minus-out['energy']))
    force=float(np.max(abs(out['gradient'])))
    return dict(minimum_eigenvalue=float(values[0]),eigenvector=mode,
                hessian_method=method,
                eigenpair_residual=residual,energy_curvature_checks=checks,
                maximum_free_force=force,force_tolerance=force_tolerance,
                positive_hessian_on_tested_fixed_boundary=bool(values[0]>residual),
                stable_on_tested_fixed_boundary=bool(values[0]>residual and force<=force_tolerance))
