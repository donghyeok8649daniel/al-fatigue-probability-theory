"""Analytic Newton stationary-point audit; never a physical time evolution.

An indefinite Hessian is allowed because the sought point may be a saddle.
Residual line search changes only numerical steps, not energy/forces. Final
force and the full finite-domain Morse index are reported independently.
"""
import numpy as np
from scipy.linalg import solve, eigh


def stationary_core(core, initial, *, force_tolerance=5e-8, max_iterations=30,
                    max_step=.2, callback=None):
    if min(force_tolerance,max_step)<=0 or not np.all(np.isfinite([force_tolerance,max_step])) or max_iterations<1:
        raise ValueError('finite positive stationary-solver controls required')
    field = np.asarray(initial,float).copy()
    history = []
    for iteration in range(max_iterations+1):
        value, action = core.linearize(field)
        g = value['gradient'].ravel()
        matrix = action.explicit_matrix()
        residual = float(np.max(abs(g)))
        asymmetry = float(np.max(abs(matrix-matrix.T)))
        if asymmetry>2e-11*max(1.,float(np.max(abs(matrix)))):
            raise ArithmeticError('analytic core matrix failed symmetry')
        history.append(dict(iteration=iteration,energy=value['energy'],maximum_force=residual))
        if callback is not None:
            callback(iteration,field,history[-1])
        if residual<=force_tolerance or iteration==max_iterations:
            break
        direction = solve(matrix,-g,assume_a='sym',check_finite=True).reshape(field.shape)
        scale = min(1.,max_step/float(np.max(np.linalg.norm(direction,axis=1))))
        accepted = False
        for step in scale*2.**-np.arange(16):
            trial = field+step*direction
            evaluated = core.evaluate(trial)
            if np.linalg.norm(evaluated['gradient']) < np.linalg.norm(g):
                field = trial; accepted = True; break
        if not accepted:
            break
    value, action = core.linearize(field)
    matrix = action.explicit_matrix()
    eigenvalues,eigenvectors = eigh(matrix,check_finite=True)
    vector = eigenvectors[:,0]
    vector *= 1 if vector[np.argmax(abs(vector))]>=0 else -1
    error = float(np.linalg.norm(action(vector.reshape(field.shape)).ravel()-eigenvalues[0]*vector))
    roundoff = len(matrix)*np.finfo(float).eps*np.linalg.norm(matrix,ord=np.inf)
    floor = max(error,roundoff)
    force = float(np.max(abs(value['gradient'])))
    return dict(field=field,energy=value['energy'],maximum_force=force,
        force_converged=force<=force_tolerance,history=history,
        eigenvalues=eigenvalues,minimum_mode=vector.reshape(field.shape),
        eigenpair_residual=error,hessian_resolution_floor=floor,
        negative_eigenvalues=int(np.sum(eigenvalues < -floor)),
        unresolved_eigenvalues=int(np.sum(abs(eigenvalues)<=floor)),
        physical_dynamics=False)
