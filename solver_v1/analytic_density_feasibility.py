"""Deterministic representability probes, not a production energy replacement.

A positive mixture of two exponential density kernels retains the same
Poisson sums. Densities are added BEFORE the per-atom embedding. This module
first tests observable compatibility; it does not admit a failed density
extension to the probability PDE or change the single-exponential default.
"""
from __future__ import annotations

from functools import lru_cache
import math

import numpy as np
from scipy.optimize import linprog, minimize, nnls

from .aluminum_full_fcc_calibration import FullFCCParameters, build_full_fcc_calibration_model
from .fcc111_active_interface import FCC111ActiveInterface
from .full_fcc_calibration_audit import IDEAL_H, bulk_basis
from .joint_fcc_interface_calibration import INTERFACE_STATES, interface_pair_basis, joint_basis


@lru_cache(maxsize=512)
def normalized_density_features(decay):
    """Bulk strain derivatives and per-depth interface density changes.

    The normalization is fixed at the reference bulk, NOT recomputed under
    strain or at each interface state. The resulting positive component
    kernel is exp(-k r)/rho_bulk(k); amplitude gauge cancels exactly.
    """
    basis = bulk_basis(float(decay))
    # Normal strain, hydrostatic strain, normal curvature, engineering shear.
    derivatives = np.array([basis[0, 3], 3*basis[0, 3], basis[2, 3],
                            basis[3, 3], basis[4, 3]])
    model = build_full_fcc_calibration_model(FullFCCParameters(.1, 1., float(decay), 1.))
    model.a0 = IDEAL_H
    rho_bulk = model.full_environment_density(IDEAL_H, 0).value
    changes = []
    for _, opening, slip, path in INTERFACE_STATES:
        interface = FCC111ActiveInterface(model, path_id=path, tolerance=2e-11)
        a = 400*IDEAL_H if math.isinf(opening) else opening*IDEAL_H
        depths, _, _, _ = interface._density_changes(a, slip)
        changes.append(np.array([row[0]/rho_bulk for row in depths]))
    return derivatives, tuple(changes)


def mixture_basis(k1, k2, weight):
    """x=(1-w)x_k1+w*x_k2, followed by F(x); never sum F componentwise."""
    if min(k1, k2) <= 0 or not 0 <= weight <= 1:
        raise ValueError("positive decay and mixture weight in [0,1] required")
    return _combined_basis((float(k1),float(k2)),(1-weight,weight))


def _combined_basis(decays,weights):
    features=[normalized_density_features(float(k)) for k in decays]
    xa,xeta,xetaeta,xaa,xgg=sum(w*feature[0] for w,feature in zip(weights,features))
    # At pristine cubic symmetry x_gamma=0. The normal variable here is
    # alpha, NOT a; factors h have already been applied by bulk_basis.
    sqrt_column = [-.5*xa, 1., .25*xeta*xeta-.5*xetaeta,
                   .25*xa*xa-.5*xaa, -.5*xgg]
    linear_column = [xa, -1., xetaeta, xaa, xgg]
    convex_column = [0., 1., 2*xeta*xeta, 2*xa*xa, 0.]
    rows = []
    for depths in zip(*(feature[1] for feature in features)):
        n=max(len(x) for x in depths)
        dx=sum(w*np.pad(x,(0,n-len(x))) for w,x in zip(weights,depths))
        if np.any(1+dx <= 0):
            raise ValueError("nonpositive mixed site density")
        rows.append([-2*np.sum(np.sqrt(1+dx)-1), 2*np.sum(dx), 2*np.sum(dx*dx)])
    top = np.column_stack((bulk_basis(float(decays[0]))[:, :2], sqrt_column,
                           linear_column, convex_column))
    return np.vstack((top, np.column_stack((interface_pair_basis(), rows))))


@lru_cache(maxsize=512)
def _bulk_density(decay):
    model=build_full_fcc_calibration_model(FullFCCParameters(.1,1.,float(decay),1.))
    return model.full_environment_density(IDEAL_H,0).value


def squared_envelope_basis(decay,shape):
    """Small non-monotone RESEARCH alternative after positive mixtures fail.

    q=exp[-k(r/L0-1)], f=q(1-shape*q)^2 >= 0 at every physical radius.
    The exact finite expansion q-2*shape*q^2+shape^2*q^3 uses the SAME infinite
    exponential Poisson sums at k,2k,3k. Signed coefficients must be added as
    DENSITY before F, not as embedding energies. Overall scale is gauge-fixed.
    This is a compatibility probe, not a claim of an Al electron orbital.
    """
    if not np.isfinite(decay) or decay<=0 or not np.isfinite(shape) or shape<0:
        raise ValueError("positive decay and nonnegative finite shape required")
    decays=(decay,2*decay,3*decay)
    weighted=np.array([1.,-2*shape,shape*shape])*np.array([_bulk_density(k) for k in decays])
    normalization=float(np.sum(weighted))
    if normalization<=0:
        raise ValueError("positive-kernel bulk density not resolved")
    return _combined_basis(decays,weighted/normalization)


def squared_envelope_profile(target,scales,*,grid=None,shapes=None):
    grid=np.geomspace(.5,4.,13) if grid is None else grid
    shapes=np.r_[0.,np.geomspace(.03,3.,17)] if shapes is None else shapes
    def fit(x):
        k,shape=x
        matrix=squared_envelope_basis(k,shape)
        coeff,_=nnls(matrix/scales[:,None],target/scales,maxiter=500)
        prediction=matrix@coeff; residual=(prediction-target)/scales
        return dict(decay=float(k),shape=float(shape),coefficients=coeff,
                    predictions=prediction,residuals=residual,
                    squared_loss=float(residual@residual))
    rows=[fit((k,shape)) for k in grid for shape in shapes]
    logs=[]
    for row in sorted(rows,key=lambda row:row["squared_loss"])[:3]:
        start=(row["decay"],row["shape"])
        result=minimize(lambda x:fit(x)["squared_loss"],start,method="Powell",
                        bounds=((min(grid),max(grid)),(min(shapes),max(shapes))),
                        options={"xtol":2e-5,"ftol":1e-8,"maxiter":60})
        rows.append(fit(result.x))
        logs.append(dict(start=start,success=bool(result.success),nfev=int(result.nfev),
                         message=str(result.message),squared_loss=float(result.fun)))
    return min(rows,key=lambda row:row["squared_loss"]),rows,logs


def mixture_fit(k1, k2, weight, target, scales, *, quadratic=True):
    count = 5 if quadratic else 4
    matrix = mixture_basis(k1, k2, weight)[:, :count]
    coeff, _ = nnls(matrix/scales[:, None], target/scales, maxiter=500)
    prediction = matrix@coeff
    residual = (prediction-target)/scales
    return dict(k1=float(k1), k2=float(k2), weight=float(weight),
                coefficients=np.r_[coeff, np.zeros(5-count)],
                squared_loss=float(residual@residual), residuals=residual,
                predictions=prediction)


def deterministic_mixture_search(target, scales, *, grid=None, weights=None):
    """Fixed grid followed by bounded deterministic local refinement.

    k1<=k2 removes the label-exchange gauge in the saved grid. The coincident
    rates and endpoint weights are the original family, not independent extra
    identifiable parameters. Bounds are research ranges, not material priors.
    """
    grid = np.geomspace(.35, 8., 13) if grid is None else np.asarray(grid)
    weights = (.1, .25, .5, .75, .9) if weights is None else weights
    rows = [mixture_fit(float(k1), float(k2), float(w), target, scales)
            for i, k1 in enumerate(grid) for k2 in grid[i:] for w in weights]
    best_starts = sorted(rows, key=lambda row: row["squared_loss"])[:3]
    local_records = []
    for start in best_starts:
        result = minimize(lambda x: mixture_fit(*x, target, scales)["squared_loss"],
                          [start["k1"], start["k2"], start["weight"]],
                          method="Powell", bounds=[(grid[0], grid[-1])]*2+[(0., 1.)],
                          options={"xtol": 3e-5, "ftol": 1e-8, "maxiter": 60})
        k1, k2, weight = result.x
        if k1 > k2:
            k1, k2, weight = k2, k1, 1-weight
        row = mixture_fit(k1, k2, weight, target, scales)
        rows.append(row)
        local_records.append(dict(start=[start["k1"], start["k2"], start["weight"]],
                                  success=bool(result.success), message=str(result.message),
                                  nfev=int(result.nfev), squared_loss=row["squared_loss"]))
    return min(rows, key=lambda row: row["squared_loss"]), rows, local_records


def intrinsic_fault_upper_bound(decay, target, scales, *, allowance=3., quadratic=True):
    """A linear-program certificate at ONE fixed decay, not a global theorem.

    Maximize the ISF energy subject to all other observations within a declared
    number of discrepancy scales. Nonnegative LJ/embedding coefficients are
    even allowed at zero: a relaxed feasible set makes the upper bound more
    generous than a strictly positive model. No residual or energy is clipped.
    """
    matrix = joint_basis(float(decay))[:, :5 if quadratic else 4]
    indices = np.array([i for i in range(len(target)) if i != 5])
    scaled = matrix[indices]/scales[indices, None]
    rhs = target[indices]/scales[indices]
    result = linprog(-matrix[5], A_ub=np.vstack((scaled, -scaled)),
                     b_ub=np.r_[rhs+allowance, -rhs+allowance],
                     bounds=(0, None), method="highs")
    upper = float(-result.fun) if result.success else None
    return dict(decay=float(decay), allowance=float(allowance), status=int(result.status),
                success=bool(result.success), message=result.message,
                intrinsic_upper_ev_cell=upper,
                intrinsic_reference_ev_cell=float(target[5]),
                coefficients=result.x.tolist() if result.success else None)
