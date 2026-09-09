"""Tail-aware static material profiling; no change to the energy or kinetics.

The perfect FCC Voronoi cells have volume b**3/sqrt(2) and covering radius
b/sqrt(2). For a decreasing radial envelope f, omitted lattice sites obey

 sum_(|R|>r) f(|R|) <= 4pi/V integral_(r-2c)^infinity (u+c)^2 f(u) du.

This integrates over a superset of the omitted Voronoi cells, after using
|R| >= |x|-c. It bounds the actual discrete infinite tail, not just its last
observed difference. Power and exponential-polynomial integrals are analytic.
The independently evaluated static Hessian is coefficient-linear. Its spectral
tail is bounded column by column, including squares of collective variations.
H_R(q)-sum_j |c_j| error_j(q) I >= 0 is sufficient at the tested q. It is
NOT a proof over untested wavevectors or away from perfect bulk equilibrium.
"""
from functools import lru_cache
from itertools import product
import math

import numpy as np
from scipy.linalg import null_space
from scipy.optimize import LinearConstraint, minimize, nnls
from scipy.special import gamma, gammaincc

from .angular_environment_reference import traceless_second, traceless_third
from .range_resolved_material import build_range_surface
from .static_bulk_stability import neighbors_in_plane_frame


class FCCTailEnvelope:
    def __init__(self, radius, b=1.):
        self.radius, self.b = float(radius), float(b)
        self.cover = b/math.sqrt(2)
        self.lower = radius-2*self.cover
        self.factor = 4*math.pi/(b**3/math.sqrt(2))
        if not np.isfinite(radius+b) or b <= 0 or self.lower <= 0:
            raise ValueError('positive FCC spacing and radius > two covering radii required')

    def power(self, exponent):
        if exponent <= 3:
            raise ValueError('integrable decreasing three-dimensional power required')
        u, c = self.lower, self.cover
        return self.factor*(u**(3-exponent)/(exponent-3)
            + 2*c*u**(2-exponent)/(exponent-2)+c*c*u**(1-exponent)/(exponent-1))

    def exponential(self, k, degree):
        """Bound sum r**degree exp(-k r); integer degree >=0."""
        if k <= 0 or int(degree) != degree or degree < 0 or k*self.lower <= degree:
            raise ValueError('envelope must be decreasing on the enlarged tail domain')
        def integral(n):
            return gamma(n+1)*gammaincc(n+1, k*self.lower)/k**(n+1)
        return self.factor*(integral(degree+2)+2*self.cover*integral(degree+1)
                            +self.cover**2*integral(degree))


@lru_cache(maxsize=128)
def normalized_amplitude(k):
    probe = build_range_surface(k, k, k, np.ones(8))
    return math.exp(k)/probe.face.bulk.embedding.rho_ref


class TailBulkCoefficientBasis:
    """Direct coefficient operator plus analytic omitted-neighbor envelope.

    No 9-times repeated model subtraction; columns are the exact differentiated
    LJ/embedding/STF expressions. Independent legacy-operator tests are required.
    The canonical state energy remains the existing infinite Poisson/Bessel sum.
    """
    def __init__(self, decays, *, radius=12.):
        self.decays = tuple(map(float, decays))
        if len(self.decays) != 3 or not np.all(np.isfinite(self.decays)) or min(self.decays) <= 0:
            raise ValueError('three positive finite radial decays required')
        self.radius = float(radius)
        bulk = build_range_surface(*decays, np.ones(8)).face.bulk
        self.geometry = bulk.geometry
        self.R = neighbors_in_plane_frame(bulk, radius)
        self.r = np.linalg.norm(self.R, axis=1)
        self.e = self.R/self.r[:, None]
        self.rr = self.e[:, :, None]*self.e[:, None, :]
        self.tail = FCCTailEnvelope(radius, bulk.geometry.b)
        self.pair = [self.r[:, None, None]**-14*(168*self.rr-12*np.eye(3)),
                     self.r[:, None, None]**-8*(6*np.eye(3)-48*self.rr)]
        k = self.decays[0]
        self.rho_weight = normalized_amplitude(k)*np.exp(-k*self.r)
        self.rho_gradient = -k*self.rho_weight[:, None]*self.e
        self.rho_hessian = self.rho_weight[:, None, None]*(k*k*self.rr
            -k/self.r[:, None, None]*(np.eye(3)-self.rr))
        self.moments = [(3, self.decays[1]), (1, self.decays[1]), (2, self.decays[2])]
        self.gradients = [self._moment_gradient(rank, k) for rank, k in self.moments]

    def _moment_gradient(self, rank, k):
        r = self.R; n = len(r); weight = normalized_amplitude(k)*np.exp(-k*self.r)
        raw = np.empty((n, 3)+ (3,)*rank)
        for b in range(3):
            for indices in product(range(3), repeat=rank):
                polynomial = np.prod(r[:, indices], axis=1)
                derivative = -k*self.e[:, b]*polynomial
                for position, index in enumerate(indices):
                    if index == b:
                        rest = indices[:position]+indices[position+1:]
                        derivative += np.prod(r[:, rest], axis=1) if rest else 1.
                raw[(slice(None), b)+indices] = weight*derivative
        projected = {1: lambda x:x, 2:traceless_second, 3:traceless_third}[rank](raw)
        return projected.reshape(n, 3, -1)

    def wavevector(self, fractional_cubic):
        return (self.geometry.plane_basis_in_stacked_cubic_axes()
            @(2*math.pi/self.geometry.lattice_constant*np.asarray(fractional_cubic, float)))

    def evaluate(self, fractional_cubic):
        q = self.wavevector(fractional_cubic); qn = np.linalg.norm(q)
        phase = self.R@q; sine = np.sin(phase)
        # Accurate acoustic limit: 1-cos(x)=2sin(x/2)^2.
        one_minus_cos = 2*np.sin(phase/2)**2
        columns = [np.einsum('n,nij->ij', one_minus_cos, h) for h in self.pair]
        rh = np.einsum('n,nij->ij', one_minus_cos, self.rho_hessian)
        rg = sine@self.rho_gradient
        columns.extend([-rh+.25*np.outer(rg, rg), 2*rh, 2*np.outer(rg, rg)])
        L = []
        for (rank, _), gradient in zip(self.moments, self.gradients):
            variation = np.einsum('n,nij->ij', sine if rank % 2 == 0 else -one_minus_cos, gradient)
            L.append(variation)
            columns.append(2*variation@variation.T)
        t = self.tail
        def cosine_power(p):
            return min(2*t.power(p), .5*qn*qn*t.power(p-2))
        errors = [156*cosine_power(14), 42*cosine_power(8)]
        k = self.decays[0]; amplitude = normalized_amplitude(k)
        # k/r <= k/lower on the enlarged integration domain.
        dh = amplitude*(k*k+k/t.lower)*min(2*t.exponential(k, 0),
                                                  .5*qn*qn*t.exponential(k, 2))
        dg = amplitude*k*min(t.exponential(k, 0), qn*t.exponential(k, 1))
        square_error = 2*np.linalg.norm(rg)*dg+dg*dg
        errors.extend([dh+.25*square_error, 2*dh, 2*square_error])
        for (rank, k), variation in zip(self.moments, L):
            amplitude = normalized_amplitude(k)
            def envelope(power):
                return amplitude*(rank*t.exponential(k, rank-1+power)
                                   +k*t.exponential(k, rank+power))
            dl = (min(envelope(0), qn*envelope(1)) if rank % 2 == 0 else
                  min(2*envelope(0), .5*qn*qn*envelope(2)))
            errors.append(2*(2*np.linalg.norm(variation, 2)*dl+dl*dl))
        return np.array(columns), np.array(errors)


def project_affine_face(matrix, rhs, point):
    """Euclidean projection onto a consistent affine face without A A^T.

    Nearly parallel spectral polarizations can be independently resolvable
    in A but vanish in the squared-condition Gram matrix. Solve for a
    particular point and its orthogonal null space directly. The caller
    independently checks primal/dual/KKT, including an inconsistent face.
    """
    A=np.asarray(matrix,float);h=np.asarray(rhs,float);x=np.asarray(point,float)
    if (A.ndim!=2 or h.shape!=(len(A),) or x.shape!=(A.shape[1],)
            or any(np.any(~np.isfinite(v)) for v in (A,h,x))):
        raise ValueError('finite affine constraint matrix, right side and point required')
    particular=np.linalg.lstsq(A,h,rcond=None)[0]
    Z=null_space(A)
    return particular+Z@(Z.T@(x-particular))


def convex_profile(matrix, observations, inequalities=None, *, nonnegative=(0,1,2,4,5)):
    """Equality-eliminated, whitened convex LS with verified linear inequalities.

    SVD turns the positive-definite objective into Euclidean projection. Active
    face polishing then supplies independent primal/dual/KKT residuals. This is
    not a radial global optimum claim. A failed numerical profile raises.
    """
    matrix = np.asarray(matrix, float); p = matrix.shape[1]
    if matrix.shape[0] != len(observations) or np.any(~np.isfinite(matrix)):
        raise ValueError('finite matching observation matrix required')
    exact = [i for i,o in enumerate(observations) if o.role == 'exact']
    selected = [i for i,o in enumerate(observations) if o.role == 'fit']
    target = np.array([o.target for o in observations]); scales = np.array([o.scale for o in observations])
    eq = matrix[exact]/scales[exact, None]; rhs = target[exact]/scales[exact]
    offset = np.linalg.lstsq(eq, rhs, rcond=None)[0]; T = null_space(eq)
    design = matrix[selected]@T/scales[selected, None]
    U, sv, Vt = np.linalg.svd(design, full_matrices=False)
    if sv[-1] <= np.finfo(float).eps*max(design.shape)*sv[0]:
        raise ValueError('unidentifiable coefficient directions: fix gauge before profiling')
    transform = T@(Vt.T/sv)
    x0 = U.T@((target[selected]-matrix[selected]@offset)/scales[selected])
    extra = np.empty((0,p)) if inequalities is None else np.asarray(inequalities,float).reshape(-1,p)
    G = np.vstack([np.eye(p)[list(nonnegative)], extra])
    A = G@transform; h = -G@offset; norms = np.linalg.norm(A,axis=1)
    row_ids=np.arange(len(G))
    if np.any(norms <= 1e-15):
        redundant = norms <= 1e-15
        if np.any(h[redundant] > 1e-10): raise ArithmeticError('infeasible constant constraint')
        A,h,norms,row_ids = A[~redundant],h[~redundant],norms[~redundant],row_ids[~redundant]
    A /= norms[:,None]; h /= norms
    run = minimize(lambda x:.5*np.sum((x-x0)**2),x0,
        jac=lambda x:x-x0, method='SLSQP', constraints=[LinearConstraint(A,h,np.inf)],
        options=dict(ftol=2e-12,maxiter=350))
    x = run.x; slack = A@x-h
    active = np.flatnonzero(slack <= 1e-7*max(1.,np.linalg.norm(x)))
    if len(active):
        # Redundant sign/polarization cuts can have nonunique duals. An
        # unrestricted minimum-norm dual may be negative even though a valid
        # nonnegative dual exists; solve the actual dual-cone membership.
        multipliers = nnls(A[active].T,x-x0,maxiter=max(500,20*len(active)))[0]
        support=active[multipliers>1e-10*max(1.,np.max(multipliers))]
        if len(support):
            polished=project_affine_face(A[support],h[support],x0)
            dual=nnls(A[support].T,polished-x0,maxiter=max(500,20*len(support)))[0]
            # NNLS returning nonnegative duals alone is insufficient: an
            # over-screened equality face can have a nonzero dual residual.
            # Never replace the raw point with such a feasible non-optimum.
            face_kkt=np.linalg.norm(polished-x0-A[support].T@dual,np.inf)
            if np.min(A@polished-h)>=-1e-8 and face_kkt<=1e-6:
                x=polished;active=support;multipliers=dual
    else:
        multipliers = np.empty(0)
    kkt = float(np.linalg.norm(x-x0-A[active].T@multipliers,np.inf))
    violation = max(0.,float(-np.min(A@x-h)))
    if violation > 1e-7 or kkt > 1e-6 or (len(multipliers) and np.min(multipliers)<-1e-7):
        raise ArithmeticError(f'QP not verified: {run.message}; primal={violation}; KKT={kkt}')
    c = offset+transform@x
    # Re-solve the accepted active face in physical coefficient coordinates.
    # In particular active nonnegative amplitudes are coordinates fixed to
    # zero, not a cancellation of two large whitened numbers followed by clip.
    # A near-active inequality with ZERO dual weight is not an equality of
    # the optimal face. Forcing all slack-tolerance neighbors to equality
    # can produce a feasible but nonstationary, much worse coefficient point.
    active_raw=row_ids[active[np.asarray(multipliers)>0]]
    zero=[nonnegative[i] for i in active_raw if i<len(nonnegative)]
    free=[i for i in range(p) if i not in zero]
    other=[i for i in active_raw if i>=len(nonnegative)]
    # Coefficient units can differ by many orders (quartic moment amplitudes).
    # Equilibrate columns for this *numerical solve*, preserving the exact
    # physical coefficients, objective and constraints on reconstruction.
    column_norm=np.linalg.norm(np.vstack([eq[:,free],matrix[np.ix_(selected,free)]/scales[selected,None]]),axis=0)
    if np.any(column_norm==0):raise ValueError('unidentified free coefficient in active face')
    column_scale=1/column_norm
    face_eq=np.vstack([eq[:,free],G[np.ix_(other,free)]])*column_scale
    face_rhs=np.r_[rhs,np.zeros(len(other))]
    eqnorm=np.linalg.norm(face_eq,axis=1)
    face_eq=face_eq/eqnorm[:,None];face_rhs=face_rhs/eqnorm
    particular=np.linalg.lstsq(face_eq,face_rhs,rcond=None)[0];Z=null_space(face_eq)
    scaled_design=matrix[np.ix_(selected,free)]*column_scale
    design_face=scaled_design@Z/scales[selected,None]
    z=np.linalg.lstsq(design_face,(target[selected]-scaled_design@particular)/scales[selected],rcond=None)[0]
    candidate=np.zeros(p);candidate[free]=column_scale*(particular+Z@z)
    # Certify the returned coefficients, not the pre-polish optimizer point.
    # The equality-eliminated map is full column rank. Its inverse recovers
    # the objective's whitened coordinates without changing the coefficients.
    # transform=T V diag(1/sv) is ALREADY factorized, with orthonormal T,V.
    # A second least-squares SVD of the physical transform can discard or
    # contaminate directions merely because coefficient units differ greatly
    # (e.g. a quartic invariant amplitude). Use its exact factored left inverse.
    # This changes no objective/constraint/tolerance or material coefficient.
    def certificate(coefficients):
        final_x=sv*(Vt@(T.T@(coefficients-offset)))
        final_slack=A@final_x-h
        final_active=np.flatnonzero(final_slack<=1e-7*max(1.,np.linalg.norm(final_x)))
        dual=(nnls(A[final_active].T,final_x-x0,
                   maxiter=max(500,20*len(final_active)))[0]
              if len(final_active) else np.empty(0))
        stationarity=float(np.linalg.norm(final_x-x0-A[final_active].T@dual,np.inf))
        primal=max(0.,float(-np.min(final_slack)))
        complementarity=float(np.max(abs(dual*final_slack[final_active]),initial=0.))
        return stationarity,primal,complementarity
    polish_accepted=False
    if np.min(G@candidate)>=-1e-9 and np.max(abs(eq@candidate-rhs))<1e-8:
        candidate_certificate=certificate(candidate)
        if candidate_certificate[0]<=1e-6 and candidate_certificate[1]<=1e-7:
            c=candidate;polish_accepted=True
    # Never replace a verified primal point merely because an unverified
    # physical-coordinate polish happened to be feasible. Certify the actual
    # returned vector with the SAME existing tolerances in either case.
    final_kkt,final_violation,complementarity=certificate(c)
    if final_kkt>1e-6 or final_violation>1e-7:
        raise ArithmeticError(f'returned profile not verified: primal={final_violation}; KKT={final_kkt}; pre-polish KKT={kkt}')
    residual = (matrix@c-target)/scales
    return dict(coefficients=c,residuals=residual,predictions=matrix@c,
        squared_loss=float(residual[selected]@residual[selected]), selected_rows=selected,
        exact_residual=float(np.max(abs(residual[exact]))),kkt_residual=final_kkt,
        pre_polish_kkt_residual=kkt,complementarity_residual=complementarity,
        active_face_polish_accepted=polish_accepted,
        scaled_primal_violation=final_violation,optimizer_success=bool(run.success),
        optimizer_message=str(run.message),strictly_positive_LJ=bool(np.all(c[:2]>0)))


def spectral_profile(matrix, observations, columns, tails, *, max_cuts=32, nonnegative=(0,1,2,4,5), inequalities=None):
    """Adaptive polarization separation at EVERY supplied q, including tails."""
    columns,tails=np.asarray(columns),np.asarray(tails)
    p=matrix.shape[1]; signed=[i for i in range(p) if i not in nonnegative]
    signs=np.ones((2**len(signed),p))
    signs[:,signed]=np.array(list(product((-1.,1.),repeat=len(signed))))
    cuts=[] if inequalities is None else list(np.asarray(inequalities,float).reshape(-1,p)); history=[]
    fit=convex_profile(matrix,observations,cuts,nonnegative=nonnegative)
    for iteration in range(max_cuts+1):
        c=fit['coefficients']; operators=np.einsum('c,qcij->qij',c,columns)
        eigen,vec=np.linalg.eigh(operators); uncertainty=tails@abs(c)
        margin=eigen[:,0]-uncertainty
        worst=int(np.argmin(margin))
        history.append(dict(loss=fit['squared_loss'],minimum_robust_margin=float(margin[worst]),worst_index=worst))
        roundoff=1024*np.finfo(float).eps*max(1.,np.max(np.einsum('c,qcij->qij',abs(c),abs(columns))))
        if margin[worst]>=-roundoff:
            fit.update(robust_sampled=True,cut_iterations=iteration,spectral_history=history,
                minimum_robust_margin=float(np.min(margin)),minimum_sampled_H=float(np.min(eigen[:,0])),
                tail_at_worst=float(uncertainty[worst]),whole_zone_proved=False)
            return fit
        if iteration==max_cuts: break
        v=vec[worst,:,0]; row=np.einsum('i,cij,j->c',v,columns[worst],v)
        cuts.extend(row-signs*tails[worst])
        fit=convex_profile(matrix,observations,cuts,nonnegative=nonnegative)
    raise ArithmeticError(f'spectral separation did not converge: {history[-1]}')


def declared_wavepoints(step=.25):
    """FCC irreducible wedge plus the previously offending Gamma-K point."""
    if not np.isfinite(step) or not 0<step<=1:
        raise ValueError('finite wavevector grid step in (0,1] required')
    axis=np.arange(0.,1.+step/2,step)
    points={(float(x),float(y),float(z)) for x in axis for y in axis for z in axis
            if 0<=z<=y<=x<=1 and 0<x+y+z<=1.5+1e-12}
    points.update({(.375,.375,0.),(.01,0.,0.),(.01,.01,0.),(.01,.01,.01)})
    return np.array(sorted(points))


def augmented_wavepoints(step=.25, additional=()):
    """Deterministic stability-constraint refinement, never target fitting.

    Independent validation counterexamples may be added explicitly. The
    resulting finite set is STILL not a proof over the whole Brillouin zone.
    Coordinates use the existing corrected cubic reciprocal convention.
    """
    base=declared_wavepoints(step)
    extra=np.asarray(additional,float)
    if extra.size==0:return base
    if (extra.ndim!=2 or extra.shape[1]!=3 or np.any(~np.isfinite(extra))
            or np.any(extra[:,2]<0) or np.any(extra[:,0]>1)
            or np.any(extra[:,0]<extra[:,1]) or np.any(extra[:,1]<extra[:,2])
            or np.any(extra.sum(axis=1)<=0) or np.any(extra.sum(axis=1)>1.5+1e-12)):
        raise ValueError('nonzero finite wavepoints in the declared FCC irreducible wedge required')
    return np.unique(np.vstack([base,extra]),axis=0)
