"""Long-wave nonlocal registry ENERGY reference, never a fatigue solver.

x and s use the fixed microscopic L0. gamma is in eV/L0^2, elastic moduli
in eV/L0^3. Energies below are eV/L0 of DISLOCATION LINE, not an activation
energy in eV. No arbitrary coherent area, correlation area, time, or mobility.
The periodic misfit series is sampled from the unchanged analytic interface;
sampling and derivatives must be independently checked against that surface.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.sparse.linalg import LinearOperator, cg


@dataclass(frozen=True)
class PeriodicMisfit:
    period: float
    coefficients: np.ndarray

    @classmethod
    def from_samples(cls, samples, *, period):
        """Equispaced periodic samples including s=0, excluding s=period.

        Real part of sum c_m exp(i 2 pi m s/period). No new empirical fit.
        No polynomial spline is inserted into any production energy model.
        """
        values = np.asarray(samples, dtype=float)
        if values.ndim != 1 or len(values) < 8 or len(values) % 2 or not np.all(np.isfinite(values)):
            raise ValueError("finite even-length periodic sample vector required")
        if not np.isfinite(period) or period <= 0:
            raise ValueError("finite positive registry period required")
        coefficients = 2*np.fft.rfft(values)/len(values)
        coefficients[[0, -1]] *= .5
        return cls(float(period), coefficients)

    def evaluate(self, s, order=0):
        if order not in (0, 1, 2):
            raise ValueError("energy, analytic first/second derivative supported")
        wave = 2*np.pi*np.arange(len(self.coefficients))/self.period
        c = self.coefficients*(1j*wave)**order if order else self.coefficients
        z = np.exp(2j*np.pi*np.asarray(s)/self.period)
        return np.polynomial.polynomial.polyval(z, c).real


class ScrewRegistryFunctional:
    r"""Periodic straight-line screw geometry, variation x perpendicular to slip.

    E/line = 1/2 int s (mu_bar/2)|D|s dx + int gamma(s) dx - int tau(x)s dx.
    mu_bar is derived sqrt(A D-B^2), NOT a fitted gradient-penalty parameter.
    The local normal gap is clamped. No opening, thermal probabilities, defect
    sources, Peierls lattice pinning, or full two-component partial splitting.
    """
    def __init__(self, misfit, *, elastic_factor, domain_length, cells):
        if not np.isfinite(elastic_factor) or elastic_factor <= 0:
            raise ValueError("positive derived elastic factor required")
        if not np.isfinite(domain_length) or domain_length <= 0 or cells < 32 or cells % 2:
            raise ValueError("positive periodic domain and even cells>=32 required")
        self.misfit = misfit
        self.mu = float(elastic_factor)
        self.length = float(domain_length)
        self.cells = int(cells)
        self.dx = self.length/self.cells
        self.x = np.arange(self.cells)*self.dx-self.length/2
        self.q = 2*np.pi*np.fft.fftfreq(self.cells, d=self.dx)
        self.kernel = self.mu/2*abs(self.q)

    def elastic_force(self, slip):
        return np.fft.ifft(self.kernel*np.fft.fft(slip)).real

    def evaluate(self, slip, traction=0.):
        s = np.asarray(slip, dtype=float)
        tau = np.broadcast_to(traction, s.shape)
        if s.shape != (self.cells,) or not np.all(np.isfinite(s)) or not np.all(np.isfinite(tau)):
            raise ValueError("finite slip/traction arrays on declared grid required")
        force = self.elastic_force(s)
        elastic = self.dx*np.dot(s, force)/2
        misfit = self.dx*np.sum(self.misfit.evaluate(s))
        work = -self.dx*np.dot(s, tau)
        return dict(energy=float(elastic+misfit+work), elastic=float(elastic),
            misfit=float(misfit), external_work=float(work),
            force=force+self.misfit.evaluate(s, 1)-tau)

    def hessian_action(self, slip, variation):
        return self.elastic_force(variation)+self.misfit.evaluate(slip, 2)*variation

    def linear_response(self, traction, *, reference=0.):
        """Within a stable well; s_hat=tau_hat/(gamma''+mu_bar |q|/2)."""
        curvature = float(self.misfit.evaluate(reference, 2))
        if curvature <= 0:
            raise ValueError("linear reference is not a stable registry minimum")
        tau = np.broadcast_to(traction, (self.cells,))-self.misfit.evaluate(reference, 1)
        return reference+np.fft.ifft(np.fft.fft(tau)/(curvature+self.kernel)).real

    def solve_intrawell_equilibrium(self, traction, *, initial=None, force_tolerance=1e-11):
        """Newton-CG static solution while the local misfit curvature is positive.

        A rejected nonconvex branch is reported as unsupported, never clipped
        into an elastic law. Iteration count is NOT time. No residual plastic
        deformation is established by this quasistatic unload calculation.
        """
        if not np.isfinite(force_tolerance) or force_tolerance <= 0:
            raise ValueError("positive force tolerance required")
        slip=self.linear_response(traction) if initial is None else np.asarray(initial,dtype=float).copy()
        for iteration in range(30):
            value=self.evaluate(slip,traction)
            residual=float(max(abs(value["force"])))
            curvature=self.misfit.evaluate(slip,2)
            if min(curvature) <= 0:
                raise ValueError("nonconvex branch: this intrawell solver does not certify it")
            if residual <= force_tolerance:
                return dict(slip=slip,force_residual=residual,iterations=iteration,energy=value["energy"],
                            minimum_misfit_curvature=float(min(curvature)))
            operator=LinearOperator((self.cells,self.cells),matvec=lambda v:self.elastic_force(v)+curvature*v)
            preconditioner=LinearOperator((self.cells,self.cells),
                matvec=lambda v:np.fft.ifft(np.fft.fft(v)/(self.kernel+np.mean(curvature))).real)
            step,info=cg(operator,-value["force"],M=preconditioner,rtol=1e-11,atol=0.,maxiter=300)
            if info != 0:
                raise RuntimeError("static Newton linear solve did not converge")
            for reduction in range(25):
                trial=slip+2.**(-reduction)*step
                if min(self.misfit.evaluate(trial,2))>0 and max(abs(self.evaluate(trial,traction)["force"]))<residual:
                    slip=trial; break
            else:
                raise RuntimeError("static intrawell continuation failed; no state was repaired")
        raise RuntimeError("static intrawell force tolerance not reached")

    def dipole_trial(self, separation, width):
        """Periodized arctangent dipole; specified pre-existing opposite defects.

        Neither separation nor width is a specimen correlation area. Width is
        a variational coordinate. This is not a proven nucleation saddle.
        Returns s, partial_d s, partial_width s, with mean s=b*d/L exactly.
        """
        if not 0 < separation < self.length/2 or width <= 0 or not np.isfinite(width):
            raise ValueError("0<separation<L/2 and positive width required")
        q = self.q; b = self.misfit.period
        coeff = np.zeros(len(q)); dc = coeff.copy(); wc = coeff.copy()
        nonzero = q != 0
        damp = np.exp(-abs(q[nonzero])*width)
        coeff[nonzero] = 2*b/(self.length*q[nonzero])*np.sin(q[nonzero]*separation/2)*damp
        dc[nonzero] = b/self.length*np.cos(q[nonzero]*separation/2)*damp
        wc[nonzero] = -abs(q[nonzero])*coeff[nonzero]
        coeff[0] = b*separation/self.length; dc[0] = b/self.length
        # Shift center to x=0 in the returned [-L/2,L/2) grid.
        return tuple(np.fft.fftshift(np.fft.ifft(c*self.cells).real) for c in (coeff, dc, wc))

    def dipole_elastic_closed_form(self, separation, width):
        r"""Infinite Fourier-series value: mu b^2/(4pi) log(1+sin^2/sinh^2)."""
        theta = np.pi*separation/self.length
        alpha = 2*np.pi*width/self.length
        return self.mu*self.misfit.period**2/(4*np.pi)*np.log1p((np.sin(theta)/np.sinh(alpha))**2)

    def optimize_trial_width(self, separation, *, width_bounds):
        """Stationary width within a declared restricted arctangent family.

        NOT a full saddle/core optimization or a material-calibrated defect.
        Independent mesh/domain/bounds checks are required by the caller.
        """
        lower, upper = width_bounds
        if not 0 < lower < upper:
            raise ValueError("positive ordered variational-width bracket required")
        def objective(log_width):
            return self.evaluate(self.dipole_trial(separation, np.exp(log_width))[0])["energy"]
        result = minimize_scalar(objective, bounds=np.log([lower, upper]), method="bounded",
                                 options={"xatol":2e-10, "maxiter":150})
        width = float(np.exp(result.x))
        slip, ds_dd, ds_dw = self.dipole_trial(separation, width)
        values = self.evaluate(slip)
        gradient_d = self.dx*np.dot(values["force"], ds_dd)
        gradient_w = self.dx*np.dot(values["force"], ds_dw)
        spectrum = abs(np.fft.fft(slip))**2*self.kernel
        order = np.argsort(abs(self.q)); cumulative = np.cumsum(spectrum[order])
        q90 = abs(self.q[order[np.searchsorted(cumulative, .9*cumulative[-1])]])
        return dict(separation=float(separation), width=width, energy=values["energy"],
            elastic=values["elastic"], misfit=values["misfit"],
            balance_traction=float(gradient_d/self.misfit.period),
            width_force=float(gradient_w), full_euler_residual=float(max(abs(values["force"]))),
            spectral_q90_times_period=float(q90*self.misfit.period),
            elastic_series_error=float(values["elastic"]-self.dipole_elastic_closed_form(separation,width)),
            optimizer_success=bool(result.success),
            width_interior=bool(lower*1.001 < width < upper/1.001),
            variational_width_bounds=[float(lower),float(upper)])

    def solve_fixed_registry_content(self, initial, *, force_tolerance=1e-8, max_iterations=1800):
        """Relax ALL profile modes at fixed mean slip (static PN reference).

        Fixed int s dx specifies a pre-existing dipole's registry content,
        NOT a fabricated pinning potential. The mean Euler force is the
        Lagrange traction holding that content. Projected force must vanish.
        The profile is not restricted to an arctangent/one-width family.
        No time, nucleation rate, or atomistic core certification is implied.

        Fourier square-root preconditioning changes optimizer coordinates,
        not the energy. Failure to reach the force tolerance is explicit.
        """
        initial = np.asarray(initial, float)
        if initial.shape != (self.cells,) or force_tolerance <= 0 or max_iterations < 1:
            raise ValueError("finite full initial profile and positive solver controls required")
        mean = float(np.mean(initial))
        curvature = float(self.misfit.evaluate(0.,2))
        if curvature <= 0:
            raise ValueError("stable far-field registry required")
        root = np.sqrt(self.kernel+curvature)
        def transform(v, multiplier):
            out = np.fft.ifft(multiplier*np.fft.fft(v)).real
            return out-np.mean(out)
        def profile(y):
            return mean+transform(y,1/root)
        def objective(y):
            s = profile(y)
            value = self.evaluate(s)
            # Division by dx is a common objective/gradient normalization;
            # it does not modify the actual physical energy being minimized.
            return value["energy"]/self.dx, transform(value["force"],1/root)
        result = minimize(objective,transform(initial-mean,root),jac=True,method="L-BFGS-B",
            options=dict(gtol=force_tolerance/20,ftol=1e-15,maxiter=max_iterations,maxls=50,maxcor=20))
        slip = profile(result.x)
        # Objective-relative termination alone need not satisfy force balance.
        # Polish the projected Euler equation; reject failed linear steps.
        # The constant mode is the fixed-mean constraint, not a spring.
        polish_iterations = 0
        for _ in range(8):
            force = self.evaluate(slip)["force"]
            projected = force-np.mean(force)
            residual = float(max(abs(projected)))
            if residual <= force_tolerance:
                break
            local_curvature = self.misfit.evaluate(slip,2)
            def action(v):
                p = v-np.mean(v)
                w = self.elastic_force(p)+local_curvature*p
                return w-np.mean(w)+np.mean(v)
            operator = LinearOperator((self.cells,self.cells),matvec=action)
            preconditioner = LinearOperator((self.cells,self.cells),
                matvec=lambda v:np.fft.ifft(np.fft.fft(v)/root**2).real)
            step,info = cg(operator,-projected,M=preconditioner,rtol=1e-9,atol=0.,maxiter=600)
            if info != 0:
                break
            step -= np.mean(step)
            for reduction in range(20):
                trial = slip+2.**(-reduction)*step
                trial += mean-np.mean(trial)
                tf = self.evaluate(trial)["force"]
                if max(abs(tf-np.mean(tf))) < residual:
                    slip = trial; polish_iterations += 1; break
            else:
                break
        value = self.evaluate(slip)
        traction = float(np.mean(value["force"]))
        residual = float(max(abs(value["force"]-traction)))
        spectrum = abs(np.fft.fft(slip))**2*self.kernel
        order = np.argsort(abs(self.q)); cumulative = np.cumsum(spectrum[order])
        q90 = (float(abs(self.q[order[np.searchsorted(cumulative,.9*cumulative[-1])]]))
               if cumulative[-1]>0 else 0.)
        return dict(slip=slip,mean_slip=mean,registry_content=float(self.length*mean),
            mean_constraint_residual=float(abs(np.mean(slip)-mean)),
            holding_traction=traction,projected_force_residual=residual,
            force_tolerance=force_tolerance,force_converged=bool(residual<=force_tolerance),
            energy=value["energy"],elastic=value["elastic"],misfit=value["misfit"],
            iterations=int(result.nit),newton_polish_iterations=polish_iterations,
            optimizer_success=bool(result.success),
            optimizer_message=str(result.message),spectral_q90_times_period=q90*self.misfit.period,
            status="all-profile stationary check of continuum PN reference, NOT atomistic core")
