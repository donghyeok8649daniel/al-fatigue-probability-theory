from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.interpolate import RegularGridInterpolator

@dataclass
class ModelParams:
    n_cells: int = 3
    b: float = 1.0
    epsilon: float = 1.0
    sigma_lj: float = 0.82
    lower_images: int = 14
    mobility_a: float = 1.0
    mobility_s: float = 0.15
    kT: float = 0.009
    chi_axial_projection: float = 0.40
    a_min: float = 0.65
    a_max: float = 3.0

class TwoRowLJ:
    """Dimensionless interacting two-row LJ model.

    Lower row: x_j = j*b, y=0.
    Upper cells: x_i = i*b + b/2 + s_i, y=a_i.

    s_i is an unwrapped configurational coordinate. Its periodic part controls
    lower-row registry, while differences s_j-s_i enter upper-row interactions.
    """

    def __init__(self, p: ModelParams):
        self.p = p
        self.lower_n = np.arange(-p.lower_images, p.lower_images + 1, dtype=float)
        self.a0 = self._find_reference_a()
        self._opening_table_ready = False

    def phi(self, r):
        p = self.p
        x = p.sigma_lj / np.asarray(r)
        return 4.0 * p.epsilon * (x**12 - x**6)

    def dphi(self, r):
        p = self.p
        r = np.asarray(r)
        return 4.0 * p.epsilon * (
            -12.0 * p.sigma_lj**12 * r**(-13)
            + 6.0 * p.sigma_lj**6 * r**(-7)
        )

    def local_energy(self, a: float, s: float) -> float:
        p = self.p
        dx = (self.lower_n + 0.5) * p.b - np.mod(s + 0.5*p.b, p.b) + 0.5*p.b
        r = np.sqrt(dx*dx + a*a)
        return float(np.sum(self.phi(r)))

    def local_deda(self, a: float, s: float) -> float:
        p = self.p
        dx = (self.lower_n + 0.5) * p.b - np.mod(s + 0.5*p.b, p.b) + 0.5*p.b
        r = np.sqrt(dx*dx + a*a)
        return float(np.sum(self.dphi(r) * a / r))

    def local_deds(self, a: float, s: float) -> float:
        p = self.p
        smod = np.mod(s + 0.5*p.b, p.b) - 0.5*p.b
        dx = (self.lower_n + 0.5) * p.b - smod
        r = np.sqrt(dx*dx + a*a)
        return float(np.sum(-self.dphi(r) * dx / r))

    def local_deda_array(self, a, s: float):
        p = self.p
        aa = np.asarray(a, dtype=float)[:, None]
        smod = np.mod(s + 0.5*p.b, p.b) - 0.5*p.b
        dx = ((self.lower_n + 0.5) * p.b - smod)[None, :]
        r = np.sqrt(dx*dx + aa*aa)
        return np.sum(self.dphi(r) * aa / r, axis=1)

    def _build_opening_table(self, force_max: float = 6.0):
        p = self.p
        self._sgrid = np.linspace(-0.5*p.b, 0.5*p.b, 81)
        self._fgrid = np.linspace(0.0, force_max, 121)
        amin_tab = np.full((len(self._sgrid), len(self._fgrid)), np.nan)
        asad_tab = np.full_like(amin_tab, np.nan)
        barrier_tab = np.full_like(amin_tab, np.nan)
        fc_tab = np.zeros(len(self._sgrid))

        agrid = np.linspace(p.a_min, p.a_max, 1800)
        for is_, sval in enumerate(self._sgrid):
            tr = self.local_deda_array(agrid, float(sval))
            imax = int(np.argmax(tr))
            fc = float(tr[imax])
            fc_tab[is_] = fc

            left_a = agrid[:imax+1]
            left_f = tr[:imax+1]
            imin = int(np.argmin(np.abs(left_f)))
            left_a = left_a[imin:]
            left_f = left_f[imin:]
            order = np.argsort(left_f)
            lf, idx = np.unique(left_f[order], return_index=True)
            la = left_a[order][idx]

            right_a = agrid[imax:]
            right_f = tr[imax:]
            order_r = np.argsort(right_f)
            rf, idxr = np.unique(right_f[order_r], return_index=True)
            ra = right_a[order_r][idxr]

            for jf, force in enumerate(self._fgrid):
                if force >= fc:
                    continue
                if force <= max(0.0, float(np.min(lf))):
                    amin = self._find_reference_a() if abs(sval) < 1e-12 else float(
                        minimize_scalar(lambda x: self.local_energy(float(x), float(sval)),
                                        bounds=(p.a_min, 1.8), method="bounded").x
                    )
                else:
                    amin = float(np.interp(force, lf, la))
                amin_tab[is_, jf] = amin

                if force < 1e-8:
                    asad = 10.0*p.a_max
                    gmin = self.local_energy(amin, float(sval))
                    barrier = max(0.0, -gmin)
                elif force < float(np.max(rf)) and force >= float(np.min(rf)):
                    asad = float(np.interp(force, rf, ra))
                    gmin = self.local_energy(amin, float(sval)) - force*(amin-self.a0)
                    gsad = self.local_energy(asad, float(sval)) - force*(asad-self.a0)
                    barrier = max(0.0, gsad-gmin)
                else:
                    asad = 10.0*p.a_max
                    barrier = np.inf
                asad_tab[is_, jf] = asad
                barrier_tab[is_, jf] = barrier

        self._amin_interp = RegularGridInterpolator(
            (self._sgrid, self._fgrid), amin_tab, bounds_error=False, fill_value=np.nan
        )
        self._asad_interp = RegularGridInterpolator(
            (self._sgrid, self._fgrid), asad_tab, bounds_error=False, fill_value=np.nan
        )
        finite_barrier = np.where(np.isfinite(barrier_tab), barrier_tab, 1e6)
        self._barrier_interp = RegularGridInterpolator(
            (self._sgrid, self._fgrid), finite_barrier, bounds_error=False, fill_value=np.nan
        )
        self._fc_grid = fc_tab
        self._fc_interp = lambda sval: float(np.interp(sval, self._sgrid, self._fc_grid))
        self._opening_table_ready = True

    def _wrap_s(self, s: float) -> float:
        p = self.p
        return float(np.mod(s + 0.5*p.b, p.b) - 0.5*p.b)

    def _find_reference_a(self) -> float:
        out = minimize_scalar(
            lambda a: self.local_energy(float(a), 0.0),
            bounds=(self.p.a_min, 1.8),
            method="bounded",
        )
        return float(out.x)

    def energy_gradient(self, a: np.ndarray, s: np.ndarray, force: float):
        p = self.p
        a = np.asarray(a, dtype=float)
        s = np.asarray(s, dtype=float)
        U = 0.0
        ga = np.zeros_like(a)
        gs = np.zeros_like(s)

        for i in range(p.n_cells):
            smod = np.mod(s[i] + 0.5*p.b, p.b) - 0.5*p.b
            dx = (self.lower_n + 0.5) * p.b - smod
            r = np.sqrt(dx*dx + a[i]*a[i])
            dp = self.dphi(r)
            U += float(np.sum(self.phi(r)))
            ga[i] += float(np.sum(dp * a[i] / r))
            gs[i] += float(np.sum(-dp * dx / r))

        for i in range(p.n_cells):
            for j in range(i + 1, p.n_cells):
                dx = (j-i)*p.b + s[j] - s[i]
                dy = a[j] - a[i]
                r = float(np.sqrt(dx*dx + dy*dy))
                if r < 0.35*p.b:
                    return 1e18, np.full_like(a, np.nan), np.full_like(s, np.nan)
                dp = float(self.dphi(r))
                U += float(self.phi(r))
                ga[i] += dp * (-dy/r)
                ga[j] += dp * ( dy/r)
                gs[i] += dp * (-dx/r)
                gs[j] += dp * ( dx/r)

        U -= force * float(np.sum((a-self.a0) + p.chi_axial_projection*s))
        ga -= force
        gs -= force * p.chi_axial_projection
        return U, ga, gs

    def energy_gradient_batch(self, a: np.ndarray, s: np.ndarray, force: float):
        p = self.p
        a = np.asarray(a, dtype=float)
        s = np.asarray(s, dtype=float)
        B, N = a.shape
        ga = np.zeros_like(a)
        gs = np.zeros_like(s)
        U = np.zeros(B, dtype=float)

        smod = np.mod(s + 0.5*p.b, p.b) - 0.5*p.b
        dx = (self.lower_n[None,None,:] + 0.5) * p.b - smod[:,:,None]
        rr = np.sqrt(dx*dx + a[:,:,None]**2)
        dp = self.dphi(rr)
        U += np.sum(self.phi(rr), axis=(1,2))
        ga += np.sum(dp * a[:,:,None] / rr, axis=2)
        gs += np.sum(-dp * dx / rr, axis=2)

        for i in range(N):
            for j in range(i+1, N):
                dxp = (j-i)*p.b + s[:,j] - s[:,i]
                dyp = a[:,j] - a[:,i]
                r = np.sqrt(dxp*dxp + dyp*dyp)
                bad = r < 0.35*p.b
                rsafe = np.where(bad, 0.35*p.b, r)
                dpp = self.dphi(rsafe)
                U += self.phi(rsafe)
                ga[:,i] += dpp * (-dyp/rsafe)
                ga[:,j] += dpp * ( dyp/rsafe)
                gs[:,i] += dpp * (-dxp/rsafe)
                gs[:,j] += dpp * ( dxp/rsafe)
                ga[bad] = np.nan
                gs[bad] = np.nan

        U -= force*np.sum((a-self.a0) + p.chi_axial_projection*s, axis=1)
        ga -= force
        gs -= force*p.chi_axial_projection
        return U, ga, gs

    def strain(self, a: np.ndarray, s: np.ndarray) -> float:
        p = self.p
        return float(np.mean((a-self.a0)/self.a0 + p.chi_axial_projection*s/self.a0))

    def well_index(self, s: np.ndarray) -> np.ndarray:
        p = self.p
        return np.floor((np.asarray(s) + 0.5*p.b)/p.b).astype(int)

    def opening_saddle_batch(self, s: np.ndarray, force: float):
        if not self._opening_table_ready:
            self._build_opening_table()
        p = self.p
        ss = np.asarray(s, dtype=float)
        sw = np.mod(ss + 0.5*p.b, p.b) - 0.5*p.b
        fc = np.interp(sw.ravel(), self._sgrid, self._fc_grid).reshape(sw.shape)
        bound = force < fc
        fq = np.full(sw.size, float(np.clip(force, self._fgrid[0], self._fgrid[-1])))
        pts = np.column_stack((sw.ravel(), fq))
        amin = self._amin_interp(pts).reshape(sw.shape)
        asad = self._asad_interp(pts).reshape(sw.shape)
        return amin, asad, bound

    def opening_stationary_points(self, s: float, force: float):
        if not self._opening_table_ready:
            self._build_opening_table()
        sval = self._wrap_s(s)
        fc = self._fc_interp(sval)
        if force >= fc:
            return None, None
        fquery = float(np.clip(force, self._fgrid[0], self._fgrid[-1]))
        amin = float(self._amin_interp([[sval, fquery]])[0])
        asad = float(self._asad_interp([[sval, fquery]])[0])
        if not np.isfinite(amin):
            return None, None
        if not np.isfinite(asad):
            return amin, None
        return amin, asad

    def opening_barrier(self, s: float, force: float) -> float:
        if not self._opening_table_ready:
            self._build_opening_table()
        sval = self._wrap_s(s)
        fc = self._fc_interp(sval)
        if force >= fc:
            return 0.0
        fquery = float(np.clip(force, self._fgrid[0], self._fgrid[-1]))
        val = float(self._barrier_interp([[sval, fquery]])[0])
        return val
