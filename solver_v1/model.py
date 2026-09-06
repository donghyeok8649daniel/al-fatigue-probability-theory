from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.interpolate import RegularGridInterpolator

from .lattice_bessel import (
    FourierLatticeConfig,
    two_row_lj_infinite_energy_gradient,
)


@dataclass
class ModelParams:
    n_cells: int = 3
    b: float = 1.0
    epsilon: float = 1.0
    sigma_lj: float = 0.82
    # Deprecated compatibility field. Production lower-row interactions are
    # evaluated by the analytic Poisson/Bessel infinite-lattice kernel.
    lower_images: int = 14
    lattice_fourier_tol: float = 1.0e-13
    lattice_fourier_max_modes: int = 64
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

    The interaction between each upper cell and the infinite lower row is the
    canonical Poisson-summed Bessel lattice kernel W(a_i,s_i). Upper-cell
    interactions remain explicit pair interactions so cell-cell correlation is
    retained without product closure.

    s_i is an unwrapped configurational coordinate. Its periodic part controls
    lower-row registry, while differences s_j-s_i enter upper-row interactions.
    """

    def __init__(self, p: ModelParams):
        self.p = p
        self._lattice_config = FourierLatticeConfig(
            tol=p.lattice_fourier_tol,
            max_modes=p.lattice_fourier_max_modes,
        )
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

    def _lower_lattice_energy_gradient(self, a, s):
        p = self.p
        return two_row_lj_infinite_energy_gradient(
            a,
            s,
            epsilon=p.epsilon,
            sigma_lj=p.sigma_lj,
            b=p.b,
            config=self._lattice_config,
        )

    def local_energy(self, a: float, s: float) -> float:
        energy, _, _, _ = self._lower_lattice_energy_gradient(float(a), float(s))
        return float(energy)

    def local_deda(self, a: float, s: float) -> float:
        _, deda, _, _ = self._lower_lattice_energy_gradient(float(a), float(s))
        return float(deda)

    def local_deds(self, a: float, s: float) -> float:
        _, _, deds, _ = self._lower_lattice_energy_gradient(float(a), float(s))
        return float(deds)

    def local_deda_array(self, a, s: float):
        aa = np.asarray(a, dtype=float)
        ss = np.full_like(aa, float(s), dtype=float)
        _, deda, _, _ = self._lower_lattice_energy_gradient(aa, ss)
        return np.asarray(deda, dtype=float)

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

        local_u, local_ga, local_gs, _ = self._lower_lattice_energy_gradient(a, s)
        U += float(np.sum(local_u))
        ga += local_ga
        gs += local_gs

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

        local_u, local_ga, local_gs, _ = self._lower_lattice_energy_gradient(a, s)
        U += np.sum(local_u, axis=1)
        ga += local_ga
        gs += local_gs

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
