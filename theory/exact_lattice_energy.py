"""Exact one-dimensional generalized-LJ lattice and crack-gap energies.

The functions here keep three objects deliberately separate: one pair energy,
the total energy of a finite chain, and the thermodynamic-limit energy per
atom.  No spatial cut-off or Taylor approximation enters the lattice formulas.
"""
from __future__ import annotations

import math


def _check(m: float, n: float) -> None:
    if not (m > n > 1.0):
        raise ValueError("require m > n > 1")


def hurwitz_zeta(s: float, q: float, terms: int = 64) -> float:
    """Hurwitz zeta by an Euler--Maclaurin tail (double precision).

    This is a numerical evaluator of the exact special function, not a
    truncated physical interaction model.  Four Bernoulli corrections make
    ``terms=64`` ample for the exponents used by this project.
    """
    if not (s > 1.0 and q > 0.0 and terms >= 8):
        raise ValueError("require s>1, q>0 and terms>=8")
    total = math.fsum((q + k) ** (-s) for k in range(terms))
    x = q + terms
    total += x ** (1.0 - s) / (s - 1.0) + 0.5 * x ** (-s)
    # B_2/(2!), B_4/(4!), B_6/(6!), B_8/(8!).
    coefficients = (1.0 / 12.0, -1.0 / 720.0, 1.0 / 30240.0, -1.0 / 1209600.0)
    rising = s
    for order, coefficient in enumerate(coefficients, start=1):
        if order > 1:
            rising *= (s + 2 * order - 3) * (s + 2 * order - 2)
        total += coefficient * rising * x ** (-s - 2 * order + 1)
    return total


def riemann_zeta(s: float) -> float:
    return hurwitz_zeta(s, 1.0)


def pair_potential(r: float, epsilon: float, sigma_lj: float, m: float, n: float) -> float:
    _check(m, n)
    if min(r, epsilon, sigma_lj) <= 0.0:
        raise ValueError("r, epsilon and sigma_lj must be positive")
    return epsilon * ((sigma_lj / r) ** m - (sigma_lj / r) ** n)


def finite_chain_energy(n_atoms: int, a: float, epsilon: float, sigma_lj: float,
                        m: float = 12.19, n: float = 6.0) -> float:
    """Exact total pair energy sum_k (N-k) v(k a)."""
    if n_atoms < 2:
        raise ValueError("n_atoms must be at least two")
    return math.fsum((n_atoms - k) * pair_potential(k * a, epsilon, sigma_lj, m, n)
                     for k in range(1, n_atoms))


def bulk_energy_per_atom(a: float, epsilon: float, sigma_lj: float,
                         m: float = 12.19, n: float = 6.0) -> float:
    """Exact thermodynamic-limit energy per atom (or representative cell)."""
    _check(m, n)
    if min(a, epsilon, sigma_lj) <= 0.0:
        raise ValueError("a, epsilon and sigma_lj must be positive")
    return epsilon * (riemann_zeta(m) * (sigma_lj / a) ** m
                      - riemann_zeta(n) * (sigma_lj / a) ** n)


def equilibrium_spacing(sigma_lj: float, m: float = 12.19, n: float = 6.0) -> float:
    _check(m, n)
    return sigma_lj * (m * riemann_zeta(m) / (n * riemann_zeta(n))) ** (1.0 / (m - n))


def critical_spacing(sigma_lj: float, m: float = 12.19, n: float = 6.0) -> float:
    _check(m, n)
    return sigma_lj * (m * (m + 1) * riemann_zeta(m)
                       / (n * (n + 1) * riemann_zeta(n))) ** (1.0 / (m - n))


def critical_stretch(m: float = 12.19, n: float = 6.0) -> float:
    _check(m, n)
    return ((m + 1.0) / (n + 1.0)) ** (1.0 / (m - n))


def energy_scale(epsilon: float, sigma_lj: float, m: float = 12.19,
                 n: float = 6.0) -> float:
    """E0 such that U_inf(a)=E0 phi(a/a0), with no additive constant."""
    a0 = equilibrium_spacing(sigma_lj, m, n)
    return epsilon * m * (m - n) * riemann_zeta(m) * (sigma_lj / a0) ** m


def epsilon_from_energy_scale(e0: float, sigma_lj: float, m: float = 12.19,
                              n: float = 6.0) -> float:
    if e0 <= 0.0:
        raise ValueError("e0 must be positive")
    a0 = equilibrium_spacing(sigma_lj, m, n)
    return e0 / (m * (m - n) * riemann_zeta(m) * (sigma_lj / a0) ** m)


def normalized_phi(lam: float, m: float = 12.19, n: float = 6.0) -> float:
    _check(m, n)
    if lam <= 0.0:
        raise ValueError("lambda must be positive")
    return lam ** (-m) / (m * (m - n)) - lam ** (-n) / (n * (m - n))


def gap_sum(s: float, q: float) -> float:
    """S_s(q)=zeta(s-1,q)+(1-q)zeta(s,q)."""
    if s <= 2.0:
        raise ValueError("the half-chain cross-gap sum requires s>2")
    return hurwitz_zeta(s - 1.0, q) + (1.0 - q) * hurwitz_zeta(s, q)


def crack_gap_energy(g: float, a0: float, epsilon: float, sigma_lj: float,
                     m: float = 12.19, n: float = 6.0) -> float:
    """All cross-gap pairs between two semi-infinite chains."""
    _check(m, n)
    if min(g, a0, epsilon, sigma_lj) <= 0.0:
        raise ValueError("all dimensional inputs must be positive")
    q = g / a0
    return epsilon * ((sigma_lj / a0) ** m * gap_sum(m, q)
                      - (sigma_lj / a0) ** n * gap_sum(n, q))


def bulk_force(a: float, epsilon: float, sigma_lj: float,
               m: float = 12.19, n: float = 6.0) -> float:
    """dU_inf/da, the tensile force conjugate to homogeneous spacing."""
    return epsilon / a * (-m * riemann_zeta(m) * (sigma_lj / a) ** m
                          + n * riemann_zeta(n) * (sigma_lj / a) ** n)


def barrier_spacing(force: float, epsilon: float, sigma_lj: float,
                    m: float = 12.19, n: float = 6.0) -> tuple[float, float]:
    """Stable and unstable roots of U'(a)=F for 0<F<F_c."""
    if force <= 0.0:
        raise ValueError("force must be positive")
    a0 = equilibrium_spacing(sigma_lj, m, n)
    ac = critical_spacing(sigma_lj, m, n)
    fc = bulk_force(ac, epsilon, sigma_lj, m, n)
    if force >= fc:
        raise ValueError("force must be below the tangent-instability force")

    def bisect(lo: float, hi: float) -> float:
        flo = bulk_force(lo, epsilon, sigma_lj, m, n) - force
        for _ in range(90):
            mid = 0.5 * (lo + hi)
            fm = bulk_force(mid, epsilon, sigma_lj, m, n) - force
            if flo * fm <= 0.0:
                hi = mid
            else:
                lo, flo = mid, fm
        return 0.5 * (lo + hi)

    stable = bisect(a0, ac)
    hi = 2.0 * ac
    while bulk_force(hi, epsilon, sigma_lj, m, n) > force:
        hi *= 2.0
    return stable, bisect(ac, hi)
