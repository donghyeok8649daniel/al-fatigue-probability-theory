import math

import numpy as np

from theory.exact_lattice_energy import (
    barrier_spacing,
    bulk_energy_per_atom,
    bulk_force,
    crack_gap_energy,
    critical_spacing,
    critical_stretch,
    energy_scale,
    equilibrium_spacing,
    finite_chain_energy,
    gap_sum,
    normalized_phi,
    pair_potential,
    riemann_zeta,
)


M, N = 12.19, 6.0


def test_finite_pair_count_identity_and_bidirectional_half_count():
    atoms, a = 13, 1.07
    direct = sum(pair_potential(abs(i-j)*a, 2.3, 0.91, M, N)
                 for i in range(atoms) for j in range(i+1, atoms))
    grouped = finite_chain_energy(atoms, a, 2.3, 0.91, M, N)
    assert math.isclose(direct, grouped, rel_tol=2e-14)
    one_sided = sum(pair_potential(k*a, 2.3, 0.91, M, N) for k in range(1, 20000))
    bidirectional_with_half = 0.5 * sum(
        pair_potential(abs(k)*a, 2.3, 0.91, M, N)
        for k in range(-19999, 20000) if k)
    assert math.isclose(one_sided, bidirectional_with_half, rel_tol=2e-14)


def test_thermodynamic_limit_converges_to_zeta_energy():
    exact = bulk_energy_per_atom(1.03, 1.7, 0.96, M, N)
    errors = [abs(finite_chain_energy(k, 1.03, 1.7, 0.96, M, N)/k-exact)
              for k in (50, 200, 1000)]
    assert errors[2] < errors[1] < errors[0]
    # The leading finite-size surface correction is O(1/N), not a zeta-tail
    # truncation error.
    assert errors[-1] < 5e-4


def test_equilibrium_critical_and_phi_equivalence():
    eps, sig = 1.8, 0.94
    a0 = equilibrium_spacing(sig, M, N)
    ac = critical_spacing(sig, M, N)
    h = 2e-5 * a0
    first = (bulk_energy_per_atom(a0+h, eps, sig, M, N)
             - bulk_energy_per_atom(a0-h, eps, sig, M, N))/(2*h)
    second = (bulk_energy_per_atom(ac+h, eps, sig, M, N)
              - 2*bulk_energy_per_atom(ac, eps, sig, M, N)
              + bulk_energy_per_atom(ac-h, eps, sig, M, N))/h**2
    assert abs(first) < 1e-7
    assert abs(second) < 2e-5
    assert math.isclose(ac/a0, critical_stretch(M, N), rel_tol=2e-14)
    assert math.isclose(critical_stretch(M, N), 1.1077715386, rel_tol=2e-10)
    e0 = energy_scale(eps, sig, M, N)
    for lam in (0.91, 1.0, 1.08, 1.3):
        assert math.isclose(bulk_energy_per_atom(lam*a0, eps, sig, M, N),
                            e0*normalized_phi(lam, M, N), rel_tol=3e-14)


def test_gap_hurwitz_identity_and_direct_cross_gap_sum():
    q, p = 1.13, 6.0
    direct = sum((ell+1)*(q+ell)**(-p) for ell in range(200000))
    assert math.isclose(gap_sum(p, q), direct, rel_tol=2e-13)
    g, a0, eps, sig = q, 1.0, 2.0, 0.9
    grouped = sum((ell+1)*pair_potential(g+ell*a0, eps, sig, M, N)
                  for ell in range(100000))
    assert math.isclose(crack_gap_energy(g, a0, eps, sig, M, N), grouped,
                        rel_tol=3e-12)


def test_stable_and_barrier_roots_merge_at_critical_force():
    eps, sig = 1.0, 1.0
    ac = critical_spacing(sig, M, N)
    fc = bulk_force(ac, eps, sig, M, N)
    stable, barrier = barrier_spacing(0.999*fc, eps, sig, M, N)
    assert stable < ac < barrier
    assert math.isclose(bulk_force(stable, eps, sig, M, N), 0.999*fc, rel_tol=2e-13)
    assert math.isclose(bulk_force(barrier, eps, sig, M, N), 0.999*fc, rel_tol=2e-13)
