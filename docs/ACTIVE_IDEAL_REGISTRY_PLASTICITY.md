# Active one-dimensional ideal-registry plasticity branch

## Status

The exact two-row Poisson--Bessel lattice energy is now an **optional active
mechanism branch**. It is not merely an archive. The primary crack-initiation
model remains the separate normal-chain model. The two energies are not added.

This distinction is geometric, not bureaucratic:

- the normal branch uses a collinear chain and homogeneous local spacing `a`;
- the registry branch uses two parallel one-index rows, prescribed separation
  `a`, row repeat `b`, and one scalar unwrapped registry `s`;
- a future jointly evolving `(a,s)` density would have two state coordinates,
  although the row lattice and the slip coordinate would remain reduced 1D;
- no 2D/3D continuum shear criterion is introduced.

## Exact energy used by the solver

For the repository coefficient convention

```text
v(r) = epsilon_c [(sigma_LJ/r)^m - (sigma_LJ/r)^n],  m > n > 1,
```

the configuration-dependent cross-row energy per upper atom, equivalently per
commensurate row repeat, is

```text
W(a,s) = sum over p in Z of v(sqrt(a^2+(p b+s)^2)).
```

Put `delta=s/b`, `eta=a/b`, and

```text
Z_nu(delta,eta) = sum over p in Z of [(p+delta)^2+eta^2]^(-nu).
```

Mellin transformation followed by Poisson summation gives the identity

```text
Z_nu = sqrt(pi)/Gamma(nu) [
  Gamma(nu-1/2) eta^(1-2 nu)
  + 4 sum over ell>=1 of
    cos(2 pi ell delta)
    (pi ell/eta)^(nu-1/2)
    K_(nu-1/2)(2 pi ell eta)
].
```

Thus

```text
W/epsilon_c = (sigma_LJ/b)^m Z_(m/2)
              - (sigma_LJ/b)^n Z_(n/2).
```

The Bessel-mode count in the implementation controls an exponentially
convergent numerical evaluation of this exact identity. It is not a physical
pair cutoff. Representative values are checked independently against the
direct real-space sum. For `q=2`, where direct truncation converges only as
`1/N`, the Bessel result is checked against the independent hyperbolic closed
form.

The derivative of the same series is used for the registry force. No fitted
sinusoidal Peierls potential and no Taylor polynomial replace the energy.

## Why the normal and registry energies are not summed

`U_infinity(a)` counts pairs in one collinear homogeneously dilated chain.
`W(a,s)` counts cross pairs between two parallel rows. Adding them without one
common atomistic cell and a disjoint pair partition would mix geometries and
can double count interactions. Therefore:

- the normal crack-initiation solver continues to use `U_infinity`;
- the ideal-registry solver uses `W` as its configuration-dependent energy;
- same-row contributions are constant when `b` is prescribed;
- two-way normal--registry coupling is deferred until it is derived from one
  common half-space or interface Hamiltonian.

## Resolved uniaxial loading

For unit loading axis `e`, slip-plane normal `n`, and in-plane slip direction
`d`, the signed Schmid projection is

```text
M = (e dot n)(e dot d),       n dot d = 0,
tau(t) = M sigma(t).
```

The sign is retained because forward and reverse registry translations must
not be merged. The reduced generalized force used in the numerical model is

```text
g(t) = tau(t) A_rep b / epsilon_c.
```

`A_rep` is a crystallographically defined interface area per repeat. It is not
the normal mechanical area `A0`, the correlation area `Ac`, or a FEM element
area unless a separate derivation proves an equality.

## Probability dynamics and nondimensionalization

Let `u=s/b` be unwrapped rather than reduced modulo one. With constant
long-time mobility `M_s`, the physical current is

```text
J_s = -M_s [P partial_s W - tau A_rep P + k_B T partial_s P].
```

The active numerical implementation scales energy by `epsilon_c`, registry by
`b`, and time by

```text
t_s = b^2/(M_s epsilon_c).
```

It then solves

```text
partial_t p = -partial_u j,
j = -[(partial_u w - g)p + beta^(-1) partial_u p],
beta = epsilon_c/(k_B T).
```

The finite-volume face flux uses Chang--Cooper exponential fitting and
backward Euler. This gives nonnegative density, exact reflecting mass balance,
and the correct discrete zero-current Gibbs ratio. The unwrapped domain is a
numerical truncation only; the reported edge probability must be negligible
and the domain must be refined.

The initial state is a metastable Boltzmann distribution conditioned on one
registry basin. A global equilibrium density on the infinite periodic line
would not be normalizable.

## Operational reduced plasticity

Write

```text
u = delta_0 + z + u_tilde,
z in Z,  -1/2 <= u_tilde < 1/2,
```

where `delta_0` is the numerically verified minimum in one period. A well
crossing alone is not declared plastic because reverse crossings are allowed.
The model reports residual reduced plasticity only when, after the applied
resolved force is removed and a relaxation interval is supplied,

```text
mean(u_tilde) -> 0,     but mean(z) != 0.
```

For a separately declared homogenization thickness `h_slip`, the kinematic
mapping is

```text
gamma_p = (b/h_slip) mean(z),
epsilon_p = M gamma_p.
```

This mapping expresses a permanent relative translation. It does not create
dislocation density, hardening, backstress, or multiple-slip interaction.

## Dissipation and symmetry

With probability chemical potential

```text
mu = w - g u + beta^(-1)(ln p + 1),
```

the closed-domain free-energy balance is

```text
dF/dt = -mean(u) dg/dt - integral j^2/p du.
```

The second term is nonpositive; dimensionally it is multiplied by
`epsilon_c/t_s`. Energy is transferred to the eliminated isothermal bath. The
periodic lattice energy itself is conservative.

Because `W(u)` is inversion symmetric about a verified minimum, a symmetric
initial density and exactly antisymmetric zero-mean load have no preferred
translation direction. The tests confirm that six complete symmetric cycles
leave `mean(z)=0.001004` at the stated discretization. A biased tensile-resolved
shear history can instead produce a directed population transfer.

## Current reproducible result

For the explicitly dimensionless demonstration

```text
m=12.19, n=6, a/b=1, sigma_LJ/b=1,
beta=20, Bessel modes=20,
peak g=0.55, ideal maximum |d(W/epsilon_c)/du|=1.0582731,
```

the load is subcritical with respect to deterministic registry instability.
After ramp, hold, unloading, and 18 reduced time units at zero force:

```text
mean(z)                    = 0.4913473
mean intrawell registry    = 1.59e-10
accumulated reduced work   = 0.2591766
maximum edge probability  = 2.93e-10
```

The result demonstrates thermally assisted, finite-time transfer between
lattice wells and residual translation. It is not a prediction of aluminum
yield strain or fatigue life.

Run:

```powershell
py -3 -m simulations.run_registry_plasticity
py -3 -m pytest tests/test_registry_plasticity.py -q
```

Outputs:

- `results/data/registry_plasticity/summary.json`
- `results/data/registry_plasticity/resolved_shear_pulse.csv`
- `results/data/registry_plasticity/symmetric_cycle.csv`
- `results/figures/registry_plasticity/active_registry_plasticity.png`

## Corrected source and remaining source-level corrections

The owner-supplied corrected PDF replaces the earlier repository PDF at
`research/source/slip_lattice_energy_mn_K_derivation_KR_v3_23pages.pdf`; its
SHA-256 is `42C3D5086CA203C76F3DC8213A1718B5121AA1273067738C5B478BCBF12D999D`.
The supplied English symbol index is preserved in `research/source/` and used
by the corrected slip TeX under `libraries/shear/docs/`.

Two project-level corrections remain intentionally stronger than the supplied
index:

1. The project coefficient `epsilon_c` and pair-well depth `epsilon_w` retain
   separate symbols. They differ by `C_mn`.
2. The index describes `A0` as a correlation patch area. This project instead
   keeps normal mechanical area `A0`, statistical correlation area `Ac`, and
   registry interface area `A_rep` distinct unless derived otherwise.

## Physical work still required

- choose the active FCC slip plane and direction for each crystal orientation;
- replace or validate the two-row central-force landscape against a complete
  Al EAM or first-principles generalized stacking-fault surface;
- determine `A_rep`, `M_s`, memory time, and `h_slip` independently;
- test the overdamped reduction against atomistic time scales;
- introduce dislocation storage/hardening only after a microscopic or
  independently calibrated derivation;
- derive a common `(a,s)` Hamiltonian before coupling registry plasticity back
  into normal crack initiation.

