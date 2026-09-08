# Statistical correlation area: what can and cannot be derived now

This is a theory note, not a numerical calibration or a change to the UI/PDE.
No physical value is assigned to A_c. The existing stable specimen aggregation
and unresolved-local-signal safeguard remain unchanged.

## 1. Three distinct objects

- A_atomic_cell: microscopic FCC cell area for eV/cell to J/m2 conversion.
- A_c: effective statistical area of an independent-equivalent initiation region.
- A_stressed: specimen surface under the modeled loading; mesh face areas may
  quadrature this surface but do not determine statistical independence.

In particular, changing the mesh must not change the physical count of
independent initiation regions. Neither b^2 nor the size of a mesh face is an
automatic calibration of A_c. A_c does not enter LJ, Bessel sums, environmental
density, strain, traction/work conversion, or mobility.

## 2. Exact meaning of the current independence ansatz

For equivalent independent local regions with survival S_l=1-p_l,

    N_eff = A_stressed/A_c,
    log S_spec = (A_stressed/A_c) log(1-p_l),
    p_spec = -expm1[(A_stressed/A_c) log1p(-p_l)].

A continuous effective count is an exponential-survival interpolation, not
proof that a noninteger number of microscopic atoms exists. Define the areal
cumulative hazard

    H_area = -log(1-p_l)/A_c.

Then S_spec=exp(-A_stressed H_area). Data from a single specimen survival
measurement identify H_area, NOT p_l and A_c independently. This is a genuine
identifiability issue; a chosen A_c cannot repair an incorrect local energy
surface or unresolved first-passage probability.

## 3. A correlation area derived from fluctuations

Suppose an external spatial model or atomistic dataset provides a stationary
initiation indicator I(x,t) on the stressed surface, with mean p(t) and
covariance C(r,t). For the surface-average indicator over domain D of area A,

    Var[I_bar] = A^-2 integral_D integral_D C(x-y,t) dx dy.

If correlations are integrable and the domain is much larger than their
range, neglecting boundary corrections gives

    Var[I_bar] ~= A^-1 integral_R2 C(r,t) d^2r
               = p(1-p) A_corr(t)/A,
    A_corr(t) = integral_R2 C(r,t)/[p(1-p)] d^2r.

This matches an effective count A/A_corr for the VARIANCE of independent
Bernoulli variables. It is undefined when p=0 or p=1 without a limiting
procedure. Finite windows and directional correlations must be retained when
the large-domain approximation is invalid. Negative covariance lobes must not
be discarded merely to produce a positive-looking area.

Crucially, matching this variance does NOT prove

    Pr[no initiation anywhere in D] = (1-p)^(A/A_corr).

The probability of no events depends on the full joint distribution, not
only its first two moments. Thus A_corr from a covariance integral is a
candidate constraint on A_c, not a general equality or a calibrated void
probability formula. Higher-order correlations matter for clustered rare
events.

## 4. Where spatial information would have to enter our theory

The current local P(a,s,t) contains no separation vector x-y and therefore
cannot identify a spatial covariance length or area. The atomistic infinite
lattice energy sums neighbors in a reference environment, but this alone
does not specify correlations among distinct stochastic active interfaces.

A future mathematically explicit route would require joint coordinates
q_i=(a_i,s_i), a derived coupling energy among active interfaces, and their
joint deterministic density. Marginalizing this joint model could supply
spatial indicator covariances and test whether independent-region survival
is justified. This is a separate theory extension, not a substitution of
atomic cell area into specimen aggregation. No new length or volume parameter
is introduced by this note.

Another route is mapped atomistic spatial-relaxation/initiation data or
spatially resolved experiments under known loading, with an explicit
definition of the local region. The required observations are joint/spatial
activity, not an S-N lifetime fit. Current local stress/material calibration
and missing collective-coordinate kinetics must be treated independently.

## 5. Resolution and uncertainty must propagate first

For fixed N_eff and a resolved local probability interval [p_low,p_high],
the independence map is monotone, so apply the stable formula to both bounds.
If the interval includes the local numerical floor, the apparent amplification
does not establish physical specimen initiation. A_c uncertainty additionally
changes N_eff; it cannot certify the local signal.

With heterogeneous loading a possible independent-region CONTINUUM ansatz is

    log S_spec = integral_D log[1-p_local(x)]/A_c(x) dA.

Its mesh discretization is a quadrature of physical area and local converged
probabilities. Mesh refinement should converge this integral without changing
A_c. It is NOT currently a validated spatial solver; one global probability
must not be painted onto a mesh as a computed heterogeneous fatigue field.

## 6. Present conclusion

The new finite-wavevector static Hessian offers one possible future connection
to spatial information. For a harmonically stable, explicitly defined crystal
and a specified equilibrium temperature,

    <u(q) u(-q)^T> = k_B T K_phys(q)^(-1)

on the non-translational subspace with boundary/zero modes treated explicitly.
For reduced-coordinate Hessians stored in eV, K_phys=E0 K_reduced/L0^2.
This covariance is of atomic DISPLACEMENTS, not of crack-initiation indicators.
A collective-coordinate projection and a justified joint first-passage theory
are still needed before using it to constrain A_c. It cannot supply a physical
mobility, physical Hz, or an independent-region survival law by itself.

A_c remains an external, statistically identifiable only with additional
spatial information, calibration parameter. The correct next step is to
derive or measure correlations and test the independence approximation, not
to select an attractive numerical area. Physical seconds/Hz are also still
unavailable and are a separate missing kinetic calibration.
