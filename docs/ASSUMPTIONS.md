# Assumptions and Approximations

This file must be updated whenever a new modeling assumption is introduced.

## Current working assumptions

1. The material of interest is high-purity / single-crystal aluminum.
2. The primary reduced structural coordinate is an interatomic-spacing descriptor `a`.
3. The state density `P(a,t)` is interpreted as a thermodynamic-limit population density, not merely a finite histogram.
4. A generalized Lennard-Jones-type interaction may be used as a baseline microscopic pair potential when explicitly stated.
5. Macroscopic affine stretch may be separated from internal/non-affine structural evolution by writing `a = λ x`.

## Controlled approximations that may be tested

- Independent adjacent spacings: `P_k ≈ P^{*k}`. This is not exact and must be validated against correlated simulations.
- Markov closure in spacing space. Allowed only if memory is shown to be negligible on the scale of interest.
- Fokker–Planck truncation. Allowed only after a small-jump/Kramers–Moyal argument.
- Moment closure. Last-resort reduced model, not a starting axiom.

## Forbidden shortcuts unless explicitly justified

- Fitting a Weibull distribution to `P(a,t)` merely because fatigue data often look Weibull-like.
- Inserting an empirical hysteresis loop law and then claiming it was derived from mechanics.
- Introducing damping, damage variables, transition rates, barriers, or kernels only to obtain the desired fatigue curve.
- Treating a single reversible LJ coordinate as a complete fatigue model.
- Confusing instantaneous unstable-tail occupancy with first-passage crack initiation.
