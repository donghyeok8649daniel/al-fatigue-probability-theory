# Open Problems

## Milestone 1 — Mechanics-derived hysteresis

Given a prescribed cyclic stress

\[
\sigma(t)=\sigma_m+\sigma_a\sin\omega t,
\]

derive an internal evolution law from microscopic mechanics such that loading and unloading do not simply retrace the same structural state, while conserving probability and satisfying the correct energy balance.

Success condition:

\[
A_H=\oint\sigma\,d\epsilon>0
\]

without inserting an empirical hysteresis law.

## Milestone 2 — Secular fatigue accumulation

A periodic hysteresis loop alone is internal friction, not fatigue. Derive conditions for

\[
P_{N+1}\neq P_N
\]

or for an equivalent slow internal state to evolve cycle by cycle.

## Milestone 3 — Crack initiation

Formulate initiation as a mechanical stability loss or first-passage event. Candidate formulations include an absorbing boundary in an enlarged state space or a distribution-level stability criterion.

## Central closure problem

Determine the minimum state required for a mechanically closed model:

- `P(a,t)` only?
- phase-space density `F(a,c,t)`?
- spacing-correlation hierarchy?
- joint structural state such as `P(a,s,t)` where `s` is a non-affine/slip coordinate?

The preferred answer is the smallest state that can be derived from mechanics without hiding essential memory or irreversibility in fitted constitutive terms.

## Numerical falsification tests

Any candidate model should pass at least:

1. zero loading -> zero hysteresis;
2. perfectly reversible conservative limit -> zero hysteresis;
3. probability normalization preserved;
4. non-negative density;
5. dimensional consistency;
6. energy balance;
7. uniform-lattice limit recovers the baseline lattice energy;
8. nonzero fatigue accumulation is not created by numerical diffusion alone.
