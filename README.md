# Al Fatigue Probability Theory

Mechanics-first research framework for fatigue crack initiation in high-purity / single-crystal aluminum.

## Research goal

Derive cyclic hysteresis, irreversible fatigue accumulation, and crack initiation from microscopic mechanics and distribution-valued state variables, while minimizing empirical fatigue laws and uncontrolled fitting.

## Core state

The primary reduced state is the thermodynamic-limit spacing density

\[
P(a,t)=\lim_{N\to\infty}\frac{1}{N}\sum_{i=1}^{N}\delta(a-a_i(t)),
\]

where `a_i(t)` is a local interatomic-spacing descriptor. The long-term goal is to derive the evolution of `P(a,t)` from mechanics rather than prescribe a probability family.

## Current backbone

- Microscopic mechanics: Newton/Hamiltonian dynamics
- Deterministic lattice baseline: generalized Lennard-Jones interaction
- Distribution energy: pair-distance distributions and correlation hierarchy
- Exact density conservation: `∂P/∂t + ∂(Pv)/∂a = 0`
- Hysteresis target: obtain `∮σ dε > 0` without inserting an empirical hysteresis law
- Fatigue target: obtain secular cycle-to-cycle evolution `P_{N+1} != P_N`
- Crack initiation: formulate as loss of mechanical stability / first-passage into an unstable state

## Repository structure

- `docs/` — theory notes, assumptions, derivations, failed approaches
- `theory/` — analytic and symbolic model code
- `simulations/` — numerical experiments
- `tests/` — conservation, limiting-case, and falsification tests
- `firmware/` — eventual C/C++ reduced model and fatigue-tester integration
- `results/` — generated summaries and figures; raw experimental data should remain outside normal Git history

## Research rule

Every important result should be classified as one of:

- **EXACT / IDENTITY**
- **DEFINITION**
- **ASSUMPTION**
- **CONTROLLED APPROXIMATION**
- **EMPIRICAL INPUT**

A model that only reproduces a known fatigue curve by fitting is not considered a successful derivation.
