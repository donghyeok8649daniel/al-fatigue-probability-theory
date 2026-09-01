# theory-core

Physical and mathematical source of truth for the aluminum fatigue probability model.

Owned here:
- `theory/`: generalized LJ lattice energy, probability evolution, slip/plasticity, crack-initiation mathematics.
- `docs/`: derivations, assumptions, governing-equation notes, symbol/variable definitions, open theory problems.
- `libraries/`: theory-side lattice/reference libraries.

Not owned here: FEM/UI simulations, generated numerical results, desktop packaging, manuscript sources, or fatigue-tester hardware/firmware.

Validated theory flows downstream to `numerical-fem`, then through `integration`.
