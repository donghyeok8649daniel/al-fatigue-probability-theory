# paper-manuscript

Manuscript-only branch for publication artifacts.

## Current manuscript

- Canonical entry point: `paper/main.tex`
- Current paper body: `paper/fatigue_probability_main.tex`
- Current direction: 1D deterministic normal generalized-LJ chain -> mechanically generated spacing probability -> exact P-u-Theta transport -> tail flux / first passage -> G1, G2 and G4.
- Initial-state supplement: `paper/initial_state_probability_note.tex` derives exact `P0`/`F0` from the microscopic initial state and records the non-uniqueness of future `P` from one-point `P0` or `F0` alone when neighbor ordering is discarded.
- G3 nonzero irreversible dissipation remains physically open.
- Registry `s` is retained only as a non-mainline plasticity extension. No FCC reconstruction is part of the current manuscript.

## Historical manuscript material

`paper/slip_lattice_energy_derivation.tex` and its old symbol indexes are retained as historical/source derivation material. They are no longer input by `paper/main.tex`; their Smoluchowski/direct-registry-drive assumptions are not the active governing theory.

## Branch ownership

Owned here:
- `paper/`: LaTeX manuscript, notation/index, manuscript assets.
- `output/`: manuscript-oriented compiled outputs.
- `research/source/`: source derivation material used to prepare the manuscript.

Not owned here: theory implementation, FEM/simulation code, numerical result generation, desktop packaging, or fatigue-tester implementation.

Validated equations/results are imported conceptually from `theory-core` and `numerical-fem`; publication integration is reviewed through `integration`.

The current tail/first-passage validation is generated on `numerical-fem` by `simulations/verify_probability_tail_first_passage.py`, with results under `results/data/probability_tail_first_passage/`.
The current initial-data sufficiency audit is generated on `numerical-fem` by `simulations/verify_initial_data_sufficiency.py`, with results under `results/data/initial_data_sufficiency/`.
