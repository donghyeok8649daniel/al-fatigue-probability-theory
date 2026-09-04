# paper-manuscript

Manuscript-only branch for publication artifacts.

## Current manuscript

- Canonical entry point: `paper/main.tex`.
- Canonical paper body: `paper/fatigue_probability_main_v2.tex`.
- Previous `paper/fatigue_probability_main.tex` is preserved unchanged as manuscript history and is not the governing body.
- Active manuscript structure is two-layered:
  1. exact finite generalized-LJ chain -> mechanically generated phase-space/spacing probability -> exact $P$-$u$-$\Theta$ projection;
  2. reduced laboratory closure -> structural/prestress $P_0$ -> quasistatic stable-branch transport -> finite-temperature survivor first passage -> $P_b$, $S$, and $F_{\rm ci}$.
- The exact finite-chain projection is the microscopic reference and is not claimed to close autonomously from one-point $P_0$.
- The reduced lab-scale closure is a controlled 1D normal-instability hypothesis. It uses a high-barrier, harmonic-well, fast-intrawell-equilibration positive-flux transition-state approximation away from the immediate spinodal neighbourhood. If the stable mechanical branch reaches $\lambda_c$, absorption is deterministic. A transmission/recrossing correction remains a validation target.
- The model contains no named lifetime distribution, scalar fatigue-damage law, fitted diffusion coefficient, or imposed permanent normalized-$P$ drift.
- The characteristic cohesive area $A_c$ remains explicit and uncalibrated; specimen correlation area/volume scaling is deliberately deferred until local calibration.
- In the strict quasistatic/fast-equilibration regime, integrated hazard per cycle scales as $1/f$; temperature dependence is activated through the mechanically derived effective-potential climb. These are falsification signatures, not fitted rules.
- Nonzero irreversible thermodynamic G3 remains a separate physical problem.
- Registry $s$ is retained only as a non-mainline plasticity/defect extension. No FCC reconstruction is part of the active manuscript.

## Symbols

- `paper/symbol_index_active_v2.tex`: canonical English symbol/computation index included by the v2 manuscript.
- `paper/symbol_index_active_v2_ko.tex`: Korean companion index for the same active model.

## Supporting manuscript notes

- `paper/initial_state_probability_note.tex`: exact $P_0/F_0$ construction from the microscopic initial state and loss of future uniqueness when neighbour ordering is discarded.
- `paper/local_correlation_hierarchy_note.tex`: why finite local correlation levels do not provide an exact autonomous projected closure in general.
- `paper/final_reduced_survival_closure_note.tex`: modular derivation of the reduced survivor equation and periodic-cycle survival law.
- `paper/peak_hazard_asymptotic_note.tex`: peak-dominated rare-event cycle-hazard approximation and temperature-slope identifiability relation for $A_c$.

## Historical manuscript material

`paper/slip_lattice_energy_derivation.tex`, `paper/fatigue_probability_main.tex`, and older symbol indexes are retained as historical/source derivation material. They are not the current canonical governing manuscript.

## Branch ownership

Owned here:
- `paper/`: LaTeX manuscript, notation/index, manuscript assets.
- `output/`: manuscript-oriented compiled outputs.
- `research/source/`: source derivation material used to prepare the manuscript.

Not owned here: theory implementation, FEM/simulation code, numerical result generation, desktop packaging, or fatigue-tester implementation.

Validated equations/results are imported conceptually from `theory-core` and `numerical-fem`; publication integration is reviewed through `integration`.

Key numerical provenance on `numerical-fem` includes:
- `simulations/verify_initial_data_sufficiency.py`;
- `simulations/verify_pair_state_insufficiency.py`;
- `simulations/verify_probability_tail_first_passage.py`;
- `simulations/audit_local_traction_lab_timescale.py`;
- `simulations/audit_collective_mode_lab_timescale.py`;
- `simulations/audit_periodic_p_first_passage_identifiability.py`;
- `simulations/audit_final_reduced_closure_sensitivity.py`;
- `simulations/verify_final_reduced_general_p0.py`;
- `simulations/verify_peak_hazard_asymptotic.py`.
