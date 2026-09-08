# theory-core-v1-solver

Integration branch for the finalized Theory Core v1 manuscript and its executable proof-of-principle solver.

## Canonical manuscript

- Entry point: `paper/main.tex`
- Canonical body: `paper/fatigue_probability_main_v4.tex`
- Historical v1/v2/v3 manuscript bodies are preserved.
- v3 established the no-go result for strict pure-normal cycle accumulation.
- v4 adds the required slow correlated configurational state and promotes the theory to:
  1. interacting LJ/configurational energy;
  2. many-body correlated Smoluchowski dynamics without product closure;
  3. plastic memory through `s_i = b n_i + xi_i`;
  4. separate plastic and normal-opening stability criteria;
  5. mechanically derived opening saddle and first-passage crack initiation.

The macroscopic target remains stress/strain hysteresis, residual plastic response, irreversible evolution, and crack-initiation probability. Microscopic coordinates are hidden constitutive states rather than proposed experimental observables.

## Solver

`solver_v1/` contains the executable dimensionless mechanism solver. The
canonical probability calculation is the direct conservative
Smoluchowski/Fokker--Planck finite-volume PDE; Monte Carlo trajectory counts are
not the production probability source.

- `model.py`: Poisson/Bessel infinite-lattice energy, relaxed axial calibration,
  strain bridge, well index, and opening mechanics.
- `probability_pde_2d.py`: direct N=1 probability-PDE reference used by the UI.
- `probability_pde_4d.py`: dense correlated N=2 reference.
- `probability_tt_6d.py`: N=3 initial-Gibbs TT compression prototype only; it
  is not yet a time-dependent 6D solver.
- `solver.py`: historical stochastic cross-check, not the production estimator.

## Desktop application

The desktop entry point is:

```bash
python -m pip install -r solver_v1/requirements.txt
python -m app.desktop_ui
```

The application always runs the probability theory. FVM/FEM is a separate
spatial-display selection, not an alternative to the theory. Its testable data
path is `app/solver_adapter.py`:

`physical stress -> sigma/E -> relaxed axial kappa -> N=1 probability PDE -> result fields`.

The desktop setup now identifies the active energy surface explicitly. It can
run the original `TwoRowLJ` reference, a clearly labeled hypothetical analytic
LJ--EAM sensitivity surface, or the stored Al-target best-feasible hybrid. Each
choice uses its own equilibrium and relaxed axial kappa, and the exact model ID,
parameter source, and calibration status remain attached to results. The
default remains the `TwoRowLJ` reference; it is never presented as calibrated
aluminum.

The conditional Gibbs initial density is prepared at the applied load at
`t=0`, which is the entered mean stress for the current sinusoid. The UI plots
normal-opening, intrawell-registry, well-index-plastic, and total axial strain
directly from the PDE result. It also exposes survival, cumulative initiation,
first-passage flux, intact/absorbed mass, mass-balance residual, and positivity
diagnostics without reconstructing them from a deterministic trajectory count.

Solver time defaults to dimensionless model time. Without a validated kinetic
calibration, the UI frequency is cycles per model-time unit, not hertz, and is
not mapped to display seconds. A physical conversion requires independently
justified reduced-coordinate mobilities or friction. Current finite-mobility
and fast-normal-coordinate diagnostics are in
[`solver_v1/DYNAMICS_FAST_A_DIAGNOSTICS.md`](solver_v1/DYNAMICS_FAST_A_DIAGNOSTICS.md).
The exact conversion framework and kinetic-calibration gate are documented in
[`solver_v1/PHYSICAL_TIME_AND_MOBILITY.md`](solver_v1/PHYSICAL_TIME_AND_MOBILITY.md).
The subsequent finite-temperature fast-$a$ derivation and its restricted
slow-$s$ validation are documented in
[`solver_v1/FAST_A_REDUCTION.md`](solver_v1/FAST_A_REDUCTION.md).
The exact surviving-ensemble strain decomposition, registry-flux bookkeeping,
current refinement study, unload/hold result, and active-UI-model audit are in
[`solver_v1/PLASTICITY_AND_STRAIN_AUDIT.md`](solver_v1/PLASTICITY_AND_STRAIN_AUDIT.md).
The separate full FCC(111) infinite-plane-stack static reference and its
direct-sum validation are documented in
[`solver_v1/FCC111_FULL_STACK_DERIVATION.md`](solver_v1/FCC111_FULL_STACK_DERIVATION.md).
It remains a research reference and is not the desktop default.
The subsequent full-FCC static refit and localized active-interface audit are
in [`solver_v1/ALUMINUM_FULL_FCC_CALIBRATION.md`](solver_v1/ALUMINUM_FULL_FCC_CALIBRATION.md)
and [`solver_v1/FCC111_ACTIVE_INTERFACE_DERIVATION.md`](solver_v1/FCC111_ACTIVE_INTERFACE_DERIVATION.md).
Their negative held-out validation result is retained explicitly; the active
interface is not connected to the production PDE.
The independent-elastic-mode audit and actual deterministic optimization run
with `python -m solver_v1.run_full_fcc_calibration_audit`; current evidence is
saved under the two `audited_v2/` result directories. See
[`solver_v1/SOLVER_VALIDATION_GATES.md`](solver_v1/SOLVER_VALIDATION_GATES.md)
for physical readiness gates and the deferred specimen/meshing UI workflow.

Run from repository root:

```bash
python -m pip install -r solver_v1/requirements.txt
python -m solver_v1.run_demo
```

## Historical mechanism-screening result

One dimensionless parameter set produces the required ordering:

| Fmax | first passage | max mean abs well index | min mean opening barrier |
|---:|---:|---:|---:|
| 2.5 | 0 | 0 | 0.3549 |
| 3.2 | 0 | 0.3021 | 0.1062 |
| 3.4 | 0.09375 | 3.1839 | 0.0609 |
| 3.6 | 0.84375 | 6.8000 | 0.0271 |

Interpretation: configurational well crossing can occur before crack initiation; further rearrangement lowers the normal opening barrier and increases first-passage probability.

## Scientific scope

This branch does **not** claim calibrated pure-Al fatigue life. The current LJ parameters, mobility values, thermal scale, and axial projection coefficient are mechanism-screening parameters. Quantitative aluminum work still requires an Al-specific EAM/MEAM or validated energy landscape, mobility calibration, a derived/configurational geometry bridge, and later characteristic correlation length/area for specimen-scale aggregation.

The fixed v1 architecture is therefore:

`interaction energy -> correlated state dynamics -> plastic/configurational memory -> opening stability/barrier -> first passage`.
