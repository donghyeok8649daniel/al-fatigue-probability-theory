# v18: source-informed line kinetics and fatigue validation boundaries

Executed research results, **not a production physical-time calibration**.
Source values, line reference and actual specimen endpoints are separated.
No LJ/static parameters, PDE mobilities, A_c or fatigue-life targets were fitted
to obtain a desired result. The v17 material candidate remains rejected.

## Files

- `line_calibration/velocity_points.csv`: all10 native-image isolated markers,
  exact CGS-to-SI conversion, reading bands,7fit/3held-out labels and residuals.
- `line_calibration/calibration.json`: slope, scatter, leave-one-out values,
  original/data hashes and explicit scope. Original image replay was performed.
- `line_dynamics/stress_scenarios.csv`:0.01/0.1/1/4MPa shear, three hypothetical
  pin spans, linear versus nonlinear existing anisotropic line graph.
- `line_dynamics/frequency_response.csv`: physical LINE-reference frequencies,
  declared b/drag/stiffness/geometry and infinite-series truncation bounds.
- `line_dynamics/spatial_refinement.csv`: independent16/32/64/128-segment
  tridiagonal solutions and work/dissipation balance.
- `line_dynamics/time_refinement.csv`:64/128/256/512steps per cycle plus a
  zero-stress hold of at least20 line relaxation times.
- `line_dynamics/cyclic_unload_hold.csv`: finest actually integrated time
  history; no permanent slip or opening probability is represented.
- `line_dynamics/fatigue_reference.csv`: seven source specimen groups, ranges,
  amplitudes, loading rates and final-fracture cycles. This is not an S-N fit.
- `line_dynamics/scope.json`: unchanged material provenance, mismatched
  temperature/line character, hypothetical pins, blocked production/fatigue gates.

Each result-producing module requires a fresh output directory and refuses
overwrite. See [derivation and commands](../../solver_v1/LINE_KINETICS_AND_FATIGUE_VALIDATION.md).
Raw source PDFs/images and test logs remain local cache artifacts; the compact
numeric source transcriptions retain DOI, units, conditions and hashes.

The diagram's MHz/GHz rolloff points validate the declared constant-drag
overdamped line equation, not real Al kinetics at those frequencies. The
0.1--100Hz comparison predicts near-equilibrium reversible bowing under the
stated geometric/elastic assumptions. It must never relabel a production
P(a,s,t) axis in seconds or Hz. Fatigue Nf, AE growth markers and atomic opening
absorption are different observables; no prediction error between them is
reported as a valid calibration.
