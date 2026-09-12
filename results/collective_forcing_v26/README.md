# v26 matched conjugate-force study

Eight completed400ps reference Al99 runs. This is not the production LJ/PDE
backend, an Al fatigue prediction or a physical-clock calibration.

- `responses.csv`: signed runs, amplitudes, temporal blocks; no clipping.
- `paired_responses.csv`: signed-pair mean and disagreement.
- `cross_responses.csv`: all measured polarizations for both applied axes.
- `fdt_cutoffs.csv`: covariance prediction and finite-endpoint correction.
- `fdt_spectral_loss.csv`: alternate spectral FDT loss and explicitly
  assumption-dependent precision planning; not a confidence interval.
- `work_checks.csv`: raw work/energy/COM diagnostics and trajectory hashes.
- `summary.json`: source hash, predeclared protocol and all summarized data.
- `response_comparison.png`: visually inspected; ranges are sensitivity.

In-phase agreement0.04–1.19% supports this coordinate/force normalization.
ALL paired quadrature block ranges cross zero; dissipative mobility remains
unresolved. Work residuals can exceed the predicted small dissipated work.
The separate `collective_forcing_phase_benchmark_v26` is a dimensionless
synthetic ENGINE/analysis test; its parameters must never be called Al data.

Full derivation and interpretation: `solver_v1/CONJUGATE_RESPONSE_VALIDATION_V26.md`.
Tests:676 solver,34 app(no skips),18 targeted passed; desktop smoke passed.
