# Completed v29 low-frequency reference MD

Twelve actual new driven trajectories, total10 ns. Raw trajectories remain
in the local ignored cache; summary.json binds source, potential, restart
and each trajectory by hashes. See solver_v1/LOW_FREQUENCY_FORCING_V29.md
for protocol, derivation, interpretation and limitations.

- responses.csv: all signs, full and half analysis windows, raw harmonics.
- paired_responses.csv: signed-pair response, parity harmonics, null envelopes,
  inverse-response mobility diagnostics (null upper limit means unbounded).
- null_windows.csv: all matched unforced windows; not independent replicas.
- fdt_spectral.csv: equilibrium spectrum sensitivity, not confidence intervals.
- work_checks.csv: clocks, force binding, whole/steady work, endpoints,
  energy residuals, temperature drift and COM diagnostics.
- control_comparison.csv/comparison.json: amplitude/dt/block comparisons.
- loss_comparison.png: blue circles dt2.5 fs, green squares dt1.25 fs;
  bars observed null envelopes, orange equilibrium FDT sensitivity range.

Frequency is reference MD0.2 cycles/ps, NOT a production fatigue-Hz setting.
Strong-probe nonlinear harmonics are detected relative to observed nulls;
combined mobility uncertainty still includes zero drag. Production clock
calibration remains false. Do not treat work-quadrature residual as dissipation.
No kinetic/static production parameters, LJ/Bessel energy or UI were changed.
