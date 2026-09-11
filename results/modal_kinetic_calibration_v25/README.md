# v25 kinetic data: scope and authority

These are reference-coordinate estimates, not a calibrated production PDE
clock. No result authorizes copying a modal or energy-normalized plane
coefficient into the local a/s model at unchanged thermal noise.

- `first_half_fit`: actual public8192-frame record, deterministic DHO fit on
  first half and prediction of second half; second-half refits are sensitivity.
- `mobility_projection`: all source DHO windows/halves in SI units.
- `experimental_linewidth`: original experimental-only vector-marker extraction.
- `phonon_reference_with_dft`: same four experimental markers plus two separate
  DFT MD markers from the same primary paper. DFT is NOT experimental data.
- `nve_*` / `nvt_*`: completed independent reference MD modal analyses.
- `interim_control_comparison`: historical first-three-case comparison.
- `control_comparison`: completed FIVE1ns control cases; supersedes the interim
  comparison, not the later long-record study.
- `spectrum_*`: direct finite-band spectra; `lowlimit` is a declared follow-up
  to the initial three-band study. A positive PSD is not a plateau certificate.
- `mode_spectrum_*`: phase-pooled mode-dependent response, not a local constant.
- `leakage_nve_1000ps_dt5`: undamped-mode taper leakage control.
  This uses the FIVE-band `spectrum_lowlimit_nve_1000ps_dt5`, not the initial
  three-band spectrum. A final hash-bound rerun exactly reproduced all ten
  band/taper records.
- `spectral_comparison`: full-record initial-band/DHO comparison, all five1ns
  cases. Plotted ranges are sensitivity, not confidence intervals.
- `thermal_channel_audit`: first kinetic-energy analysis, retained as history.
- **`thermal_channel_mass_checked` supersedes the thermal energy normalization**:
  EAM pair_coeff resets the actual engine mass to26.982, whereas the original
  postprocessing used the earlier input26.98. The original trajectory and
  stored observables are preserved; analysis explicitly applies26.982/26.98
  to the energy channels. Correlation fractions are unchanged by this common
  scaling. This never changed the MD forces/timestamps or inferred M=C0/(kBT tau).
- `sampling_interval_check`: same225ps record at5fs versus25fs sampling;
  largest diagonal finite-band mobility change0.0241%, exact old-record prefix.
- `long_record*_prefix_*`: nested1/2/5ns studies; not independent replicas.
  Keep invalid/unresolved lowest bands in the JSON instead of silently omitting
  failed estimates. No output is a zero-frequency certificate.

The engine mass was directly probed: restart26.982; after `mass 1 26.98`,26.98;
after the hash-bound Al99 `pair_coeff`,26.982. The remaining6.3e-8 energy-unit
normalization difference is reported rather than hidden. Future MD runner
queries/stores the engine mass and uses it in kinetic-energy postprocessing.

Both5ns long-record studies and their analyses completed.
`long_record_comparison` combines all prefixes/bands/time steps; its plot was
visually inspected. `all_control_comparison` combines all seven main records
and the experimental/DFT linewidth comparison. Main output sensitivity ranges
are not confidence intervals or production mobility certificates.

Final regression:12 targeted passed1.42s;673 solver passed906.20s;
34 app passed97.53s with no skips; desktop smoke passed.
The original LJ/Bessel, static calibration and production kinetic JSON were
not changed. Research MD is reference data, not a new production backend.

The detailed methods and current decision are in
`solver_v1/MODAL_KINETIC_CALIBRATION_V25.md` and `CURRENT_WORK_HANDOFF.md`.
