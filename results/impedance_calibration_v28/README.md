# v28 finite-band impedance mobility calibration

Executed new deterministic WLS on completed v25 spectra, compared against
actual completed v26/v27 driven MD. No new MD trajectory in this stage.
Full derivation: `solver_v1/IMPEDANCE_MOBILITY_CALIBRATION_V28.md`.

- Normal scalar candidate: **7.246372352e10 m²/(J s)**.
- Direct110 scalar candidate: **2.180655810e11 m²/(J s)**.
- Ratio3.0093; these are periodic plane-gap candidates, not production Ma/Ms.
- `calibration_candidate.json`: sources, normalized fits and closed clock gates.
- `low_band_residuals.csv`: all784 observations; four lower bands excluded from loss.
- `unavailable_source_bands.csv`:196 original bandwidth/bin-resolution failures with reasons.
- `fit_control_sensitivity.csv`:168 deterministic control fits, not independent samples.
- `forced_impedance.csv`:60 inverse-disk sensitivity records, including unresolved loss.
- `overdamped_predictions.csv`: actual forced-response validation, including failures.
- `raw_verification_dt*.json`: both5 ns raw records independently reconstructed.
- `mobility_controls.png`: observed spans, **not confidence intervals**.

Blank inverse upper bounds mean unavailable/unbounded, not zero mobility.
Null-plus-control uncertainty permits zero dissipative drag in all four
high-frequency full-force cases. Low-frequency forced loss remains unresolved.
High-frequency storage enhancement is inconsistent with the same-static-PMF
pure overdamped response; this does not disprove a lower-frequency limit.
No physical production seconds/Hz or local-cell mobility is certified.

Executed checks: targeted20 passed; full solver693 passed(769.26s);
app34 passed(143.03s), desktop smoke passed. No existing test removed/skipped.
