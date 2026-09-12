# Completed v31 weak-response reference MD

Ten completed 2 ns NVE records: two new initializations, each with null,
normal +/- and direct110 registry +/- forcing. The source Al99 potential is
target-only; the production LJ/Bessel energy and physical-time gate are unchanged.

- `responses.csv`: each sign, full record and halves.
- `pairs.csv`: signed partners, phase, null comparison, parity harmonics.
- `null.csv`: actual matching unforced lock-in windows.
- `checks.csv`: raw-source hashes, actual temperatures, native work residuals.
  Unforced work diagnostics are unavailable/blank, not zero measured loss.
- `fdt.csv`: independent unforced loss estimates near the prescribed drive.
- `aggregate.csv`: two-initialization response and empirical inverse-error sets.
  A null upper mobility is unbounded/unresolved, not zero mobility.
- `generator_closure`: predictions from every unforced full C/K band, no refit.

Both coordinate aggregates admit zero drag within the empirical complex-error
sets. These sets are not confidence intervals. The finite-frequency point
mobilities must not be copied into production kinetic calibration. Neither
physical local a/s Hz nor actual Al yield/fatigue is certified by this report.

The CSV/schema and console-JSON export defects found in the first postprocessing
attempts were repaired without changing any MD data. A successful fresh replay
reproduced the seven top-level result CSV/JSON files byte-for-byte. All ten
trajectories remain hash-bound in the summary; raw trajectories and incomplete
export attempts are preserved locally, not included in Git.

Derivation, numerical results and limitations:
`solver_v1/WEAK_KINETICS_AND_MATERIAL_V31.md`.
