# v31 full-plane generator diagnostics

These are finite-band, matched reference-MD coordinate calculations, not
production a/s mobility calibrations. See
`solver_v1/WEAK_KINETICS_AND_MATERIAL_V31.md` for the derivation.

- `N6_rank_checked`: existing long N6 unforced record, full zero-sum inverse.
- `N12_condition_checked`: existing shorter N12 record, with rank/conditioning
  rejection. Use this report, not the superseded raw-inverse experiments.
- `seed35461_new`, `seed49277_new`: the two new v31 2 ns unforced records.
- Other N6/N12 folders are preserved intermediate diagnostics. Particularly,
  very large numbers in early N12 inversions are not accepted mobilities.

Every final report binds its input trajectory by SHA256. Raw trajectories are
local cache, not included in Git. C/K matrices retain all three components of
the independent zero-sum plane space. Scalar projection and state elimination
are not interchangeable. Finite-band sensitivity ranges are not confidence
intervals; a numerically invertible matrix need not be statistically converged.

`results/weak_replica_v31/generator_closure` compares predictions from these
unforced matrices with the separate signed-forcing runs, without response
refitting or selecting a preferred frequency band.
