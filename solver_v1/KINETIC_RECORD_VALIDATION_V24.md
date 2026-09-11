# Actual Al kinetic validation: longer record, not an invented clock

## Scope and executed evidence

The user's request is actual yield, fatigue and production physical-time
validation. This step advances the kinetic branch. It does **not** complete
the yield/fatigue branch or change the LJ/Bessel energy, production PDE,
material parameters, mobility, UI or specimen correlation area.

Source: E. Fransson and P. Erhart, *Molecular dynamics trajectory for
crystalline aluminum*, 2023, https://doi.org/10.5281/zenodo.10014454.
The source is 6912-atom periodic Al99 NVT at 300 K, not an experimental
opening trajectory. Actual source input supplies 5 fs steps, 25 fs records,
4.065 Angstrom lattice and 1 ps Nose-Hoover damping. That damping is NOT
adopted as physical collective-coordinate friction.

The new bounded download actually completed 8192 frames / 204.775 ps,
versus 1024 frames / 25.575 ps previously. No full 3.5 GB checksum is
claimed from a prefix. Published input/potential MD5s were checked. The
old 1024-frame coordinate AND time arrays match the new prefix exactly.
Projection SHA256:
`0d4eead55eff67a4ec075152501db6c7cf39032cca04e34b32ce6c13ab7f3cb3`.
The ignored cache is `.cache/kinetic_validation_v24/zenodo_8192/`;
machine-readable results and provenance are `results/kinetic_validation_v24/`.

## Coordinates, units and two independent convergence axes

The coordinate is the adjacent difference of periodic (111) plane means,
with 12 plane classes, 576 atoms per class, three displacement components.
It is neither a rigid-half interface nor the production reduced cell.
Plane classes and conjugate Fourier modes are not independent specimens.

For C(t)=<dq(t)dq(0)^T> [m^2] and K(T)=integral_0^T C(t)dt,

    Q(T)=C0^(-1/2) sym[K(T)] C0^(-1/2) [seconds].

Under the previously derived equilibrium harmonic GLE assumptions,

    Gamma0=kBT C0^(-1) K(infinity) C0^(-1),
    M0=C0 K(infinity)^(-1) C0/(kBT) [m^2/(J s)].

A positive finite cutoff is not a converged integral. We predeclared
prefixes 1024/2048/4096/8192 and lag counts 126/254/510/1022 wherever
each block contains at least two lag windows. Four disjoint quarters
provide another time-window check. Record length and cutoff are separate
axes. All cutoffs, not only favorable ones, are saved. Block differences
are sensitivity measures, not confidence intervals.

At the SAME 3.15 ps cutoff:

| Frames | Smallest Q eigenvalue [ps] | Empirical floor [ps] |
|---:|---:|---:|
|1024|-.0107767|.0382522|
|2048|-.0112812|.0342136|
|4096|-.00849708|.0158149|
|8192|-.00719449|.00796563|

The longer record reduces this sensitivity, but does not yield a positive
resolved integral. At the SAME 8192-frame record:

| Cutoff [ps] | Three Q eigenvalues [ps] | Floor [ps] |
|---:|---|---:|
|3.15|-.00719449, -.00621997, -.00443718|.00796563|
|6.35|.00667574, .00740416, .00903552|.00920581|
|12.75|.000258314, .0136126, .0202151|.0188382|
|25.55|-.00241885, .00124789, .00499399|.0204499|

Some isolated cutoffs become positive above the local floor (e.g. 6.25%
of late cutoffs in the shortest-lag study), but no robust positive plateau
is established. Choosing those cutoffs would manufacture a mobility.
Negative finite integrals are not evidence of physical negative friction.

## Short-time rejection versus static normalization

At 0.15 ps the normalized correlation eigenvalues are approximately
(-.569963, -.500413, -.174055); block-plus-antisymmetry sensitivity is
.182309. The most-negative magnitude is 3.126 times that sensitivity,
NOT a statistical sigma or p-value. This reproduces a resolved oscillatory
response, incompatible with a direct reversible harmonic overdamped
description at that lag. It does not disprove a sufficiently coarse-time
Smoluchowski limit or justify changing the canonical PDE to inertial dynamics.

An independent equal-time check uses the SAME Al99 potential and MD box:

    C_difference,k = 4 sin^2(pi k/N) kBT H(k)^(-1)/N_atoms_per_plane,
    C_local = sum_(k=1)^(N-1) C_difference,k / N.

No fitted atom-count or area factor is used. Actual MD variance relative
to this static harmonic source prediction differs by -1.816% normal,
+5.348% direct110, +9.568% transverse. Each local difference is within
the reported four-block variance sensitivity; this is not every-mode
agreement or a proof that anharmonicity is negligible. In particular,
static normalization success cannot repair temporal oscillations through
a time-unit relabel. Al99 here is only a target/control, never the LJ base.

## What genuinely has a scoped experimental kinetic calibration

Gorman, Wood and Vreeland (1969), https://doi.org/10.1063/1.1657472,
measured edge/mixed line displacement under timed stress pulses. The
original isolated markers were replayed, not replaced by guessed numbers.
Seven calibration/three held-out markers at 23 C reproduce

    velocity/shear = 1.1627233044e-5 m/(Pa s),
    held-out velocity RMS error = 3.501185 m/s.

This is a restricted line-translation fit, not complete population or
line-character calibration and not a/s mobility. The existing Rayleigh
projection in `line_drag_reference.py` additionally requires the matching
validated core profile and local friction assumptions; none is supplied
by simply dividing a line drag by an atomic size.

## Actual yield and fatigue: still separate required validation

No new specimen-yield or fatigue prediction was executed in this step.
The previous source-to-specimen audit showed why an MPa first-source fold
cannot be called proof yield: one declared dilute 1 um source population's
recoverable bow strain is about 2901 times below a .002 axial demand.
That is a mechanism/strain-budget failure, not a unit conversion to fix.
Finite emission, interactions, source population and persistent strain
remain necessary; do not fit a population multiplier to conceal them.

Deschanel et al. (2017), https://doi.org/10.1038/s41598-017-13226-1,
provides actual fatigue validation protocols. Stress full ranges 50/62 MPa
mean amplitudes 25/31 MPa, not 50/62 MPa; reported fracture cycles cannot
be equated with local atomic opening absorption. This step does not claim
to reproduce these lives, PSBs, or detected cracks. Those observable and
microstructure mappings remain uncompleted, not a measured prediction error.

## Decision and next discriminating work

1. Short-lag direct overdamped fit: rejected for these source coordinates.
2. Low-frequency positive converged mobility: unresolved after the new
   eightfold record; do not just repeat positive-window fitting.
3. Production a/s seconds and Hz: remain disabled. The source has physical
   timestamps, but coordinate/PMF equivalence, memory and thermostat checks
   are not thereby validated.
4. Next kinetic work should resolve mode-dependent memory/thermostat
   sensitivity and the actual coordinate projection, rather than another
   arbitrary scalar t0. Yield work separately requires finite-source
   emission and matched specimen strain, not more ideal traction rescaling.

Reproduction uses `run_kinetic_record_refinement` and the existing
`run_plane_covariance_validation`. No new energy law was introduced.
Targeted42 passed (2.27s); solver661 passed (1509.98s); app34 passed,
zero skipped (128.91s); desktop smoke passed (1.483s). The six new tests
verify protocol bookkeeping, not Al physics. The actual physical analyses
above are separate saved results. Full source/projection/test ledger:
`results/kinetic_validation_v24/verification_ledger.json`.
