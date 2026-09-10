# Normal-environment response: v17 scoped research

The prior v16 static energy improvement is preserved. No production model,
probability equation, mobility, clock, specimen area or historical fit changes.
This study is about the remaining normal force/curvature error, not inducing
slip, fitting yield, or claiming an experimental aluminum fatigue result.

Source provenance is unchanged from [matched-condition calibration](MATCHED_INTERFACE_CALIBRATION.md):
[Mishin et al., Phys. Rev. B59,3393 (1999)](https://doi.org/10.1103/PhysRevB.59.3393),
NIST Al99 source SHA256
`60c8a085be79d273324ab421f5b1447578fef55c1acfc6492c0999f15ee8a284`.
These are like-for-like computed 0K source-potential targets, NOT newly measured
experimental material constants. The separate previously processed MD source
is [Fransson/Erhart dataset](https://doi.org/10.5281/zenodo.10014454):1024 actually
processed frames,300K,576atoms/plane,12plane classes. Its source metadata retains
the compressed-prefix checksum and explicitly does not claim a checksum of the
entire unprocessed trajectory. Neither source is a new production potential.

## 1. Diagnose before adding a parameter

An actual fixed-shape coefficient audit is in
`results/fcc111_active_interface/normal_response_v17/initial_diagnosis/`.
At the pristine interface the last v16 candidate has Haa=27.9164 eV/L0²
versus source 19.6609. Its rank1 term alone contributes 13.3095, while the
rank2/Eg terms contribute 3.4634/3.5247. At a=1.62h the rank1 curvature is
-1.62486, total -1.55444 versus source -0.155829 eV/L0².

All eight actual coefficient controls are preserved. Forcing perfect Haa
exactly on the last fixed-pair section increases original fit loss from
596.14 to 1373.60 and reduces Hxx from 4.7624 to 2.2565 (source 4.3936).
Setting D1=0 increases loss to 5679.94. Omitting the finite-q constraints
does not change the original optimum: this particular tradeoff is not caused
by an active spectral constraint. These are fixed-shape convex results, not
a proof that the entire old nonlinear family is impossible. New same-family
shape refinement is required as a control for any extension.

## 2. Smallest range hypothesis

The current radial constraint ties ranks 1 and 3 to k_odd. They describe
different tensor environments; sharing their decay was a modeling restriction,
not a crystallographic identity. Test ONE additional microscopic radial shape:

    Q1_i(k1) = sum_j w_k1(r_ij) r_vec_ij,
    w_k(r) = C_k exp(-k r),
    E1 = D1 sum_i |Q1_i(k1)|².

C_k is the existing normalized infinite-FCC density convention at the fixed
reference geometry. Its amplitude is fixed during deformation. D1 and k1
are independent material parameters, never functions of applied stress.
k1=k_odd EXACTLY recovers the previous family. No new amplitude, empirical
stress law, cutoff, spatial correlation area or arbitrary specimen length
is introduced. Rank3, both rank2 tensors, their saturation and scalar EAM
retain their original definitions. In particular, k1 does NOT replace k_odd
inside the radial Eg combination.

The canonical neighbor sums and all coordinate jets remain the existing
infinite plane Poisson/Bessel expressions. Each atom sums its complete
environment before the norm is taken. A rigid active interface contributes

    W1_int = 2 D1 sum_depth |sum_cross_layers Delta Q1_plane|²,

not the norm of a global crystal sum, nor a sum of separate plane norms.
Per-interface-cell energy is eV, q is in the fixed L0=4.05/sqrt(2) Angstrom,
and A_atomic_cell=sqrt(3)L0²/2 is used for J/m² and traction only. A_c is absent.

## 3. Analytic derivatives and harmonic implications

For a bond R and w=C exp(-k1 r),

    d_beta(w R_alpha) = w[delta_alpha,beta-k1 R_alpha R_beta/r].

For I=Q1:Q1, the exact coordinate derivatives remain

    I_p = 2 Q1:Q1_p,
    I_pq = 2 Q1_p:Q1_q+2 Q1:Q1_pq.

An affine deformation preserves inversion partners in a perfect Bravais
lattice; hence Q1=0 for every homogeneous strain. Its energy contributes
neither cohesion nor homogeneous elastic constants, for ANY k1. This does
NOT make k1 redundant in a nonuniform displacement or active interface.

At a cubic Bloch wavevector q, define

    L1(q) = sum_R [cos(q.R)-1] d_R(w R).
    H1(q) = 2 D1 L1(q) L1(q)^T.

The sine sum vanishes by inversion. Since L1=O(|q|²), H1=O(|q|⁴); changing
k1 changes finite-q response even though the homogeneous elastic tangent
is unchanged. The complete total H(q), including all signed coefficients,
must be recertified. Its omitted-neighbor bound uses the same analytic FCC
Voronoi integral with k1 substituted ONLY in the rank1 column. Cubic
dilations keep material kernels and gauges fixed.

Independent validation must use displaced per-atom direct energies with
the NEW k1, not an old tied-range direct evaluator. Finite radius is only
an independent validation/admissibility method, never the canonical energy.

## 4. Calibration experiment and decision rules

The inspected v16 validation observations become development data. New
off-grid/off-path states are declared before v17 optimization and excluded
from the loss. Existing 10% energy/curvature discrepancy scales and absolute
250 MPa force scale remain visible; they are not measured error bars.
Near-perfect normal jets add information about the observed stiffness failure.

Compare tied and separate k1 on the SAME dataset and exact bulk constraints.
If a fixed positive LJ pair is retained, label it an inherited CONTROL, not
two new Al measurements. Repeat the exact-constraint sensitivity/SVD with
the new log(k1) column, and keep every stopping reason and failed profile.

Required tests: nested limit, coefficient replay, analytic gradients/Hessians,
registry periodicity, per-site direct sums, independent nonzero-k1 Bloch
energies, radius/tail and dilation checks, stationary Morse indices, opening
extrema, static normal/shear/mixed loading and independent source residuals.
Matching one curvature or lowering training loss is not material adoption.
Finite q is a static Hessian, NOT Hz. Model-time PDE/finite-source/core,
actual yield, physical mobilities and specimen UI gates remain separate.

The vector environment also has an exact radial-gradient interpretation.
For psi_k(r)=C_k exp(-kr)(r/k+1/k²),

    grad_R psi_k(r) = -C_k exp(-kr) R,
    Q1_i = -sum_j grad_R psi_k(r_ij).

Thus the invariant measures a coherent local inversion/environment imbalance.
The identity does NOT add psi as a pair potential: the pair is still LJ,
and the energy is the stated per-site norm, including cross-neighbor terms.
A newly exposed interface naturally produces a large normal imbalance; this
explains why its coefficient contributes strongly to Haa. It does not justify
tuning that coefficient to a desired macroscopic yield stress.

Executed results and their limited scope are recorded below. A completed run
means the declared calculation finished, not that material acceptance passed.

## 5. Declared numerical controls (not calibration data)

Both new fits retain the positive u/v from v16 `saturated`, with a SHA256
binding to that run. This is a matched embedding/range experiment, not a
claim that u/v were independently measured. The tied-range start uses the
last v16 positive-pair shape. Independent k1 starts are 3 and 9; both and
all shape bounds are recorded before the new objective is evaluated.
k1 is a reduced inverse microscopic distance, with physical inverse range
k1/L0. Its [2,12] search bounds match the prior odd-range numerical window;
they are not measurements. A bound optimum is not an identified material
parameter. The dimensionless alpha normalization is fixed under strain.

There are 104 fitted interface observations, 5 exact bulk rows, 2 exact
inherited-pair CONTROL rows, and 48 new held-out interface observations.
Each fit uses a deterministic bounded simplex, maximum96 profile evaluations
per start. Every accepted/rejected coefficient profile and termination is
retained. Finite-q constraints include undeformed and .99, 1.0037037037037038,
1.01 fixed-material dilations; no density gauge is reset by deformation.

A third, explicitly nested start k1=k_odd uses bounded two-point trust-region
least squares (maximum12 objective iterations, numerical shape probes recorded
separately). It was declared while the simplex searches were in progress,
before inspecting the new validation errors. This is deterministic numerical
continuation, not a new independent physical target or an equal-budget claim.

After completing a shape search, also remove the inherited u/v CONTROL and
re-solve coefficients at that exact shape. This separates a restrictive pair
control from a true same-shape material tradeoff. It is NOT a new nonlinear
optimization. If the least-squares optimum is a zero-LJ closure, retain it as
a diagnostic and refuse material instantiation; do not add a tiny positive
pair floor. A fixed-pair experiment alone cannot establish failure of the
full LJ-parameter family.

The first completed nested run reached the previously declared k_even=12
boundary. Two subsequent fixed-shape CONTROL profiles at k_even15 and18
check whether that numerical search boundary alone explains failure. Other
shapes and inherited positive pair are held at the nested result; no extra
material term is introduced. These are not completed nonlinear searches,
not new measurements, and not proof of the optimal decay. Their outputs
are separate `even_bound_probe`, with the enlarged bound explicitly saved.

The independent `run_even_environment_validation` deliberately reuses several
historical validation locations for derivative/direct-sum reproducibility.
Their legacy labels containing `heldout` do NOT make them blind in v17.
Only the 48 `v17_...` rows marked heldout in the new definition are omitted
from the current fit. Once inspected, these are no longer prospective data
for a future v18 fit. The existing 300 K MD covariance is a separate,
already-inspected validation source, not a source of real a/s mobility.

## 6. Independent source-normal audit (completed)

`run_source_normal_audit` reconstructs the geometric neighbor-cutoff crossings
of the external source and checks its energy/force/Haa stencils. No source or
target is changed. Source nominal cutoff6.28721Å, last tabulated point
6.286581279Å, spacing0.000628721Å. The evaluator uses the latter edge (no
extrapolation of the table); pair force/curvature there are 2.55e-8eV/Å and
-8.12e-5eV/Å². The density and its derivatives are essentially zero there.

Five layer-shell crossings in h<=a<=4h occur at a/h=1.28656257, 1.59390987,
1.93089834, 2.28656257, 2.59390987. All 31 source normal-state FD stencils
avoid them. At step1e-5L0, the maximum source Haa FD discrepancy is9.07e-5
eV/L0²; at the disputed1.62h it is6.21e-9. Source Haa there is-0.155829,
not the candidate's approximately-1.55. Near geometric crossings the two-sided
Haa differences at epsilon1e-5h are at most0.00262eV/L0², which include normal
smooth variation as well as the small table-edge discontinuity.

Therefore these source numerical effects do NOT explain the observed order-one
candidate curvature discrepancy. A finite-cutoff source is not the infinite
LJ model, but that difference cannot be used to excuse the measured near-field
error. This check neither replaces the source by a tabulated production
potential nor silently changes discrepancy weights.

## 7. What the inverse calculation can actually certify

At fixed radial/shape parameters theta the observation map is EXACTLY
coefficient-linear: y=M(theta)c. Let E(theta)c=e contain the declared exact
bulk and optional inherited-pair controls. Choose a particular c_p and a
null-space basis Z with E Z=0. If S is the diagonal discrepancy scale,

    B = S^(-1) M_fit Z,
    r_p = S^(-1) (M_fit c_p-y_fit),
    min_z ||B z+r_p||² = ||(I-B B^+) r_p||².

This last expression is the unconstrained fixed-shape lower bound, NOT a
global bound over theta. The actual coefficient solve also enforces declared
amplitude signs and spectral PSD cuts. Its primal/dual/KKT certificate is
checked in the original full coordinate space after numerical whitening.
An unexplained small whitened residual cannot excuse a negative physical
amplitude. Conversely a fixed-shape tradeoff cannot prove all analytic
environment families insufficient.

For local sensitivity with exact constraints maintained, differentiate

    E dc/dtheta_j = -(dE/dtheta_j)c,
    dy/dtheta_j = (dM/dtheta_j)c+M dc/dtheta_j.

The equality tangent contains 10-7=3 coefficient directions plus4 tied or5
separate-range shape directions. Shape derivatives use log coordinates;
coefficient scales max(1,abs(c)) are numerical normalization, not physics.
The actual singular values and derivative-step changes are saved. Active
inequality/bound cones are NOT confidence intervals. A full local numerical
rank with large residual still denotes a misspecified/inadequately fitted
candidate, not a quantitatively validated material.

## 8. The source total curvature contains cancellation

Direct source-only component evaluation (forces eV/L0, Haa eV/L0²) gives:

| a/h | source pair Haa | source embedding Haa | total Haa |
| --- | ---: | ---: | ---: |
| 1.0 | 17.59821923 | 2.06269364 | 19.66091286 |
| 1.1 | 7.55952389 | 3.31542756 | 10.87495145 |
| 1.42 | -3.59108457 | -1.95695028 | -5.54803485 |
| 1.62 | 1.04462163 | -1.20045022 | -0.15582859 |
| 1.855 | 0.68855519 | -0.35097346 | 0.33758173 |

The source total has multiple traction extrema, not merely one cohesive peak.
The first nested candidate has one resolved maximum over [h,5h], whereas the
source has five Haa=0 crossings there. Both129/257 mixed geometric/uniform
bracket studies resolve this distinction. A small energy-curve error can hide
large derivative errors and a missing shoulder. At1.62h the source's small
total curvature is a genuine cancellation, not floating-point uncertainty.

Pair/embedding decomposition is gauge-dependent: the table does NOT demand
that our LJ pair match the source's tabulated pair separately. It instead
locates the total-response feature that the tested analytic environment must
recover while keeping the LJ base. Do not replace LJ by that table, declare
its nonmonotonicity a cutoff bug, or fit a mobility/yield threshold to mask it.
Any subsequent nonlinear-environment hypothesis needs an explicit analytic
derivation and another same-source force/Hessian validation, not only energy.

## 9. Completed nested comparison and numerical controls

All comparisons here replay the SAME v17 dataset:104 interface fit rows and
48 held-out rows. Do not compare the old596.14 loss directly with a new v17
loss. Its replayed value is1335.56 because the observation set changed.

| calculation | v17 loss | held-out normalized RMS | perfect Haa [eV/L0²] |
| --- | ---: | ---: | ---: |
| old v16 positive-pair candidate, replay only | 1335.5595 | 4.8765 | 27.9164 |
| v17 tied-range,96 shape profiles | 1256.9468 | 4.9991 | 27.7419 |
| v17 nested independent range,72 profiles | 1256.8219 | 5.0089 | 27.6537 |
| v17 independent k1=3/9 starts,192 profiles | 1263.2421 | 4.9570 | 27.8577 |
| same nested shape,remove pair CONTROL,1 profile | 1254.6575 | 5.0045 | 27.5465 |
| external source target | not a fit | not a fit | 19.6609 |

The nested shape is
`[2.3172001610,6.2719556784,12.0,436.6576096,6.2801140650]` in the order
`[k_scalar,k_odd,k_even,alpha_even,k_rank1]`. Its coefficients, in the existing
`[u,v,A,B,C,D3,D1,D2,D_E,K3]` convention, are

    [0.01719039249,0.10069343124,1.42724786649,1.60649560795,
     2.91581618017,37.61096334505,4.76541234049,1.98877158778,
     1.58123232938,32276.17846077].

This is a research vector, not a default Al parameterization. The pair
CONTROL corresponds to epsilon_LJ=0.1474539790eV and
sigma_LJ=0.7448133324L0; these are inherited, not new measurements. The extra
k1 returned close to k_odd. Its training-loss improvement over the tied fit
is0.00994%, with slightly worse held-out error. Removing the pair CONTROL
at that shape improves loss only0.1722%; both pair coefficients remain
positive. The analogous tied-shape free-pair profile gives1255.4174. These
controls do not rule out a different nonlinear joint optimum.

The k_even=15/18 fixed-shape checks give1293.0220/1401.78, worse than the
boundary12 result. This rejects neither all larger decays nor the whole
analytic family. Both main searches ended at their predeclared evaluation
budgets; optimizer convergence/global optimality is NOT claimed.

The exact equality-tangent singular values for tied and nested are:

    tied:   [934.5781,107.0231,28.8494,23.0552,20.2949,11.5054,7.42669]
    nested: [902.8181,181.9299,75.8991,27.0304,22.0660,20.2090,10.7450,7.54423]

Numerical ranks7/8 and condition numbers125.84/119.67 use the stated
coefficient/log-shape scaling. The nested equality-derivative residual is
3.31e-14 and Jacobian step change2.47e-4. Bounds/inequality cones and target
misspecification remain; these are NOT statistical confidence intervals or
proof of unique microscopic identification.

## 10. What is and is not physically improved

For the nested candidate the perfect Haa remains40.65% too large. Its
normal-relaxed fault/saddle energies are0.175144/0.198270J/m² versus source
0.150479/0.172002. The relaxed saddle Haa is24.3787 versus12.7173eV/L0².
Each saddle Hessian is evaluated at its OWN stationary point; the separate
component tables also provide like-for-like fixed-coordinate comparisons.

Held-out normalized RMS, in the order energy/normal-force/Haa/Hxx, is
`[0.8817,4.6549,5.8300,6.6274]`. Only8/12 energy observations are within the
declared scale: RMS<1 does not mean every point passed. Only2/12 force,
5/12 Haa and5/12 Hxx points are within scale. The new fit is NOT accepted
as a quantitatively validated Al active-interface surface.

The nested opening maximum is10.2741GPa at a/h=1.185781; source12.9696GPa
at1.185846. These are ideal rigid coherent-interface traction extrema, NOT
experimental yield, NOT finite-loop nucleation stresses, NOT stochastic
first-passage thresholds. The source's additional traction extrema and
small negative lobe remain visible rather than being clipped away.

The actual static load checker uses eleven declared MPa stress tensors and
nine signed loading/unloading states for each surface. For source plus one
candidate this is198 solved states,99 each, not198 independent experiments.
The nested solutions all passed their root checks. At50MPa pure normal
traction, source/candidate opening increments are
0.0009253659/0.0006592599Angstrom (error-28.76%). At4MPa pure shear they are
0.0003309979/0.0003242181Angstrom (error-2.05%). Maximum static unload return
is3.62e-15L0. This is reversible branch/root bookkeeping, NOT a dynamic
zero-stress hold or evidence of residual plasticity. General tensor traction
is still a research calculation, not a newly connected production tensor PDE.

Existing source MD plane variance comparisons at300K give nested errors
(-12.41%,-22.26%,-14.84%) for normal/slip/transverse. The model uses the
actual fixed-material thermal box stretch; no density gauge is reset, no
atom count is fitted, no mobility is inferred. A0K harmonic approximation
and an atomistic finite-temperature covariance are distinct objects.

## 11. Independent verification and regression

Each completed validation uses243 independent reciprocal q points, direct
bulk radii12/16L0 and continuous local searches; six material dilations each
use254 points at both radii. These sampled modes are positive above their
analytic omitted-neighbor bounds. They are not a whole-Brillouin-zone proof.
At the nested reference radius16 the minimum margin is0.00653958eV/L0².

Nested direct/reciprocal unit-channel energy disagreement at18layers/radius18
is at most3.06e-16. Coordinate-Hessian FD errors at step1e-5L0 are at most
8.63e-8eV/L0² after4e-5/2e-5 refinement. Reciprocal tolerance2e-11->2e-13
changes H by at most4.41e-13. Direct displaced per-site finite-q energies
also agree, with derivative-step/roundoff changes explicitly retained; the
smallest FD amplitude is not automatically the most accurate one.

The final code actually ran:

- targeted5files:29passed,129.63s;
- full `solver_v1`:528passed,1102.12s;
- full `app`:31passed,214.92s,0skipped;
- desktop smoke:exit0,2.01s, then final recheck exit0,0.98s;
- `git diff --check`:clean apart from Git's informational LF/CRLF warnings.

The earlier526-test full run also finished, but the528-test run is the final
regression after adding the last two tests. XML logs are local under
`.cache/normal_response_v17/`; aggregate verified counts belong in the saved
manifest. Production smoke retains a0=0.7713438268704838 and
kappa=86.29296488740997. No app/PDE/historical calibration files changed.

## 12. Reproducibility and retained evidence

Run from the repository root with Python3 and existing dependencies. All
output directories must be fresh; runners refuse to overwrite prior results.
`run_tail_constrained_calibration` actually optimizes; `report_*` only replay
stored candidates. Historical profiles are not newly executed fits.

For example, the independent-range simplex calculation is:

```text
python -m solver_v1.run_tail_constrained_calibration
  --normal-development --quadrupole-saturation --rank-one-range
  --symmetry-channel --quartic-angular --exact-bulk
  --starts results/fcc111_active_interface/normal_response_v17/range_starts.json
  --fixed-pair-from results/fcc111_active_interface/even_environment_v16/saturated
  --additional-wavepoints results/fcc111_active_interface/interface_development_v15/strain_counterexamples.json
  --wavepoint-step 0.1 --stability-stretches .99 1.0037037037037038 1.01
  --optimizer feasible-simplex --max-nfev 96 --out FRESH_RUN_DIRECTORY
```

Join the lines into one shell command or use the shell's continuation syntax.
For the tied control omit `--rank-one-range` and use `tied_starts.json`.
For the nested run use `nested_start.json`, `--optimizer least-squares` and
`--max-nfev 12`. A free-pair coefficient control uses the saved final shape,
omits `--fixed-pair-from`, and uses `--profile-only --max-nfev 1`.
Neither a profile-only successful stop nor a validation replay is a successful
nonlinear shape optimization.

The two independent simplex starts also completed their96 evaluations each.
Their best verified combined profile has shape
`[2.4841497658,6.3099824750,11.9823531066,485.2569646,6.3310019895]` and
loss1263.2421. The last verified trial and the optimizer's last accepted
simplex point can differ at the evaluation limit; both are retained. Neither
start stopped on convergence. Removing the pair CONTROL at this exact shape
gives1258.2939/held-out4.9824, not a resolution of the normal mismatch.
Its equality tangent is rank8,condition122.82; the same practical-uncertainty
warnings apply. Seven run directories contain365 actual verified coefficient
profiles:360 during four nonlinear starts and5 coefficient-only controls.
The eight initial v16 fixed-shape audits are separate, not v17 fit iterations.

Independent validators for each admissible completed run are:

```text
python -m solver_v1.validate_tail_calibration RUN --grid-step .08 --out FRESH_VALIDATION
python -m solver_v1.run_even_environment_validation --models RUN --out FRESH_JET_CHECKS
python -m solver_v1.run_tensor_traction_audit --models RUN --out FRESH_STATIC_LOADS
python -m solver_v1.run_material_plane_covariance MD_SOURCE_DIRECTORY
  --models RUN --additional-stretches .995 1.007 --wavepoint-step .08 --out FRESH_COVARIANCE
```

`MD_SOURCE_DIRECTORY` contains the previously acquired source-bound projection,
metadata and source potential; do not replace it by fabricated coordinates.
The source-only cutoff check is `run_source_normal_audit --out FRESH_SOURCE_CHECK`.
The eight original coefficient diagnoses use `run_normal_response_audit` with
the two named v16 models and a fresh output directory. The versioned result
directory retains start JSON, definitions, all optimizer histories/checkpoints,
independent CSV tables and source/model checksums.

The v17 result folder's local `.gitattributes` disables JSON newline conversion:
these records bind one another by raw-byte SHA256. The actual Windows output
uses CRLF; Git's default text normalization would otherwise change those bytes
in the repository. This scoped setting preserves the verified records on any
checkout without altering historical results or numerical values. Raw/index
blob identity is checked separately from the numerical regression. New runs
may have different byte hashes on another platform and must record their own
bindings; a byte hash is not a platform-independent semantic JSON hash.

The new source-only check does not independently validate DFT/experiment.
Its purpose is narrower: distinguish a source-evaluator error from a material
response mismatch in a like-for-like comparison.

## 13. Decision and next gate

**Do not adopt the independent rank1 range as a validated Al correction.**
The nested mathematical implementation and independent derivative checks
are useful, but the added freedom gives negligible local fit improvement
and does not cure the independently observed normal-force/Hessian errors.
Keep it as an explicitly named research comparison; do not add it to the
production selector or silently change the previous research candidate.
The numerical evidence excludes neither an unexplored shape optimum nor
all possible analytic many-body families. There is no global infeasibility
certificate, and budget stops are disclosed.

The narrowed diagnosis is:

1. Reference derivative/cutoff uncertainty is far too small to explain the
   observed mismatch at1.62h.
2. Finite-q stability constraints are inactive in the inspected coefficient
   optima, so they are not the cause of those particular residuals.
3. A scalar removal/reweighting of the rank1 contribution trades normal
   stiffness against registry response; an independent exponential range
   alone does not resolve that tradeoff in the executed searches.
4. The missing feature concerns the total nonlinear normal environment
   response, including a shoulder/cancellation in source traction. Matching
   energies, bulk moduli or one tangent alone cannot certify that function.

Before another material term is proposed, use the saved complete force/jet
budgets and equality-tangent residuals to specify which extra functional
dependence is required. Preserve the positive LJ base, infinite analytical
sums, per-atom counting, homogeneous constraints and fixed units. If further
same-family continuation is attempted, label its budget and use the v17
observations as already-inspected data, not fresh blind validation. Do not
attach an arbitrary stress-dependent coefficient or empirical yield clamp.

There is also a SEPARATE ideal-versus-actual-strength issue. Even the external
reference's first coherent normal traction is12.97GPa. Reproducing it exactly
would validate an ideal rigid interface curve, not ordinary specimen yield.
Real-strength progress requires the subsequent spatial/core/source problem
already identified in `DISCRETE_FCC_SCREW_DERIVATION.md`: a resolved nonlinear
core/vector/normal-relaxed configuration and a finite-source/loop energy,
with independently calibrated defect state and appropriate kinetics. A
uniform interface-cell barrier cannot become a finite activation barrier by
multiplying it by A_c or an invented line length. Those calculations were
NOT performed by this static v17 audit and are NOT newly validated here.

Production Smoluchowski dynamics, strain/registry/opening bookkeeping,
reduced-model status and specimen probability safeguards are unchanged.
M_a,phys, M_s,phys and t0 remain unavailable; physical seconds/Hz remain
disabled. No dynamic PDE, fatigue lifetime, experimental yield or residual
plasticity claim follows from the static tensor branch returns. The solver
is not complete in the user's material-validated sense; the geometry/mesh
UI redesign remains behind its existing validation/user-confirmation gate.
