# v15: full-environment calibration development, not production promotion

This study follows actual v14 interface-shape failures. It preserves the LJ
pair, all infinite Poisson/Bessel kernels, physical atomistic cell units and
per-atom environment counting. No yield stress, fatigue lifetime, mobility,
desired slip/opening order or specimen correlation area appears in the loss.
The matched source is 0 K Mishin Al99, SHA256
`60c8a085be79d273324ab421f5b1447578fef55c1acfc6492c0999f15ee8a284`.
It is source/validation data, not a replacement production pair potential.

## 1. Honest development / validation split

The nine previously inspected v14 held-out failures become development data.
They cannot remain called blind validation after influencing this task.
`interface_development_targets.py` adds the measured saddle's a/slip coupling,
normal opening forces at 1.1/1.3/1.5/2/2.5/3h and three Shockley path energies.
There are five exact staged bulk constraints, thirty interface fit rows and
thirty NEW off-grid/off-path held-out rows. The holdout locations are fixed
before fitting. They never enter the loss or selection of deterministic starts.
These are tests against the same source surface, not independent experiments.

Force scale is 250 MPa converted by the microscopic cell's conjugate work.
New energy scales are max(10% of source magnitude,0.01 J/m²), Hessian scales
max(10%,0.1 eV/L0²). Inherited scales are retained. These explicitly specified
model-discrepancy tolerances are not reported statistical uncertainties.

## 2. Executed baseline diagnosis

The unchanged quartic angular family ran 71 actual profiles/246.71 s and
reached training squared loss1251.068874, held-out normalized RMS78.29824.
The three source cubic constants still equal114/62/32 GPa, cohesion3.36 eV
and the nominal FCC force is zero. Interface validation fails: its relaxed
fault/saddle energies are0.218260/0.218464 J/m², versus0.150479/0.172002.
The tiny reverse barrier is not a successful Al result.

Four deterministic continuation starts ran154 profiles/479.72 s, reaching
the same loss1251.068873. Three starts stopped by ftol and one was numerically
rejected. Neither ftol nor a finite budget proves a global minimum. At the
best radial shape, unconstrained coefficient least squares and sign-constrained
least squares both give1251.068874: active sign constraints do not explain the
remaining residual at that shape. This is not a lower bound over all shapes.

A predeclared100-shape grid, including long/short radial ranges, recorded85
verified sign-constrained profiles and15 numerical failures. Two selected
low-loss *training-only* seeds were retried with spectral constraints; both
stopped on numerical certificate failures, not accepted stationary optima.
These failures are preserved. They do not prove that those basins are bad
physics, or that the whole analytic family is impossible.

## 3. Positive exponential-density mixture: tested, NOT adopted

Let the fixed reference-normalized full environments be x_0 and x_1, with
the already-defined scalar/odd exponential decay scales. Test exactly one
extra shape parameter w in [0,1]:

    x_i=(1-w) x_0,i+w x_1,i,
    F_i=-A sqrt(x_i)+B(x_i-1)+C(x_i-1)^2.

F acts AFTER the mixed site density. It is not a mixture of embedding energies.
All first/second x derivatives combine linearly before applying F',F''.
The two normalizations remain fixed during strain. w=0 recovers the existing
family; coincident decays make w unidentifiable. No new length or area is fitted.
The coefficient operator, finite-q density-gradient square and tail envelopes
are changed consistently; endpoint, derivative, periodicity and tail tests pass.

Two actual starts w=.2/.7 ran274 profiles/761.17 s. The best weight tends to
the declared logit lower bound w=.0001234129. Loss1251.096986 and heldout78.30371
do not improve the unchanged model. This extra parameter is **not adopted**.
Independent stationary/finite-q/refinement validation was executed and is saved
under `mixed_density/validation/`. Local numerical rank is not material confidence.

## 4. Spatial stiffness information, not phonon Hz

Finite-q positivity alone cannot calibrate the magnitude of the atomic
restoring force. The independently differentiated SOURCE EAM Bloch matrix is

    H(q)=sum_R [1-cos(q.R)] Hess(phi+2F'_bulk rho)(R)
         +F''_bulk g(q) g(q)^T,
    g(q)=sum_R sin(q.R) grad rho(R).

This is eV/Angstrom² (then explicitly eV/L0²), not a frequency. No mass or
mobility is introduced. It was independently checked by evaluating actual
per-atom EAM energies of sinusoidal plane displacements, then refining their
amplitude. Nonlinear F is still per atom after summing all its neighbors.

Seven predeclared longitudinal/transverse projections at q=(.25,.25,.25),
(.5,0,0),(.4,.4,0), in cubic2pi/a_lat units, are development targets. Seven
projections at (.5,.5,.5),(.9,0,0),(.6,.6,0) are excluded from fitting.
Symmetry-equivalent transverse modes are not counted twice. A10% stiffness
discrepancy scale is declared; target data are static0K at4.05Angstrom,
NOT the300K MD box at4.065Angstrom. Actual fits and independent radius12/16
tests are stored separately in `spatial_stiffness_fit/`.

## 5. Next minimal mixed-invariant hypothesis

After the preceding failures, consider a differentiable per-site environment
energy F(x,I), with I=||Q3||² the existing full STF rank3 invariant. The first
mixed Taylor coefficient absent from the current separable family is

    Delta F_i=E_xI (x_i-1) I_i.

This adds ONE energy coefficient, no new radial range, orientation, empirical
yield threshold or pair function. It is an analytic many-body research
extension, not a claim of standard EAM/MEAM or new Al calibration.
The density and tensor are summed over infinite neighbors before forming
this site product, and the two half-crystal site energies are then summed.

Put y=x-1. Its analytic derivatives are

    (y I)_j=y_j I+y I_j,
    (y I)_jk=y_jk I+y_j I_k+y_k I_j+y I_jk,
    I_j=2 Q:Q_j,
    I_jk=2 Q_j:Q_k+2 Q:Q_jk.

The per-site baseline x=1,Q3=0 gives no energy, first or second derivative.
For centrosymmetric affine bulk Q3=0 exactly. A general small displacement
has y=O(u),I=O(u²), so the new term also has zero harmonic finite-q column.
At separated atoms I=0, so it does not change the isolated-atom reference.

Do not let a signed mixed term introduce unbounded negative nonlinear energy:

    C y²+E_xI y I+K3 I² = [y,I] M [y,I]^T,
    M=[[C,E_xI/2],[E_xI/2,K3]] >=0,
    C>=0, K3>=0, E_xI²<=4 C K3.

This convex positive-semidefinite coefficient constraint is imposed during
profiling with separating polarizations, independently of the finite-q
stability cuts. No force or coefficient is clipped after the calculation.
The guarantee applies only to this nonlinear sector: it is NOT proof of
global crystal stability or validity of the other signed angular terms.
Zero E_xI recovers the verified old implementation. The source derivative,
per-atom product, Hessian, periodicity, known convex-projection solution and
actual observation-matrix replay must pass before any fit is launched.
The implementation remains research-only regardless of fit improvement.

## 6. Executed tradeoffs and family ablations

The first completed runs, followed by independent root/curve and155-q
radius12/16 verification, give the following. Larger loss definitions are
NOT compared as though their targets were identical.

| Run | Training squared loss | Interface heldout normalized RMS | Other result |
|---|---:|---:|---|
| unchanged exact-bulk quartic |1251.068874|78.29824|C=114/62/32GPa|
| positive density mixture |1251.096986|78.30371|weight approaches lower bound|
| extra static Bloch targets |2300.289575|79.10517|seven additional Bloch heldouts RMS9.19104|
| one mixed density/angular invariant |1128.234326|71.20138|PSD sector near boundary|
| rational angular ablation |1146.171120|79.79577|one stopped start, two numerical failures|
| elasticity included with5% scales rather than exact |311.221239|14.82148|C=83.42508/70.94344/34.24797GPa|

The last experiment preserves exact force/cohesion but makes the three
elastic constants ordinary5%-scaled fit observations. No target/force is
clipped. C'=(C11-C12)/2 falls to6.24082GPa versus26GPa: the lower interface
loss is bought at the expense of physically independent bulk elasticity.
It is rejected, not a replacement calibrated material.

The rational form is the already derived v14 alternative

    F_angular(I)=D3 I + K3 I²/(1+alpha I), alpha>=0.

It preserves analytic derivatives and the I->0 harmonic limit. Three fixed
starts alpha=1e3/1e5/1e7 ran96 profiles; the best alpha~934.875. The relaxed
fault/saddle energies0.184491/0.184578J/m² still have an almost vanished
reverse barrier. This is not an improvement in the required full topology.
Numerical failures are retried separately after the precision repair below;
old run directories are historical records, not overwritten optima.

All parameter vectors, target-by-target physical-unit residuals, optimizer
terminations and heldout counts are in `completed_report/`. The script is a
report/replay, not an optimization disguised as one. Future report directories
include any subsequently completed retrials without erasing these results.

## 7. Exact-stage identifiability, not the wrong larger tangent space

With observations y=M(k)c and exact calibration rows M_E(k)c=y_E,
use coefficient normalization D=diag(max(1,abs(c))) and

    A=M_E D / scale_E, Z=null(A),
    dc_coefficient=D Z dz,
    A z_j=-(M_E,j c)/scale_E,
    dy_j=M_j c+M D z_j.

The least-norm particular correction z_j keeps ALL exact targets fixed.
Held-out observations are absent from this Jacobian. For exact-bulk stages,
force, cohesion and all three cubic constants must be eliminated, not just
the two force/cohesion rows used in the older, larger-space sensitivity.
Both diagnostics are clearly distinguished in `validation_exact_tangent/`.
Finite-difference shape steps are independently halved.

| Stage | Tangent rank | Normalized condition number |
|---|---:|---:|
| quartic |8|311.58|
| mixture |9|6.06e6|
| static-Bloch constrained |8|280.36|
| density/angular cross |9|921.72|
| rational |9|3006.47|
| non-exact-elastic joint fit |11|3497.86|

Exact-constraint derivative residuals are below1.1e-12. These numbers depend
on the declared normalized coordinates. Active inequality **cones** are not
replaced by equality constraints; this is neither statistical confidence nor
a global uniqueness proof. The tiny density-mixture weight is practically
unidentifiable despite numerical full rank.

## 8. Numerical defects repaired without relaxing physical constraints

Two reproducible numerical mechanisms were separated from material failure.

1. A zero-dual, merely near-active inequality was incorrectly promoted to an
   equality in physical-coordinate face polishing. This changed already
   verified KKT~1e-14 solutions into KKT.003/33.9. Only positive-dual face
   constraints are now imposed; columns are numerically equilibrated and the
   exact factored inverse is used. A polish is accepted only after the SAME
   primal/KKT certificate. Independent scale/zero-dual regressions cover this.
2. Projecting almost parallel spectral cuts through A A^T squared their
   condition number. The fixed probe
   k=(3.08175732286109,3.73956565371123,6.84896287148422)
   failed32 cuts at margin -5.93447e-9. A direct affine SVD/null-space
   projection of the SAME face instead needs one cut: margin -1.59333e-12,
   KKT1.86810e-12, primal1.42464e-12. The negative margin is within the
   unchanged separately computed arithmetic cancellation envelope; no
   material stiffness or spectral-tail criterion was lowered. A mathematical
   near-parallel regression independently exposes the Gram failure.

The corrected projection is x_particular+Z Z^T(x0-x_particular), A Z=0.
It still must pass primal/dual/KKT checks. It is not an unconstrained fit,
force clipping, or tolerance inflation. Actual before/after probes are saved
in `precision_probe_before/` and `precision_probe_direct_svd/`; additional
deterministic optimizer retrials have separate directories.

## 9. Interpretation gates

A stationary bulk and smaller training residual do not validate local GSF,
cleavage, atomic cores, experimental yield, collective-coordinate activation
normalization or kinetics. Dynamic PDE/UI promotion is not part of this study.
The exact bulk matching is a staged inverse constraint, not zero uncertainty
in real Al. Numerical profile failures, small singular values, parameter bounds
and independent heldout errors must accompany any reported vector.

Primary reference: [Mishin et al., Phys. Rev. B59 (1999)3393](https://doi.org/10.1103/PhysRevB.59.3393).
Actual public kinetic/strength evidence is documented separately, with0K/293K/
300K and site/plane/line/specimen normalizations kept distinct.

## 10. Newly accessible radial basin and independent stability counterexample

After repairing the numerical face projection, two previously blocked starts
were actually retried. `grid_svd_certified` reached loss1084.893445 and
heldoutRMS63.43026 in77 profiles. Independent central shape derivatives
(3-point, relative step5e-5) give `central_shape_refinement`:loss1084.852151,
35profiles,ftol termination,first-order optimality.26176. This is not a global
minimum or proof of satisfactory stationarity. Failed profile coordinates and
their numerical error messages are now retained explicitly.

The lower loss was NOT accepted. The independent.1-step grid and continuous
search found a large negative finite-q curvature near

    q_cubic=(.144635212,.144635212,.144635212),
    min eig H=-5.0821493 eV/L0², tail bound=.00205967 eV/L0².

Matching all three elastic constants and a sparse q grid missed a real
finite-wavelength instability, not a sign-roundoff effect. Reproduction files
remain in the rejected run directories. A new recorded counterexample is
added as a stability constraint, not as an Al target or an energy modification.

`spectral_refined_fit` uses the same family/loss with a.1 constraint grid plus
that q. Two starts ran38profiles/129.31s:loss1117.294014,RMS65.36096.
Independent.08 and.04 grids (243 and1641 points at each radius12/16), plus
continuous local searches, are positive above tails. The continuous small-q
minimum is~.00822224 with tail~.00205054. This is sampled/local-search evidence,
not a whole-zone proof. The new mixed-invariant continuation
`spectral_cross_fit` givesloss1077.905674,RMS64.83695 and also passes its243q
and continuous-search tests. All retainexact114/62/32GPa,3.36eV/atom.

Both still fail the actual interface data. For example the refined quartic
predicts heldout direct-pathHaa312.75 vs36.3168eV/L0², and normal force
-12.52845 vs-1.96851eV/L0. Stability alone does not repair an inaccurate
force/curvature surface. An84-case actual static MPa normal/shear sweep,
including this refined candidate, is saved separately. No rate is inferred.

## 11. A single stable reference geometry was also insufficient

The same fixed material, with unchanged L0/rho_ref/coefficients, was tested
at the actual MD dilation4.065/4.05=1.0037037. Both new spectral candidates
have unstable periodic plane modes there. `spectral_material_modes/` refuses
harmonic covariance rather than inverting negative Hessians. Reference-point
elasticity does not establish stability in a neighborhood or at the300K box.

An explicit subsequent constrained study is therefore run at stretches
.99,1,1.0037037037,1.01 with the already verified fixed-material harmonic
formulas. The calibration observations remain0K and unchanged; positive
nearby stiffness is an admissibility check, not a fit of finite-temperature
free energy. The cross invariant contributes (x-1)H_D3 away from x=1 and
must not be dropped. Results are separate `strain_stable_*` directories;
their actual completion and independent checks, not this proposal alone,
determine their status. Neither PDE nor UI is automatically changed.

### Completed neighborhood-stability recalibrations

All three subsequent runs finished and were independently replayed. The
same exact bulk targets, coefficient gauges and numerical certificates were
kept. Static Bloch magnitudes add seven training observations only in the
third row, so its loss is not directly comparable with the other two.

| Run | Profiles / wall seconds | Training loss | Interface heldout RMS | Bloch heldout RMS | Exact-tangent condition |
|---|---:|---:|---:|---:|---:|
| strain_stable_quartic |62 /396.01|1251.068874|78.29824|not fitted|311.58|
| strain_stable_cross |112 /720.95|1128.234352|71.20608|not fitted|920.60|
| strain_cross_bloch |95 /503.31|1957.053393|78.90034|8.70063|207.26|

The first two return essentially to the earlier, less accurate but more
stable basins. This is not evidence that a different global optimum cannot
exist. Independent tests use six stretches .99/.995/1/4.065*4.05^-1/1.005/
1.01, 254 sampled wavepoints and radii12/16. All three pass those sampled
tail-bounded tests, plus reference-geometry continuous local searches.
This is still not a full-zone or finite-temperature stability proof.

The extra Bloch targets improve that dataset's conditioning but do not solve
the independent interface discrepancy. At the actual MD box its harmonic
normal variance is30.94% above the measured estimate. The cross-only fit
gives normal/slip/transverse variance errors +7.94%/-22.47%/-15.07%; the
quartic gives -8.89%/-11.81%/-3.40%. Sampling, thermostat and anharmonic
differences are separate from these fixed-material predictions. None of the
MD variances was fitted, and no mobility follows from matching variance.

`final_report/` supersedes earlier report snapshots for this task: it contains
all20 actually executed optimization runs, exact parameter vectors, every
optimizer termination, independently replayed target residuals, separated
interface/Bloch heldouts and exact-stage singular values. Numerical failures
remain visible. All `full_material_accepted` flags are false. The conclusion
is **current tested low-complexity candidates are inadequate**, not a theorem
that every possible analytic LJ/environment family is impossible. A lower
loss purchased with unstable phonon curvature or incorrect cubic elasticity
is not a successful Al calibration.

## 12. The large heldout curvature is not a unit-conversion error

A direct analytic observation-matrix replay separates the terms at the
previously excluded direct-path state (a=h,s1=.39L0,s2=0). All models and the
source use the SAME physical normalization. The normal Hessian budget is:

| Candidate | LJ pair | Scalar embedding | Quadratic angular | Quartic angular | Cross invariant | Total Haa |
|---|---:|---:|---:|---:|---:|---:|
| strain_stable_quartic |91.689|68.539|230.133|13.215|0|403.575|
| strain_stable_cross |97.426|75.620|182.976|10.841|2.338|369.201|
| strain_cross_bloch |78.542|80.521|245.705|20.205|-.559|424.413|

Units are eV/L0² throughout; source Haa=36.316785 in those units. Exact
values and reconstruction residuals are in
`final_report/heldout_curvature_term_budget.csv`. This is an evaluation of
already fitted parameters, NOT a fit to those excluded points. The analytic
rows are produced by `development_observations` and the model's observation
cache; each coefficient multiplies its unweighted jet column.

The discrepancy is distributed across the paired scalar cancellation and
angular curvatures. It cannot honestly be repaired by changing MPa/eV units,
normalizing the strain plot, or simply deleting K3. The density/angular
extension improves some opening forces but does not suppress this large
off-calibration registry curvature while preserving the bulk constraints.
The task therefore leaves the explicit interface-adoption gate closed.

## 13. Frozen-code reoptimization and completed verification

The final same two-start cross optimization was independently executed,
not reconstructed from its saved best vector:112 profiles,532.31s. Shape,
coefficients, predictions, normalized residuals, cubic constants and loss
1128.2343520341506 reproduce exactly in this environment. Timing fields are
retained and naturally differ. The independent243-q/radius12/16 plus
stationary/curve/identifiability checks were rerun as well.

`release_report/` now contains21 runs /1770 coefficient profiles, including
that reproducibility run. The earlier20-run `final_report/` remains a
historical snapshot and contains the additional curvature term budget.
The frozen code passed73 targeted tests,504 full solver tests,31 app tests
with no skips, and the desktop smoke. A software-test pass or deterministic
fit reproduction does not change the failed material-adoption decision.
