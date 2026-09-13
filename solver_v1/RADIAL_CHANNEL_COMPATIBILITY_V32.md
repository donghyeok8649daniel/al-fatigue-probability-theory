# v32 — diagnose the remaining static obstruction before kinetic transfer

Starting branch `probability-pde-solver-v1`, local and freshly fetched remote
`639248b67319a01d9b80009a7a657e00c1fb52e3`, clean worktree. This work does not
change the LJ pair, infinite Poisson/Bessel geometry, production PDE, material
calibration files, mobilities, physical clock, or specimen correlation area.

## Why another calculation is needed

The v31 shape continuation retained seven exact bulk/pristine-interface
anchors but its best maximum normalized discrepancy was 12.50738498, not the
required <=1. Merely reaching an optimizer budget or a search boundary is not
a physical impossibility proof. We test distinct, minimal environmental
extensions **individually**, then a larger dictionary as a capacity diagnostic.
The latter is not a proposed 22-parameter material.

All source targets remain matched 0 K Al99 calculations with the previously
recorded hash and conditions, not experimental aluminum yield measurements.
The source potential is target-only. Previously inspected interface jets are
development data. Four excluded Bloch matrices are not used in selection;
their previously inspected provenance is not represented as blind validation.
Four additional off-grid states are declared in each definition before the
study, but have not yet been used to approve any candidate.

## Fixed energy, units, and added invariants

For each site, use the existing normalized exponential infinite sums

    x_i(k) = sum_j C_k exp(-k r_ij*) ,  x_bulk(k)=1,
    Q_li(k) = sum_j C_k exp(-k r_ij*) STF_l(r_ij*).

The amplitude gauge C_k is fixed at the reference FCC geometry; it is not
reset during deformation. Distances are in L0=4.05/sqrt(2) Angstrom. Each
invariant is dimensionless and its amplitude is in eV per atom. The interface
energy difference sums sites in both half-crystals, with factor two only for
the actual symmetric environment. Energy/gradient/Hessian units are
eV/interface cell, eV/L0, eV/L0^2. Bulk source Hessians in eV/Angstrom^2 are
multiplied by (L0/Angstrom)^2. A source binding inherited from the core loader
also describes row-repeat units; that metadata does not change these explicitly
specified interface/Bloch observable units.

One added coefficient at a time multiplies one of:

    x_i(k_scalar)^z ||Q_1i(k_new)||^2,
    ||Q_3i(k_new)||^2,
    ||Q_2i(k_new)||^2,
    x_i(k_new)-1,
    [x_i(k_new)-1]^2.

The first three and the quadratic density amplitude are nonnegative; the
linear-density amplitude is signed, as the existing B coefficient is. The
chosen inverse ranges 3,5,8,11 are deterministic **research choices**, not
measured aluminum lengths. The extra rank2 term is unsaturated; the original
saturated term remains untouched. Its homogeneous noncubic elastic response
is retained in the bulk anchor matrix, not incorrectly set to zero.

The two scalar densities define a per-site multivariate embedding hypothesis,
not a sum of plane embeddings. A linear density term is gauge-equivalent to
an exponential pair redistribution: B sum_i sum_j w_ij = sum_pairs 2B w_ij.
Consequently its fitted coefficient cannot be called independently identified
many-body physics. The explicitly retained base pair remains LJ. This is a
diagnostic family expansion, not a silent replacement of the production pair.

Every added coefficient is exactly zero in the nested original-model limit.
A repeated same-range column creates a gauge/null direction, not a new
independent measurement. No fitted interpolated force or numerical derivative
is used as the energy generator.

### Coherent vector-channel follow-up

After the separate positive radial invariants fail, test the actual per-site
coherent field, not the sum of its separate energies:

    Lambda x_i^z ||Q1_i(k0)+t Q1_i(k1)||^2, Lambda>=0.

k0 is the old vector range; k1=3 or5 and t=-3,-1,-.3,1 are declared before
evaluation. A signed radial angular weight does not make scalar electron
density negative; x is unchanged. The cross term is exactly
2t Q1(k0):Q1(k1). At t=0 this duplicates the old rank1 column; at k0=k1,
t=-1 it vanishes. These are tested gauge limits, not fitted exceptions.

The actual site jets are padded in depth, combined linearly, and then passed
through the same analytic product/norm chain rule. In pristine cubic bulk:

    L_mix(q)=L_0(q)+t L_1(q), H_mix(q)=2 L_mix L_mix^T,
    dL_mix<=dL_0+|t|dL_1,
    dH_mix<=2[2||L_mix|| dL_mix+(dL_mix)^2].

The bulk operator uses the existing analytic FCC omitted-neighbor envelope.
An independent six-phase direct sinusoidal per-atom energy, with amplitude
refinement, verifies its Hessian. The infinite reciprocal site jet remains
canonical; the finite neighbor sum is only this independent test.

All10 profiles completed in49.66s. Best single eta12.439883; the full
eight-column capacity diagnostic reaches12.071181, with excluded-q error
239.4069%. Its conditional rank is7, not11: coherent weights contain
redundant radial quadratic directions. This follow-up also fails acceptance.

## Analytic chain rule and positive-power test

For I=Q:Q and coordinate indices p,q,

    I_p=2 Q:Q_p,
    I_pq=2 Q_p:Q_q+2 Q:Q_pq,
    (g I)_p=g I_p+g_x x_p I,
    (g I)_pq=g I_pq+g_x(x_p I_q+x_q I_p+x_pq I)+g_xx x_p x_q I.

Each per-site expression is summed only after its complete environmental
moments have been formed. This is implemented independently in
`coordination_power_research.py` for g=x^p, p>=0:

    g_x=p x^(p-1),  g_xx=p(p-1)x^(p-2).

The old p<=1 limit was a search convention, not a universal physical bound.
We explicitly test p=0,1,1.5,2,3,4,6,8 as **alternative** single-exponent
models at fixed other shapes, not their sum. Positive exponents preserve the
isolated-neighbor decay: x^p I1 decays like a polynomial times
exp[-(p k_scalar+2 k1)r]. At a pristine centrosymmetric site x=1,Q1=0,
all these choices have the same harmonic column. They change finite-state
response, not the fixed-coefficient pristine finite-q response.

## Positive nonlinear cross-sector control

The existing verified mixed invariant y I3, y=x-1, is also tested without an
unchecked indefinite polynomial sector. Add one nonnegative amplitude Lambda:

    Delta F_i = Lambda [y_i+t I3_i]^2,
    Delta column = column_C+2t column_cross+t^2 column_K3.

Its sector matrix is Lambda (1,t)^T(1,t), hence PSD. The original C,K3
nonnegative terms remain. The deterministic t choices are -100,-10,-1,-.1,
.1,1,10,100. Their full cone is a diagnostic only: the eight columns add
**one** independent linear direction beyond C and K3, not eight physical
parameters. The SVD detects the redundant directions. Separated-atom reference
energy is included through column_C; it must not be dropped from cohesion.

## Compatibility and identifiability

At fixed shape, all observable jets are coefficient-linear: y=M c. Solve

    min eta,  E c=e,
    |M_i c-y_i| <= eta scale_i,
    c_j>=0 for explicitly declared indices.

Keep the same seven exact anchors and 115 inspected interface curvatures,
plus four prescribed fit-q symmetric matrices. Off-diagonal components use
sqrt(2) weights, so Euclidean feature norm equals Frobenius matrix norm.
The q scale .05 ||H_source||_F is a discrepancy convention, not experimental
uncertainty. No fatigue, yield, phase, or crack probability enters the loss.

Independent primal/dual/KKT checks accompany every LP. With column scaling D
and N spanning null(E D), singular values of diag(1/scale) M D N report
conditional identifiability. They do not provide confidence intervals or
full shape-parameter identifiability. No spectral constraint was included in
these necessary tests; a failed necessary target box is sufficient to reject
the trial, while a passed box would still require spectral/fresh validation.

## Actual fixed-shape results

| Study | Profiles including baseline | Best eta | Excluded-q maximum relative error at reported best |
|---|---:|---:|---:|
| Baseline | — | 12.507385 | 210.0802% |
| Best single angular radial channel, rank3 k=5 | 12 single trials | 12.271819 | 189.7212% |
| All 12 angular channels, capacity only | 14 total | 11.248055 | 200.6077% |
| Alternative positive powers | 9 | 12.507385 | 210.0802% |
| Best single scalar-density channel, quadratic k=3 | 8 single trials | 12.475869 | 210.8576% |
| All 8 scalar-density channels, capacity only | 10 total | 11.255225 | 218.1568% |
| Positive cross-square cone, capacity only | 10 total | 12.502764 | 209.9249% |

These four studies actually completed in 226.19,36.76,78.34,25.22 seconds
respectively (wall times depend on concurrent work). They contain 43 LP
profiles including four baseline replays. No trial passes eta<=1. The scalar
dictionary also drives a base LJ coefficient to zero and is ineligible even
apart from its large residual. Several p>=2 alternatives have the same issue.

Including the separately declared coherent-vector and vector-density-shape
follow-ups, there are 61 actual fixed-shape profiles including six baseline
replays. None passes. The latter tests the positive site expression
Lambda [1+t(x-1)]^2 ||Q1||^2 for six declared t values. Its eight-profile
study reaches eta=12.462730, excluded-q error=209.0168%, conditional rank5.

The angular dictionary has conditional rank15, singular values from1.79350
to.00206962 (condition about867). The scalar dictionary rank11 spans1.02248
to2.33239e-5 (about43838). The cross dictionary rank4 has seven additional
singular values around1e-16: genuine algebraic redundancy, not new physics.

For the best angular dictionary the opening curvature at the inspected
v17_opening2.455 state is -0.189856 versus source +1.521219 eV/L0^2.
The v20_new_11 opening curvature is -0.357621 versus +2.865426. The
fit-q K(.375,.375,0) H22 is66.0350 versus37.1732. These simultaneously active
errors explain why adding more simple radial channels did not resolve the
static obstruction. A curvature component alone is not a full Morse-index
or global barrier-topology statement.

An independent deterministic same-family shape start (2.4,6.1,8,30,6.3,0)
is separately declared to test the previous near-saturation-bound basin.
It completed 180 profiles in 1526.05 seconds, ending at the evaluation budget
(optimizer_success=false). Best eligible eta=12.899739, worse than the v31
baseline. Independent full replay agrees; excluded-q error=217.2508%, minimum
tested eigenvalue=1.448810. This is not optimizer convergence or a no-go proof.

### Radial primitive follow-up

`QUADRATIC_EXPONENTIAL_ENVIRONMENT_V32.md` derives the precise obstruction of
positive exponential mixtures: their logarithmic radial curvature cannot be
negative. The actual target auxiliary density has negative curvature at the
nearest-neighbor distance. A positive quadratic-envelope exponential retains
the exact Poisson structure through derivatives of the exponential transform.
Its source-density-only fit improves, but the fixed-shape material test gives
eta=19.349784 and excluded-q error=230.8101%, hence fails. This result is saved
separately from the 61 channel profiles. A subsequent direct material-shape
search and fresh validation have their own status/summary files. Neither
auxiliary density matching nor a passing derivative test is Al calibration.

### Is imposed convexity the actual obstruction?

The source auxiliary embedding has F_xx(1.2)=-3.182568 eV in its normalized
density gauge, whereas -A sqrt(x)+B(x-1)+C x-shift-squared has
F_xx=A/(4 x^(3/2))+2C>=0 for A,C>=0. Convexity of one embedding contribution
is not a universal requirement for stability of the **total** EAM energy.
Nevertheless a sign constraint cannot be removed and the result assumed stable.

`embedding_convexity_audit` actually solves four fixed-v31-shape diagnostic LPs:
inherited signs eta12.507385; signed C eta12.507385; signed rank3 quartic K3
eta12.404692; both signed eta12.404692. All retain positive LJ in this audit,
but none passes. Thus inherited scalar convexity alone is not responsible for
the mismatch at this shape. Negative quartic coefficients may introduce
large-amplitude instabilities; these ablations are not proposed material laws.
They are four additional profiles, not part of the 61 channel profiles.

## Reproduced zero-amplitude implementation bug

Constructing the screened research model with D1=0 previously raised
StopIteration because the base energy intentionally omitted that zero term.
The constructor now builds a unit geometric moment for optional diagnostics
without adding it to the energy list. Tests compare energy, gradient and
Hessian exactly against the unscreened zero-D1 limit, for tied/separate ranges.
This is a genuine boundary-constructor repair, not a mechanism for changing
the nonzero-D1 v31 material results or creating plasticity.

## What these calculations do NOT finish

1. A usable Al interface calibration still requires a joint material pass;
   fixed-shape failures do not prove the entire LJ/Bessel analytic family
   impossible. These extensions are not adopted.
2. No new kinetic trajectory is claimed here. v31's actual20ns data retain
   their unresolved low-frequency dissipation and local-PMF mapping limits.
   Neither static fit amplitudes nor finite-frequency point mobilities can
   supply a validated local a/s production clock.
3. No newly failed material is passed into a dislocation/source/yield or
   fatigue simulation to manufacture a result. Core, finite-source geometry,
   specimen slip budget and kinetic validation remain distinct requirements.

Production physical seconds/Hz remain disabled. Existing theory and historical
results are preserved. A negative compatibility result is a completed
diagnostic, **not** completion of the user's three physical validation goals.
