# v20: initial interface tangents versus finite-deformation calibration

This is a continuation of the completed v19 negative material assessment.
No new energy term is introduced at this stage. Neither the production PDE
nor its LJ/Bessel calibration, mobility, time units or area semantics changes.

## Question and exact algebra, before running fits

v19 leaves Haa about41% too high even though the five bulk quantities are
enforced exactly. For an intact stable interface and small applied traction:

    H_int delta_q = A_atomic_cell L0 T / E0,
    delta_q_phys = L0 delta_q.

The tensor projection, actual atomic area, E0=1eV and L0=4.05/sqrt(2)Å
are unchanged. Thus a wrong H_int gives a wrong static compliance even with
correct stress units. It must not be repaired by tuning mobility or kappa.

At fixed microscopic shape k, the existing analytic energy/force/Hessian
observations are exactly linear in the ten coefficients:

    y = M(k)c, c=(u,v,A,B,C,D3,D1,D2,D_E,K3),
    E_bulk(k)c = e_bulk,
    h_a(k)c = Haa_target, h_x(k)c = Hxx_target.

All infinite scalar/vector/tensor environment sums precede the per-atom
functions. The five exact bulk rows constrain equilibrium, cohesive energy,
C11,C12,C44. The initial normal and registry Hessians are distinct LOCAL
interface observables, not extra cubic elastic constants. We now test whether
the same family can satisfy them simultaneously, and at what cost elsewhere.

Solve deterministic constrained least squares with the existing spectral/tail
inequalities and coefficient signs. At each supplied wavevector/dilation:

    H_R(q;c) - sum_j |c_j| tail_j(q) I >= 0.

This condition is only sampled stability; it is not an all-zone theorem.
Use the existing polarization-exchange QP and its primal/dual/KKT checks.
A QP numerical failure is not a proof of physical infeasibility. A separate
scaled equality/sign linear-program check distinguishes that simpler question.
Resolved positive u,v are required for an eligible LJ material candidate.

For two exact tangent rows the remaining fit has102 observations, the SAME
other observations/scales as v19. Also report the original104-row source loss
for every control so changing which rows are exact never hides a tradeoff.
Source discrepancy scales are not experimental standard errors. The 0K Al99
potential supplies targets only; it does not replace the analytic LJ pair.

If t is a controlled normal tangent, the fixed-shape value function is:

    V(t) = min_c sum_other102 [(M_i c-y_i)/scale_i]^2
    subject to bulk constraints, h_a c=t, h_x c=Hxx_source,
               original coefficient/spectral inequalities.

This is convex in t on its feasible domain. It quantifies the compromise
between local compliance and the other energy/force/curvature targets. The
interpolated t values are DIAGNOSTIC CONTROLS, not newly measured material
data. All reports retain both imposed and original source targets. Endpoint
or finite-range failure does not prove global impossibility of the family.

To see convexity, interpolate two feasible pairs(c1,t1),(c2,t2). Exact
constraints interpolate linearly, sign constraints stay feasible, and
sum_j|c_j|tail_j is convex. Consequently the tail-subtracted operator at
the interpolated c dominates the interpolation of the two feasible operators.
The quadratic observation loss is convex too, giving

    V((1-z)t1+z t2) <= (1-z)V(t1)+z V(t2),  0<=z<=1.

This statement is for FIXED microscopic shapes, not the nonconvex outer
radial optimization. Seven independent equality rows leave three coefficient
directions from ten. Adding five (or six) shape derivatives requires solving
the differentiated seven-row constraint, not simply differentiating M c at
fixed coefficients and calling that constrained material identifiability.

## Declared staged experiment

1. Use both saved v19 joint shapes: old family p=0 and power family p=-1.
   Replay their original coefficient profiles first.
2. Keep the five bulk equalities. Fit with no extra tangent constraint,
   normal only, registry only, and both source tangents.
3. With registry tangent exact, impose normal tangent at five fractions
   0,.25,.5,.75,1 between the unconstrained prediction and source. Record
   the common102-observation loss, original104-source loss, KKT and stability.
4. Only after these controls are evaluated, consider actual radial/shape
   refinement with BOTH tangents exact; do not widen the old shape windows
   merely to lower loss. Record starts, evaluations, bound hits and stops.
5. Evaluate any eligible candidate independently under signed normal,
   two-shear and mixed MPa conditions; compare energy/force/Hessian and
   finite-q/tail refinement. Static unload is not a dynamic zero-stress hold.

Previously inspected v17/v19 excluded data stay excluded but are retrospective.
Declare new off-grid observations before calculation, exclude them from the
objective and selection, and check their state identities do not overlap.
No fitted tangent may subsequently be presented as held-out validation.

## Scientific limits

Even a successful initial-tangent fit would not establish the finite-opening
or GSF shape, defect-source strength, actual yield, fatigue or dynamics.
The rigid-interface ideal traction is not measured yield stress. Production
Ma_phys, Ms_phys and t0 remain unavailable; physical seconds/Hz remain disabled.
Keep original and new candidates separate. Do not add an empirical plasticity
law, force clipping, A_c in atomic energy, or a hidden time conversion.

Source: Mishin, Farkas, Mehl and Papaconstantopoulos (1999), PRB59,3393,
DOI10.1103/PhysRevB.59.3393; matched0K source potential from
[NIST Al99](https://www.ctcms.nist.gov/potentials/entry/1999--Mishin-Y-Farkas-D-Mehl-M-J-Papaconstantopoulos-D-A--Al/).
Its existing SHA-256 and the parent fit hashes must be saved with this study.

## Completed fixed-shape control: what failed, exactly

Sixteen profiles completed, all with verified equality, sign and sampled
spectral/KKT checks. The original target tangents are19.6609128633 and
4.39359303512eV/L0². The two initial five-equality profiles exactly reproduce
the saved v19 values/losses, before adding either tangent constraint.

| Fixed shape | Original104 loss | Both-exact102 loss | Pair after both-exact |
|---|---:|---:|---|
| v19 old p=0 |1238.611009|2831.110021|u>0, v=0|
| v19 power |1227.634888|2859.423588|u>0, v=0|

For the code's pair coefficient convention,

    phi(r*) = u r*^-12 - v r*^-6,
    sigma_LJ/L0 = (u/v)^(1/6), epsilon_LJ/E0 = v²/(4u).

Thus v=0 with u>0 is NOT an ordinary finite positive-epsilon/sigma LJ model.
It is the closure of the coefficient optimization domain, not a candidate to
instantiate or feed into the PDE. No arbitrary v floor was inserted. Strictly
positive feasible profiles may exist at greater loss; a boundary optimum of
one fixed shape is not proof that the whole analytic family is impossible.

The registry-exact-only controls remain close to the baseline losses. The
normal-exact-only controls instead reduce registry tangent to2.394/2.419,
well below4.394. With both exact, major conflicting observations include the
intermediate opening forces (a/h≈1.1–1.3) and the registry saddle curvature.
The full rows, units, source/imposed targets and residuals are exported.
Fixed-shape V(t) has positive sampled second differences for both families.

## Shape search and the distinction between its endpoint and a candidate

Two deterministic starts use the corresponding completed v19 shapes, not
random sampling or a reset to invented material parameters. The inherited
numerical windows are k_scalar∈[.7,8], k_odd,k_even,k_rank1∈[2,12],
alpha_even∈[1e-5,1e6], and p∈[-1,1] (p=0 fixed for the old family).
They are numerical search windows, not independently measured Al ranges.
Both tangent equalities and the five bulk equalities remain active at EVERY
profile. Every trial/difference evaluation and numerical stop is retained.

The closed-domain optimizer may end at a zero-LJ profile while an earlier
positive-LJ trial has higher loss. Report both, and never describe the latter
as a converged material-shape optimizer endpoint. `shape_selection_status`
checks this explicitly. An eligible positive-LJ profile still needs independent
material validation; coefficient feasibility/KKT alone does not adopt it.

## Numerical helper repair

The v20 toy equality tests reproduced a boundary-case error in the existing
QP helper: when all inequalities are constant and satisfied on the equality
manifold, SciPy received an empty LinearConstraint and raised IndexError.
The objective has already been whitened to||x-x0||²/2, so in this case x=x0
is the exact solution. The helper now returns that analytic solution only
after checking the reconstructed equality/KKT residuals and coefficient signs.
It does not drop a varying inequality, clip a force or change a material bound.

## Additional source-tail diagnosis, before its calculation

Al99 is a published finite-cutoff TARGET potential. The candidate LJ lattice
is infinite, deliberately without that cutoff. At fixed perfect registry and
large opening a*, the reciprocal-zero attractive pair term is

    E_pair,cross = -v*pi/(2 A*) sum_(k>=1) k/[a*+(k-1)h*]^4 + O(a*^-8),
    sum = h*^-4 [zeta(3,q)+(1-q)zeta(4,q)], q=a*/h*,
    W_a,pair3 = 2*pi*v/(A* h*^5) [zeta(4,q)-(q-1)zeta(5,q)],
    W_a,pair3 ~ pi*v/(6 A* h*² a*³).

Here A*=sqrt(3)/2 at b*=1 is the dimensionless ATOMIC lattice area. The
exponential environment interactions die exponentially (times polynomial
factors); they are not converted to a cutoff either. The source force is
exactly zero once every cross-interface distance exceeds its tabulated
support. Its cell-area normalization remains the same, but its far-range
functional family is different.

Compare actual per-term forces at a/h=2.5,3,4,6,10,20,40,80 and the exact
Hurwitz expression in the asymptotic region. This diagnostic does not alter
the fit, its targets, or their roles. A far-tail discrepancy alone must not
be used to pronounce the analytic theory impossible; equally, it cannot
excuse incorrect forces at a/h≈1.1–1.3. Quantify which terms dominate rather
than assuming every disagreement is caused by the infinite LJ attraction.

The actual two-model tail audit passed the exact large-a pair-mean check with
maximum absolute error5.42e-20eV/L0. Al99 declares6.28721Å cutoff; the
implemented table support ends at6.286581279Å. At a/h=2.5 the positive v20 old
trial has1242.61MPa normal traction:29.86 from LJ,1212.74 from environmental
terms, versus4.06MPa from the source. At a/h=3 the source is exactly zero;
the candidate has454.69MPa, of which440.37 is environmental. Thus the large
disagreement is NOT mainly the unavoidable algebraic LJ tail. It points to
the fitted environment response over changing coordination. Do not add a
cutoff, refit the targets, or turn this attribution into a proof that an entire
embedding family is impossible. The final three-candidate, source/hash-bound
data are in `final_tail_audit/`; `tail_diagnosis/` and `tail_decomposition/`
are earlier completed snapshots, not replacements for the final audit.

## Completed shape refinement and identifiability

The two actual bounded shape searches finished recording 82 old-family and
162 power-family coefficient profiles, in 3126.49 wall-clock seconds. Together
with the 16 fixed controls this is 260 recorded profiles, in addition to the
selected-profile replay checks. No random search, new term, enlarged bound,
force clipping or source-target change was used. Here **old family means the
v19 p=0 many-body research family**, not the original production TwoRowLJ.

| Family | Selected positive-LJ trial loss | Minimum recorded closure loss | Outer optimizer stop |
|---|---:|---:|---|
| old, p=0 | 2215.203159 | 2012.434779 | ftol, 17 nfev / 13 njev |
| power | 2223.200136 | 2014.002476 | 30 nfev / 22 njev budget stop |

Both lowest profiles have u>0,v=0 and are ineligible LJ closures. The last
accepted optimizer losses are 2012.434839 and 2014.002599; finite-difference
trial evaluations account for the slightly lower recorded minima. The
reported positive trials are NOT the optimizer endpoints. Old-family
optimality is 1.408 despite its ftol stop; power optimality is 4.716 and its
optimizer did not converge. Neither result proves a global minimum or rules
out another admissible shape basin. Numerical coefficient-profile failures: 0.

The selected positive trials, saved without replacing any existing material:

| Parameter | old-family trial | power-family trial |
|---|---:|---:|
| u | .00469129023710 | .00673566805554 |
| v | .0107489011795 | .0198388390094 |
| A | 2.164270268 | 2.325506683 |
| B | 1.115115085 | 1.259870412 |
| C | 2.261620018 | 2.191847381 |
| D3 | 25.90773763 | 26.81474088 |
| D1 | 2.323302160 | 1.878201303 |
| D2 | 2.638465959 | 2.406179115 |
| D_E | 1.625425138 | 1.609646147 |
| K3 | 9456.761143 | 9989.000862 |
| k_scalar | 2.858861194 | 2.911456258 |
| k_odd | 6.746880127 | 6.679723848 |
| k_even | 10.63009070 | 10.91797438 |
| alpha_even | 94.80299452 | 94.35212010 |
| k_rank1 | 5.694618712 | 6.165522180 |
| power p | 0 (fixed) | -.9988414764 |

The original per-site normalization and coefficient units from v19 are
unchanged. The numerical coefficients multiply the same dimensionless
channels in E0=1 eV; radial lengths are in L0, and screening shapes are
dimensionless. Raw full precision is in the two `parameter_sets.csv` files.
For example the old positive trial maps to epsilon_LJ=.006157094889 eV and
sigma_LJ/L0=.8709396023. These are research trials, not a published Al fit.

All seven equality residuals are below 4.45e-14 in their declared scaled
units; coefficient KKT residuals are below 1.44e-11. At the nominal pristine
h*=sqrt(2/3), both candidates reproduce Haa=19.6609128633 and
Hxx=4.39359303512 eV/L0². Matching these two fitted targets is not held-out
validation. Al99's actual zero-traction root shifts slightly from nominal h
because the source table/rounded lattice parameter is retained, not forced
to have a zero source gradient by subtracting a linear term.

For constrained local sensitivity, differentiate E(k)c=e:

    E dc = -(dE)c.
    dc = dc_particular + T dz,  E T = 0.
    dy = M T dz + [(dM)c + M dc_particular] dk.

Apply the documented coefficient/shape and observation scales before SVD.
The seven exact rows have rank 7, leaving three coefficient tangent
directions, plus five/six shape coordinates. Saved-vector replay error is 0.

| Family | Rank / dimension | Singular values (normalized coordinates) | Condition |
|---|---|---|---:|
| old | 8/8 | 724.985,134.987,52.812,37.858,24.227,14.216,6.769,2.409 | 300.943 |
| power | 9/9 | 736.765,140.440,52.965,37.728,24.819,15.023,7.034,2.677,.1721 | 4280.321 |

Halving the deterministic derivative step from 2e-4 to 1e-4 changes the
Jacobian by at most 1.925e-4 / 1.900e-4 respectively. Differentiated equality
residuals are below 4.49e-14. These are local equality-manifold diagnostics:
active inequality cones, a statistical likelihood/confidence interval and
whole-family identifiability are NOT inferred. The added power coordinate
costs substantial conditioning and did not improve this bounded experiment.

## Final independent static response and material rejection

The final validator recomputed the source, v19 power parent and BOTH new
positive trials, using the SAME 48 predeclared excluded jets. Each field has
12 observations. The scales are 250 MPa for force; for energy the larger of
10% of the target magnitude and .01 J/m²; for curvature the larger of 10%
and .1 eV/L0². They are declared discrepancy scales, not measured standard
deviations or experimental confidence bounds.

| Candidate | Energy RMS | Normal-force RMS | Haa RMS | Hxx RMS |
|---|---:|---:|---:|---:|
| v19 power parent | 1.0025 | 5.6153 | 6.0263 | 6.1650 |
| v20 old positive trial | .8786 | 7.8649 | 6.0596 | 5.1335 |
| v20 power positive trial | .8752 | 7.8780 | 6.0714 | 5.1376 |

Both v20 candidates have only 1/12 normal-force observations within its
declared scale, versus 2/12 for the parent. Their small energy RMS does not
validate force or curvature. **Neither candidate passes the material gate.**

The validator also solved 11 explicitly oriented normal/two-shear/mixed
tensor scenarios, with signed fractions 0,.5,1,.5,0,-.5,-1,-.5,0.
All 396 states across four models have verified stable local roots.

| Applied interface traction | Source response [Å] | v19 power | v20 old | v20 power |
|---|---:|---:|---:|---:|
| 50 MPa normal, delta a | .0009253659 | .0006562111 | .0009286880 | .0009286672 |
| 4 MPa shear1, delta x | .0003309979 | .0003272311 | .0003309946 | .0003309946 |
| 4 MPa shear2, delta y | .0003309325 | see raw table | .0003308315 | .0003308295 |

The normal primary-response error improves from -29.09% to +.3590%/+.3568%.
The primary shear1 errors are -.0009855%/-.0009839%. This is a real scoped
compliance improvement and is largely expected from the fitted initial
tangents; it is not independent Al yield validation. All secondary responses
are exported too. For example normal opening under shear1 is 8.67e-8 Å for
the old trial versus 2.57e-8 Å for the source, still an appreciable relative
mismatch in that small second-order quantity. Nothing was multiplied to
make plastic strain visible. The paths return statically; no dynamic hold,
residual plasticity, fatigue, transition rate or physical Hz is established.

| Static landscape observable | Source | v19 power | v20 old | v20 power |
|---|---:|---:|---:|---:|
| relaxed fault energy [J/m²] | .150479 | .173907 | .122424 | .123257 |
| adjacent saddle energy [J/m²] | .172002 | .197371 | .162561 | .162566 |
| W(a=40h,0) [J/m²] | 1.741285 | 1.901684 | 1.883279 | 1.887106 |
| first fixed-registry normal peak [MPa] | 12969.568 | 10194.649 | 10338.973 | 10338.437 |
| peak a/h | 1.185846 | 1.182816 | 1.125712 | 1.126205 |

The last row locates a fixed-registry traction extremum, not a coupled
spinodal or actual yield. W(40h,0) is a finite large-opening sample, not
exactly infinite separation for the LJ candidate. All stationary Morse
indices and gradients are reported. Source cutoff-shell features include a
later negative rigid-opening traction; the source is not presumed monotone.
Two root-bracketing resolutions agree closely on the FIRST peak but a later
source minimum shifts from a/h=1.631682 to 1.632270 (traction changes only
.00201 MPa). We do not call every source extremum location converged.

The image `response_report/opening_tradeoff.png`, inspected after creation,
shows the corrected initial tangent and the remaining finite-opening error,
including the new candidate's curvature shoulder over the displayed a/h<=2
domain. The full curves, including the larger-opening states, remain in CSV.

## What the error attribution permits, and the next gate

At a/h=2.5, the old trial's 1212.74 MPa environmental force contains
1093.50 MPa from scalar embedding, 91.62 from rank1, 21.83 from rank3,
.045 from rank2, 1.23 from Eg and 4.52 from the quartic term. The LJ net
force is only 29.86 MPa. The source gives 4.06 MPa. Thus "infinite LJ versus
finite cutoff" is not the main explanation of this large mismatch.

This identifies the partial-coordination scalar environmental response as
the next specific question, not an excuse to add an arbitrary polynomial.
Read the previous v14 cubic-density and v15 positive-mixture ablations before
proposing any extension. First derive their remaining derivative freedom
after bulk/tangent equalities, test rank and source force/curvature budgets,
and retain independent held-out states. A tested pair of candidates does not
prove the whole analytic family insufficient. No new term is added in v20.

## Verification, reproduction and preserved boundaries

Final independent validation took 1221.61 wall-clock seconds while other
work was running. The maximum finest-step analytic-versus-FD Hessian error
is 3.753e-7 eV/L0². Interface tolerances 2e-11/2e-13 and independent per-site
real-space checks are recorded separately from full nonlinear displaced-site
Bloch Hessian checks. Five dilations (.98,.995,1,1.007,1.02), radii12/16 and
243 independently sampled q points give minimum tail-subtracted margin
.00565074 eV/L0². This is not a proof over every q and strain. Tiny-step
energy differences also show roundoff; the refinement tables retain it.

Final executed regression: targeted 47 passed in 59.99 s; full solver
565 passed in 923.85 s; app 34 passed in 86.14 s; desktop smoke exit 0.
Failures, errors and skips: 0. JUnit logs are locally retained under
`.cache/tangent_calibration_v20/`; `verification.json` records the summary.
Execution timings are computer wall-clock timings, NOT calibrated PDE time.

Reproduce with Python 3 (use `py -3` when available), from the repository.
Every runner requires a NEW output directory and refuses to overwrite one.
The following paths illustrate a separate replay, not replacements for the
committed experiment:

```powershell
$parentFit = 'results/fcc111_active_interface/coordination_screening_v19/joint_refinement'
$newRun = 'results/fcc111_active_interface/tangent_calibration_v20_reproduction'
py -3 -m solver_v1.run_interface_tangent_calibration --parent-joint $parentFit --out "$newRun/controls"
py -3 -m solver_v1.run_tangent_shape_refinement --parent-joint $parentFit --out "$newRun/shapes" --max-nfev 30
py -3 -m solver_v1.validate_tangent_calibration --models "$parentFit/power_family" "$newRun/shapes/old_family" "$newRun/shapes/power_family" --out "$newRun/validation"
py -3 -m solver_v1.report_tangent_calibration --parent-joint $parentFit --refinement "$newRun/shapes" --validation "$newRun/validation" --out "$newRun/report"
py -3 -m solver_v1.audit_tangent_source_tail --models "$newRun/shapes/old_family" "$newRun/shapes/power_family" "$parentFit/power_family" --out "$newRun/tail"
```

The first two commands run optimization; the report replays saved vectors
and computes sensitivity, not a new fit. Intermediate `checkpoint.json`
and `progress.json` deliberately retain completed=false snapshots. Final
`completion.json`/`calibration.json` and the per-optimizer stop fields govern
completion; scientific acceptance is a separate false flag.

Source binding is Al99 SHA-256
`60c8a085be79d273324ab421f5b1447578fef55c1acfc6492c0999f15ee8a284`.
The old production energy registry, static material files, conservative PDE,
UI, specimen aggregation and kinetic JSON are unchanged. Physical a/s
mobility and t0 remain unavailable. This is completed scoped calibration
research with a negative full-material assessment, not a finished Al fatigue
solver. No production, physical-time, actual-yield or specimen UI gate opens.
