# Stable-branch domain continuation and material-conflict audit

Research stage `stable_core_v13`, continuing `range_core_v12` without changing
the LJ/Bessel potential, production PDE/UI, historical parameters or kinetics.

## Declared core comparison before execution

A centered force-converged stationary saddle is not an admissible reference
for stress-induced residual deformation. Start from the independently checked
lower zero-load core. For free disks I_R contained in I_Rprime, initialize

    u_Rprime^0(i) = u_R(i) + u_affine(delta_tau,i),  i in I_R
                 = u_Volterra(i) + u_affine(tau_new,i), otherwise.

The material hash, physical length, logical lattice indices and prescribed
Volterra center must match. This only initializes a new calculation: every
affected site's environment, per-site embedding, force and full vector
Hessian is recalculated on the target disk/neighborhood. It does not certify
the omitted infinite transverse changes. The x-row sum remains infinite and
analytic. Copying an existing free displacement is not interpolation, force
clipping or an empirical constitutive correction.

Refine the free radius and transverse environment ring independently. Require
force and Morse checks on each final state. Compare the SAME inner rows, both
actual displacement and analytic forces under the refined generator. Compare
the declared Taylor-remainder radial energy partition with the same-potential
outer logarithm. Never use a larger unstable state's plausible energy plateau
as a substitute for stable-branch convergence.

Selected sequence starts from historical R6/r6 stable+, then R8/r6; follow
with ring and larger-disk refinement subject to measured computational cost.
Keep stable small-disk static load/unload separate from infinite-domain motion.
No core minimum alone supplies a finite-source geometry or 0.2% yield strain.

### Optional analytic Newton continuation

The previous L-BFGS iteration is retained as the default. A separately selected
Newton-CG method uses the SAME analytic Hessian-vector operator and energy,
with the same independent final force/Morse criteria. It is an optimization
comparison, not a harmonic core approximation or accuracy reduction. Results
must agree on the same stationary branch before using it for larger domains.

### Smooth radial matching, with its exact subtraction convention

A hard radial sum jumps whenever a discrete row enters the integration disk.
Keep that raw diagnostic, and additionally use a declared numerical window
w_alpha(t): one for t<=alpha, zero for t>=1, and 1-3x^2+2x^3 between, where
x=(t-alpha)/(1-alpha). This window is NOT a fitted material length or core size.

For the already defined Taylor-remainder site partition,

    E_c,b^w(R) = sum_i w_alpha(r_i/R) e_i^remainder
                - k_row [ln(R/b) + c_w],
    c_w = integral_0^1 [w_alpha(t)-1]/t dt.

Writing delta=1-alpha, direct polynomial integration gives equivalently

    c_w = -6 delta sum_(n>=0) delta^n / [(n+1)(n+3)(n+4)].

The omitted positive absolute terms have ratio below delta, so the geometric
tail bound is explicit. For alpha=1/2, c_w=19/6-5 ln(2). Independent quadrature
tests check both forms. The continuum log coefficient is never fitted.
For alpha R>=b the subtraction is exactly k integral_b^R w_alpha(r/R) dr/r.
Compare alpha=1/2 and 3/4, radial extent, free domain and environment ring.
Window smoothing cannot turn an unstable core or boundary bias into validated
matching. It changes no force, core state or total fixed-boundary energy.

## Material audit, separate gate

The previous candidate matches neither all independent cubic constants nor
held-out interface curvature. Decompose its error into hydrostatic, tetragonal
and shear symmetry channels before adding another analytic term or changing
a loss metric. An odd STF moment is zero in centrosymmetric affine bulk; it
cannot independently repair a bulk elastic constant. A perfect-cubic rank2
moment is zero but its deviatoric derivatives contribute to elasticity.

At the actual v12 candidate B=75.8374 GPa, C'=(C11-C12)/2=8.75937 GPa,
C44=32.8059 GPa, versus source 79.3333/26/32. Its rank2 column contributes
only 0.306852 GPa to C' but 29.7808 GPa to C44. Uniformly increasing D2
therefore cannot repair the tetragonal channel without changing others.
These are an exact decomposition of the total strain jets, not independently
equilibrated per-term materials.

Before adding any term, reprofile the SAME eight coefficient columns and
three ranges. At fixed ranges this is a convex quadratic least-squares problem
with two exact linear constraints (reference force/cohesion) and the unchanged
nonnegative sector u,v,A,C,D3. Enumerate the 32 active faces: solve equalities,
minimize on each null space, check bounds, then choose the lowest feasible
residual. Active coefficients set to zero are the face parameterization, not
post-fit clipping. Retain equality/roundoff and KKT residuals.

The outer log-range least-squares runs three declared deterministic starts,
using the same bounds and target scales as v12. Held-out values are excluded.
This separates an optimization-stopping problem from an energy-family conflict.
No radial global-optimality or family-impossibility theorem is claimed.

Any proposed extension must have a zero-extension limit, exact energy jets,
gauge/rank analysis, actual deterministic fit and independent stability checks.
An empirical yield target or desired slip ordering never enters the loss.
The source remains a matched 0 K atomistic comparator, not measured fatigue.

### Spectral admissibility inside the coefficient problem

Positive C11-C12 and C44 only check the long-wavelength limit. At fixed
radial ranges the independent finite-q bulk Hessian is also coefficient-linear:

    H(q;c) = sum_j c_j H_j(q).

For a normalized negative-eigenvalue polarization v at a sampled wavevector,
positive-semidefinite H requires the necessary linear halfspace

    sum_j c_j [v^T H_j(q) v] >= 0.

The research runner adds the most violated such halfspace, re-solves the
coefficient quadratic problem, and rechecks all declared q. No force, energy
or target is clipped; there is no fitted instability cutoff. On each active
face E c = r, write c=c_particular+Z z with E Z=0 and minimize the same
normalized target residual over z. Enumerating independent active faces
selects the global convex minimum at THOSE fixed ranges and those finite
constraints. Bounds u,v,A,C,D3>=0 and the force/cohesion equalities remain;
the LJ pair must independently remain strictly positive. B,D1,D2 remain the
same signed research sector. KKT/feasibility residuals are saved.

The static operator uses the independently checked direct lattice evaluator,
not a replacement finite-cutoff energy. For each of the two declared radial
sets, probe Gamma-X/L/K at fractions .02,.1,.25,.5,.75,1 in the corrected +ABC
cubic frame. Add at most six separating planes, saving a nonconverged result
if the budget does not remove all violations. This is NOT radial optimization
under spectral constraints and NOT a proof over the whole Brillouin zone.
Re-evaluate final coefficients at independent radii 6,8,10 L0. A numerically
zero eigenvalue at one radius is not a positive stability margin. The matrix
norm change between radii is reported separately from roundoff.

`run_spectral_range_continuation.py` then performs a DISTINCT actual radial
optimization, starting at both completed radial basins, at most 30 nonlinear
function evaluations per start. It retains the offending q/polarization as
one necessary halfspace and recomputes its H_j at every new radial vector.
It does NOT freeze a halfspace coefficient vector from an old energy. As the
polarization is fixed, satisfying this condition is still not full PSD. All
declared paths and radius refinements are independently reapplied at the end.
The energy/target family and bounds remain unchanged. Optimizer exhaustion
and nonzero optimality are explicitly retained, not relabeled convergence.

### Units and scope of the load test

L0=4.05/sqrt(2) Angstrom, b*=1. For the historical unchanged candidate,
k_row=0.277580619166 eV per straight-row repeat multiplies ln(R/b). Applied
resolved shear is converted as tau*=tau_MPa 10^6 L0^3/(1.602176634e-19 J/eV).
The affine exterior solves C_antiplane grad(u_x)=(0,tau*), i.e. sigma_xy=0,
sigma_xz=tau. Its tensor comes from the SAME potential, not source Al moduli.
This is a prescribed static shear experiment with a straight pre-existing
screw and a fixed Dirichlet exterior. It is not tensile 0.2% proof stress,
measured source pinning, an atomistic loop activation energy, or a physical
zero-stress time hold. There is no mobility, time or statistical area in it.

## Completed numerical evidence and remaining scope

The actual results are stored under
`results/fcc111_active_interface/stable_core_v13/`; v12 evidence is unchanged.
The stable-core material is the HISTORICAL `angular_monotone_opening` candidate,
not either new fit. Its parameter-file SHA256 remains
`9d00fbf54831c134fe9961e17f0603defdc7958bc3590743b2094264a31b6655`.
It is still not a quantitatively accepted Al material. Never silently apply a
new fit to an old core state/far field.

| Free radius / ring | Force residual [eV/L0] | Minimum H [eV/L0^2] |
|---|---:|---:|
| 6 / 6, independent Newton replay | 5.02e-10 | 0.4045673 |
| 6 / 8 | 1.58e-7 | 0.4058644 |
| 8 / 6 | 3.35e-10 | 0.2613280 |
| 8 / 8 | 1.60e-9 | 0.2615344 |
| 10 / 6 | 3.00e-10 | 0.1723491 |
| 10 / 8 | 2.80e-10 | 0.1723558 |

All are verified fixed-boundary minima, winding one, not centered saddles.
The same R6/r6 Newton/L-BFGS comparison differs by 1.8e-13 eV in energy and
1.48e-7 L0 on the 18 common inner rows. R6->8->10 at fixed ring6 instead
changes those rows by .0047665 then .0027026 L0: optimizer error is much
smaller than remaining free-domain error. Ring6->8 changes are .0001982
and .0002693 L0 at radii6 and8. These are convergence differences, not
statistical confidence intervals or measured core-displacement tolerances.

The smooth finite part at each full free radius, ring8, is:

| Radius/L0 | alpha=0.5 [eV/row] | alpha=0.75 [eV/row] |
|---|---:|---:|
| 6 | .36928403 | .37703345 |
| 8 | .36908278 | .37263122 |
| 10 | .36736394 | .36951753 |

The much smoother curve does not establish its infinite-domain limit.
At the SAME integration radius6, changing free radius8->10 changes the
alpha=.5 value from .36014162 to .35605935 eV. Do not hide that boundary
effect by comparing only each disk's outer endpoint. Raw/hard partitions,
common-radius ring/domain changes and both window constants are also saved.
Some runs include large elapsed wall-clock gaps of unverified origin;
iteration/evaluation counts are retained, and wall time is not a CPU benchmark.

At radius8/ring8, 0->50 MPa changes the common inner rows by .00809328 L0,
with force4.87e-10 and minimum H+.2609368. Newton reports a precision-loss
warning after seven iterations/40 evaluations; the independently evaluated
force and energy-curvature tests nevertheless pass. Both facts are retained.
Static unloading is evaluated separately; it must not be called a kinetic hold.
The completed unload returns to the zero-load field within7.10e-10 L0 on
the same inner rows; final force4.20e-10, minimum H+.2615344, energy difference
1.48e-14 eV. The old smaller stable disk likewise recovered within1.01e-7 L0.
There is no resolved residual registry change in these fixed-boundary tests.
The smaller case used the opposite symmetry-related stable branch, so its
loaded increment is a legacy comparison, not a same-branch domain certificate.
In particular none of these numbers is a macroscopic epsilon_p or yield stress.
Radius8 and10 correspond to free radii only2.291 and2.864 nm, with the
exterior Volterra center fixed. Image/restoring forces from this boundary can
inhibit translation. Recovery at50MPa in this setup does NOT establish that
real Al stays elastic at50MPa. A translating-core barrier needs a compatible
outer response and a declared work/translation coordinate, followed by actual
finite-source geometry. Likewise a smaller minimum fixed-disk Hessian value
in a larger domain can be a longer-wavelength elastic mode; it is not by
itself a reduced Peierls/yield stress.

### Optimization failure versus physical failure

The same-family unconstrained range study made 588 distinct profile calls
(the saved array is authoritative if this count changes on a rerun).
Its three actual outer stops are:

| Start | Final loss | nfev/njev | Actual stop |
|---|---:|---:|---|
| previous v12 C-fit | 38.90653745 | 11/8 | relative loss tolerance |
| previous v12 exact-C radial start | 23.91482825 | 14/12 | relative loss tolerance |
| declared [2.7,6.4,10] | 49.19626938 | 65/63 | evaluation budget exhausted |

Neither a relative-loss stop with nonzero optimality nor a budget stop is a
proof of global convergence. The best second basin has decays
(5.5058225511,4.8746400515,3.4803990843), coefficient order
(u,v,A,B,C,D3,D1,D2) and values
(.77912067894,7.74231038585,0,48.2370475974,.36977545858,
117.088000214,-14.3765010903,-145.673598913).
It satisfies the two exact constraints to1.78e-13 normalized residual and
free KKT residual1.21e-11. This makes it a good mathematical profile, NOT a
good material: C=(99.6492,66.0161,31.1875)GPa, held-out normalized RMS19.56,
and a large resolved negative finite-q eigenvalue near -76.96 eV/L0^2.

Its local nine-direction scaled sensitivity singular values are approximately
(1288.16,177.38,120.78,44.49,24.92,18.67,9.00,2.63,.3504), condition3676.
The finite radial-difference step check changes entries by1.04e-4. The active
A=0 face, coefficient scaling, model discrepancy and physical instability
preclude confidence/unique-material claims from this algebraic rank.

At those fixed ranges the necessary spectral cut changes the optimum loss to
103.622311, C=(95.6664,69.2851,32.1251)GPa. Its sampled minimum at direct
radius8 is about zero; at radius10 it is -.0018605 with radius8->10 matrix
change .0041120 eV/L0^2. Classify this as marginal/unresolved, not a stable
material. A separate actual radial continuation is performed with that
necessary condition retained; its outputs must not be confused with this
fixed-range failure. Failure of these particular ranges is not a theorem
that every admissible analytic LJ family must fail.

The actual necessary-spectral radial continuation completes99 profile calls.
The two outer starts stop at losses38.90654191 and31.62468250 with16/12 and
17/10 nfev/njev, respectively. The second has decays
(5.3707449113,5.6320174464,3.2556291352), coefficients
(.41306528163,3.49172990503,4.09109696972,23.4600618913,0,
65.6450471980,-2.11694765208,-58.3030883841).
Its C=(98.3205,68.0336,30.6618)GPa and C'=15.1434GPa still miss the source
114/62/32GPa and26GPa. Held-out normalized RMS17.3383 is worse than v12's
9.4344. KKT residual3.87e-13 checks the inner problem, not global nonlinear
optimality. Local nine-direction condition1269 is not physical identifiability
or a confidence interval; radial-step change3.15e-5 is saved independently.

The optimized own-interface stationary fault/saddle energies are .1415213
and .1705550 J/m2, versus matched source .150479/.172002. Those plausible
energies alone do not rescue the inaccurate elastic/saddle-shape response.
All three local perfect/fault/saddle roots have the expected Morse index.
That is still not bulk finite-q stability.

**The apparent marginal spectral sign is resolved by a further check.**
At the actually offending Gamma-K halfway wavevector, extend the independent
direct radius10->12->16->20. Fixed-range and reoptimized candidates have
minimum H=-.00281457 and-.00354150 eV/L0^2 at radius20, while the last matrix
changes are6.36e-5 and8.00e-5. The negative values are about44 times these
last changes and persist under refinement. Reject both candidates, rather
than reading the radius8 constrained zero as stability. This is empirical
tail convergence, not a proved absolute tail bound. The refinement table
keeps the earlier radius10 unresolved classification separately from the
final negative-mode finding.

Thus this step repairs the optimization/admissibility diagnostic, but does
not adopt a new Al material. Future fitting must impose a tail-controlled
stability margin over an adaptively verified q set DURING radial optimization,
not merely at one finite-radius operator. Neither arbitrary positive margins
nor force clipping are appropriate. Further expressivity must be justified
only after this robust optimization and independent channel/gauge analysis.
The stronger radius study does not prove that all radial starts or all
minimal analytic embedding extensions fail.

The 5%-C and10%-interface scales are declared MODEL DISCREPANCY scales, not
experimental confidence intervals. In particular the small source saddle Hxx
has a small normalization denominator; always inspect its dimensional error
and the full saddle Hessian as well as normalized RMS. All held-out values
remain excluded from the fitted loss. The source is the same matched 0 K
Mishin Al99 comparator with SHA256
`60c8a085be79d273324ab421f5b1447578fef55c1acfc6492c0999f15ee8a284`.

### Reproducibility and interpretation

Actual optimizations: `run_profiled_range_calibration`,
`run_spectral_material_profiles`, `run_spectral_range_continuation`.
Existing completed paths are protected against overwrite. Intermediate
duplicated optimizer snapshots live in ignored `.cache/stable_core_v13/`;
complete profiles/stops are retained in each final calibration JSON.
Actual core minimization: `run_isolated_screw_core` with explicit independent
source/output roots; each case records the source-state and material hashes.
`report_stable_core_diagnostics` and `report_stable_material_diagnostics` only
REPLAY stored outcomes under their actual analytic generators and regenerate
tables/plots. They must never be described as rerunning the optimization.

Material acceptance, full-domain/all-character core matching, a measured
finite source and macroscopic0.2% deformation remain separate gates. There
are no new empirical energy terms, stress clipping, fitted pin lengths,
activation areas or mobility adjustments in this step. The LJ base/Bessel
math, production PDE/UI, historical calibration and A_c semantics are intact.
Physical mobility, seconds/Hz and actual Al yield validation remain unavailable.

## Executed verification

- Combined targeted new/old core and profile tests:34 passed,42.53s. The21
  new tests cover continuation/provenance, method agreement, window integrals,
  exact coefficient profiling/held-out isolation and finite-q linearity/cuts.
- Full `pytest solver_v1 -q`:415 passed,1810.82s, including all prior394 tests.
- `pytest app -q`:31 passed,200.25s; no skips.
- `app.desktop_ui --smoke`:passed,2.51s; historical LJ a0=.7713438268704838
  and kappa=86.29296488740997 unchanged.
- `git diff --check` and staged diff check:passed.
- Generated evidence:36 JSON,45 CSV and1 SVG parsed/checked; actual plot inspected.
- Python3.13 was used because the `py` launcher is unavailable. BLAS/OMP
  threads were one in these processes. Some tests ran concurrently with
  research calculations; their wall times are not isolated performance benchmarks.

No ongoing calculation is represented as completed. Raw optimizer warnings,
budget stops and the rejected negative eigenvalues remain in the evidence.
