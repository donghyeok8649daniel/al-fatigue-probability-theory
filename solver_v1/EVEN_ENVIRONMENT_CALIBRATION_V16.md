# v16: finite-distortion even-environment response

Research only. This follows v15's measured interface failures, not a desired
plasticity or fatigue result. Original LJ, all infinite Poisson/Bessel sums,
per-atom counting, bulk targets, kinetic status and production defaults stay.
The source remains the matched 0 K Mishin Al99 target-only reference.

## 1. Mechanism diagnosis before extending

At the previously inspected direct-registry point (a=h,s1=.39 L0,s2=0), the
v15 stable quartic candidate's Haa [eV/L0^2] decomposes as

    LJ repulsive 499.139; LJ attractive -407.450; scalar embedding 68.539;
    rank3 20.489; rank1 25.093; rank2 102.195; combined Eg 82.356;
    quartic rank3 13.215; total403.575 versus source36.317.

Thus saturating only the rank3 quartic term, already tested in v15, does not
directly address the largest even-environment contributions. This budget is
not proof that the pair/other terms are correct or that one extension will fit.

## 2. One nested shape with a fixed amplitude gauge

For each atom let Q2 and QE be the existing full rank2 tensors; QE is the
previously derived radial combination with zero T2g response at the reference.
The site environment is summed over all neighbors BEFORE nonlinear functions.
Define Ij=Qj:Qj and the fixed reference normalization

    Nj^2 = 2 ||partial_gamma Qj|reference||^2,
    N_E^2=1 in the existing Eg gauge.

gamma is the existing [111] engineering-shear reference mode. It specifies
an amplitude convention, not a fitted specimen orientation or length. The
same Nj is held fixed under ALL subsequent states/strains. Rescaling Qj and
Nj by c and Dj by 1/c^2 preserves its physical site energy. The dimensionless
normalized distortion Ij/Nj^2 is invariant under that amplitude gauge.

Test a single common alpha>=0:

    f_j(I)= I/(1+alpha I/Nj^2),
    E_even,i=D2 f_2(I2,i)+D_E f_E(IE,i).

All other energy terms, including K3 ||Q3||^4, remain the original family.
This is a Padé-type analytic research hypothesis for the missing nonlinear
angular response, NOT a known Al parameterization, standard MEAM, or a yield
cutoff. No tabulated pair potential or empirical plastic-flow law is added.
The broad use of nonlinear invariant functions is discussed in
[Drautz2019, PRB99,014104](https://doi.org/10.1103/PhysRevB.99.014104);
that paper does NOT supply or validate the present specific rational law.

For k=alpha/Nj^2, f'=1/(1+kI)^2, f''=-2k/(1+kI)^3. With coordinate indices p,q,

    I_p=2 Q:Q_p,
    I_pq=2 Q_p:Q_q+2 Q:Q_pq,
    E_p=D f' I_p,
    E_pq=D(f'' I_p I_q+f' I_pq).

The neighbor derivatives are unchanged analytical reciprocal derivatives.
There are no poles for I,alpha>=0. f is nonnegative/increasing/bounded when
alpha>0, but concave in I. A bounded scalar term does NOT guarantee a positive
crystal Hessian. Signed D2 remains subject to total-energy stability tests.

At a cubic reference Q2=QE=0, f=I+O(I^2). The difference starts at fourth
order in displacement. It changes neither cohesive reference, equilibrium,
cubic elastic constants, nor ANY harmonic finite-q column. This also holds
at isotropic cubic dilation. It does NOT hold at an arbitrary sheared/noncubic
state with Q_bulk!=0. No noncubic Hessian may silently reuse the cubic formula.
The pre-existing analytic bulk tail certificates therefore remain applicable
to the cubic harmonic operator, not automatically to every nonlinear state.

More explicitly, for a small atomic displacement field u, Q_j=L_j u+O(u^2)
and the added change is -D_j alpha (Q_j:Q_j)^2/N_j^2+O(u^6) as a function
of the invariant. Its displacement expansion starts at order4 (and may have
higher odd orders through Q). The second variation is identically zero for
every displacement pattern, not only a tested homogeneous strain. Exponential
moment kernels justify termwise differentiation for bounded sufficiently
small displacements away from atomic collisions. No finite-q frequency or
kinetic information follows from this static statement.

**Large-alpha warning:** the preserved tangent is a nonuniform limit.
For a path deltaq with I~||Q_q||^2 deltaq^2, the saturation crossover is

    |deltaq| ~ N/(sqrt(alpha) ||Q_q||).

As alpha increases, a perfectly matched harmonic curvature may act only in
an extremely narrow neighborhood. Taking alpha->infinity before deltaq->0
then differs from the reverse order. This can create a narrow traction peak;
it is NOT evidence of a successful correction from ideal to real yield.
The crossover derives from the fitted tensor jet, not an invented source
length. Independent opening checks combine logarithmically refined offsets
near h with the ordinary far-branch grid, at two resolutions. A large-alpha
fit is not accepted from its reference Hessian alone.

## 3. Geometry and normalization

For the active interface, the site tensors are the convergent cumulative
cross-layer changes about the zero-tensor cubic reference. Lower and upper
half-crystals contribute equally by inversion of the neighbor environment:

    W_even = 2 sum_depth [D2 f_2(I2,depth)+D_E f_E(IE,depth)].

The factor2 is a site multiplicity, not pair double counting. It is not
f(sum_depth I). Units remain eV/interface cell with L0 fixed; no A_c enters.
The large-a limit and full registry period must be checked on the total model.

## 4. Predeclared experiment

All previously inspected v15 holdouts become DEVELOPMENT data. The new target
generation adds36 off-grid/off-path holdout jets, never used in the loss. The
same source condition, numerical target scales and exact five bulk constraints
are retained. Compare the old family AND the new family on this SAME dataset;
the numerical value of a v15 loss is not comparable to this larger loss.

Run deterministic coefficient profiles, then bounded shape optimization. At
fixed decays/alpha all energy coefficients remain linear; the same QP with
primal/KKT certificates and analytic finite-q tail inequalities is used.
No tolerance is relaxed. alpha=0 exactly reproduces the old family; positive
alpha is sampled in log coordinates, with the zero boundary separately tested.
SVD must use the exact five-bulk-constraint tangent and include shape columns.

Independent tests: analytic derivatives, per-site vs global counting, gauge
and rotation invariance, alpha=0, direct/reciprocal neighborhood convergence,
nonzero-alpha direct sinusoidal site energies, coefficient replay, stationary
Morse indices, full registry/opening curves, reciprocal/layer tolerance,
finite-q radius/grid refinement and small cubic dilation. A good training
loss or matching Haa at one state cannot authorize production adoption.

Physical Ma/Ms/t0/seconds/Hz and actual specimen yield remain uncalibrated.
The next core/source/PDE step is conditional on these static checks, not on
the desire to make plasticity visible.

## 5. Open-boundary LJ and feasibility controls

The coefficient QP uses the CLOSED nonnegative cone, whereas a real LJ pair
requires u>0 and v>0 in u r^-12-v r^-6. The corresponding constants are
sigma/L0=(u/v)^(1/6) and epsilon=v^2/(4u). A profile with either coefficient
zero cannot define this finite positive LJ pair. Its low residual is only a
boundary diagnostic, not a fitted LJ material. The loader rejects it. If all
profiles lie on that boundary, the calculation saves that negative result
explicitly instead of failing at min(empty admissible set).

After the first free-saturation run approached that boundary, a separate
embedding-only CONTROL was declared: retain the two positive LJ coefficients
from v15's strain_stable_quartic. Both the old and new angular families are
tested with that same pair. These two constraints are not new measured Al
targets and do not establish that the inherited pair is uniquely calibrated.
They are included in the exact-stage tangent/SVD and excluded from the loss.
No arbitrary tiny positivity bound or different pair functional form is used.

Some exact-bulk/fixed-pair shape trials have no coefficient vector satisfying
the declared constraints. A trial-domain rejection is not a proof that the
whole material family is impossible. A separately recorded continuation uses
bounded deterministic Nelder-Mead in the same log-shape coordinates, with a
.025 initial simplex step and the same objective/certificates. An undefined
profile receives the extended-real +infinity indicator (rejected), never a
finite artificial penalty or a repaired unstable model. Feasible optima,
infeasible trials, unverified numerical profiles and budget termination remain
distinguishable; no global-minimum proof is claimed.

## 6. Completed experiments and quantitative outcome

Seven actual runs produced420 certified coefficient profiles, not420 material
fits. All attempts, rejected trials, derivative probes and optimizer stopping
reasons are saved in `results/fcc111_active_interface/even_environment_v16/`.
The data generation has5 exact bulk rows,60 fitted interface rows and36
non-fitted validation rows (12 states x energy/normal force/normal curvature).
Pair-controlled runs add2 exact CONTROL rows, never two fitted observations.
The new validation states were declared before fitting. They were subsequently
inspected; they are not repeatedly rebranded as fresh blind tests.

| Run | Certified profiles | Fitted loss | Validation normalized RMS | Interpretation |
|---|---:|---:|---:|---|
| baseline |75|31834.67959|19.72518|Old quartic family, new common dataset |
| saturated |192|616.15700|4.00688|Best strictly positive LJ profile; not final unconstrained minimum |
| fixed_baseline |17|112542.88734|41.00151|Old v15 pair retained; trial left admissible coefficient domain |
| fixed_baseline_continuation |33|105727.41455|39.18280|48-evaluation bounded simplex,15 rejected profiles |
| fixed_saturated |54|55543.47099|28.01651|Archived pre-fix tiny negative K3; constructor correctly refused |
| fixed_saturated_reprofile |1|55543.47099|28.01651|Same shape, exact-boundary coefficient re-solve only |
| positive_pair_refinement |48|596.14414|4.13966|Retain positive pair selected by training loss; local simplex section |

The192-profile free search contains172 zero-pair closure profiles and20
positive-LJ profiles. Its lower closure loss553.90588 is NOT an accepted LJ
result. No positivity threshold was invented to disguise that boundary. Both
free-shape starts used their24-evaluation budget, not a proved optimum. The
last positive-pair continuation lowers training loss but slightly worsens
validation: report both rather than select by the validation score. Neither
surface is accepted as a complete Al material or a production PDE generator.

### 6.1 Actual positive coefficients and units

The last positive-pair section has shape

    (k_scalar,k_odd,k_even,alpha)
      = (2.612878645739,6.427100269117,11.957043235971,425.060826361)

and coefficients in the existing normalized site-basis order

    (u,v,A,B,C,D3,D1,D2,D_E,K3)
      = (0.01719039248874,0.10069343124175,2.024529775081,
         1.868099716373,2.580138379998,37.782668436283,
         4.695053645070,1.995355709895,1.581882666800,18140.5041598844).

Energy remains eV/interface cell, coordinate scale L0=2.863782463805517 Angstrom,
and A_atomic_cell=sqrt(3)L0^2/2. The pristine spacing is sqrt(2/3)L0, not the
historical two-row a0. Per-atom cohesion 3.36 eV and cubic C=(114,62,32) GPa remain
exact to numerical precision. During strains neither material coefficients,
reference density, L0 nor the moment gauge is recalibrated. A_c is absent.
The angular amplitudes' units refer to their explicitly normalized site bases;
their numerical sizes are not directly comparable to an ordinary pair epsilon.

Exact-constraint tangent ranks/condition numbers are 8/177.97 (baseline),
9/146.44 (free positive saturated) and 7/122.53 (last pair-controlled section).
The last two exact pair controls reduce the manifold dimension; the improved
condition number is not two new experimental data. The old-pair saturated
section is very poorly conditioned (~8.77e5). Even decay remains close to its
declared upper bound12. These are local sensitivity diagnostics, not confidence
intervals, unique global identification or an analysis of active inequality cones.

### 6.2 What now agrees, and what does not

All 12 non-fitted ENERGY observations of the positive free saturated candidate
fall within the declared 10%/0.01 J/m² scales (maximum normalized error 0.93468,
RMS 0.55221). This is real energy-level progress, not force/kinetic validation.

| Static quantity | Matched Al99 source | Old-family baseline | Saturated | Positive-pair refinement |
|---|---:|---:|---:|---:|
| Normal-relaxed fault energy [J/m²] |0.150479|0.373858|0.167822|0.163016|
| Full-vector saddle energy [J/m²] |0.172002|0.421187|0.190537|0.186420|
| Perfect Haa [eV/L0²] |19.660913|41.038072|28.010901|27.916359|
| Perfect Hxx [eV/L0²] |4.393593|9.043722|4.414200|4.762392|
| Saddle Haa [eV/L0²] |12.717338|36.120639|24.802905|24.112704|

The fault/saddle values are energies of the actual verified Morse-index0/1
stationary states, not values at unchanged source coordinates. The registry
curves in the figure instead explicitly hold a=h; do not mix the conventions.
At40h the work-of-separation approximation changes from0.59175 to1.87166J/m²
for the saturated candidate, versus source1.74129J/m². The analytic lattice
still has a nonzero infinite tail at40h; it is not an exactly separated limit.

The energy fit does NOT remove the normal-force/curvature failure. Free
saturated validation RMS is5.07158 for normal force and4.70525 for Haa, using
the original scales. At a=1.62h sourceHaa=-0.15583 whereas candidate=-1.51371
eV/L0². At the relaxed saddle the normal curvature is still about90–95% too
large. Thus the next material problem is specifically the normal response and
its nonlinear coupling, not arbitrary visibility of plastic strain.

Opening extrema were bracketed with129/257 logarithmic-near-h PLUS uniform
far-branch samples (258/514 total points). The baseline's large negative
traction minimum is -4.284GPa at3.151h. The two positive saturated candidates
show only one traction maximum on the checked [h,5h] branch and do not show
that negative-force pocket. Their peak is10.57/10.49GPa versus source12.97GPa.
These are ideal uniform-interface normal tractions, NOT experimental yield.
The actual Al99 rigid source itself has small later extrema, including a
-0.786GPa minimum around2.387h; those features are retained in source tables,
not censored to make a monotonic candidate look exact. One EAM source is not
experimental or DFT truth. Its force-shape discrepancy remains in the loss.
The high-alpha OLD-pair control still has a -3.807GPa pocket and is rejected.
No additional narrow near-h traction extremum was found for the reported
positive candidates by these refined brackets; the general nonuniform-limit
warning above is not a claim that such an extremum occurred in these data.

### 6.3 Independent numerical and loading checks

- Reference-state checks:243 independent wavevectors at radii12/16 plus
  continuous local searches are positive above the analytic tail envelopes
  for the baseline and both positive saturated candidates. This is not a
  whole-Brillouin-zone proof.
- Fixed-material dilation check:6 stretches x254wavevectors x2radii per model.
  All36 model/stretch/radius records are positive above the tail; worst margins
  across radii/stretches are0.0043054,0.0060222,0.0060245eV/L0², respectively.
- Direct/reciprocal even-environment energy comparison at radius18 has maximum
  absolute unit-amplitude error1.04e-17 for the first positive saturated candidate.
  Tightening reciprocal tolerance2e-11->2e-13 changes tested Hessians by at most
  3.61e-12eV/L0² across that candidate and the baseline.
- Three central-difference levels converge by the expected factor~4; final
  Hessian errors are1.45e-6 (baseline) and8.00e-8 (saturated)eV/L0². These are
  much smaller than the actual force/curvature material discrepancies.
- Nonzero-alpha finite-q harmonic response is independently checked from
  displaced per-atom nonlinear energies, with all six amplitude levels saved.
  Large nonlinear truncation errors and eventual roundoff are not deleted.
- Actual static tensor scenarios:297 states (source/baseline/saturated) plus
  198 states (source/last refinement), covering normal axes, shear components,
  mixed loading and static unloading. All states are stable roots; maximum force
  residuals are 5.51e-14 and 1.69e-14 eV/L0 in the two suites, respectively.
  Maximum static return error is 3.62e-15 L0. The source cases are repeated in
  both suites, so 495 evaluations are not 495 independent loading protocols.
  At4MPa interface shear, source displacement0.000330998Angstrom versus
  saturated0.000329450Angstrom (about0.47%). At50MPa normal traction it is
  0.000925366 versus0.000650806Angstrom (about-29.7%). This is NOT a cyclic PDE,
  dynamic hold, physical yield, or residual-plasticity result.

The independent300K source MD covariance comparison is also rerun with fixed
geometry/material units,576atoms per plane and12periodic classes. Harmonic
variance errors change from (+44.4%,-57.7%,-53.7%) to
(-14.4%,-20.9%,-13.3%) for normal/slip/transverse components at the source box.
The last positive-pair section gives(-14.3%,-26.9%,-19.9%). This is an additional
validation discrepancy, not a kinetic fit. Finite-temperature anharmonic
renormalization is absent, and these values are not included in the0K loss.

### 6.4 Verified QP repair, not coefficient clipping

`boundary_reproduction/actual_qp.npz` preserves a bit-for-bit reproduction of
the pre-fix issue: K3=-1.7148860e-9 despite whitened KKT3.55e-15. In the heavily
scaled basis a good whitened certificate did not ensure an exactly nonnegative
physical coefficient. The material constructor correctly refused that vector.

The fix resolves the constrained least-squares problem on the exact active
nonnegative face and then tests the ORIGINAL full-space KKT, equality and
inequality certificates. It does not independently clip K3 while retaining
the other coefficients. If the selected boundary is wrong the original KKT
test rejects it. Actual fixed-shape re-profiling gives K3=0 exactly; maximum
prediction change5.90e-13, exact residual2.04e-13, original KKT3.55e-15. It
does not fix the physical curve mismatch of that old-pair candidate.

## 7. Reproduction and release status

Run from the branch worktree, with a Python3 interpreter. `py -3` is not
installed in the validation environment; the actual tests used Python3.13.
Every result output path must be fresh. The inherited source target is
Al99 SHA25660c8a085be79d273324ab421f5b1447578fef55c1acfc6492c0999f15ee8a284.

```text
python -m solver_v1.run_tail_constrained_calibration --symmetry-channel --quartic-angular --exact-bulk --even-development --quadrupole-saturation --additional-wavepoints results/fcc111_active_interface/interface_development_v15/strain_counterexamples.json --wavepoint-step .1 --stability-stretches .99 1.0037037037037038 1.01 --max-nfev 24 --out <fresh-free-run>
python -m solver_v1.validate_tail_calibration <completed-run> --grid-step .08 --out <fresh-validation>
python -m solver_v1.run_even_environment_validation --models <completed-positive-runs> --out <fresh-checks>
python -m solver_v1.run_tensor_traction_audit --models <completed-positive-runs> --out <fresh-static-scenarios>
```

`definition.json` records each run's exact starts, bounds, optimizer, constraints,
source binding and fixed-pair control. For the last section add
`--fixed-pair-from .../saturated --starts .../positive_pair_refinement_starts.json
--optimizer feasible-simplex --max-nfev 48`. The old-pair control instead uses
v15/strain_stable_quartic. `--profile-only` re-solves coefficients without shape
optimization; it must not be advertised as an additional optimized fit.

Final verification: targeted40PASS39.02s; full solver518PASS585.55s;
app31PASS63.25s,0skips; actual desktop smoke exit0,1.01s. The earlier full
test process was intentionally stopped after adding the narrow-extremum
regression and is NOT counted as a pass. Only the later complete run counts.
The smoke retains two-row a0=0.7713438268704838, kappa=86.29296488740997.

**Decision:** a material-energy improvement is demonstrated within the
LJ/infinite-series framework. Bulk and finite-distortion energy agreement no
longer justify saying "nothing has been calibrated." But derivative, thermal,
identifiability and finite-defect/kinetic issues prevent full Al adoption.
Do not promote a new PDE/UI default. Physical Ma/Ms/t0/seconds/Hz, actual
specimen yield, finite-source/core activation and A_c remain uncalibrated.
