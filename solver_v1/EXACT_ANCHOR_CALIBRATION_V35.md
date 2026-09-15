# v35: exact-anchor coordinates for material calibration

## Scope

This study improves optimization of the existing coefficient-linear
LJ/Poisson/Bessel/environment family. It retains its six shape coordinates,
ten energy coefficients, seven exact anchors, coefficient sign conditions,
source targets, physical units and discrepancy scales. No new energy term,
production parameter, mobility, temperature or specimen area is fitted.

The146 observation rows contain seven exact anchors,115 tested interface
Hessian components and24 tested finite-q components at four wavevectors.
The inherited LP flag `spectral_constraints_included=False` refers to omitted
additional positive-semidefinite stability constraints; it does not mean the
finite-q discrepancy rows were omitted. No full-Brillouin-zone stability
certificate is claimed.

The starting positive-LJ v34 result has eta=11.17775779075355, versus the
necessary target-box criterion eta<=1. Its 12-iteration joint optimization
stopped at its iteration limit. During intermediate steps its normalized
exact-anchor residual reached about 4.1; only the final separate coefficient
LP restored those anchors. Eliminating the exact equalities should prevent
that drift, but is not by itself evidence for an adequate material family.

## Exact elimination and derivative

At shape x, write the scaled observation matrix as A(x)=M(x)D/s and the
scaled targets as b=y/s. D is a fixed positive numerical column scale.
Partition the seven exact rows into dependent and free coefficient columns:

    E = A_exact = [E_B, E_F]
    c = D*z
    z_F = u
    z_B = E_B^(-1) * (b_exact-E_F*u).

Seven pivot columns B are chosen by QR with column pivoting at the starting
matrix, then held fixed. The remaining three coefficients u, the six shape
coordinates x and the epigraph variable eta are optimized. These are ten
numerical unknowns; no material degree of freedom is added or removed.

Matrix inversion notation above denotes a linear solve in the implementation.
A pivot whose condition estimate times `7*machine_epsilon` exceeds 1e-6 is
rejected. This is a numerical arithmetic budget, not a physical tolerance or
a statistical uncertainty. Rank loss is not repaired by a pseudoinverse.

Define N_F=I and N_B=-E_B^(-1)*E_F. At fixed u,

    partial_k z_F = 0
    partial_k z_B = -E_B^(-1) * (partial_k E)*z
    r = A*z-b
    partial_k r = (partial_k A)*z + A*(partial_k z)
    partial_u r = A*N.

The minimax constraints are

    eta-r_test >= 0,  eta+r_test >= 0,
    z_i >= 0 for the original nonnegative coefficient subset.

The latter constraints include dependent coefficients; retaining signs only
on the free coefficients would change the energy admissibility conditions.
The optimizer uses the full implicit Jacobian above. It does not differentiate
the nonsmooth coefficient-LP optimum using one arbitrarily selected dual.

Equality elimination is a standard optimization coordinate construction; see
Boyd and Vandenberghe, *Convex Optimization*, section 10.1.2,
[book](https://web.stanford.edu/~boyd/cvxbook/bv_cvxbook.pdf).
The shape-dependent matrix and nonlinear feasible domain here are **not**
asserted to be convex. The implementation uses [SciPy SLSQP](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-slsqp.html)
with supplied constraint Jacobians and reports its actual termination status.

## Actual energy derivatives and tests

Energy gradients and Hessians remain the established analytic lattice jets.
Only differentiation with respect to calibration shape parameters uses the
existing bounded second-order difference stencil at log-coordinate step1e-4
in the primary runs, and5e-5 in the declared finer continuation.
The complete implicit constraint Jacobian is also checked at the actual Al
starting point against two fresh one-sided directional differences, including
the inherited exponent boundary. In the local run, relative maximum errors
are 4.298e-7 and 4.056e-7; the declared numerical check tolerance is 2e-5.
The starting anchor pivot condition is about159.35 and the exact residual is
2.84e-14.

Synthetic tests check a shape-dependent anchor solve against finite differences,
physical-coefficient reconstruction under numerical rescaling, singular-pivot
rejection, and agreement with an independently posed full coefficient LP.
A separate feasible-cone LP detects possible first-order descent at a final
point and respects active shape boundaries and nonsmooth epigraph cusps.
It is a necessary local diagnostic with finite numerical derivatives, not a
proof of a global optimum or material validation.
The audit compares the executed box with infinitesimal directions inside
the original shape domain. The isolated-neighbour decay constraint remains
included. It does not evaluate inadmissible finite steps or certify a global
optimum merely because a local box prevents a descent direction.

## Search protocol

The two joint runs and the separate scan start from the same completed v34
best profile.

1. `local`: each of the first five log-shape coordinates is restricted to
   +/-0.35 about the start, intersected with the previously declared bounds.
   The existing screening exponent keeps [-1,1]. The budget is30 iterations.
2. `full_saturation`: retain the same local decay ranges and exponent bounds,
   but give the **existing** saturation parameter its original entire positive
   range [1e-5,1e6]. The budget is40 iterations. This control checks sensitivity
   to the narrower local saturation bounds.
3. `saturation_scan`: hold the other five shapes fixed and independently
   profile coefficients at the initial alpha and at alpha values
   1e-5, 1e-2, 1, 1e2, 1e4, 1e5 and1e6. These points test distant parts of
   the same positive range without relying on a local shape derivative.
   This is a deterministic scan, not an optimizer or an exhaustive search.
4. `refined_continuation`: after the full-saturation run terminates, start
   from its best independently certified positive-LJ profile. Reset only the
   numerical coefficient scales at that point, retain the original full alpha
   range, and halve the matrix shape-difference step to5e-5. The budget is12
   iterations. This follows the active-cone tolerance sensitivity described
   below; it does not use excluded-state results.

The local boxes preserve the original isolated-neighbour decay condition
throughout: `2*k1 + min(p,0)*k_scalar > 0`. The saturation parameter does not
enter that condition. No target box or material discrepancy is loosened.
The actual run bounds and fixed numerical scales are in
`anchored_definition.json`; the inherited `definition.json` records the
original ShapeProfile data-loading context.

At every optimization callback, a separate fixed-shape coefficient LP is
solved and its primal/dual certificate retained. The final shape is selected
only by the smallest training eta among completed positive-LJ profiles,
including the initial profile. A worse or ineligible optimizer endpoint is
not substituted for the best eligible training point. All callback iterates,
including any inequality violation, remain recorded.

An optimizer success flag, exact anchors or an improved training eta does not
certify the remaining interface surface, spectral stability, material response,
physical clock, actual yield or fatigue. A boundary optimum in a local box is
not a global-family optimum.

## Independent validation plan

Before each fit begins, `validation_plan.json` declares ten interface states
and six finite-q points, excluded from that fit and its selection. The validator
compares the baseline and the training-selected candidate on the same source
values, and rejects any interface-state overlap with the training states.
Full-history blindness is not claimed; the older four excluded finite-q points
remain explicitly labelled as previously inspected development validation.

Validation includes full 3x3 interface Hessians, energy and force differences,
analytic tolerance refinement from2e-11 to2e-13, two-step independent
energy/gradient differentiation at two off-axis states, bulk radius12/16
comparison and analytic finite-q tail bounds. Positive tested eigenvalues
would cover only those wavevectors, not the entire Brillouin zone.
The reference is source Al99 at the matched static geometry, not experimental
yield data or the finite-temperature plane PMF from the kinetic campaign.

## Reproduction and actual completion

    python -m solver_v1.run_anchored_calibration \
      --parent results/rank_four_environment_v33/fixed_shape \
      --resume results/profile_shape_v34/joint_angular_domain24_restored/joint_summary.json \
      --out <fresh fit directory> --maxiter 30

For the saturation control add `--saturation-global --maxiter 40` with a fresh
output directory. The rank-four parent above supplies the unchanged ten-column
**baseline** profile and target provenance; a rank-four energy term is not used.
For the separate fixed-other-shapes scan use `--scan-saturation` instead.

`report_anchored_calibration` selects the smallest certified positive-LJ
training eta across completed runs, including the baseline, **before** the
excluded validator executes. It stores all training residuals and source
report hashes. Its diagonal-curvature sign diagnostic counts only opposite
signs for which both magnitudes exceed the already declared discrepancy
scale; that scale is not reinterpreted as statistical uncertainty.

    python -m solver_v1.validate_anchored_calibration \
      --candidate <completed fit directory> \
      --baseline results/profile_shape_v34/joint_angular_domain24_restored/joint_summary.json \
      --out <fresh validation directory>
    python -m solver_v1.audit_anchored_stationarity \
      --study <completed fit directory> --out <fresh local audit directory>
    python -m solver_v1.report_anchored_calibration \
      --fits <completed fit directories> \
      --baseline results/profile_shape_v34/joint_angular_domain24_restored/joint_summary.json \
      --out <fresh training selection directory>
    python -m solver_v1.plot_anchored_calibration \
      --report <training selection directory> --validation <excluded validation directory>
    python -m pytest solver_v1/test_anchored_shape_calibration.py -q

Actual optimization, validation and regression completion belongs in
`CURRENT_WORK_HANDOFF.md` and each saved summary. Running calculations are not
reported as completed calibrations.

### Completed local control

The `local` run actually terminated successfully after23 optimizer iterations
and3376.664 s. Its independently restored coefficient LP gives
eta=11.077728587023318, versus11.17775779075355 at the start (0.8949% reduction).
The maximum callback exact-anchor residual is4.263e-14. The selected shape is
the optimizer endpoint; the coefficients are separately LP-certified there.
The saturation parameter reaches its local upper limit929742.386 and the
screening exponent remains at-1. This is still far outside eta<=1.

The active-cone diagnostic gives minimum d_eta about3.88e-10 in the executed
box, and-2.50e-6 when only the local box is removed. Both active tolerances
(2e-6 and1e-5) give the same result. The latter very small numerical direction
is not promoted to resolved physical improvement or global shape stationarity;
the shape derivatives remain finite differences. The completed full-saturation
control and the new excluded-state validation are recorded separately below.

### Completed full-saturation control

The second run actually completed40 iterations in5072.299 s and stopped at
its iteration limit (`optimizer_success=False`). Its best eligible training
eta is10.99216621215597, selected from an earlier callback, rather than the
optimizer endpoint. The maximum callback exact residual is5.684e-14.
The selected even decay is23.99999923, at the declared24 upper bound, and
p=-0.99999999997 remains at its lower bound. Alpha is1662.17327. Scanning
alpha alone at the starting other shapes did not find this improvement:
its best value was still the starting11.17775779.

The final active-cone minimum d_eta changes from-0.12347 to-0.00014003 when
the active tolerance changes from2e-6 to1e-5. Removing only the local box
gives the same results because the relevant shape limits are original ones.
Near-active discrepancy rows therefore matter. This is not a robust local
stationarity certificate, and motivates the declared shorter continuation
with fresh numerical coefficient scaling and a finer shape derivative.

    python -m solver_v1.run_anchored_calibration \
      --parent results/rank_four_environment_v33/fixed_shape \
      --resume results/anchored_calibration_v35/full_saturation/joint_summary.json \
      --out <fresh continuation directory> --maxiter 12 \
      --saturation-global --shape-step 5e-5

The continuation's actual starting Jacobian checks have relative errors
8.406e-8 and8.868e-7, with absolute errors8.292e-6 and8.747e-5. Both pass
the declared numerical check, but the smaller independent check step has
the larger error: monotone derivative convergence is **not** established.
The new pivot condition is119.07. Neither this condition estimate nor the
finer nominal step is a material-identifiability certificate.

### Completed continuation and excluded validation

The continuation terminated successfully after6 iterations and1008.120 s.
Its selected shape is its optimizer endpoint, with an independent coefficient
LP eta=10.99216570097229. This is a1.66037% training reduction from v34.
The additional5.112e-7 change from the full-saturation control is below the
LP's2.198e-6 certificate tolerance; that final tiny change is not resolved
material improvement. The saved raw minimum determines the selection under
the previously declared rule. Neither excluded interface nor finite-q errors
were available to this selection. The physical shape and all ten coefficients
are stored in `refined_continuation/joint_summary.json`.

The maximum normalized exact-anchor residual over continuation callbacks is
3.197e-14; the final chart value is1.776e-14. The separately solved coefficient
LP's unscaled exact-row residual is7.105e-14. Final active-cone minima are
-1.566e-7 in both domains at both declared active tolerances. This small
numerical value, finite shape derivatives and active original k2=24/p=-1 bounds
do not certify a global optimum. All intermediate trial violations remain in
the checkpoint records; exact anchors do not imply every trial satisfies the
other inequalities.

The selected candidate and v34 were independently evaluated at the ten
predeclared new interface states and six new wavevectors, with no refitting:

| Metric | v34 baseline | v35 selected |
|---|---:|---:|
| Training eta (necessary limit1) | 11.177758 | 10.992166 |
| Maximum new-interface full-Hessian relative error | 118.2114% | 118.0056% |
| Maximum previously inspected excluded-q error | 176.7758% | 166.6661% |
| Maximum new excluded-q error | 170.1685% | 165.1135% |

The norms are Frobenius matrix norms. A smaller maximum is not uniform
improvement: for example, interface state6 worsens from about29% to55%, and
state9 from about25% to45%. Four training diagonal curvatures still have
opposite signs with both magnitudes exceeding the original discrepancy scale.
All individual targets, residuals and full excluded matrices are retained.
The common material gate therefore **fails** despite the improved objective.

For the selected candidate, the maximum Hessian change when tightening analytic
tolerance is1.785e-13. Independent two-step energy/gradient differences give
maximum gradient error4.991e-10 and Hessian error1.330e-9. Bulk radius12-to16
changes at most2.135e-4 in Frobenius norm (1.466e-4 in operator norm); the
radius16 tail bound is at most5.147e-4 in operator norm. These numerical errors
are separate from the large source discrepancies. The minimum tested bulk
eigenvalue is1.4920; positivity at these14 wavevectors is not all-q stability.

The validator reports interface energy per primitive interface cell and bulk
stiffness per atom, using fixed L0=2.863782463805517e-10 m. Source parameter
provenance is copied without the unrelated straight-row energy/tensor labels
from the shared source loader. Source Al99 is a static reference here. No new
MD, local PMF, mobility, production seconds/Hz, yield or fatigue-life calibration
was performed. The comparison PNG/SVG contains both the enlarged progress axis
and the full eta scale with its required limit.

Final actual test results and the separately requested STL UI fix are recorded
in `results/anchored_calibration_v35/COMPLETED_SUMMARY.md` and the handoff.
