# Validation and research roadmap after PART I–III

20 October 2026. This is an ordered set of falsifiable gates, not a statement
that atomistic fatigue has been validated. Full directive PART IV–VII remains
deferred. No new material constants are fabricated to close these gates.

## Completed in this delivery

The isolated `research/first_principles_20261020/verify_foundations.py` ran
to exit code 0, with **8 passed, 0 failed**, and saved the full finite-valued
JSON in `results/run_001/verification.json`. Python 3.12.14, NumPy 2.5.3 and
SciPy 1.18.1 were used. It completed in a reaped owned process; the private
execution receipt records start/end, PID and timeout status. No production
file was edited and no full production solver/application suite was run.

| Check | Evidence | Limit |
|---|---|---|
| Row Poisson identity | Independent direct value/gradient/Hessian, mode/tail refinements; scaled error $<9.84\,10^{-15}$ | Relative/scaled channel error, not an absolute error claim for near-singular huge energies |
| 2D plane transform | Direct triangular lattice versus reciprocal $K_{p-1}$, multiple real/reciprocal cutoffs | Energy analytic tail bound only; not a claimed Hessian tail certificate |
| Actual cross-cut counting | Enumerated finite $4\times5$ row pairs versus $c_k$ Bessel sum; error $<9.44\,10^{-16}$ | Explicit zigzag geometry, not every orientation |
| Nonlinear cut derivatives and endpoint | Stress-free stack, positive local Hessian, finite differences and approach to separation limit | Dimensionless illustrative stack; refinement equality at machine precision is not an exact omitted-tail bound |
| PMF counterexample | Exact transverse Gaussian integral versus quadrature with same full Bessel potential | Hypothetical bath, not a measured Al entropy correction |
| Hamiltonian bath / exact GLE | Nonlinear force identity, external-work balance, analytic FDT, recurrence and initial slip | Finite selected bath; no claim of Brownian Al forcing |
| Two distinct stochastic limits | Local Bessel-Hessian OU covariance and frequency-response refinements | Explicit local quadratic control, not nonlinear fatigue trajectories |
| Multiplicative-noise convention | Correct Gibbs flux versus omission of Itô drift | Fixed-opening periodic registry quotient only |

## Gate 1: fix the microscopic contract

Declare the potential and its validation range, finite patch and boundary
conditions, crystal plane/path, atom groups, $C$ or nonlinear $\Xi$,
area/length/energy/mass normalization, temperature/ensemble, imposed normal
or projected loading, and preparation. Keep the row reference, full plane
and counted cut distinct. Reject a comparison if the MD observable cannot
be pushed forward to the paper's declared coordinate and energy measure.

## Gate 2: test the thermodynamic bridge at the same coordinate

Estimate a constrained PMF or thermodynamic-integration mean force with
convergence in sampling/window overlap, system size, boundary distance and
temperature. Check the coarea/metric convention explicitly. Compare
$\nabla A(q)$ and free-energy differences with the rigid $U_\Gamma$ and
relaxed zero-temperature curve over the actually occupied and escape region.
Measure uncertainty and check overlap; agreement only near one minimum does
not validate an opening barrier. If hidden preparations give different mean
forces, record history variables or a nonequilibrium reference rather than
calling a single equilibrium PMF complete.

**Pass condition:** the equality $A=U_\Gamma+C$, or a quantitatively bounded
correction, is supported on the declared operating region. Otherwise use the
measured/derived $A$ and retain the Bessel energy as the rigid baseline.

## Gate 3: identify and test reduced dynamics

Use same-coordinate atomic trajectories to estimate conditional acceleration,
velocity statistics and memory/response across states and preparations.
Compare held-out multi-time correlations, responses and first-passage laws,
not only equilibrium histograms. Test off-diagonal friction rather than
imposing it to vanish everywhere. Separate bath decorrelation, momentum
relaxation, resolved relaxation and loading periods with uncertainty bounds.
Test for recurrent/slow tails and non-Gaussian or nonstationary residuals.
Reconstruct the residual force with a causal kernel or the appropriate
nonlinear memory, and check its preparation dependence. A finite fitting
interval is not a certificate of the zero-frequency friction integral.

**Decision:** retain the GLE/inertial law when memory or inertia is resolved;
adopt Markov and then overdamped closure only if their errors are controlled
on the target loads and event-time range. Existing v31 Al measurements are
the starting evidence, with their zero-drag and coordinate-matching limits
carried forward; do not silently reuse their central mobility as calibrated.

## Gate 4: bind physical time and boundaries

Establish a finite compatible mobility/friction interval and thermal
normalization before reporting seconds or cycles. Resolve continuous domain,
reflection/absorption and initial-layer contracts. Check system-size and
patch-coherence scaling; dividing a plane energy by $N_p$ requires both
$M\mapsto N_pM$ and $k_BT\mapsto k_BT/N_p$. Validate inertial/colored-noise
boundary layers if they affect the event. A numerical no-flux outer box is
not evidence of physical confinement.

## Gates beyond this first delivery

Only after the foundations are tested should one assess irreversible defect
evolution, spatial interactions/coalescence and the map from local opening
to observed crack initiation. Use specimen geometry, signed load histories,
synchronized negative/positive inspections and detection metadata. Evaluate
what a fixed two-coordinate state can predict under repeated loading before
adding any history variable. Do not introduce an arbitrary damage variable,
barrier drift or fitted clock merely to produce a fatigue curve.

The immediate review decision is concrete: preserve Rev.11 as a conditional
reference and publish this reconstruction alongside it. The next scientific
work is same-coordinate PMF plus dynamical validation, not algebraic
promotion of $W_\Gamma$ into a material probability current.
