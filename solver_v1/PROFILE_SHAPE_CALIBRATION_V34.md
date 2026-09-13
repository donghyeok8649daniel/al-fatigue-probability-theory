# v34 — nonsmooth profile audit before further material extensions

Research only. The user-facing UI was launched from the verified branch;
its entered values and results are not modified by this study. The LJ pair,
Poisson/Bessel energies, target units/scales, production PDE and kinetic
calibration remain unchanged.

The previous outer shape optimizations stopped at evaluation budgets. A
certified coefficient LP at one shape does not certify shape stationarity.
Before adding more physics, audit the shape objective itself.

Let x be five logarithmic existing radial parameters and the existing signed
screening exponent. At fixed x, the current problem is

    v(x) = min_(c,eta) eta
    E(x)c=e,
    J(x)c-eta <= y,
    -J(x)c-eta <= -y,
    c_i>=0 for the declared subset; eta>=0.

E,J,e,y include the unchanged observable discrepancy scales. With SciPy's
RHS sensitivity duals lambda,mu_plus,mu_minus, the envelope expression is

    partial_k v = -lambda^T (partial_k E)c
                  -(mu_plus-mu_minus)^T (partial_k J)c.

The numerical column scaling inside the LP does not change this physical-c
formula: its derivative cancels by stationarity and complementarity. Tests
independently solve perturbed LPs with positive/negative residuals and a
parameter-dependent column gauge. Calibration matrix derivatives may use
deterministic finite differences of the analytic energy observables. This
does not replace the energy gradient/Hessian by finite differences.

At active-set changes, one returned dual is not necessarily a unique
derivative. The actual study compares two bounded derivative steps against
freshly reprofiled LPs, retaining asymmetric boundary stencils. A failed
gradient comparison is a nonsmooth/numerical diagnostic, not a material
impossibility certificate. No arbitrary drag, temperature, stress or target
tolerance adjustment follows from it.

If the profiled objective is nonsmooth, the equivalent explicit epigraph
formulation can retain x,c,eta jointly. Its constraints M(x)c are smooth
where the same analytic lattice series is resolved. This avoids treating a
selected LP dual as the gradient of a globally smooth function. Independent
reprofiling and exact-anchor checks remain necessary before accepting even
an optimizer result. Material and excluded-state validation remain separate.

Actual running/completed status belongs in CURRENT_WORK_HANDOFF.md and
results/profile_shape_v34; seven initial algebra/stencil tests passed before
the actual material sensitivity run. No calibration success is asserted here.
# Joint epigraph follow-up

## Completed numerical results (not accepted Al calibration)

| Study | SQP exit | seconds | certified final coefficient LP eta | excluded-q max relative error |
|---|---|---:|---:|---:|
| Original bounds | success,4 iterations |252.005|12.2042548163|203.7672%|
| Angular bounds24 | iteration limit,8 |652.576|11.5036122926|178.5307%|
| Feasible-LP restart, same bounds24 | iteration limit,12 |1691.183|11.1777577908|176.7758%|

Last shape=[4.1937068058,7.0313346624,19.5130752125,655178.38566,
6.4437319068,-1]. Final LP exact residual1.865e-13 and positive LJ. The exponent
is at its lower bound. The 12-iteration SQP is NOT converged, and neither the
parameter vector nor its modest improvement is a material calibration. No local
mobility, seconds/Hz, actual yield or fatigue validation follows from this.

At the first locally converged shape the active constraints include sign errors:
v19_power_new_2_Haa prediction -0.33337 vs reference +1.51240;
saddle_Hxx prediction +0.041093 vs reference -0.186425. Thus the discrepancy is
not remedied by one global energy/unit multiplier. This is evidence for ongoing
landscape-shape mismatch, NOT a proof that every allowed analytic family fails.

Actual full regression:779 solver tests/702.93s;50 app tests/62.27s. All source
EAM data remain target-only, LJ/Poisson/Bessel production and kinetic gates stay
unchanged. Completed excluded-q checks are development validation, not blind.

First actual run: SLSQP success in 4 iterations / 252.005 s, independent final
LP eta=12.204254816344415 versus initial12.507384978807881. Joint final exact
residual2.52e-9; the fresh coefficient LP restores/certifies the exact anchors.
All LJ coefficients remain positive. This is a local bounded optimization result,
not a global-family optimum. Excluded finite-q maximum relative Hessian error is
2.03767208965 (203.767%), so material adoption fails despite optimizer success.
These are previously inspected excluded points, not new blind validation.

The final rank2 radial decay reaches12, screening saturation reaches1e6, and
existing exponent reaches1. To distinguish a radial search boundary from a
family limitation, joint_angular_domain24 continues from this completed point
with ONLY the three angular decay upper bounds expanded12->24. Scalar decay,
saturation, exponent bounds, signs, exact targets and target tolerances remain
unchanged. This is explicitly a search-domain sensitivity, not physical tuning
or an accepted parameter set. Its actual completion is recorded separately.

The completed initial audit used 24 profiles / 607.087 seconds. Its derivative
gate failed: dual-envelope versus reprofiled finite differences disagreed by
9.9928 and 9.2599, despite envelope step refinement changing only 8.99e-7.
No dual-guided profile optimization was executed. This does not establish a
global optimum, nor failure of the entire nonlinear family.

The follow-up keeps exactly the same observable matrix M(x), targets y, scales
s, coefficient signs, and shape bounds. Instead of differentiating min_c, solve

    min_(x,u,eta) eta,     c = D u
    r_exact = 0
    eta - r_test >= 0,    eta + r_test >= 0
    r = (M(x) D u - y)/s.

D is fixed positive numerical scaling computed at the starting certified LP.
The Jacobian is dr/dx_k=(dM/dx_k)c/s, dr/du=M D/s, dr/deta=0.
The two inequality Jacobians have opposite residual signs and +1 in eta.
The energy jets remain analytic; only the observable-matrix shape dependence
is differenced with the existing bounded second-order stencil.

Synthetic full-Jacobian and nonsmooth-profile cusp tests pass. The actual
study writes every iterate's exact residual and inequality violation; an SQP
trial iterate can be infeasible and must not be adopted. The final proposed
shape is reprofiled with the independently certified coefficient LP. A finite
iteration budget is not convergence and a successful optimizer is not material
validation. Excluded q points remain outside the loss, and fresh validation is
still required before adoption. Actual results are under joint_epigraph.
