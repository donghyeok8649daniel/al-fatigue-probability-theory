# Tasks 1–3: completion gates, measured results, and remaining obstructions

This is a research checkpoint, **not completion of the three physical goals**.
The study does not change production LJ/Bessel, material parameters, mobilities,
probability bookkeeping, specimen aggregation, or the UI energy selector.

## 1. One credible Al bulk/interface energy

Bulk equilibrium, cohesion and three elastic constants have been fitted in
earlier work. Calling everything uncalibrated would conceal that achievement.
But those same candidates do not reproduce the independent interface and
finite-wavelength restoring response. A correct modulus at q=0 is not a correct
core or interface. A production-ready common material law is still missing.

This v32 investigation actually tested:

- 61 coefficient profiles across six radial/environment-channel studies;
- an independent 180-profile search in the previous radial family;
- a positive polynomial-times-exponential kernel, with an exact differentiated
  Poisson transform and analytic gradient/Hessian;
- a 120-profile material search in its three scalar-shape coordinates;
- a completed 320-evaluation joint scalar/angular shape search (4714.573 s,
  evaluation-budget termination, not optimizer convergence);
- four embedding/angle-sign controls and three auxiliary pair/gauge controls;
- the separately derived true spatial density gradient and its one-amplitude
  material compatibility test;
- the v33 STF rank-four environment with its nonzero cubic bulk background,
  five fixed-decay studies and independently evaluated new interface states.

The best large angular dictionary gives eta=11.24805 for a necessary box
criterion eta<=1. It is a capacity diagnostic, not an adopted many-parameter
model. The true gradient extension gives eta=12.50710, versus12.50738 before;
excluded-q error remains210.08%. Analytic/direct derivative checks pass, so
the mismatch cannot be called a missing finite-difference or unit factor.
The new scalar-only material trial also fails, and some smaller objective
values occur at zero LJ attraction: these are explicitly rejected.

The completed eight-shape search gives best eligible eta=15.27014, with
100.563% maximum fresh-interface Hessian error and249.667% excluded-q error.
The independently derived rank-four term retains the required background
contribution H4=2LL^T+4Q_bulk:K. Its best training-selected signed control has
eta=12.45013; fresh-interface and excluded-q errors are99.202% and208.557%.
These are actual material discrepancies, while analytic derivative errors
are below5e-10 in the tested reduced-coordinate jets. Thus these trials do
not repair the common-material gate, even though their pristine forces and
tangents are satisfied. A negative rank-four coefficient was a diagnostic,
not an adopted stable material extension or an all-wavevector stability proof.

The final joint five-rank4-range capacity control actually completed two LPs
after exactly replaying all15 single-range predictions. Its signed eta is
11.55613, still above1, and it sets the LJ repulsive coefficient to0. The
smaller31.731% excluded-q error of that ineligible solution is not evidence
of a successful retained-LJ calibration. With all five new coefficients
nonnegative, eta remains12.50730. This closes that particular fixed-shape
multi-range control, not the question of all analytic environment families.

Finally,13 actual profiles restore the pure-exponential endpoint omitted by
the earlier finite envelope-width domain. Mixing the two normalized densities
at the same decay is still the same quadratic family for nonzero mixture;
the zero-mixture matrix is exactly the baseline. No sampled mixture improves
eta12.50738 (density identity error2e-15). This removes that particular
comparison loophole without claiming a global nonlinear-family no-go result.

The exact convex-profile certificates apply at their fixed radial shapes.
Finite search budgets and several unsuccessful analytic extensions do **not**
prove that every LJ-plus-environment family is impossible. Conversely, taking
a candidate merely because it is the best found does not make it Al-calibrated.
No additional empirical yield or desired-event-ordering term has been fitted.

The next material hypothesis must explain the simultaneous interface-curvature
and finite-q conflict. Repeating mobility fits, increasing temperature, loosening
the material discrepancy box until it passes, or changing the plot cannot do so.
Before a new term is adopted it needs an independent mathematical definition,
the same per-site density counting, analytic infinite sums, nested limits,
identifiability, and fresh material validation.

## 2. Local mobility, seconds and Hz

There are real kinetic measurements in this repository. The v31 study ran
10 weak-forcing/unforced Al99 reference records totaling20 ns, after separate
initial equilibration. These are reference MD, not the production LJ potential.
The point estimates for the plane-coordinate finite-band mobility are
approximately1.36214e11 (normal) and2.06114e11 (slip) m^2/(J s). They are not
the production local-coordinate M_a and M_s.

Two independent checks prevent physical-clock promotion:

1. The observed inverse-response error sets include zero drag, so the mobility
   upper bound is unbounded. The error sets are not statistical confidence
   intervals. Replacing a signed noisy drag by its absolute value would invent
   a calibration.
2. A plane-average PMF is not the local interface energy used by the PDE.
   Even exact per-plane-atom normalization changes the thermal coefficient:

       F_bar=F_plane/Np,
       M_bar=Np M_plane,
       theta_bar=kBT/Np.

   This preserves the same plane-coordinate generator, not a new local law
   at the original numerical kBT. Eliminating other coordinates can also
   produce memory; marginal and projected mobility are generally different.

The dimensional identity is complete:

    M_i* = t0 E0 M_i,phys/L0^2,
    t0 = M_i* L0^2/(M_i,phys E0),
    t_seconds=t0 t_model,       f_Hz=f_model/t0.

This identity determines no numerical t0 without a physically validated M.
The remaining experiment is a matching-coordinate, matching-PMF kinetic
measurement with a resolved dissipative response and a demonstrated Markov
window (or an explicitly derived reduction). A generic bulk self-diffusivity,
atomic-mass time, or dislocation-line mobility does not supply that datum.
Current production M_a,phys/M_s,phys/t0 remain unavailable and Hz disabled.

## 3. Actual yield and fatigue, rather than ideal traction

There are retrieved condition-specific experimental benchmarks, not an absence
of strength data. For example, Pigato et al. report different room-temperature
YS for three ultra-pure Al microstructures, and Krebs et al. measure a wire
CRSS at0.002 plastic shear. These are not interchangeable definitions or
identical initial specimens. The v21/v11 datasets preserve those distinctions.

[Pigato et al. (2026)](https://doi.org/10.3390/ma19061195) reports tensile YS
16.4,74.6,34.1 MPa for its recrystallized6N, fragmented5N5, recovered5N states
respectively at25C. The paper's retrieved text still does not specify the
numerical proof offset; it is not silently relabeled Rp0.2. The actual source
populations and pin coordinates are not supplied. The publisher text was
rechecked during v32; no missing value was guessed from a generic Al table.

[Krebs et al. (2017)](https://doi.org/10.1038/nmat4911) concerns particular
cast single-crystal wires, not an arbitrary pure-Al specimen. The retained
CRSS/G points cannot be converted with an unverified normalization modulus.
Single-arm assumptions in its supplement are not measured pin lengths and
are not the same as our double-ended source reference.

The exact geometry already implemented is

    epsilon_slip = sym[(1/V) sum_k A_k b_k tensor n_k],
    V sigma:d epsilon_slip = sum_k (b_k.sigma.n_k) dA_k.

It exposes an identifiability obstruction: two specimens with the same atomic
energy and kinetic law but different source populations/pins can have different
macroscopic strain histories and yield. Static atomic targets alone cannot
select the unknown initial source state. An assumed micrometre pin length that
produces MPa bow-out is therefore not validation of a specimen's proof stress.

Required remaining work is a validated common material/core, finite activation
geometry, independently specified initial source/microstructure state, kinetic
evolution and unload/persistence, followed by a like-for-like experimental
comparison. Source nucleation, slip transfer, line interactions, irreversible
strain and crack first passage must not be collapsed into one damage number.
Fatigue-life data must not be used to tune a mobility or energy barrier.

## What is preserved

No Monte Carlo trajectory-counting solver, empirical yield clipping, fatigue
life fit, characteristic-volume scaling or hidden force rescaling was added.
The statistical A_c remains external probability aggregation only. The new
gradient coefficient is not a characteristic area. Research calculations and
negative findings are saved; production remains behind the validation gates.

## Where to continue without repeating failed work

- `RADIAL_CHANNEL_COMPATIBILITY_V32.md`: six completed channel controls.
- `QUADRATIC_EXPONENTIAL_ENVIRONMENT_V32.md`: density-shape derivation and
  completed material/gauge/all-ranges diagnostics and independent validation.
- `SPATIAL_DENSITY_GRADIENT_V32.md`: true spatial gradient, exact derivatives,
  Bloch consequence, test results and failed fixed-shape material check.
- `RANK_FOUR_ENVIRONMENT_V33.md`: nonzero cubic STF background, full Bloch
  Hessian, five completed decay controls and failed independent validation.
- `WEAK_KINETICS_AND_MATERIAL_V31.md`: actual kinetic records and uncertainty.
- `SPECIMEN_YIELD_AND_SLIP_BUDGET_V21.md`: source-to-specimen mechanics.
- `CURRENT_WORK_HANDOFF.md`: actual active sessions, final tests and Git status.

An implementation/test completion must never be reported as completion of
physical goals1–3 while any of these gates remains unmet.
