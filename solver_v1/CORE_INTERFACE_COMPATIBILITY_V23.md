# v23 — source-core and interface-jet compatibility

Continuation of v22, starting clean local/fresh remote
`d1ac6695518d3738b8303687d8c6a7b58877203e`, probability-pde-solver-v1.
This is static calibration research, not a new production energy or fatigue law.
No LJ/Bessel kernel, default parameters, probability PDE, mobility or A_c changes.

## The question before further fitting

The v22 wider candidate reduced frozen-source atomic force RMS by47.23% but
still predicted the wrong sign of saddle Hxx and of a previously inspected
opening Haa. Does the old objective merely prefer this tradeoff, or do the
exact initial tangents and declared coefficient signs exclude the source jets?
Test that question separately from a local nonlinear optimizer's termination.

At a FIXED existing shape theta the whole current law is coefficient-linear:

    y(theta,c)=M(theta)c,  E(theta)c=e,
    c=(u,v,A,B,C,D3,D1,D2,D_Eg,K3).

The five bulk and two pristine-interface equalities are retained in the main
study. A separately labeled five-equality control removes the two tangent
equalities only to diagnose their role; it is not adopted as a new calibration.
Coefficient signs are those of the existing family. Scalar shapes, per-atom
environment summation, source normalization and all derivative formulas stay
unchanged. A fixed-shape LP cannot prove impossibility of the entire family.

For a selected set of observed jets J with source values y and previously
declared positive discrepancy scales s, solve

    eta_min=min_(c,eta>=0) eta
    subject to E c=e, c_j>=0 (declared indices),
               -eta s_i <= J_i c-y_i <= eta s_i.

This is a linear program, not an empirical strength cutoff. eta_min>1 means
that no coefficient vector under THOSE constraints matches every selected jet
within its declared discrepancy scale. Removing sign or spectral restrictions
enlarges the feasible set; its optimum is a lower bound for the restricted
problem, not an admissible material. A failed numerical LP is not a proof.

The implementation independently checks primal feasibility, dual signs,
stationarity, complementarity and the primal/dual objective gap. Column/row
scaling changes numerical units only. Save the multipliers so a claimed
incompatibility can be traced to a weighted combination of exact anchors,
target jets and coefficient signs. No clipping of coefficients or forces.
These are numerical certificates for the recorded coefficient matrix with
reported tolerances, not interval-arithmetic proofs covering every possible
amplification of infinite-series truncation in unbounded coefficient limits.

## Testing a feasible target box

Where intervals are feasible, minimize the SAME v22 two-block objective while
enforcing the specified source jets inside their original discrepancy scales:

    L=L_interface/L_interface,baseline + L_core/L_core,baseline.

These weights are a transparent diagnostic tradeoff, not reported experimental
uncertainties. Previously inspected excluded rows used for a new constraint
are development data for this step and cannot be called independent validation.
Unselected states remain excluded; separately declare new states before fitting.

An affine bound l<=Jc<=u is passed to the existing verified spectral/QP helper
using c_const=1 as an exact bookkeeping coordinate:

    [J,-l](c,1)>=0, [-J,u](c,1)>=0.

The appended coordinate is not a material coefficient or energy term. Its
bulk spectral operator and truncation tail are identically zero. The original
positive-LJ requirements, spectral/tail checks and returned-vector KKT tests
remain. A good frozen-force fit still requires actual atomic re-relaxation,
Morse, neighborhood/domain and signed static load-return tests before any
conclusion about the core. Passing these would still not establish physical
yield, finite-loop activation or production a/s seconds/Hz.

## Results

### Completed fixed-shape audit

48 actual LPs completed with independent primal/dual certificates. For each
of the original and v22 wider shapes, both named offending curvatures can be
matched together (eta about2e-14) with the seven anchors and original signs.
Thus their two sign errors alone do NOT prove algebraic impossibility.

Larger previously inspected sets expose a different result:

| Shape / constraints |27 fit curvatures eta_min|115 inspected curvatures eta_min|102 old fit observations eta_min|
|---|---:|---:|---:|
| Original, exact7 / declared signs |9.245760|26.309501|16.367467|
| Wider, exact7 / declared signs |7.971412|18.464447|15.570254|
| Wider, exact7 / all signs relaxed |7.950753|18.338403|13.728641|
| Wider, bulk5 only / declared signs |7.375726|13.940478|13.497389|

These are MINIMUM worst normalized errors over all coefficient vectors in
the stated LP, not simply the residuals of the previous least-squares fit.
The declared scales are unchanged and are not experimental uncertainties.
Spectral conditions are omitted from this necessary-condition audit; adding
them cannot improve the bound. None of these fixed-shape results excludes
other radial shapes or proves whole-family failure.

In the wider115-curvature LP, the active source-jet multipliers are
saddle_Hxx .03635737, v19_new_4_Haa .69798448 and
v19_power_new_2_Haa .26565816 (signed upper/lower forms saved separately).
Their sum is1. The sole active coefficient lower-bound multiplier is K3's;
the exact-anchor multipliers and scaled column map are saved too. The latter
two opening states have a*=1.39212667048 and .953668006524, respectively.
Allowing negative K3 does NOT solve the whole set: the all-signs-relaxed
lower bound is still18.3384, with other active rows. A dual support is a local
compatibility witness, not a rule to negate an environmental term by hand.

### Completed same-shape joint fits with source boxes

The original v22 objective was independently replayed before the new boxes.
One original discrepancy scale is enforced around each selected source value.
The formerly inspected v19 opening state is now explicitly a development
constraint, not an independent validation observation.

| Protocol |Other-interface squared loss|Core-force RMS eV/L0|saddle Hxx|opening Haa|
|---|---:|---:|---:|---:|
| Unchanged replay |2276.471|.152417|+.182166|-1.234287|
| Saddle box only |4376.632|.201940|-.167782|-.772407|
| Opening box only |39502.901|.219379|-.509881|+.803279|
| Both source boxes |41943.919|.175221|-.205067|+.803279|

Curvatures are eV/L0^2. The saddle-only optimum has a zero LJ coefficient
and is an INELIGIBLE closure, not a new LJ potential. The other two constrained
fits retain positive LJ and pass the sampled bulk-tail conditions, but their
large other-interface failures forbid adoption. Both-box exact residual is
5.33e-14, KKT1.33e-15 and sampled robust margin.005045: it is not merely a
failed coefficient optimizer. Its actual atomic-core calculation is separate.

### Completed bounded radial check, not a converged material fit

Before introducing another energy term, the existing five radial/saturation
shapes are varied deterministically within the inherited numerical windows.
The minimax objective contains all27 old fit curvatures plus the two declared
opening witnesses above, not a desired core width, yield stress or crack rate.
All seven anchors and coefficient signs remain. The zero-screening exponent
stays fixed. These windows are not independently measured Al parameter bounds.
The selected-state Bessel matrix is required to match the corresponding rows
of the FULL existing matrix before optimization. It is not a surrogate.
The separate spectral and full-data/core checks remain necessary afterward.

160 actual profiles ran in908.009s;154 have strictly positive LJ amplitudes.
Powell stopped at maxfev, not at a verified shape optimum. The retained trial
has shape

    (k_scalar,k_odd,k_even,alpha_even,k_rank1,power)
      =(1.6563525852,6.1294145190,11.9987627137,
        616.94046646,7.7272930862,0).

Its selected29-jet minimax error is9.983579, down from18.464447. However, the
full115-jet error of THAT vector is40.233385 at v17_direct023_Hxx. Improvement
on selected targets is not improvement everywhere. An independent LP with
ALL115 curvatures at this new shape gives13.923143. This improves the fixed-
shape necessary lower bound but still cannot reach the unit target box.
There is no global radial-convergence or whole-family impossibility claim.

The entire285-observation Bessel matrix, sampled bulk stability/tail operators
and two independently restored source-core configurations were then evaluated.
Subset-to-full prediction replay error is0; actual model jets at critical
states agree with the coefficient matrix to8.89e-16 eV/L0^2. The existing
seven exact constraints were NOT replaced by the two troublesome curvatures.

| New-shape coefficient protocol |Other-interface squared loss|Maximum115-jet error|R8 force RMS eV/L0|Excluded R16 configuration RMS|
|---|---:|---:|---:|---:|
| Retained selected minimax trial |9956.044|40.2334|.623487|.592920|
| Original static LS |4318.881|25.4908|.187314|.179693|
| Same v22 joint LS |4520.295|30.5027|.099569|.096993|
| Joint LS with both source boxes |21284.734|34.1069|.504378|.482114|

All four returned profiles have positive LJ and positive tested tail-subtracted
bulk margins(.005858 to.006189). This is sampled bulk stability, NOT proof of
every interface configuration or the whole Brillouin zone. The R16 core state
is excluded from these coefficient objectives; it is an already inspected
source configuration, not a newly acquired blind Al validation dataset.
Frozen gradients are not actual re-relaxations of these candidate cores.

The new joint force RMS is65.53% below the original v20 candidate and34.67%
below v22 wider, but other-interface loss roughly doubles relative to v22
wider. It predicts Hxx=+.382221 at the SOURCE saddle and Haa=-.272193 at the
inspected opening point. It is saved separately for static tests, not adopted.

### One explicit minimal analytic ablation, rejected

The previously derived v14/v15 rational quartic was NOT adopted then. Here it
is explicitly retested with the later full environmental family at the fixed
v22 wider shapes; it is not silently merged into the default model. See
`TAIL_CONTROLLED_AL_CALIBRATION.md` and
`INTERFACE_CALIBRATION_DEVELOPMENT_V15.md` for those earlier negative trials.

For EACH atom, after its infinite Bessel environment sum, let

    I3_i = Q3_i:Q3_i,
    E3,quartic = K3 sum_i H_alpha(I3_i),
    H_alpha(I)=I^2/(1+alpha I), alpha>=0,
    H_alpha'=I(2+alpha I)/(1+alpha I)^2,
    H_alpha''=2/(1+alpha I)^3.

All first and second coordinate derivatives still follow the exact chain
rule. At an affine centrosymmetric bulk state Q3=0. Around the pristine
interface I3=O(delta q^2), so H=O(delta q^4): no initial Hessian, equilibrium,
cohesion or harmonic bulk column is added by changing alpha. Large alpha
gives H approximately I/alpha at nonzero I, so the K3/alpha amplitude becomes
nearly indistinguishable from D3. It is not independent physics by fiat.

Use beta=alpha I_ref with I_ref=.001090396280746736, the actual maximum site
invariant at the pre-existing source saddle. This is a dimensionless moment
gauge, NOT a physical length, area or fitted strength. Five declared beta
values0,.1,1,10,100 were actually profiled with the same seven equalities,
coefficient signs and bulk-tail checks. Static losses were2115.902,2141.428,
2233.879,2234.272,2225.861. The first two optima are zero-LJ closures and are
ineligible. These are static-only losses, not the v22 two-block joint loss.

For EVERY beta the115-jet minimax lower bound remains18.464447, with K3=0
at its optimum. The proposed extra shape therefore does not fix the measured
compatibility witness. D3/K3 normalized-column correlation rises from.904417
to.999929, and the two-column singular-value ratio falls from.224032 to
.005950. This adds poor identifiability without the required improvement.
It is rejected for this tested configuration, not ruled out at every shape.

`RationalQuarticStaticInterface` is a separate research class. The existing
22-channel core bridge now explicitly refuses nonzero quartic saturation;
otherwise it would silently omit this unimplemented law. Default saturation0
and all existing energy behavior are unchanged. Tests verify the zero limit,
finite-state analytic derivatives and this refusal. No new canonical term.

### Topology must be tested at actual stationary points

A wrong Hxx sign evaluated at the SOURCE's stationary coordinate does not
by itself prove that the CANDIDATE lacks an index-one saddle: its stationary
location can move. Conversely a root with a small force residual can still
be a minimum or a higher-index saddle. `run_v23_static_topology.py` therefore
re-solves pristine/fault/saddle states with the full3x3 Hessian, checks downhill
endpoints, scans the first fixed-registry opening traction extrema with two
bracket resolutions and executes a signed0,+50/-50MPa normal,+4/-4MPa shear
static return protocol. These are static tests; no time hold or fatigue claim.
The normal-only traction maximum is not called a coupled spinodal or yield.

The source Al99 curvatures are also checked by step refinement. They remain
matched atomistic-model observations, not experimentally exact Al stiffness
measurements at every registry/opening coordinate. A difference from them
must be reported quantitatively, not used as proof that Bessel sums are wrong.

The ACTUAL full-vector stationary solve and downhill checks completed for
source, v22 wider and v23 radial-joint candidates. All9 stationary states have
the expected full Morse index; all3 saddles connect distinct lower minima.
In particular, v22 wider's own saddle has Hxx=+.068616 but minimum Hessian
eigenvalue=-1.387798 eV/L0^2. The full coupled Hessian, not a diagonal entry,
is the stability test. Any earlier implication "positive Hxx means no slip
saddle" is explicitly withdrawn. The missing physical validation remains
quantitative fidelity and kinetics, not the mere existence of a saddle.

| Matched rigid-interface quantity |Al99 source|v22 wider|v23 radial joint|
|---|---:|---:|---:|
| Relaxed fault energy J/m^2 |.150479|.079457|.075309|
| Adjacent saddle energy J/m^2 |.172002|.119750|.114757|
| Reverse barrier J/m^2 |.021523|.040292|.039448|
| Finite large-opening W(40h), J/m^2 |1.741285|1.840583|2.524404|
| First fixed-registry normal traction maximum MPa |12969.568|8924.474|6536.542|
| a/h at that first maximum |1.185846|1.120683|1.210490|

The near-h geometric+uniform33/65 bracket studies agree on all three FIRST
normal maxima. They do not fully resolve every later source feature: a later
minimum shifts about.000588 in a/h and the fine grid finds an additional
small late maximum. Neither the source nor the candidate traction is clipped.
The source-curvature FD errors at step1e-5 are2.475e-7(saddle) and3.549e-9
(opening) eV/L0^2, far smaller than the candidate jet mismatch. The source
interpolation was not altered to improve agreement.

All15 signed mixed-traction states(0,+50/-50MPa normal,+4/-4MPa shear in both
registry directions, static return to zero) pass local force/Hessian checks.
This is an equilibrium loading-path test, not model-time integration or
physical-time unloading/hold. Its small returns do not certify residual
plasticity. The newly reduced ideal coherent opening peak is not an Al yield
strength calibration: the interface is still uniformly translated/opened.

### Actual nonlinear core tests are separate

The wider-shape both-box material's stationary-root attempt stopped at its
600s budget. An independent checkpoint evaluation gives force.242865 eV/L0,
not a converged root. The same-operator energy-descent retry then stopped at
900s with force1.863155 eV/L0 and relative energy-48.332675 eV/straight repeat.
These are unfinished atomic relaxations, NOT verified minima or proof of an
unbounded energy. They must not be hidden behind the successful coefficient
KKT certificates. Both checkpoints and exact failure causes are retained.

The new radial-joint material, starting from the SOURCE's coordinates only
with its own exterior/energy, reaches a stationary core in7 Newton iterations:
force1.7161e-8 eV/L0. Its smallest Hessian eigenvalue is-.380258829, with
independent energy differences approaching that negative curvature. Thus a
small force residual does NOT make it a stable core. The L-BFGS mode-descent
trial encountered an unconverged reciprocal series and stopped without a
repaired force. This is a numerical trial failure, not a demonstrated physical
spinodal; a separate same-energy Newton-CG descent is recorded independently.

The Newton-CG retry actually converged in17 iterations/487.82s, force
4.343e-10 eV/L0 and minimum Hessian eigenvalue+.172653874 eV/L0^2. The winding
is1 and direct energy differences confirm positive curvature. Thus the
failed L-BFGS trial was not proof that no finite-boundary minimum exists.
This validates one R8/ring5 static minimum of the rejected material, not an
infinite-domain core or a calibrated Al energy. Separate neighborhood and
domain comparisons are required and saved separately.

The same minimum was then actually re-relaxed, not merely re-evaluated:

| Free radius / environment ring |Free atom rows|Maximum force eV/L0|Minimum Hessian eV/L0^2|Elapsed seconds|
|---|---:|---:|---:|---:|
| R8 /5 |286|4.343e-10|.172653874|487.823|
| R8 /7 |286|3.297e-8|.172236076|259.084|
| R10 /7 |444|2.771e-12|.130868746|461.814|

All three have winding1, independently checked positive energy curvature and
zero saved-energy replay discrepancy. Within R2 the largest vector change is
2.1912e-5 L0 for ring5->7 and1.9930e-3 L0 for free radius8->10. The latter is
not zero and a two-radius comparison does NOT certify the infinite domain.
The quarter-crossing sampled registry spans are.885664,.885690,.880949 L0;
the same source R8 diagnostic is2.760327 L0. These discrete interpolated spans
are not continuum core widths, partial separations, fitted targets or proof
of a unique minimum. The candidate still remains atomically narrow.

Independent frozen-source force RMS values are.0995813535,.0995688748 and
.0995676686 eV/L0 at environment rings5,7,9. The ring7->9 difference RMS is
5.8643e-6 eV/L0, whereas the remaining source-force error is.09957. Tightening
the ring9 reciprocal tolerance from2e-12 to2e-14 changes no displayed value.
The mismatch is much larger than those tested numerical changes; untested
nonlinear tails and infinite free-domain errors are not declared absent.

### Short-wave restoring forces: a concrete remaining material defect

Static cubic paths Gamma-L and Gamma-K perpendicular to the straight screw
line were tested at11 wavevectors, two controlled neighbor radii12/16 and
three candidate materials:66 full3x3 matrices. The existing analytic harmonic
tail is included; this is not a finite-cutoff replacement of LJ/Bessel energy.
Source matrices in eV/Angstrom^2 are multiplied by(L0/Angstrom)^2. The corrected
+ABC cubic basis is used on both sides. No mass, frequency or time is inferred.
These are post-fit audits, NOT blind held-out certification: overlap with some
earlier selected bulk projections can exist. The initial overly broad
`used_in_fit=False` metadata was corrected and the study actually rerun as
`bulk_dispersion_audited_partition` without changing any numerical matrix.

For the radial-joint candidate, the matrix-relative errors near Gamma at
q=(.005,.005,.005) and(.005,.005,0) are.6325% and.6353%. At q=(.5,.5,0) the
FULL matrix relative error is470.812%; its screw projection alone differs
by-13.1325%. Do not confuse those two metrics. In the largest-discrepancy
polarization the source and candidate curvatures are60.086257 and375.872020
eV/L0^2, respectively. At this q the radius12->16 matrix change is8.4149e-6
eV/L0^2 and the radius16 harmonic tail bound is3.6665e-5. This error is not
explained by the tested truncation. Positive bulk Hessians do not establish
correct short-wave material stiffness or a correct nonlinear defect core.

An independent actual energy calculation checks the matrix itself. At this
q the FCC sites occupy four exact phase classes theta_i. For a displacement

    u_i=epsilon e cos(theta_i),
    r_ij(epsilon)=R_ij+epsilon e[cos(theta_i+q.R_ij)-cos(theta_i)],
    mean Delta E=epsilon^2(e.H(q).e)/4+O(epsilon^4).

Each site's LJ pairs, scalar density, Q1/Q2/Q3 and nonlinear invariants are
recomputed before averaging over the four sites. Source EAM uses its own
published finite support only as a target comparator. Candidate direct
neighbor radii12/16 are independent validation, never the canonical generator.
The two displacement polarizations and two amplitude studies give60 actual
curvature evaluations. In the large-error direction, candidate energy-FD
curvature at epsilon=.0005,.000125,.00003125,.0000078125 L0 is
368.593978,375.407957,375.842980,375.870193 eV/L0^2. It approaches the analytic
375.872020, with finest relative discrepancy4.861e-6. The source approaches
60.086257 (finest discrepancy about1.0e-6 relative, with eventual roundoff).
In the screw direction the finest candidate discrepancy is2.68e-8 relative.
The large source/candidate difference therefore cannot simply be dismissed
as a Bloch-matrix prefactor or coordinate-conversion bug in this check.

The large even-invariant saturation shape also yields a very narrow harmonic
amplitude window: epsilon=.002 L0 underestimates that candidate curvature by
23.53% in energy FD, whereas smaller steps recover it. This is the behavior
of the same smooth energy, not a force repair. It is an additional warning
against interpreting a fitted tangent or a narrow core as a credible Al law.
This does NOT by itself prove that one short-wave polarization causes the
entire nonlinear core discrepancy. No extra damping or scalar strength factor
is used to hide it.

### Interpretation and next scoped question

This round identifies a real coefficient/shape tradeoff and corrects an overly
strong diagonal-curvature interpretation. It does NOT prove the entire
analytic family incapable of fitting Al. It does rule out adopting the
actually inspected vectors as full Al material calibrations merely because
one force loss improves. No parameter vector replaces the preserved default.

The next derivation should identify which partial-coordination scalar or
angular derivative freedom is missing AFTER the seven exact anchors, using
these saved dual witnesses, actual stationary topology AND short-wave vector
stiffness. A next fit should declare full short-wave validation components
before fitting, test rank/tradeoffs, and not improve core forces by making
other polarizations artificially stiff. The v14 cubic-
density, v15 positive-mixture/cross/rational and this v23 quartic ablation
are negative controls, not license to add all terms together. Independent
material/source uncertainty is still distinct from optimizer discrepancy
scales. Any new term requires an explicit per-atom energy, analytic derivatives,
gauge/rank audit and repeat full-interface/core validation. A core force fit
alone does not justify a finite-loop activation law, yield or fatigue claim.

Current a/s mobilities M_a,phys and M_s,phys, t0, physical seconds/Hz and
specimen correlation area remain unavailable/uncalibrated. Straight-line
energy per repeat is not a finite defect activation energy. No production
PDE, strain/absorption bookkeeping, stress map, A_c path or UI default changes.

The first launcher failed at constructing a Python result dictionary after
its first LP, not due to physical infeasibility. That record and an explicit
corrected retry remain separate. Two core-seed preflight failures (historical
missing geometry metadata and mismatched free domains) were likewise retained
without treating them as atomic instabilities or changing the input guards.

Source provenance remains Mishin et al.(1999), PRB59,3393,
DOI10.1103/PhysRevB.59.3393, with the unchanged target-only Al99 source from
[NIST](https://www.ctcms.nist.gov/potentials/entry/1999--Mishin-Y-Farkas-D-Mehl-M-J-Papaconstantopoulos-D-A--Al/).
It supplies matched0K atomistic comparisons, not calibrated macroscopic yield.

## Executed verification and reproduction status

Targeted48 tests passed in99.99s, followed by2 topology tests in1.40s.
The FINAL complete solver run passed655 tests in1945.86s; app passed34 in
167.68s with zero skips. The actual desktop smoke passed in1.523s and printed
the preserved a0=.7713438268704838, kappa=86.29296488740997. The earlier full
attempt was intentionally interrupted after test collection changed and is
NOT counted as a pass. Completed JUnit hashes, protected input hashes, nine
explicit failed-attempt records and all result bytes are inventoried in
`results/core_interface_compatibility_v23/validation_ledger`. No atomic run
remains live or unaccounted in that ledger. Tests certify the stated numerical
contracts, not Al-material adoption or actual yield/time.

The first ledger runner encountered the two historical `error`/`message`
failure-record schemas. Both schemas are now read explicitly; no failed case
is ignored. The complete ledger and desktop smoke were actually rerun. No
scientific quantity was changed by this report-only correction.
