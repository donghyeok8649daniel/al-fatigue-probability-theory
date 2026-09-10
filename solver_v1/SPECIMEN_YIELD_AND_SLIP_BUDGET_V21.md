# From source motion to specimen yield: signed slip and the strain budget

Status: completed scoped research calculation, **not calibrated specimen yield**.
The latest user request is to move beyond ideal-strength agreement toward actual
Al yield. This study addresses the missing source-to-specimen strain connection;
it does not restart the static material fit, lower an ideal barrier, or add a
phenomenological yield/hardening law. Production energy/PDE/UI and physical-time
calibration are unchanged. The work starts at `a0b611b1`.

## 1. Compare the same mechanical observable

The previously reported10.33897GPa is a **fixed-registry ideal normal opening
traction maximum** of the v20 research candidate. It is neither its ideal shear
fold nor a tensile proof/yield stress. Dividing it by a specimen's reported YS
does not measure a yield-prediction error. This distinction does not excuse a
model that fails to produce specimen plastic strain at an appropriate stress.

The refreshed primary experimental source is Pigato et al., *Microstructural and
Mechanical Characterization of Ultra-Pure Aluminum for Low-Amplitude-Vibration
Cryogenic Applications*, Materials19(6),1195(2026),
[doi:10.3390/ma19061195](https://doi.org/10.3390/ma19061195).
The original EuropePMC XML is downloaded, parsed and hashed; no graph-digitized
or guessed number is substituted. Table2 supplies12 YS values at25,-50,-100,-150C.
The three starting microstructures and temperatures remain separate.

| Material at25C | Source-reported YS [MPa] | Initial state |
|---|---:|---|
| 6N,99.9999wt%Al |16.4 +/-1.3|recrystallized/equiaxed|
| 5N5,99.9995wt%Al |74.6 +/-3.7|deformation-fragmented substructure|
| 5N,99.999wt%Al |34.1 +/-1.1|recovered/subgrain structure|

Tension follows the rolling direction at0.00025/s (+/-20%). The paper describes
95% confidence for its reported uncertainties; rows without uncertainty are
not assigned one. Although ISO6892 is cited, the retrieved text does not state
the numerical offset used for YS. We therefore leave the plastic-strain criterion
null: **reported YS is not silently relabeled Rp0.2**. Individual slip systems,
source pin spacings and mobile-source populations are not supplied. The figures
show texture/microstructure; they do not supply these missing values.

These are condition-specific validation records, not loss targets for LJ,
embedding, mobility, or an adjustable source length. The existing Krebs2017
0.002-*plastic-shear* CRSS/G benchmark is preserved as a different observable.
No room-temperature specimen is declared matched to a0K material comparator.

## 2. Exact kinematics of resolved slip surfaces

For an oriented slip surface S_k inside an explicitly defined specimen volume V,
let b_k be its Burgers vector [m], n_k its unit normal, and A_k its **signed swept
area** [m2] relative to the initial state. A negative A means reversed motion in
the fixed sign convention. At small strain, the distributional displacement jump
gives the volume-average distortion and strain

\[
 \overline{\beta}^{slip}={1\over V}\sum_k A_k\,\mathbf b_k\otimes\mathbf n_k,
 \qquad
 \overline{\epsilon}^{slip}={\overline{\beta}^{slip}+
                  \overline{\beta}^{slip,T}\over2}.
\]

Tangential b implies trace(epsilon_slip)=0. For b_k=b m_k and axial unit vector e,

\[
 \epsilon^{slip}_{ee}={b\over V}\sum_k A_k
             (\mathbf e\cdot\mathbf m_k)(\mathbf e\cdot\mathbf n_k).
\]

This is a geometric identity, not a law for when A grows. It applies to arbitrary
supplied slip systems in one common physical frame; it does not assume source
independence or replace the existing local surviving-ensemble strain fields.
Simultaneously reversing b and A, or n and A, preserves the physical result.
Doubling the physical averaging volume AND all resolved swept surfaces preserves
strain. Enlarging volume alone reduces the contribution of a single event.

With uniform symmetric stress, independent tensor and resolved-traction forms
must give the same virtual work:

\[
 V\boldsymbol\sigma:d\overline{\boldsymbol\epsilon}^{slip}
   =\sum_k (\mathbf b_k\cdot\boldsymbol\sigma\mathbf n_k)\,dA_k.
\]

Nonuniform stress requires the local surface integral; applying this uniform
identity to a spatially varying field is not implemented. The tests include
normal/shear mixed tensors, signs, rotations, superposition, engineering-shear
factor2 and independent work evaluation. No stress conversion is fitted.

If a resolved dislocation line sweeps its surface at normal speed v_perp,

\[
 \dot A_k=\int_{\mathcal L_k}v_\perp\,d\ell,\quad
 \dot\gamma={b\over V}\sum_k\int_{\mathcal L_k}v_\perp\,d\ell
            =b\rho_{mobile}\langle v_\perp\rangle_{line}.
\]

The final equality is the Orowan **kinematic identity**, not a supplied velocity,
mobile-density evolution, source-emission rule, or fitted constitutive law.
Those quantities must still follow from the energy/dynamics and matched initial
microstructure. Ordinary total diffraction dislocation density is not
automatically rho_mobile. Similarly, a wire diameter is not a pin separation.
In the CSV names, `line_density_m2` denotes line length per volume [m^-2],
and `source_number_density_m3` denotes count per volume [m^-3], not areas or
volumes. Swept areas themselves use `area_m2` [m^2].

## 3. Recovery, net transport, and gross area

`transport_history` computes changes relative to the initial signed area, not
relative to an assumed defect-free state. For resolved snapshots,

\[
 A^+_j=\sum_{i\le j}\max(\Delta A_i,0),\quad
 A^-_j=\sum_{i\le j}\max(-\Delta A_i,0),\quad
 A_{net}=A^+-A^-=A_j-A_0,\quad A_{gross}=A^++A^-.
\]

The max operations define positive/negative measures; they do not clip a strain
or repair an energy. Gross area is never used in place of signed A in strain.
Unobserved reversals between saved snapshots make sampled gross area a lower
resolution description. This quantity is not the divergent nearest-neighbor SG
recrossing traffic and is not asserted to count committed atomistic transitions.

A pinned dislocation can bow and retract, producing nonzero slip strain during
loading but zero net change after relaxation. Calling all such strain permanently
plastic would overcount recovery. The API therefore supplies `signed_slip_strain`
and does not automatically certify residual plasticity from any nonzero value.
Similarly, a loading-branch offset crossing is not a separate proof of permanent
registry change after hold. No scalar threshold is used to erase a residual.

## 4. Connection to deterministic probability evolution

This mapping does not require Monte Carlo or empirical plasticity. If a future
resolved defect configuration Gamma has deterministic density P(Gamma,t), apply
the same observable B(Gamma)=sum A_k b_k tensor n_k/V to that density. Then

\[
 {d\over dt}\int B P\,d\Gamma
 =\int\nabla_\Gamma B\cdot\mathbf J\,d\Gamma
  -\int_{\partial\Omega}B\,\mathbf J\cdot\mathbf n\,dS.
\]

The first term transports slip; an absorbing boundary removes its carried
moment. For survival S and opening outflux j_open, conditional differentiation
adds the normalization term, giving the selective-opening contribution

\[
 {1\over S}\left[-\int_{opening} B j_{open}\,dS
                  +\langle B\rangle_{surv}\int_{opening}j_{open}\,dS\right].
\]

It is not plastic flow. This is the same integration-by-parts structure already
used by the registry-moment audit. The current two-coordinate production density
does not contain a spatial population of emitted curved lines; this derivation
does not silently claim that it does. No new PDE generator is registered here.

## 5. Finite-source calculation from unchanged LJ/Bessel elasticity

Use the saved v20 old-family **research** candidate with its already fitted
C11/C12/C44=114/62/32GPa at0K. Its parameter file is byte-hash bound and untouched;
the full interface material gate still fails. The line coefficient follows

\[
 k(\theta)={\mathbf b^T\operatorname{Re}K_0\mathbf b\over2\pi},\qquad
 \gamma(\theta)=k(\theta)\ln(R/r_{core}),\quad
 T=\gamma+\gamma''.
\]

This is the SAME infinite-energy bulk Hessian -> half-crystal Schur -> long-wave
line-energy construction. Core energy and finite nonlocal parts remain omitted,
not measured to be zero. The current candidate's initial interface tangents do
not independently determine these terms.

Use the declared edge-oriented double-ended reference, HYPOTHETICAL spans
L=0.2/1/5um, R=L and r_core=2b. Nothing is selected to match experimental YS.
For the positive-T branch, define Q=gamma sin(theta)+gamma' cos(theta),
r=-gamma cos(theta)+gamma' sin(theta), p=tau b. At endpoint theta_e,

\[
 p={2Q(\theta_e)\over L},\quad x={Q(\theta)\over p},\quad
 y={r(\theta_e)-r(\theta)\over p},\quad
 A=2\int_0^{\theta_e}y(\theta){T(\theta)\cos\theta\over p}\,d\theta.
\]

Adaptive dimensionless quadrature computes A; an independently minimized
piecewise-linear line gives a second value. For constant tension,
A=R_curve^2(theta_e-sin(theta_e)cos(theta_e)), and its fold gives pi L^2/8.
For the present anisotropic edge reference, A_fold/L^2=0.5897092570.

At theta_e=pi/2, tau_c=2gamma(pi/2)/(bL). Explicit maximum-{111}<110> Schmid
projections give the following **conditional first-source stresses**:

| Hypothetical L [um] | Resolved shear [MPa] | Axial[100]/[110] [MPa] | Axial[111] [MPa] |
|---:|---:|---:|---:|
|.2|38.4914|94.2844|141.4266|
|1|9.81419|24.0398|36.0596|
|5|2.38602|5.84453|8.76679|

These use favorably oriented equivalent sources, not measured texture or an
observed active system. They are not experimental yield predictions.72 actual
2–50MPa/zero stress cases record subcritical areas or absence of a subcritical
branch; above-fold rows remain undefined rather than inventing an emitted loop.

## 6. The newly quantified obstruction: first source is not proof strain

For identical same-sense sources with count N in V, define initial pinned-line
density rho_line=NL/V and geometric overlap parameter eta=NL^3/V=rho_line L^2.
This rho counts initial straight pinned length, not the current curved length.
For a selected Schmid factor m,

\[
 \epsilon_{bow}=m b\rho_{line}{A\over L}
               =m{b\over L}{A\over L^2}\eta,
 \qquad
 \eta_{required}={\epsilon_{criterion}L\over m b(A/L^2)}.
\]

The last expression is a **strain demand**, not a calibration of density.
We explicitly test epsilon_criterion=.002 as a conventional axial demand,
without asserting that it is the unstated criterion in Pigato's YS values.

For the illustrative dilute-spacing scenario eta=.01, L=1um and [100], even
the fold has epsilon_bow=6.89449e-7, only1/2900.87 of that demand. Reaching
.002 by simply adding identical instantaneous bows requires eta=29.0087.
Across the nine cases eta_required ranges5.8017–217.5649. The implied typical
3D source separation is smaller than the source span; the **declared dilute
independent-source approximation** cannot justify that extrapolation.
Eta<<1 is a spacing heuristic, not a universal theorem of zero interaction.
Special correlated arrangements require their own interaction calculation.

This prevents a false repair: multiplying the number of reversible sources
until an apparent proof strain appears cannot certify real yield. Even if a
loading offset were reached, these bows recover. Continued slip requires
resolved source operation/emission, propagation/escape and interaction, not
just a lower atomistic barrier or a fitted source-count multiplier.

Nor does one measured first-source stress identify the source geometry. With
the outer convention R=L and g=ln(L/r_core), the logarithmic sensitivities are

\[
 \frac{\partial\ln\tau_c}{\partial\ln L}=g^{-1}-1,\qquad
 \frac{\partial\ln\tau_c}{\partial\ln r_{core}}=-g^{-1}.
\]

One threshold supplies a rank-one row for two unknown geometric parameters;
source population does not enter that threshold at all. Fitting L to a desired
YS with a guessed core cutoff would neither identify the microstructure nor
validate its strain production. This is an analytic non-identifiability, not
a numerical optimization failure. Additional independent source/core/transport
observations are required. The current study does not fit along this null family.

## 7. Executed convergence, cycles and hold

The run executes297 static cycle states (three cycles per orientation/span),
36 independent spatial-refinement solves,9 time/grid hold histories,72 stress
cases and12 parsed experimental YS records. Angular32/64 coefficient change is
8.27e-25J/m; this is numerical precision, not material accuracy.

The independent line-area errors divided by L^2 at theta_e=1 decrease

    3.55951e-4 -> 8.91282e-5 -> 2.22909e-5 -> 5.57328e-6
    segments:32,64,128,256.

At the nearer-fold angle1.4 they decrease6.59562e-3 ->1.92400e-3 ->
5.10591e-4 ->1.30005e-4. A finite graph is not called exact at the fold.
Quadrature1e-9/1e-11 and force residuals are saved independently.

All quasistatic zero-load returns have zero incremental slip area by the
resolved reversible branch symmetry, while gross area increases over cycles.
This is a static return, not a physical-time hold. Separately, the existing
linear pinned-line diffusion equation is actually integrated with32/64/128
segments and64/128/256 steps per period, six cycles and >24 line relaxation
times at zero stress. At the finest pair the axial-bow peak is3.87441e-10,
the remaining value is-1.25857e-20 (3.25e-11 of the peak), and maximum slope
is.001645, consistent with the small-slope scope. This nonzero decaying
transient is **not** a claim of resolved residual plasticity.

The sampled hold endpoint is24.05282tau at64 steps/period and24.00373tau at
128/256. These are all sufficiently long decay checks, not an exactly common-
endpoint convergence table for the tiny final residual. No residual-plasticity
claim is based on their floating-point nonzero values.

This separate linear-line check uses0.01MPa resolved-shear amplitude and an
explicit ideal test projection m=1/2. It is not the[100] static stress sweep or
a prediction for the experimental specimen's unknown texture.

The time calculation reuses the explicitly conditional v18 Gorman line-drag
reference with the0K elastic comparator. We output t/tau, not new production
seconds. Temperature/geometry matching and atomistic core dynamics are not
certified by this numerical hold. There is no opening probability in this
line model. Data are not mislabeled as full P(a,s) evolution.

## 8. What is implemented, and the next physical step

Implemented here: signed multisystem tensor strain, work consistency,
reversible/gross separation, analytic minor-arc area, strain-budget audit,
actual signed stress/refinement/hold studies and a provenance-bound measured
strength reference. There is no empirical hardening, source population fit,
stress reduction factor, material refit, activation-volume input or Monte Carlo.

To obtain a defensible specimen yield, the next defect state must support
finite-source **emission and transported slip after loss of the minor branch**,
with a same-energy validated core and finite/intersegment elasticity. Matching
measured source/obstacle populations and strain protocol then permits a
like-for-like yield criterion. A just-emitted line can escape or be stored;
neither event can be prescribed to produce the desired residual strain.
The existing material/interface discrepancy remains an independent gate.

Atomic interface cell area, statistical correlation area, slipped surface and
specimen volume are distinct objects. A_c is not an argument anywhere in this
calculation. Existing local probability and specimen probability aggregation
remain unchanged. M_a,phys, M_s,phys and t0 remain unavailable; production
physical seconds/Hz stay disabled. No UI redesign or energy-model promotion.

Reproduce the actual research run (fresh result directory required):

    python -m solver_v1.run_specimen_yield_bridge --out results/specimen_yield_bridge_v21 --source-xml .cache/specimen_yield_bridge_v21/pigato2026.xml

The original article stays in the ignored cache; compact extracted facts,
source and material hashes, all hypotheses, actual results and a plot are in
`results/specimen_yield_bridge_v21/`. An independent complete replay reproduced
all six CSVs byte for byte. Targeted62 passed(1.75s), full solver591 passed
(988.27s), app34 passed(65.89s), desktop smoke exited0(1.28s); no skipped tests.
Tests, actual timings and byte hashes are recorded separately in verification.json.
