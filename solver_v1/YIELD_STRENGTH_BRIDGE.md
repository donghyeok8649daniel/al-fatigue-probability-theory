# Toward measured yield: independent elasticity and a finite-source reference

Status: `yield_bridge_v11`, **research only; measured Al yield NOT reproduced**.
The purpose is to make a condition-specific measured yield the validation
objective, rather than lowering a pristine uniform-interface ideal fold and
calling it yield. Production LJ/Bessel, probability theory, static parameter
files, mobility, specimen probability aggregation and UI are unchanged.

## 1. Which measured strength?

The current benchmark is 99.99% as-cast Al single-crystal wires at room
temperature, tensile displacement rate300nm/s. Krebs et al. define their
Figure2d CRSS at **0.002 plastic shear strain**. Seven isolated markers are
now read directly from that figure, with pixel/provenance hashes. This is a
different observable from the previous Figure2b large-strain flow points.
The measured ordinate is CRSS/G; its normalization modulus is not numerically
specified in the retrieved article/SI, so the MPa column remains null. No
guessed modulus, axial Schmid factor or source length is substituted.
[Krebs et al., Nature Materials16,730–736 (2017)](https://doi.org/10.1038/nmat4911),
[accepted article](https://eprints.whiterose.ac.uk/id/eprint/126662/1/Article%20Text.pdf).

The official [supplement](https://infoscience.epfl.ch/server/api/core/bitstreams/2cff7b76-19e7-4591-b7ca-370af1bcd090/content)
was retrieved through the repository ORIGINAL bundle after the legacy URL
failed. Its single-arm model assumes L=D/3, discusses D/2, and an orientation-
dependent coefficient .068–.148. These are **model assumptions**, not measured
source geometry. They are not inserted as our yield law. Its DD mobility is
not our collective-coordinate mobility. The paper's single-arm source also
is not the double-ended pinned line derived below.

Annealing is not an interchangeable label for these specimens. The follow-up
study reports oxide growth and higher/more variable micro-wire flow after
500C/2h treatment; its condition is kept separate, not imported as a generic
bulk-Al yield target.
[Verheyden et al., Scripta Materialia161,58–61 (2019)](https://doi.org/10.1016/j.scriptamat.2018.10.009).

The selected normalized CRSS points span approximately1.51e-4–5.22e-4, with
reading halfwidth8.57e-6. These are selected isolated points, NOT the full
population range, statistical uncertainty, or seven identical microstructures.
Specimen orientation, pin locations, core state and protocol must match
before a model/experiment prediction error can be claimed. The existing
`StrengthConditions` guard remains in force.

## 2. An elastic metric issue, not a change of units

Keep the existing exact LJ + per-site environment energy and coefficient
identity y=B(d_rho,d_Q)c, c=(u,v,A,B,C,D3,D1,D2). Define

\[
 H=Q C,\qquad C=(C_{11},C_{12},C_{44})^T,\qquad
 Q={V_{atom}\,10^9\over eV_J}
 \begin{pmatrix}3&6&0\\1/3&2/3&4/3\\1/3&-1/3&1/3\end{pmatrix}.
\]

H contains hydrostatic, [111] normal and engineering-simple-shear strain
Hessians in eV/atom; C is in GPa at zero pressure. Q has rank3. In particular,
individual C components involve subtraction of these modes. Small relative
errors in each H do **not** guarantee small relative errors in each C.

If the old diagonal discrepancy covariance is Sigma_H, a coordinate-only
change requires

\[
 \Sigma_C=Q^{-1}\Sigma_H Q^{-T},\qquad
 \delta C^T\Sigma_C^{-1}\delta C
 =\delta H^T\Sigma_H^{-1}\delta H.
\]

Dropping the off-diagonal entries changes the statistical/discrepancy model.
The new diagonal5%-per-C comparison is thus an **explicit alternative
model-discrepancy metric**, not a claim that the old unit conversion was
wrong. Neither metric is measured uncertainty. All nonelastic targets,
fixed-density gauge, held-out roles and per-atom normalization are retained.
Yield data, mobilities and desired event ordering never enter the loss.

Two actual deterministic studies are run:

1. Exact equilibrium/cohesion; least squares with independently scaled C.
2. All five bulk targets exact; minimize the remaining interface residuals.

For the second case, at fixed radial decays there are three coefficient null
coordinates. The convex quadratic problem can be solved by enumerating its
independent active faces (at most3 of5 nonnegative coefficient bounds).
Each face fixes the corresponding coefficients to zero by parameterization,
then solves the remaining exact constraints and least squares. No energy or
force is clipped. This replaces initial SLSQP roundoff-sensitive bound
classification; the initial run is retained, not called the final fit.
Every tested decay records feasibility, the selected active bounds, the number
of feasible faces, exact residuals and the local SVD. The individual rejected
active-face solutions are not saved as material candidates.

The same predeclared radial grid and bounded Powell searches are executed.
Budget exhaustion and local convergence are recorded separately; neither
proves a global optimum or impossibility of the entire analytic family.
An exact bulk match alone never passes interface/held-out validation.
The final table is `yield_bridge_v11/material_summary.csv`; no candidate is
promoted. The full basis, not a fitted gamma lookup table, is evaluated.

### Executed material result

The final active-face run took402.17s:26 predeclared radial starts, then70
Powell evaluations for the C metric (local convergence reported), and100 for
exact bulk (budget exhausted). The latter's best *evaluated* loss223.0087 is
slightly below the optimizer's returned223.0933; both are retained. No radial
global optimum or family-wide impossibility is inferred from this budget.

| Quantity | Historical candidate | Independent-C metric | Five bulk targets exact | 0K source |
|---|---:|---:|---:|---:|
| C11 [GPa] |87.8162|85.7348|114.0000|114|
| C12 [GPa] |64.8326|70.6774|62.0000|62|
| C44 [GPa] |49.2710|32.4645|32.0000|32|
| Relaxed fault energy [J/m2] |.126800|.145773|.027060|.150479|
| Index-one stationary energy [J/m2] |.173238|.178987|.042449|.172002|
| Opening at40h [J/m2] |1.762768|1.798845|1.741053|1.741285|
| Held-out normalized RMS |15.6863|12.2948|25.4113|—|

All three candidate perfect-interface Hessians are positive, with minimum
eigenvalues4.9452/4.9997/4.7104eV/L0². Bulk cohesion is3.36eV/atom and the
equilibrium is a/L0=sqrt(2/3); these constraints do not certify the rest of the
landscape. The independent-C metric improves C44 but misses C11 by24.8%.
Making C exact instead gives fault/saddle energies about82.0%/75.3% too low.
The opening endpoint alone looks excellent in that case, yet its held-out
shape is worse. This is the quantitative reason **neither new fit is adopted**.
New stationary roots have the stated Morse indices; global minimum-energy
connectivity and full-Brillouin-zone stability are not newly certified here.

The exact-bulk best coefficient vector is

    (u,v,A,B,C,D3,D1,D2) =
    (1.2186187186, 11.3791470207, 0, 71.8188898934,
     .3343117520, 77.8502469712, -33.1401737520, -118.9403833792)
    d_rho=5.3349561253, d_Q=3.5167312912.

The active A=0 boundary and signed angular terms must not be disguised as
a unique physically interpreted material fit. The five exact bulk rows have
numerical rank5, with maximum normalized equality residual6.02e-13. At fixed
decays the three-column, column-normalized interface-nullspace Jacobian has
singular values1.58840,.688437,.0550942 (condition28.83), before the active-A
tangent restriction. The separate six-coefficient comparisons in
`identifiability.json` remove only the first two bulk constraints, for common
comparison of all fits. Neither analysis includes the two radial directions
or supplies an uncertainty interval under model discrepancy/boundary activity.

The constitutive family, density gauge, and original parameter-file hash
9d00fbf54831c134fe9961e17f0603defdc7958bc3590743b2094264a31b6655
remain unchanged. New vectors are saved only as rejected research results.

## 3. Why finite sources introduce a different stress scale

The earlier uniform W_int(q) moves a complete infinite half-crystal. Even a
well-grounded Al source potential has a GPa fold in that state space. To
approach real yield, a pre-existing defect must move or a finite source must
operate; its restoring force is not the ideal traction on every atom at once.

Start with the **same potential's** bulk elastic tensor and the already
validated two-half-space Schur kernel

\[
 K(\mathbf q)=|\mathbf q|K_0(\widehat{\mathbf q}),\quad
 K_0=(Z_+^{-1}+Z_-^{-1})^{-1}.
\]

Let a straight line tangent be t(theta)=(cos(theta),sin(theta)), perpendicular
in-plane variation n(theta)=(-sin(theta),cos(theta)), and the full Burgers
vector be b_vec=(b,0,0). A step disregistry has nonzero Fourier amplitude
b_vec/(ik). Per unit line length, including both signs of k,

\[
 {E_{outer}\over L_{line}}
 ={1\over2}\int {dk\over2\pi}{\mathbf b^T K_0\mathbf b\over |k|}
 =k_{line}(\theta)\ln(R/r_{core}),\qquad
 k_{line}(\theta)={\mathbf b^T\operatorname{Re}K_0(\mathbf n)\mathbf b\over2\pi}.
\]

K0 has unitsPa, k_line has unitsPa m²=N=J/m. Its coefficient is derived from
the existing LJ/Bessel Hessian; no empirical yield, line-tension prefactor or
new potential is fitted. The isotropic limits are

\[
 k_{screw}={\mu b^2\over4\pi},\quad
 k_{edge}={\mu b^2\over4\pi(1-\nu)},\quad
 k(\theta)=k_{screw}\cos^2\theta+k_{edge}\sin^2\theta.
\]

For anisotropic FCC the existing corrected +ABC cubic-to-plane basis is
used. Pi-periodic angular Fourier interpolation and its first two analytic
series derivatives are refined32/64/128samples, including independent
midpoint evaluations. At high resolution second derivatives reach amplified
floating-point noise; further samples are not claimed to improve that floor.

This is a **long-wave leading-log expansion**, not an atomically exact curved
core. The complete line energy would also require an independently derived
core term and finite nonlocal self/intersegment contributions. Those are
missing; they are not assigned zero as a physical material fact.

## 4. Actual pinned-line equilibrium calculation

For an explicitly specified double-ended source of spanL, set
gamma(theta)=k_line(theta) ln(R/r_core) for the outer-only diagnostic. Keep
R/r_core fixed during this variational problem. The functional is

\[
 E[y]=\int_{-L/2}^{L/2}\gamma(\arctan y')\sqrt{1+y'^2}\,dx
       -p\int_{-L/2}^{L/2}y\,dx,\qquad p=\tau_{res}b.
\]

Both terms have unitsJ. This is finite line energy minus mechanical work,
NOT an infinite-line energy multiplied by a guessed activation length and
inserted into exp(-E/kT). It has no thermal rate or fatigue-lifetime meaning.

The first variation gives the anisotropic line stiffness, not gamma alone:

\[
 T(\theta)=\gamma+\gamma'',\qquad T\,\mathcal K=p,
 \qquad Q(\theta)=\gamma\sin\theta+\gamma'\cos\theta,
 \quad Q'=T\cos\theta.
\]

For mirror symmetry and T>0, the endpoint angle theta_e gives

\[
 p={2Q(\theta_e)\over L},\quad
 x(\theta)={Q(\theta)\over p},\quad
 y(\theta)={R_1(\theta_e)-R_1(\theta)\over p},\quad
 R_1=-\gamma\cos\theta+\gamma'\sin\theta.
\]

These are exact within the declared local-line functional. Its smooth branch
has a fold at theta_e=pi/2,

\[
 \tau_{c,outer}={2\gamma(\pi/2)\over b L}.
\]

The code refuses this simple branch if reflection symmetry or positive line
stiffness is unresolved. A faceted/anisotropic instability needs separate
treatment. The formula is **not added to the constitutive solver as yield**.

An independent piecewise-linear energy minimization solves every graph degree
of freedom. Its segment conjugate force is Q, Hessian coefficient is
T cos³(theta)/dx, and interior balance is Q_left-Q_right-p dx=0. Newton
iteration is a static solve, not seconds/model time. Shape, force residual,
swept area and32/64/128/256segment refinement are saved. Above the analytic
fold we report lack of this branch, not atomistically verified multiplication.

## 5. Hypothetical source sizes: an informative but conditional result

Explicit tests use L=.5,1,2,5,10micrometres, R=L and r_core/b=.5,1,2.
These are a **sensitivity grid**, not measured source lengths, a chosen best
Al radius, or parameter calibration. R is the asymptotic outer cutoff; it is
not the local radius of curvature. Order-one choices change the omitted
finite part. There is no physical default ofL orr_core in the API.

For the unchanged historical candidate and the middle r_core=b convention:

| L [micrometre] | Outer-only critical resolved shear [MPa] |
|---:|---:|
|.5|29.0441|
|1|15.8705|
|2|8.60943|
|5|3.80027|
|10|2.03498|

Thus a finite-source force balance can have an MPa scale without lowering
the atomic ideal barrier by hand. This is progress in mechanism/normalization,
**not proof of agreement with experimental yield**. Model/source elastic
comparators and2/4/10/25/50MPa static scenarios are separate rows. They are
not a fit to these experimental stresses. The material with improved C44 but
incorrect C11-C12 changes the line coefficient significantly, illustrating
why a scalar elastic match is insufficient for source-strength predictions.

The declared r_core/b=.5/1/2 sensitivity atL=1micrometre yields
17.2189/15.8705/14.5221MPa for the historical candidate. This spread is not a
complete error bar: the finite nonlocal part and real core are still missing.
The independent-C candidate's threshold changes because its long-wave tensor
changes, not because any source parameter was fitted to experimental strength.

Actual verification includes125 static stress/length/model combinations:
53 independent subcritical graph solves and72 above-fold exclusions. There
are80 separate32/64/128/256segment refinements. For the historical candidate
at endpoint angle.8, maximum-bow error divided byL decreases

    2.03325e-4 -> 5.13447e-5 -> 1.28710e-5 -> 3.21995e-6.

Near angle1.4 the256segment errors are6.33e-4–9.92e-4 ofL across the models,
so finite-grid graphs are not called exact near the fold. The analytic branch
within the declared leading-log functional remains the comparison reference.
Maximum normalized force residual is1.99e-10. Angular64-to128 coefficient
change is at most1.16e-23J/m; stiffness change is4.82e-20J/m. These tiny
numerical errors do not remove the substantially larger physical modeling
uncertainty. Runtime28.71s, plus18.43s for material stationary/FD and reporting.

Material analytic-gradient/Hessian checks use the actual infinite lattice
surface. Halving the FD step4e-5 to2e-5L0 reduces the errors by approximately4;
the worst fine gradient/Hessian differences are2.07e-8eV/L0 and
5.91e-7eV/L0². The original LJ calibration is untouched.

## 6. Bow-out, first slip, CRSS and residual plasticity differ

Below the fold this line bows and retracts on quasistatic unload. It does not
by itself establish residual plasticity or emitted-loop traffic. The volume
average from a genuinely swept slip area in a specified specimen would be

\[
 \Delta\gamma_p={b A_{swept}\over V_{specimen}},\qquad
 \Delta\epsilon_{p,axial}
 ={b A_{swept}\over V_{specimen}}(\mathbf e\cdot\mathbf m)
                                      (\mathbf e\cdot\mathbf n).
\]

Here V_specimen is actual specimen geometry, not an invented activation
volume. A single first microscopic event need not reach gamma_p=.002. Source
population, multiplication, intersections, surfaces/oxide and loading history
must be resolved to compare to that criterion. The current code deliberately
returns `experimental_yield_prediction=null`. It neither adopts a forest-
hardening law nor fits an effective source length to force a match.

Statistical correlation area, atomistic interface area and swept area are
different objects. No specimen probability parameter enters energy, stress,
elasticity, line tension, strain or the present source calculation.

## 7. Next admissible scientific step and reproduction

The material/held-out interface mismatch remains. The next atomistic step is
a jointly acceptable surface and a topologically/structurally specified
vector/partial core with finite line variation and real pin/free-surface
geometry. Its energy derivatives must reproduce the long-wave coefficient
above before source activation can be claimed. Independent microstructural
information is required; fitting ideal strength, line span or mobility to
the desired yield value is not a substitute. A minimal analytic environment
extension requires the recorded compatibility failure plus a new gauge/rank
audit, not an arbitrary high-order fit.

Run `run_yield_elastic_audit` (actual fits), `run_yield_benchmark` (cached
source pixels), `run_finite_source_reference` (actual equilibrium/refinement),
then `report_yield_bridge` (new material stationary/FD checks plus plots of
those completed runs). The initial `material_metric` record predates the
active-face QP; final calibration is `material_metric_refined`.

Physical collective-coordinate mobilities, t0, seconds and Hz remain
unavailable. No new PDE, first-passage probability, cyclic Al response, source
kinetics or mesh/UI capability is claimed. Read `verification.json` for
actually executed test counts/timings, not expected results.
