# Material compatibility before specimen strength — material_strength_v10

## Scope and present gate

This stage answers the request to move toward **actual material strength**. It
does not call the rigid-interface ideal fold an experimental yield prediction.
The existing LJ/Bessel family, per-atom environment sums, static parameter
files, production probability PDE and UI defaults are preserved. The new
research candidate is separate. There is no fatigue, mobility, specimen-area
or strength-factor fitting.

The material gate has two independent requirements: correct energetic physics
and the state space needed for the specimen mechanism. Changing coordinates
from a constrained straight registry path to a relaxed vector path helps the
first issue; it cannot by itself create a finite dislocation source.

## 1. Same analytic family, better specified observations

In the existing nondimensional geometry, with energy in eV, write

\[
 E=\frac12\sum_{i\ne j}(u r_{ij}^{-12}-v r_{ij}^{-6})
   +\sum_i\left[-A\sqrt{x_i}+B(x_i-1)+C(x_i-1)^2
          +D_1 I_{1i}+D_2 I_{2i}+D_3 I_{3i}\right],
 \quad x_i=\varrho_i/\varrho_{\rm ref}.
\]

Here \(u=4\epsilon_{LJ}\sigma_{LJ}^{12}\),
\(v=4\epsilon_{LJ}\sigma_{LJ}^{6}\). The pair functional form is LJ, not a
tabulated replacement. The per-site STF invariants \(I_l\) and their existing
radial normalization are unchanged. Each environment is summed before taking
its invariant or applying nonlinear embedding. No new energy term is added.
Isolated-atom embedding is \(F(0)=-B+C\), which must enter atomization energy.

Let \(d_\rho,d_Q\) be the dimensionless scalar/angular radial decays. At fixed
decays every energy, gradient and Hessian observation is **exactly linear** in
\(c=(u,v,A,B,C,D_3,D_1,D_2)^T\):

\[
 y_j={\cal B}_j(d_\rho,d_Q)c.
\]

`VectorCoefficientBasis` evaluates this coefficient identity from the actual
infinite reciprocal energy. It is not a surrogate interpolation, a fitted
gamma table or a finite-neighbor approximation. The common source/candidate
state is \(q=(a,u_x,u_y)\), not the old scalar straight-line constraint.

The reference conditions are the 0 K Mishin Al99 model on rigid half-crystals,
allowing only the active gap and two registry coordinates to move. The source
file checksum is
`60c8a085be79d273324ab421f5b1447578fef55c1acfc6492c0999f15ee8a284`.
See [Mishin et al., PRB 59, 3393 (1999)](https://doi.org/10.1103/PhysRevB.59.3393)
and the existing `ALUMINUM_CALIBRATION_SOURCES.md`. It is a reference potential,
not experimental truth, and is used only to generate/validate observations.

The declared observations are saved with role, units and normalization:

* exact: zero bulk normal force at 4.05 angstrom and 3.36 eV/atom cohesion;
* fitted: three independent bulk elastic combinations, perfect interface
  normal/registry stiffness, energy and normal force at the SOURCE fault and
  saddle, saddle path force, opening energies at 1.1h, 1.5h, 2h, 40h;
* held out: direct quarter/half path and oblique energy, opening at 1.25h and
  3h, fault/saddle normal and registry curvatures (nine observations).

The reverse barrier equals saddle minus fault energy. It is not an additional
independent fit target. Source stationary points are solved afresh with the
full Hessian/Morse index. Predictions evaluated at source points are not
silently treated as candidate stationary points: candidate roots are solved
independently during validation.

## 2. Units and elastic information

\[
 L_0=4.05/\sqrt2\ {\rm angstrom}=2.863782463805517\,10^{-10}{\rm m},
 \quad h/L_0=\sqrt{2/3},
 \quad A_{atomic\ cell}=\frac{\sqrt3}{2}L_0^2
          =7.102490842787125\,10^{-20}{\rm m^2}.
\]

\[
 \gamma_{surface}=W_{eV}eV_J/A_{atomic\ cell},\qquad
 T_i[Pa]=W_{q_i^*}eV_J/(A_{atomic\ cell}L_0).
\]

The atomic area converts interface energy, not a probability or specimen
force. Statistical correlation area \(A_c\) is absent from this calculation.
There is no modification of the production reduced stress map
\(f^*=\kappa_{axial}\sigma/E\).

For \(V=A_{atomic\ cell}h\), the independent strain curvatures are

\[
 H_{hydro}=3V(C_{11}+2C_{12}),\quad
 H_{normal}=V(C_{11}+2C_{12}+4C_{44})/3,\quad
 H_{shear}=V(C_{11}-C_{12}+C_{44})/3.
\]

The same fixed potential and length/density normalization are used during
strain. The rounded reference constants 114, 62, 32 GPa are 0 K model values.
The independent modes do not replace local GSF or cleavage validation.

## 3. Gauge fixing, optimization, identifiability

Density amplitude is absorbed into the reference density convention; it is
not fitted independently. `rho_ref` is fixed during displacement and strain.
The LJ/embedding split remains a declared model convention; the linear
embedding coefficient does not by itself constitute independently identified
many-body physics.

The first two linear constraints eliminate the LJ coefficients exactly:

\[
 c=c_0+Tz,\quad z=(A,B,C,D_3,D_1,D_2),\qquad
 c_{0,uv}={\cal B}_{0:2,uv}^{-1}y_{0:2},\quad
 T_{uv}=-{\cal B}_{0:2,uv}^{-1}{\cal B}_{0:2,z}.
\]

The coefficient subproblem is deterministic constrained least squares,

\[
 L=\sum_{j\in fit}[(\mathcal B_jc-y_j)/s_j]^2,
 \quad u,v>0,\quad A,C,D_3\ge0.
\]

Signed B, D1 and D2 follow the existing research sector. Five-percent bulk
and ten-percent interface energy/curvature scales, and a 0.25 GPa force scale,
are **declared model-discrepancy allowances, not experimental uncertainty**.
No weights are selected from a desired strength. Held-out rows never enter L.

The executed radial profile has 25 predeclared grid points plus the previous
decays and 160 Powell evaluations. Powell exhausted its evaluation budget;
its best recorded point is **not claimed to be a converged global optimum**.
The coefficient subproblems at selected candidates converged. Stage ablations
at fixed best decays separate bulk, registry and opening constraints; they are
not mislabeled as independent radial optimizations. Bulk-only zero columns
D1/D3 are fixed to zero, not divided by zero during column normalization.

After removing the two exact constraints, the local coefficient Jacobian has
six columns. Adding log-decay derivatives gives eight. Both are scaled and
saved with SVD/right vectors/column cosines. These are local sensitivities,
not confidence intervals. Finite-difference radial sensitivity is repeated
at steps 4e-4 and 2e-4; numerical full rank is not global uniqueness.

## 4. Continuous-path check and constraint exchange

A candidate with positive force at 11 opening sample points still had
\(T_n=-46.4765\) MPa at \(a/h=2.7963372\). This is far above the computed
numerical error and cannot be discarded as roundoff. It was found by solving
\(W_{aa}=0\) between independent curvature sign brackets (91/181 samples).

The separate `opening_exchange` stage keeps the selected radial decays and
adds the actual violating minimum as a new **linear coefficient constraint**:

\[
 \partial_a\mathcal B(a_{min},0,0)c\ge0.
\]

It refits coefficients, then locates the changed minimum again. The stopping
comparison uses actual 2e-11 to 2e-13 tolerance changes and a floating-point
dot-product error envelope, not an arbitrary MPa floor. Independent 121/241
brackets are rechecked. No force or energy is clipped. This resolves only the
tested zero-registry opening interval through 12h, not all registry directions,
all distant secondary extrema, or a rigorous global monotonicity theorem.

**Important source/prior audit:** the SOURCE itself has an opening traction
minimum of -785.6754MPa at a/h=2.386631. Two independent181/361 brackets,
energy finite differences and the separate scalar-source energy agree.
Therefore nonnegative traction (monotone opening **energy**) is a declared
shape prior used in this fit, **not a universal physical law or a property
proved for the Al reference**. A negative lobe cannot by itself establish a
coding defect. The unrelaxed source's oscillatory separation curve also is
not experimental cleavage truth. Its derivatives are retained, never clipped.
No claim about interpolation/cutoff versus inherent SOURCE-potential behavior
is made without additional validation of that distinction.

Removing this shape prior at the same decays gives L=28.4043 and
C11/C12/C44=85.16/64.07/52.70GPa. Thus the elastic mismatch is not explained by
that restriction alone. `opening_shape_prior_ablation.json` preserves this
additional coefficient fit without pretending its complete stability/shape
has been certified. A future calibration should decide this prior from
independent relaxed decohesion evidence, not enforce it to look successful.

Bulk finite-q Hessians are checked independently on corrected cubic GX/GL/GK
paths with real-space radii 8 and 12L0. These are static diagnostic sums, not a
replacement for canonical reciprocal energy, and not calibrated frequencies.

## 5. Executed material results and rejection

The exchange-corrected coefficient vector is

    decays = (2.906199643752374, 5.308553273298447)
    (u,v,A,B,C,D3,D1,D2) =
    (.20739749079909, 1.10528334351375, 8.23440762764772,
     12.7427207949929, 1.13853483400765, 89.5058742267111,
     5.00455384545968, 85.3789543582300)

This is a research fit, **not a newly adopted Al potential**.

| Observable | Source | Previous candidate | Exchange-corrected candidate |
|---|---:|---:|---:|
| C11 [GPa] | 114 | 87.8162 | 85.2587 |
| C12 [GPa] | 62 | 64.8326 | 64.2424 |
| C44 [GPa] | 32 | 49.2710 | 52.5428 |
| Connected forward barrier [J/m2] | .1720024 | .1732380 | .1668228 |
| Relaxed intrinsic fault [J/m2] | .1504795 | .1268002 | .1519913 |
| Reverse fault barrier [J/m2] | .0215229 | .0464378 | .0148315 |

The selected 14-point fit squared loss drops 87.4386 -> 28.6052, but the nine
held-out normalized RMS errors are still 15.6863 -> 9.7749. The largest new
held-out discrepancy is saddle Hxx at the SOURCE state: -0.18642 versus
+0.33977 eV/L0^2. This component is not a Hessian eigenvalue; the candidate's
own saddle still has Morse index one. Normalized errors can be large for a
small reference component; the dimensional errors are preserved in the CSV.

The initial joint fit local SVD (after two exact constraints, including two
decays) is 328.949, 103.687, 77.947, 19.850, 15.372, 4.382, 2.474, 1.592;
condition number206.60. Coefficient-only condition number141.96. A/C sensitivity
cosine is .9913 and A/B is -.9460. This is not evidence of uniquely identified
atomic parameters, especially with active constraints and model discrepancy.

The corrected negative traction is -5.72e-15 eV/L0 versus a computed local
3.77e-14 eV/L0 resolution envelope; the refined bracket gives -7.05e-15.
Direct-sum energy errors at radius/layers24/48/72 decrease
4.01e-5 -> 5.24e-6 -> 1.58e-6 eV/cell (worst of two states). The tight reciprocal
check changes H by at most6.39e-14 eV/L0^2. Fine finite-difference gradient/H
errors are at most7.24e-9 eV/L0 and1.18e-7 eV/L0^2, below the material errors.
Finite-q samples stay positive (minimum .0097165 eV/L0^2); increasing the
independent validation cutoff8->12 changes a sampled minimum branch by at
most .0011247 eV/L0^2. This is not an all-wavevector theorem.

The candidate's local all-free uniform shear fold is 2.58174GPa, versus the
old3.12402GPa, unchanged by independent31/61 brackets within2.5e-11MPa.
These numbers are **not improved predictions of measured yield strength**.
At 4.595,6.929,9.881MPa the corrected local static branch has ux/L0 roughly
.0001214,.0001831,.0002611; it returns after static unload to within8.72e-16L0.
The actual source potential also has only a reversible local response under
this rigid-perfect-interface calculation. This comparison points to the
missing finite-defect state space, not a conversion factor in stress units.

Decision: reject promotion. Bulk elasticity and held-out shape remain wrong
despite improvements in particular fault energies and opening monotonicity.
This bounded study does not prove that every parameterization of the whole
analytic family is impossible. The discrepancy normalizations also need
transparent sensitivity analysis; scalar fit loss alone is not adoption.

## 6. Experimental strength benchmark: what is and is not measured

J. Krebs, S. I. Rao, S. Verheyden, C. Miko, R. Goodall, W. A. Curtin,
A. Mortensen (2017), *Cast aluminium single crystals cross the threshold from
bulk to size-dependent stochastic plasticity*, Nature Materials 16, 730–736,
[doi:10.1038/nmat4911](https://doi.org/10.1038/nmat4911).
The [accepted manuscript](https://eprints.whiterose.ac.uk/id/eprint/126662/1/Article%20Text.pdf)
provides a room-temperature, 99.99% as-cast Al single-crystal wire tensile
reference at 300 nm/s displacement rate. Figure 2b supplies the following
**flow points**, digitized from the dark-blue 103 micrometre wire:

| Plastic shear strain | Resolved shear stress [MPa] |
|---:|---:|
| 0.10 | 4.5952 |
| 0.20 | unavailable: overlapping cyan trace |
| 0.50 | 6.9286 |
| 0.80 | 9.8810 |

These are not axial proof stresses. Pixel reading bands are +/-0.00357 strain
and +/-0.1905 MPa, not specimen scatter. The figure cannot resolve the 0.002
plastic-shear criterion, so that CRSS remains null. Exact orientation and
source geometry require further specimen data. Diameter is not source length.
No literature strength datum enters the atomic fit. PDFs remain ignored cache;
only short factual provenance and these numerical readings are saved.

`StrengthConditions` prohibits an ideal fold, resolved shear/axial mismatch,
different strain criterion, or missing specimen/protocol binding from being
reported as a matched prediction error. It does not establish physical truth
merely because two metadata strings match. `stress_at_plastic_strain` requires
a resolved monotonic loading bracket; first microscopic nonzero slip is not
the finite-strain strength criterion. The 0 K rigid-half probes at the same
MPa values are an explanatory comparison, not a reproduction of those wires.

## 7. What must produce lower specimen strength without breaking the theory

For the uniform interface, a full Hessian fold under applied traction is an
ideal homogeneous instability. A real finite source has nonuniform registry,
normal relaxation, free surfaces/obstacles and a specified defect history.
The thermodynamic work acts on an accumulated slipped region; it does not
require every atom of an infinite plane to cross the barrier simultaneously.

The safe atomistic starting point remains the *same* energy difference,

\[
 \Delta E[\{U_i\}]=\frac12\sum_{i\ne j}
 [\phi(|R_{ij}+U_j-U_i|)-\phi(|R_{ij}|)]
 +\sum_i[F(\rho_i[U])-F(\rho_i[0])]
 +\sum_{i,l}D_l[I_{li}[U]-I_{li}[0]].
\]

The perfect infinite background can still be summed with Poisson/Bessel;
defect perturbations require a controlled convergent difference, not a
finite-neighbor replacement of LJ. Nonlinear embedding stays per site.
For a finite slipped area, shear work is \(-\int\tau\,u\,dA\). Its spatial
extent must come from a resolved source configuration, not an invented
activation area or the specimen probability correlation area.

In a justified long-wave limit the Peach–Koehler driving force per line length
is \(\tau b\), whereas curved-line restoring forces depend on core energy and
curvature. The often used dimensional comparison \(\tau b\sim\Gamma/R\)
explains why source size matters; **it is not inserted here as an empirical
yield law**, nor used with a guessed R to report a strength. An infinite
straight-row energy in J/m is not a finite activation energy in eV. Finite-line
variation and source bow-out are absent from the current row-invariant code.

Next requirements: an acceptable common bulk/registry/opening surface;
validated vector/partial core and finite source with size convergence;
independently specified microstructure and loading protocol; then matched
net transport/unload behavior. Kinetics additionally requires physical
collective-coordinate mobility. Static energy fitting cannot supply seconds,
Hz, fatigue life, or spatial independence.

## Reproduction and status ledger

Run `solver_v1.run_vector_material_calibration fit --local-evaluations 160`,
then its `validate` stage. Run `solver_v1.run_material_strength_rechecks repair`,
validate with `--output .../material_strength_v10/opening_exchange`, then
`solver_v1.run_material_strength_rechecks verify`.
Optional `fetch_strength_reference article --pages 5 --images` requires pypdf
only in the ignored research tool cache; `run_strength_benchmark` records
selected pixels and hashes. No PDF package is added to the application.

Read `results/fcc111_active_interface/material_strength_v10/SCHEMA.md`, the
machine-readable residuals and `verification.json` for executed results.
The first report attempt failed at a zero sensitivity column and is preserved
as incomplete; it is not a successful calibration. A later captured SOURCE
NaN involved an exactly empty neighbor set. The source-only evaluator now
returns the exact zero jet before that empty reduction; nonfinite values for
nonempty sums remain errors. This is not clipping or a LJ cutoff change.
