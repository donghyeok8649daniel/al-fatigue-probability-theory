# From ideal interface strength toward defect-controlled strength: nonlinear rows

Status: **static scalar anti-plane research only**. This module is not a new
production PDE, a measured Al yield law, or a kinetic calibration. It implements
the nonlinear next-stage contract in `DISCRETE_FCC_SCREW_DERIVATION.md`, using
the unchanged `angular_monotone_opening` research candidate. The LJ base,
infinite-row Poisson/Bessel identities, scalar embedding, and fitted angular
coefficients are unchanged. This is a state-space extension, not a strength fit.

The user's question is how to progress from ideal crystal strength to the
strength of real defect-containing metal. Artificially dividing a GPa threshold
by a number would not answer that question. We instead introduce spatially
different row registries and an explicitly declared existing screw pair, solve
the actual nonlinear energy, apply low-MPa stress, and test what is still missing.

## 1. Exact geometry and deliberate restriction

In the verified +ABC FCC construction,

\[
 R_{njl}=(nb+\delta_{jl},y_{jl},z_l),\quad
 \delta_{jl}=b(j+l)/2,\quad y_{jl}=d(j+l/3),\quad z_l=hl,
 \quad d=\sqrt3 b/2,\quad h/b=\sqrt{2/3}.
\]

Each cross-section site represents a complete infinite atomic row parallel to
e1, not one isolated atom with a finite-neighbor cutoff. Let i=(j,l), and take
the physical x displacement divided by L0 to be

\[
 U_i=u_i+\gamma z_i,\qquad
 \eta_{iR}=\delta_R+u_{i+R}-u_i+\gamma z_R.
\]

The field u is periodic on a logical j,l parallelogram. The affine shear gamma
is an additional mechanical variable. Transverse coordinates y,z are **fixed**
in the solver. This permits a nonlinear screw core in a restricted subspace,
not normal accommodation, a two-component partial core, a curved finite line,
a Frank–Read source, or a finite nucleation event.

Changing 24x24 to 32x32 or 48x48 changes the number of actual atomic rows and
periodic image spacing. It is not a continuum mesh refinement: no fictitious
points are inserted between atoms. The same actual core separation must be
checked independently of the seed's nominal separation parameter.

## 2. Nonlinear infinite row coefficients

For p=3,6 the existing LJ power sum is

\[
 S_p(A,\eta)=C_{p0}(A)+\sum_{m\ge1} C_{pm}(A)\cos(g_m\eta),
 \quad g_m=2\pi m/b,
\]
\[
 C_{pm}=\frac{4\sqrt\pi}{b\Gamma(p)}
  \left(\frac{\pi m}{bA}\right)^{p-1/2}K_{p-1/2}(g_m A),
\quad
 C_{p0}=\frac{\sqrt\pi\Gamma(p-1/2)}{b\Gamma(p)}A^{1-2p}.
\]

The LJ channel is exactly 4 epsilon [sigma^12 S6 - sigma^6 S3]. No pair fit or
empirical pinning function is introduced. The scalar density uses

\[
 H(A,g)=\frac{2A\kappa}{q}K_1(Aq),\qquad q=\sqrt{\kappa^2+g^2},
 \qquad f(r)=C e^{-\kappa r}.
\]

Its paired positive-mode coefficient is (2C/b)H. For a raw angular component
with x,y,z powers p,v,w, its coefficient is

\[
 \frac{2C_{\rm ang}}b i^p\partial_g^p H(A,g)\,y^v z^w.
\]

An orthonormal symmetric trace-free basis compresses the full tensors to
3,5,7 components without changing their squared norms. The original angular
decay, amplitude, and rank-dependent energy coefficients are preserved. The
scalar density kernel and angular kernel are not silently identified.

Include the exact ABC phase exp(i g_m delta_R)=(-1)^{m(j+l)} in a complex
coefficient c_(m,R,c), with c indexing pair, scalar density and the 15 angular
components. Every site's change is

\[
 \Delta t_c(i)=\operatorname{Re}\sum_{m,R}c_{mRc}
   [e^{i g_m(u_{i+R}-u_i+\gamma z_R)}-1].
\]

The reciprocal zero mode cancels **only because y,z are fixed** and this is an
x-translation difference. It must be restored for transverse forces below.
Same-row atom distances do not change in this state space and cancel exactly;
the self atom is not reintroduced.

## 3. Correct nonlinear per-site many-body energy

Let rho_b be the full infinite FCC density of the unchanged cubic reference.
For one repeat b L0 of each infinite line,

\[
 \mathcal E_b[u,\gamma]=\frac12\sum_i\Delta t_{\rm pair}(i)
 +\sum_i\{F[\rho_b+\Delta t_\rho(i)]-F(\rho_b)\}
 +\sum_{i,r=1}^3D_r\|\Delta Q_r(i)\|^2.
\]

Q_(r,b)=0 for this cubic reference. A noncubic reference is rejected, not
handled by incorrectly dropping its Q2 background. Each site's **complete**
environment is summed before evaluating F or a tensor norm. Scalar F' is not
frozen at rho_b. F'' contributes to the full nonlinear Hessian even though its
perfect harmonic anti-plane contribution vanishes by symmetry. The candidate's
negative D1 is not clipped or changed to produce stability.

The half factor occurs only in pair counting. Embedding is per site. The rigid
two-cut test independently recovers the existing full-plane interface energy;
no separate gamma energy is added to a kernel already containing local energy.

## 4. FFT evaluation is an algebraic convolution, not a new approximation

For z_m(i)=exp(i g_m u_i), define the periodic convolution

\[
 (A_m(\gamma)z)_i=\sum_R c_{mR}e^{i g_m\gamma z_R}z_{i+R}.
\]

Then the current channel is Re sum_m [conj(z_m) A_m z_m] minus its perfect
reference. FFTs implement this retained convolution exactly. The slow independent
real-index implementation and large real-space atomic rows check it separately.

Let w_ic=partial E/partial t_ic. The channel weights are

\[
 w_{\rm pair}=1/2,\quad w_\rho=F'(\rho_i),\quad
 w_{r}=2D_r\Delta Q_{ri}.
\]

The analytic field gradient is

\[
 \mathcal E_{u_i}=\operatorname{Re}\sum_{m,c}i g_m
 \{z_i[A_m^T(w\bar z)]_i-w_i\bar z_i(A_mz)_i\}.
\]

T denotes the ordinary complex transpose, not an accidentally substituted
Hermitian transpose. An analytic directional differentiation evaluates the
Hessian-vector product. It includes dw_r=2D_r dQ_r and dw_rho=F'' drho.
The affine derivatives follow by multiplying each row coefficient by
i g_m z_R, and then -(g_m z_R)^2. Finite differences are independent tests,
not the solver's force implementation.

At u=gamma=0 this Hessian recovers the separately verified harmonic FCC
infinite-row symbol. A uniform u translation is the null gauge. The local
stability calculation excludes that mode; if an artificial gauge shift is
returned below every physical eigenvalue, it explicitly fetches the nongauge
eigenpair rather than reporting the shift as material stiffness.

## 5. Stress control and dimensional normalization

This research calculation has explicit engineering shear gamma=dx/dz. One
atomic cell has volume b*d*h*L0^3, so the actual periodic cell volume is

\[
 V_{\rm cell}=N b d h L_0^3.
\]

This is geometric counting of represented atoms, not an arbitrary activation
volume or a fitted length. At applied resolved shear tau, minimize

\[
 \mathcal G_b=\mathcal E_b-\tau V_{\rm cell}\gamma/\mathrm{eV_J},
 \qquad \tau_{\rm internal}
   =\frac{\mathrm{eV_J}}{V_{\rm cell}}\partial_\gamma\mathcal E_b.
\]

With physical tau in Pa, `tau_internal_units=tau*L0^3/eV_J` in eV/L0^3.
For the saved L0, 1 MPa corresponds to 0.00014659180163330851 in these units.
The perfect state's computed small-shear modulus is about 24.085 GPa for this
orientation and candidate. It is not replaced with the literature target.

**gamma=0 is strain control, not zero applied stress for a defect-containing
cell.** At zero applied stress a pre-existing dipole has nonzero total affine
slip content. The initial stress-free state is the reference for any subsequently
generated registry change. Neither the initial gamma nor its initial integer
content is a newly produced plastic strain.

The minimized E is eV per b L0 of infinite straight line. Convert to line energy
by E*eV_J/(b L0) in J/m. It is not an energy barrier for a finite thermal event.
No A_c, fitted activation volume, atomic mass or mobility enters this static
calculation. The production kappa*sigma/E convention is unchanged.

## 6. Integer relabeling, winding and shear decomposition

For every row, u_i -> u_i+b k_i with integer k_i relabels identical infinite
atomic positions. A field composed only of integer b steps is not a core.
Using phase theta_i=2 pi u_i/b, principal bond phase differences define the
plaquette circulation and Burgers winding

\[
 n_\square=\frac1{2\pi}\sum_{\partial\square}
           \operatorname{Arg}e^{i(\theta_{\rm next}-\theta_{\rm current})}.
\]

We store signed core locations, integer residual and distance from half-period
bond ambiguity. This is a discrete topological diagnostic; it does not alone
certify continuum Burgers stress fields or a dissociated physical Al core.

For layer bonds v_i=u_(i+e_l)-u_i+gamma h, set

\[
 k_i=\lfloor v_i/b+1/2\rfloor,\qquad \xi_i=v_i-bk_i.
\]

By periodic telescoping,

\[
 \gamma=\frac{b}{h}\overline{k}+\frac1h\overline{\xi}
        =\gamma_{\rm registry}+\gamma_{\rm intrabond}.
\]

The mean registry component is invariant under integer row relabeling even
though individual k_i change. This is a separate spatial shear decomposition,
**not** the production axial chi survivor bridge or a newly derived probability
flux. We track changes relative to the initial stress-free defect state.

## 7. Static solver and what unload means here

All independent row values and, under stress control, gamma are minimized. A
uniform row translation is removed. L-BFGS is followed when needed by analytic
Newton-CG with descent checks. `optimizer_success` and actual maximum atomic
force/scaled stress residual are separate fields. The tolerance is a numerical
stationarity criterion, not a plasticity cutoff.

The protocol is imposed shear 0,4,15,25,50,25,0,-25,-50,-25,0 MPa. Previous
equilibria initialize subsequent solves. This is athermal quasistatic local
branch following: optimization iterations are not time, temperature or mobility.
Any branch loss/jump must be checked separately; successful minimization is not
a global-minimum or activation-path proof.

The final zero-stress solve is static unload, **not a finite-temperature
zero-stress hold**. There is no production P(a,s,t), SG interwell traffic,
opening probability, physical frequency or fatigue lifetime in this study.

## 8. Crucial check: forces excluded by the scalar subspace

Solving E_u=0 at fixed y,z cannot establish full atomic equilibrium. The separate
`nonlinear_screw_transverse_audit.py` evaluates gradients in these omitted
directions using the same infinite-row potential. It does not implement their
relaxation or declare a vector core complete.

For H_j=partial_(g^2)^j H,

\[
 H_j=2A\kappa(-A/2)^j q^{-j-1}K_{j+1}(Aq),\quad
 \partial_A H_j=-2A\kappa(-A/2)^j q^{-j}K_j(Aq).
\]

Together with y,z polynomial derivatives and dA/dy=y/A, dA/dz=z/A, this gives
analytic transverse derivatives of scalar and tensor row sums. The LJ radial
derivatives are reused from the verified power-kernel implementation.

Here the **m=0 term must remain**. Its scalar-density gradient can couple to
different F'(rho_i) near a defect even though its slip derivative is zero. The
audit therefore includes m=0 with its single-mode, not doubled, Fourier factor.
It differentiates at fixed physical x; the reference affine x phase is not
accidentally differentiated as an extra z displacement.

If C_(m,alpha) is a row-channel derivative and A_(m,alpha) its convolution,

\[
 E_{v_{i\alpha}}=\operatorname{Re}\sum_{m\ge0,c}
 \{z_i[A_{m,\alpha}^T(w\bar z)]_i
        -w_i\bar z_i(A_{m,\alpha}z)_i\},\quad \alpha=y,z.
\]

A separate rigid-plane calculation verifies the normal force including nonlinear
embedding. Zero-force perfect bulk and total force cancellation are also tested.
The actual core-force audit extends transverse rings 6,10,14,20. It does not
use the anti-plane harmonic tail estimate for a zero-mode-containing derivative.
Report eV/L0 (or the corresponding N per atomic line repeat), not an invented
local GPa obtained by choosing a convenient volume.

## 9. Convergence and recorded results

Raw files are in `results/fcc111_active_interface/nonlinear_screw_v7/`:

- `metadata.json`: unchanged parameter hash, model status, units and mode/ring
  controls; independent FFT and derivative errors;
- `stress_history.csv`, `static_states.csv.gz`: every actual static solve and
  selected full atomic states, not synthetic extrapolated histories;
- `stability.csv`: smallest nongauge scalar/affine Hessian eigenvalues;
- `tail_refinement.csv`, `relaxed_core_tail_refinement.csv`: actual spatial-row
  and reciprocal-mode refinements, distinct from periodic domain effects;
- `domain_and_unload_summary.csv`: matched discrete core separation and
  changes from each initial stress-free defect state;
- `transverse_force_audit.csv`: omitted force, ring changes and net residual;
- `execution_summary.json`, `physical_and_numerical_status.json`: executed
  numerical results versus physical acceptance.

The anti-plane default uses 6 rings /168 neighboring infinite rows /11 positive
reciprocal modes, checked against extra rings. Bounds are accumulated per
channel with powers 1,g,g^2; they are not a dimensionally mixed eV error bar.
The unobserved infinite remainder is assessed by refinement, not falsely called
a rigorous bound. Exponentially convergent phase-dependent rows are distinct
from a finite physical neighbor cutoff.

### Actual executed results

The 44 static solves completed in 656.08 wall seconds. Every state passed the
independent scalar/affine force test. Maximum field-force residual was
2.40159e-11 eV/L0 and the shear decomposition residual was 2.50722e-18.
There was **no new registry shear and no change in signed core locations**
over the stated +/-50 MPa protocol in any of the three domains.

The actual opposite-core separation was 1.98408669 nm in all domains; the
7.3b arctangent seed parameter is not substituted for this measured distance.

| Atomic rows | Zero-stress E [eV/b-repeat] | Initial registry shear | Change in affine shear at +50 MPa | Final new registry shear |
|---|---:|---:|---:|---:|
| 24x24 | 2.22178949 | .0170103454 | .00208300997 | 0 |
| 32x32 | 2.26973530 | .00956831931 | .00208024362 | 0 |
| 48x48 | 2.30313892 | .00425258636 | .00207803167 | 0 |

Initial mean slip scales with the declared dipole content divided by the
periodic cross-sectional area. It is not a spurious newly generated plastic
strain. The final affine difference from each initial state was at most
2.12538e-13. The perfect 24x24 reference also returned to gamma=0 to numerical
precision and had no winding or integer registry change.

The remaining 32->48 energy change is .03340363 eV/b-repeat (about 1.45% of
the latter energy). Thus the energy of an **isolated** dipole is not yet
domain-converged. The incremental +50 MPa response changes by 2.21194e-6
(about .106% of the response). The absence of a hop is robust in these runs;
it does not establish an isolated dislocation Peierls threshold.

At 24x24, the smallest scalar+affine nongauge Hessian eigenvalues for
0,+50,-50 MPa were .25271401,.25283546,.25257973 in the explicitly scaled
optimization coordinates. Eigen-residuals were 1.21e-8,9.31e-9,3.19e-8.
These are static curvatures, not frequencies and not the full vector Hessian.

6->8 transverse row refinement changed relaxed-state energy by at most
1.60e-14 eV, x-gradient by 2.71e-15 eV/L0 and stress by 4.22e-13 MPa. The
FFT/real-index channel difference was 2.22e-16; a nonlinear directional
gradient finite-difference check differed by 1.24e-9 and its Hessian-vector
check by 7.66e-10. These numerical checks do not erase periodic image effects.

The crucial negative result is the **omitted transverse force**:

| Domain | Max omitted norm, zero stress [eV/L0] | At +50 MPa [eV/L0] |
|---|---:|---:|
| 24x24 | 2.2298290 | 2.2146541 |
| 32x32 | 2.2610920 | 2.2448068 |
| 48x48 | 2.2833823 | 2.2661844 |

The maximum is about 1.28 nN per line repeat. The 14->20-ring gradient change
was at most 1.69867e-10 eV/L0. A direct real-index nonlinear energy
re-evaluation after actual y,z perturbations independently checked the force
formula to 9.03e-11 at the final finite-difference step. The omitted force is
therefore **not a roundoff or row-tail ambiguity**. It is concentrated near the
cores in `static_defect_audit.png`.

This rejects **full mechanical-core equilibrium** for the computed scalar
solution. It does not prove how much transverse relaxation will reduce a
barrier or guarantee the desired real-Al strength. It identifies a necessary
missing degree of freedom quantitatively, without tuning any material parameter.

## 10. What advances toward actual strength, and what still cannot be claimed

This establishes a nonlinear, spatially nonuniform energy with correct counting,
analytical forces, stress control and explicit existing-defect topology. It
removes the assumption that an entire perfect infinite interface must slip at
once. It does **not** guarantee that this restricted candidate generates the
correct Al core, glide stress, source strength, fatigue or residual plasticity.

The material fit itself is still limited: C11/C12/C44 remain approximately
87.816/64.833/49.271 GPa versus the earlier 114/62/32 reference targets. No
parameter has been altered in response to the defect result.

The next justified step is nonlinear transverse/vector relaxation with the
same per-site energy and tested zero-mode terms, followed by explicit partial
core/path validation and matched-boundary/atomistic comparison. Only then can
an unforced migration saddle or depinning threshold be interpreted. Finite
line curvature/source geometry and measured defect populations are separate
requirements for specimen strength. None may be replaced by a fitted stress
multiplier or an arbitrary conversion of J/m into eV.

In that extension A_(iR)=sqrt[(y_R+v_(i+R)-v_i)^2+(z_R+w_(i+R)-w_i)^2]
depends on the site. The present fixed-coefficient FFT convolution cannot be
reused as though the row geometry were still constant. Preserve the analytic
infinite-line sum but evaluate the actual state-dependent radii/moments, keep
zero modes and their long-range tails, and test vector forces/Hessians against
direct atomic geometry. A cheap frozen harmonic transverse correction alone
would not validate a large nonlinear core reconstruction.

Primary context (not numerical fitting targets in this task):

- Lu, Kioussis, Bulatov and Kaxiras, *Generalized stacking fault energy surfaces
  and dislocation properties of aluminum*, Phys. Rev. B 62,3099 (2000),
  [DOI](https://doi.org/10.1103/PhysRevB.62.3099),
  [author manuscript](https://arxiv.org/abs/cond-mat/9903440): core structure and
  its detailed energy representation affect computed Peierls stress. This
  motivates validating the missing core freedoms, not assigning their values.
- Xu, Xiong, Chen and McDowell, *An analysis of key characteristics of the
  Frank–Read source process in FCC metals*, JMPS 96,460–476 (2016),
  [DOI](https://doi.org/10.1016/j.jmps.2016.08.002): finite source geometry,
  core behavior and image forces matter. The present infinite straight-line
  study is not a Frank–Read source simulation.

Reproduction:

```text
python -m pytest solver_v1/test_nonlinear_fcc_screw.py -q
python -m solver_v1.run_nonlinear_screw_study --sizes 24 32 48
python -m solver_v1.report_nonlinear_screw_study
```

All generated plots are actual static energies/fields. No global crack
probability is painted on a mesh. UI/production promotion remains blocked by
material, vector-core, finite-event and kinetic validation gates.
