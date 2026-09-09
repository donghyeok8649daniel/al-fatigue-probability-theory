# Vector infinite-row FCC core: derivation and validation contract

Research only. Same saved LJ/exponential-density/embedding/angular candidate;
no fit, kinetic calibration, finite-loop claim or production registration.
This extends nonlinear_screw_v7 by allowing the two previously frozen directions.

## Geometry and energy

For logical site i=(j,l), let U_i=(u_i,v_i,w_i), constant along the infinite e1
atomic row. Its displacement includes an affine x shear gamma*z_i. For row
offset R=(j_R,l_R), the actual three components entering each infinite sum are

\[
 x_{iR}=b(j_R+l_R)/2+u_{i+R}-u_i+\gamma h l_R,\quad
 y_{iR}=d(j_R+l_R/3)+v_{i+R}-v_i,\quad
 z_{iR}=h l_R+w_{i+R}-w_i,
 \quad r_n^2=(nb+x_{iR})^2+y_{iR}^2+z_{iR}^2.
\]

All three local displacements relax. The transverse periodic cell vectors
remain fixed bulk vectors; this is not a zero-normal-stress cell-shape ensemble.
The shear gamma can be fixed or stress-controlled. Atomic rows remain straight
and infinite: this is a three-component cross-section, not a curved 3D source.

The site channels T contain the LJ row energy, scalar environment density and
orthonormal STF moments of ranks 1,2,3. Sum ALL retained neighbor-row channel
changes first, adding the exact infinite cubic rho_bulk and Q_bulk=0 reference.
Then use

\[
 E_b=\tfrac12\sum_i\Delta T_{\phi i}
 +\sum_i\{F(\rho_b+\Delta T_{\rho i})-F(\rho_b)\}
 +\sum_{i,r} D_r\|\Delta Q_{ri}\|^2.
\]

There is no global-crystal embedding or extra local GSF energy. Unchanged
same-row distances cancel; self atoms are excluded. Cubic background moments
are not assumed for a deformed site: actual Q_i is recomputed there.

## Infinite row sum and analytic vector derivatives

Set A=sqrt(y^2+z^2). Use the existing Bessel row sum C_p0+sum C_pm cos(g_m x),
including m=0. For an exponential kernel H=2 A kappa K1(Aq)/q,
q=sqrt(kappa^2+g^2), define

\[
 H_j=\partial_{g^2}^jH=2A\kappa(-A/2)^j q^{-j-1}K_{j+1}(Aq),
\]
\[
 H_{j,A}=-2A\kappa(-A/2)^j q^{-j}K_j(Aq),\qquad
 H_{j,AA}=H_{j,A}/A+2A\kappa(-A/2)^j q^{1-j}K_{j-1}(Aq).
\]

K_-1=K1. Ordinary g derivatives follow from
H_g=2gH1, H_gg=2H1+4g^2H2, H_ggg=12gH2+8g^3H3.
Angular raw components are i^p H_(g^p) y^v z^w. Radial and polynomial product
rules yield their full vector gradients/Hessians before the STF projection.

For f=h(A)P(y,z), t=(y,z)/A,

\[
 \nabla_\perp f=h' t P+h\nabla P,\quad
 \nabla^2_\perp f=h''tt^T P+\frac{h'}A(I-tt^T)P
 +h'(t\nabla P^T+\nabla P t^T)+h\nabla^2P.
\]

The x derivative multiplies each reciprocal coefficient by ig, the xx
derivative by -g^2, and mixed x/transverse derivatives by ig. The m=0 density
and moment terms survive transverse differentiation and must not be deleted.
Actual radii depend on i, so the previous fixed-coefficient FFT is not reused.

For LJ p=3,6, half-integer K has an EXACT exponential-polynomial expression;
using it is evaluation of the same Bessel function, not a finite-neighbor model.
Its radial mode coefficient satisfies C''=g^2 C-2p C'/A. Tests must compare
these identities to the unchanged general-K implementation.

## Pair symmetry, forces and Hessian

Only one of R,-R need be evaluated because T_c(-r)=p_c T_c(r), with p=+1 for
pair/density/even-rank moments and p=-1 for odd moments. Each evaluated bond
adds T to site i and p*T to site j. For site energy weights
z_phi=1/2, z_rho=F'(rho), z_Q=2D Q, define b_c=z_ic+p_c z_jc. Then

\[
 g_{ij}=\sum_c b_c\nabla T_c(r_{ij}),\qquad
 E_{U_i}\mathrel{-}=g_{ij},\quad E_{U_j}\mathrel{+}=g_{ij}.
\]

For a direction delta U, differentiate both b and grad T, including F'' and
all angular cross-site terms. This is an analytic Hessian-vector product, not
a frozen harmonic approximation or a finite-difference production Hessian.
Uniform translations in all three components are gauges, not stiffness modes.

## Stress and scientific acceptance

Use the actual geometric V_cell=N*b*d*h*L0^3 and
G_b=E_b-tau V_cell gamma/eV_J. This preserves the previous stress-control
normalization; no characteristic activation volume or A_c enters it.
E_b is eV per infinite-line repeat b L0, not a finite activation barrier.
Optimizer iterations and static unloading are not physical time/hold tests.

Required checks: independent direct atomic row values/derivatives, analytic
vector Hessian, scalar energy/force recovery, omitted-force recovery, global
force balance, vector local stability, mode/row-tail/domain convergence,
actual loaded/unloaded state comparisons. Opening, slip, spontaneous pair
annihilation and failed continuation must not be conflated.

### Independent perfect-bulk Fourier Hessian

At the perfect state set theta_R=q_j j_R+q_l l_R. For each channel c of parity
p_c, its linear response to a vector Fourier displacement is

\[
 L_c(q)=\sum_{R>0}\nabla T_c(R)
 \{e^{i\theta_R}-1-p_c(e^{-i\theta_R}-1)\}.
\]

Since Q_bulk=0 in the cubic reference,

\[
 H(q)=2\sum_{R>0}(1-\cos\theta_R)
 [\nabla^2T_\phi+2F'(\rho_b)\nabla^2T_\rho]
 +F''(\rho_b)L_\rho^\dagger L_\rho
 +2\sum_cD_c L_c^\dagger L_c.
\]

The last sum covers angular channels only. This3x3 assembly, independent of
the real-index Hessian assembly, agrees with a finite periodic-cell
Hessian-vector product to2.84e-14. It is a perfect-bulk reference, not a
substitute for the nonlinear core Hessian.

The24x24 lowest nongauge curvature changes from0.255137823 at ring6 to
0.257829013 at ring10 and0.257935086 at rings16,20,28. However, the minimizing
polarization switches from a y/z mode to x. The **full** matrix difference
ring16 versus28 is1.18e-5 and ring20 versus28 is3.81e-6, even though their
minimum eigenvalues coincide. Do not infer spectrum convergence just from
one eigenvalue. All are static curvatures, not physical frequencies.

## Truncation and the newly active zero mode

Transverse motion revives a long-range term that cancelled in v7. For each
power the zero-mode row coefficient is

\[
 C_{p0}(A)=c_p A^{1-2p},\qquad
 c_p=\sqrt\pi\,\Gamma(p-\tfrac12)/(b\Gamma(p)).
\]

Its largest radial Hessian magnitude is bounded by
\(2p(2p-1)c_p A^{-2p-1}\). Let R denote the retained square logical-offset
ring, alpha the smallest singular value of the j,l-to-y,z geometry matrix,
and D=2 max_i |U_perp,i-mean U_perp|. For k>R there are 8k offsets and their
deformed radius is at least alpha*k-D. Put c=alpha-D/(R+1)>0. A conservative
bound on the omitted LJ zero-mode Hessian sum is

\[
 B_R=32\epsilon\sum_{p=3,6}\sigma^{2p}
       2p(2p-1)c_p\,c^{-2p-1}\zeta(2p,R+1).
\]

Both attractive and repulsive terms are bounded by absolute magnitudes.
The resulting per-site force norm bound is D B_R. The **total cell energy**
error is bounded by N D^2 B_R/4: the first-order displacement term cancels
exactly when summing a periodic cell. The attractive remainder is therefore
O(R^-5), not exponentially small. These are bounds for the LJ zero mode
ONLY. They do not certify all density/angular/nonzero-mode channels.
`pair_zero_mode_tail_bounds` reports that restriction explicitly, and an
independent two-neighborhood test checks the bound.

Positive reciprocal modes use a per-row envelope for values and derivatives
and require two successive small terms. At actual large radii fewer modes are
needed. The safety ceiling raises an error instead of returning repaired
forces. This is empirical Bessel truncation control checked against tighter
series and direct atoms, not a proof that two terms alone bound every tail.

At a saved nonzero-transverse 32x32 checkpoint (max transverse displacement
0.078884 L0), actual energy/force changes versus ring16 were:

| ring | energy difference [eV/line repeat] | max force difference [eV/L0] |
|---:|---:|---:|
|4|0.00766374|0.00400383|
|6|0.00132058|0.000657118|
|8|0.000196158|0.000110746|
|10|0.0000250004|0.0000153228|
|12|0.00000190814|0.00000111542|

Thus it would be wrong to reuse v7's near-machine-precision transverse
omission claim after releasing y,z. Ring16 is a comparison reference, not an
infinite-tail certificate. Final nearly perfect states have much smaller
differences, but a nonzero core/fault needs its own tail study.

## Independent numerical identities

Twelve fixed points include both registry signs and transverse radii from
about 0.7 to 4.7 L0. Direct atom sums with +/-96,192,384 x images are used
only as validation. At +/-384, maximum absolute discrepancies over all
channels are:

- values: 1.78e-14;
- gradients: 1.42e-13;
- Hessians: 5.12e-13.

Channel units differ (pair eV, scalar density, rank-dependent moments), so
these are arithmetic checks, not one material error bar. Per-field scaled
and non-near-zero relative errors are saved in `direct_atomic_validation.csv`.
The reciprocal calculation used at most14 positive modes; largest last-mode
envelope8.28e-15 at requested2e-12. Additional tests verify all-vector energy
gradients, analytic Hessian products, stress work, translation gauges,
per-site direct atom assembly and exact scalar/omitted-force recovery.

## Actual vector relaxation: what changed physically

The same previously saved candidate and actual v7 scalar states are used,
without optimizing any material parameter. All three local displacements are
released, as well as affine shear under imposed resolved shear stress.
The scalar equilibrium had an excluded transverse norm about2.23eV/L0.
The 24x24 vector relaxation reduces it to a full maximum component residual
6.09e-8eV/L0 (ring6). The excess energy falls from2.22178949eV/line repeat
to about1.4e-13, its opposite x windings disappear, and the remaining local
vector displacement is near numerical zero. The smallest nongauge static
curvature is0.25513782, not a frequency. Ring10 and the32x32 domain are
checked separately; the machine-readable records contain their final values.

This is an actual release of a missing mechanical degree of freedom, not a
change of mobility, LJ coefficients, stress scaling, hardening, or a plastic
cutoff. The supplied close opposite pair can remove its internal elastic
energy at zero external shear. That does NOT mean Al has zero yield stress:
it is a particular existing defect pair's static relaxation, not a creation,
source-activation, measured dislocation population or kinetic experiment.

The low-MPa continuation is performed **after** this zero-stress relaxation.
The primary24x24 sequence is0,4,25,50,0,-50,0MPa. At+50MPa the affine shear
increment is about0.002076; no new x winding is generated. Static unload
recovers the nearly perfect state within the recorded errors. There is no
physical-time zero-stress hold, hysteretic dissipation rate or fatigue-life
prediction here. Optimizer iteration plots are explicitly labeled NOT time.

### Counterexample: x winding is not vector topology

A separate8x8 close-seed test initially failed a proposed assertion that
zero x winding implies a perfect crystal. The assertion was scientifically
wrong, not a failure to adjust the force tolerance. Actual three-component
relaxation at ring4 yields E=0.575047975eV/line repeat, full force residual
7.22e-8 and positive smallest curvature1.33224. There is no x winding, yet
max transverse displacement is0.129675L0. Ring10 gives E=0.585289091eV,
residual8.17e-8 and curvature1.39727. The nonzero state must not be erased.

Inspection shows an almost j-independent **cell-spanning registry fault**,
not a surviving localized screw pair. Across its active layer the additional
in-plane shift is approximately(0.5,0.25579)L0, normal change0.01553L0;
other layer shifts compensate the fixed transverse cell shape. This is close
to extra tau=(b/2,d/3) modulo the triangular lattice, but not an independently
validated relaxed-Al partial core. The 16x16 test of the same seed protocol
already loses its pair during the initial scalar relaxation, so it is NOT
a matched-final-core domain comparison. It does demonstrate that the8x8
outcome cannot be promoted to an isolated-defect material result.

An actual ring16 re-relaxation preserves this small-cell fault:
E=0.585337628eV/line repeat, full force1.36e-9 and curvature1.39755637.
After0->50->0MPa the transverse state changes by at most1.51e-10L0 and the
global-phase-removed x coordinates by1.20e-10L0. Nevertheless the scalar
x-projected registry index can report a change of -0.07654655! At ring10
it reports +0.07654655 instead. This is the b/2 partition flipping, not a
large residual plastic deformation. The observed counterexample is saved,
not clipped or silently renamed as physical slip.

Most importantly, the fault's x shift lies essentially on b/2. The old scalar
formula floor(delta_x/b+1/2) is then ill-conditioned and is **not** a definition
of vector plastic slip. New `layer_registry` diagnostics report both components,
distance to 0,tau,2tau modulo the triangular primitive lattice, actual normal
shift, row-to-row dispersion and distance to that scalar partition. Every
scalar-registry column in raw case files is a kinematic x projection ONLY.
The consolidated table prefixes those columns with `x_projected_`.
No production s=bn+xi convention is silently changed by this separate model.

## What this does and does not solve

This advances from a coherently slipping plane or constrained scalar core to
a nonlinear three-component atomic-row equilibrium. It removes the concrete
omitted-force obstruction found in v7. It does not establish a quantitative
Al flow strength. A periodic pair can annihilate or leave a cell-spanning
registry fault; neither is a finite dislocation source or fatigue nucleation.

For a pre-existing line, the mechanical work per line length when slipped
area advances by b*delta d is tau*b*delta d. Thus an externally small stress
can act on existing defects without overcoming the whole perfect plane's
ideal traction peak. The internal pair force and image interactions must be
included, as they are here. This virtual-work statement is not an empirical
yield correction, an activation energy (a finite line length is absent), or
a calibrated rate law.

Still required: matched-Al vector/core/GSF validation, finite-line/loop/source
geometry and boundary conditions, domain-converged states and barriers,
collective-coordinate kinetic data. The saved candidate's independent
elastic errors persist. A_c remains untouched, and physical mobility,
seconds/Hz and production/UI promotion remain unavailable.

Coincident transverse rows (A=0) require a separate limiting representation;
the present reciprocal evaluator refuses that singular representation instead
of returning repaired forces. None of the accepted states approaches this
limit. Finite periodic cells and retained-neighborhood equilibria are not an
infinite isolated-core proof. There are19 actually solved static states across
7 completed cases, plus retained interrupted precursor records; these must not
be conflated with19 independent material validations.

## Execution and regression record

- New vector tests17, together with existing nonlinear tests:34 passed,
  100.54s. No old production test was removed or weakened.
- Full solver_v1:305 passed,1154.60s,0 skipped.
- App:31 passed,268.00s,0 skipped (actual Tk tests, not headless exclusions).
- Desktop smoke:PASS,2.80s; original LJ a0/kappa were reported unchanged.
- `git diff --check`:PASS including all newly added files.

The Python3 launcher was unavailable; the installed Python3.13 interpreter ran
the equivalent commands. Some regressions ran concurrently with research jobs,
so reported durations are actual wall times, not performance promises.
All results and the positive/negative scientific findings are under
`results/fcc111_active_interface/vector_core_v8/`; see SCHEMA.md first.
