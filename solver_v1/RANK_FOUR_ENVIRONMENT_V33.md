# First nonzero cubic STF background: rank-four analytic environment

Status: **separate research hypothesis**, not an adopted Al material. LJ pair,
production probability PDE, kinetic calibration, and existing STF1/2/3 models
remain unchanged. This is the next bounded test after v32's failed radial and
true-gradient controls, not a many-parameter fit to desired yield or fatigue.

## 1. Why this channel is mathematically different

The existing odd moments vanish in any affine centrosymmetric crystal.
STF rank2 vanishes in cubic FCC, although not under arbitrary affine strain.
Their positive squared energies therefore contribute a positive-semidefinite
harmonic square at pristine cubic bulk. The true scalar-density spatial
gradient also vanishes there. These properties limit how those channels can
repair the observed competition between interface curvatures and an excessively
large finite-wavelength restoring response.

A fourth-order symmetric traceless environmental tensor generally does **not**
vanish in cubic FCC. Its background produces a second harmonic contribution
that is not constrained to be positive. It is the first new STF order with this
property in the present hierarchy; simply squaring the existing rank3 invariant
is **not** a rank4 tensor and does not supply this background term.

Systematic invariant descriptions of atomic environments are discussed in
[Drautz, PRB99,014104(2019)](https://doi.org/10.1103/PhysRevB.99.014104).
The construction here is explicitly derived below, is not a claim to implement
that complete expansion, and does not adopt a neural or tabulated pair model.
All radial sums retain the exponential Poisson representation. The first trial
uses the **existing scalar decay and amplitude gauge**, hence adds only one
energy coefficient D4 and no fitted radial length.

## 2. Tensor definition, trace projection and units

With reduced relative position R_ij and the fixed scalar kernel f(r)=C exp(-kr),

    T_i,abcd = sum_(j!=i) f(r_ij) R_a R_b R_c R_d,
    Q_i = STF(T_i),
    E4 = D4 sum_i Q_i:Q_i.

The double contraction includes all81 Cartesian tensor entries. A9-component
orthonormal STF basis preserves exactly the same Frobenius norm. For symmetric T,
let U_cd=T_aacd and V=T_aabb. Then

    Q_abcd = T_abcd
      -(delta_ab U_cd+delta_ac U_bd+delta_ad U_bc
        +delta_bc U_ad+delta_bd U_ac+delta_cd U_ab)/7
      +(delta_ab delta_cd+delta_ac delta_bd+delta_ad delta_bc) V/35.

Explicit symmetrization before projection makes the implementation an
orthogonal projector for arbitrary test tensors. It has rank9, is idempotent,
and is rotationally covariant. For one neighbor this is the familiar harmonic
polynomial R_a R_b R_c R_d minus its traces. Q is constructed per atom after
the infinite neighbor sum; neither a global crystal norm nor plane-wise
embedding is substituted.

R is reduced by L0, and normalized f is dimensionless, so D4 is in eV. A
dimensional fourth moment differs by L0^4 and its squared-energy coefficient
by L0^-8. This is tensor dimensionalization, not a new physical characteristic
volume or a statistical area. A_c is absent.

## 3. Exact plane Fourier series

Use the existing triangular plane and ABC phases. Let

    t=|G|^2, Q=sqrt(k^2+t),
    Fhat(d,G)=2 pi k exp(-dQ)(Q^-3+dQ^-2).

For a tensor entry with z normal indices and4-z in-plane indices,

    Fourier[R_a R_b R_c R_d exp(-kr)]
       = d^z i^(4-z) partial_(G_i1)...partial_(G_i_(4-z)) Fhat.

The exact derivative recurrence is

    partial_Gj[P(G) Fhat^(n)(t)]
      = (partial_Gj P) Fhat^(n)(t)+2G_j P Fhat^(n+1)(t).

It requires at most four t derivatives. Normal derivatives differentiate the
whole d^z exp(-dQ) polynomial. The existing `_radial` recurrence computes them
analytically; registry derivatives of the Poisson series multiply by iG_x/iG_y.
Projection to nine components is linear and commutes with all derivatives.
Reciprocal G=0 is included. The full ten-jet contains value,a,x,y,aa,ax,ay,xx,xy,yy.

For the same plane, triangular sixfold symmetry gives

    sum x^4 f = sum y^4 f = (3/8) S4,
    sum x^2 y^2 f = (1/8) S4,
    S4=sum_(R!=0) |R|^4 C exp(-k|R|).

Other degree-four entries with odd x/y powers vanish, and normal entries are
zero. The scalar triangular shell series excludes self and uses a Voronoi
covering-radius/incomplete-Gamma tail, not an unexplained fixed neighbor cutoff.
The interplane sum remains Poisson-summed. Reciprocal/layer stopping envelopes
are empirical diagnostics, to be checked by direct/tolerance comparisons.

## 4. Nonzero bulk, active interface, and pressure identity

For the perfect FCC reference,

    Q_bulk = Q_same_plane + 2 sum_(l>=1) Q_plane(lh,l tau).

Even parity gives the factor2; Q_bulk must not be dropped. At depth j below
an active interface, define

    delta Q_j = sum_(l>=j+1) [Q_plane(a+(l-1)h,l tau+u)
                              -Q_plane(lh,l tau)],
    Q_j = Q_bulk+delta Q_j.

The interface energy difference per atomic interface cell is

    W4 = 2D4 sum_j [2 Q_bulk:delta Q_j+delta Q_j:delta Q_j].

This avoids subtracting infinite bulk constants. Bulk reference terms have
zero state derivative; the active derivatives are retained. At a=h,u=0,
W4=0 exactly, but a component's normal force need **not** vanish. The total
calibrated energy, not every additive component separately, must be stress-free.
An exact multiplicity/work check is

    partial_a W4(h,0) = (1/h) partial_eta_n E4_bulk(eta_n=0),

where eta_n is homogeneous normal engineering strain at fixed lateral lattice.
This identity is tested. The five bulk calibration rows include its normal
pressure, cohesive contribution **-D4||Q_bulk||^2** (isolated atom has zero
moment), and three independent elastic modes. Mode-to-C11/C12/C44 conversion
uses the unchanged atomic volume and eV/GPa factors. Hydrostatic curvature
currently uses the documented five-point strain stencil, independently refined;
interface gradients/Hessians are analytic.

## 5. Full Bloch Hessian, not only the positive square

Write t(R)=STF(R^4) f(r). For a cubic bulk displacement mode, let

    L_ic(q)=sum_R sin(q.R) partial_i t_c(R),
    K_ijc(q)=sum_R [1-cos(q.R)] partial_i partial_j t_c(R).

Even parity of t makes the first mode variation imaginary. The per-atom
Hessian coefficient is nevertheless real and symmetric:

    H4_ij(q)=2 sum_c L_ic L_jc+4 sum_c Q_bulk,c K_ijc.

The second term is mandatory. It can have either sign and is not removed by
using a positive D4. Therefore positivity of the complete material Hessian must
be checked; positivity of D4 alone is no certificate. H4 has an acoustic q^2
contribution in general, unlike odd-moment/true-gradient q^4 squares.

The independent direct evaluator differentiates degree-four monomials and
exp(-kr) analytically. Tensor norm majorants give

    ||t|| <= C r^4 exp(-kr),
    ||grad t|| <= C exp(-kr)(4 sqrt(3) r^3+k r^4),
    ||Hess t|| <= C exp(-kr)[36 r^2+(8 sqrt(3)+2)k r^3+k^2 r^4].

Existing FCC Voronoi exponential integrals bound omitted radial tails. With
errors dQ,dL,dK, an absolute Hessian bound is

    dH <= 2(2||L||dL+dL^2)
          +4(||Q||dK+||K||dQ+dQ dK).

The direct radius is a controlled validation/profiling approximation to the
infinite interaction, not a replacement of the canonical plane series.

## 6. Actual checks before material fitting

`test_rank_four_environment.py`: six tests actually passed in21.89s:

1. STF projector, trace-free condition and rotation-invariant norm.
2. Plane reciprocal sum versus independent direct tensor sum; derivative jets.
3. Nonzero cubic bulk background and three-mode elastic mapping, with stencil
   refinement and independent affine direct-atom energies.
4. Interface reference, normal-pressure work identity, derivative jets and
   registry periodicity.
5. Neighbor-level analytic Jacobian/Hessian versus independent differences.
6. Actual sinusoidal per-site energies versus the complete Bloch Hessian,
   amplitude extrapolation and radius/tail comparison.

The first fixed-shape study uses `run_radial_channels_v32 --rank-four` as a
shared coefficient-profile harness. It saves a primary D4>=0 trial and a
separate signed-D4 **diagnostic**. Neither bypasses full stability or material
validation. Three new states are declared in its definition before fitting and
are not training rows. Previous inspected source points are not called fresh.
Refer to `results/rank_four_environment_v33/fixed_shape` for actual completion
and fit results; an in-progress study must not be described as successful.

## 7. Completed fixed-range controls and independent validation

Five studies actually completed: the inherited scalar range plus predefined
k4=3,5,8,11. They contain15 LP profiles including five baseline replays. The
additional ranges change an atomic radial kernel, not stress, temperature,
mobility, A_c or a source pin length. Their D4 amplitude uses a fixed gauge
normalized at the undeformed bulk. The first nonnegative-amplitude trial and
the separately labeled signed control must not be conflated.

| k4 | eta, D4>=0 | eta, signed diagnostic | Excluded-q error, signed | Actual seconds |
|---|---:|---:|---:|---:|
| inherited4.946062 |12.507334|12.507334|210.0855%|110.946|
|3|12.507295|12.507295|210.0811%|128.653|
|5|12.507340|12.507340|210.0856%|129.952|
|8|12.507385|12.450127|208.5574%|92.405|
|11|12.507385|12.499613|209.8870%|84.349|

The minimum training eta selects k4=8, D4=-.225446749 eV, solely as an
independent diagnostic candidate. The eta<=1 criterion still fails by a large
factor. It is not adopted, and a signed coefficient is not a stability proof.

`independent_validation` actually evaluated the three predeclared new states,
not training data, and all eight q matrices. Maximum fresh Hessian error is
99.2024%; maximum excluded-q error208.5574%. Independent gradient/Hessian
difference errors are4.1783e-10 eV/L0 and4.7480e-10 eV/L0^2. Tightening the
new reciprocal/layer tolerance changes total H by2.7756e-17 eV/L0^2. The
largest finite-q radius change is1.2363e-4, with bound2.9809e-4 eV/L0^2.

The full reconstructed pristine energy is-1.64e-18 eV/interface cell, maximum
force3.36e-15 eV/L0 and H approximatelydiag(19.66091286,4.39359304,4.39359304)
eV/L0^2. Thus the exact tangent anchors and bulk-background pressure cancel
correctly, but matching that one equilibrium does not validate the rest of
the interface or spectrum. The large material errors are not accounted for
by measured differentiation/tail errors.

This establishes that the tested single rank4 invariant/range controls are
insufficient additions at the retained parent shape. It is not a global
impossibility theorem for all angular/radial analytic energies. No new core,
yield, fatigue or physical-time prediction is manufactured from this failed
candidate. The full source-binding JSON inherits the loader's historical
straight-row metadata; **actual observables here** use eV/atom for bulk,
eV/interface atomic cell for interface energy, GPa for cubic constants and
eV/L0^2 for displacement Hessians. They are not energies per dislocation-line
length, and none of these static Hessians is a phonon frequency.

## 8. Completed joint five-range capacity control

To avoid inferring a multi-range limitation from single-range failures,
`run_rank_four_capacity_v33` rebuilds the five existing columns from their
analytic energy and permits them simultaneously. Before the joint profiles,
all15 previous single-column prediction vectors are independently replayed;
their maximum discrepancy is exactly0 in the executed run. Targets, units,
parent shape, exact anchors and excluded-q roles are identical. No new
material term or radial range is invented for this control. Two joint LPs
completed in284.663s; results are in `joint_capacity_replay`.

| Joint capacity control | eta | Max excluded-q error | LJ base retained? |
|---|---:|---:|---|
| All five D4>=0 |12.50729550|210.0811%|positive coefficients, but target box fails|
| Five signed D4, diagnostic only |11.55612833|31.7309%|NO: repulsive LJ coefficient exactly0|

The signed solution includes cancellations among coefficients about-588.38
and+596.92eV in the fixed moment normalization. It is not a minimal or adopted
material law. In particular the reduced excluded-q error cannot be reported
as successful calibration while the LJ repulsive base is removed and eta>1.
The signed LP is a relaxed capacity bound, not an approval of those signs.
Its necessary target-box lower bound also exceeds1; restricting that same
fixed-shape problem to strictly positive LJ cannot improve the bound.

Both matrices have exact rank7 and conditional rank8, with conditional
singular values approximately[1.673621,1.028827,.949846,.613555,.129376,
.051149,.014217,.005187]. This is a scaled equality-tangent sensitivity,
not a statistical confidence interval or nonlinear radial identifiability.
Primal/dual gaps are below9e-15; maximum normalized primal violation in the
signed control is7.88e-10. No all-q spectral stability or fresh-interface
validation is claimed for this rejected capacity solution. Its excluded
points were already inspected in previous studies and are not called blind.

The first `joint_capacity` invocation stopped at the cached matrix boundary
because a JSON list was passed where a hashable tuple was required. It wrote
only a declaration, no completed calculation. That directory is preserved;
the tuple conversion was corrected and the separate successful replay above
actually executed. This was orchestration input handling, not a change to
the verified energy or derivatives.
