# Analytic angular-environment research extension

Status: a new, explicitly separated research hypothesis, not a replacement of
the LJ pair and not an accepted aluminum potential. The full probability PDE
and its selectable energy models are unchanged. Physical kinetic time is
still unavailable. Mathematical verification does not by itself establish
material validity.

## 1. Why test angular information?

At the perfect triangular plane registry, let the ABC translation be tau.
A Shockley partial is u=tau, up to a lattice-equivalent representative. For
any radial plane interaction w(d,delta), inversion and periodicity imply

    w(h,tau) = w(h,-tau) = w(h,2 tau),  3 tau in Lambda.

Consequently the first opposing plane makes EXACTLY ZERO contribution to the
intrinsic stacking-fault energy, despite changing the local stacking. In
reciprocal space theta=G.tau=2 pi (m+n)/3, so

    cos(2 theta)-cos(theta)=0.

For the active interface the remaining pair difference is

    Delta W_pair = sum_(k>=2) k sum_G B_G(kh)
                         [cos((k+1)theta)-cos(k theta)].

The reciprocal-zero term cancels at fixed opening. The first nontrivial
terms start at k=2. With the shortest reciprocal vector

    G_min = 4 pi/(sqrt(3)b),  h/b=sqrt(2/3),
    exp(-2h G_min) = exp(-8 pi sqrt(2)/3).

The Bessel coefficient has additional algebraic prefactors, so the exponential
alone is NOT a numerical error bound or an exact barrier value. It does
explain why a smooth radial scalar environment can poorly distinguish two
stackings even when its bulk stiffness is accurate. The actual scalar-family
fit and fixed-decay feasibility tests quantify this limitation; they do not
constitute a theorem for every conceivable radial potential.

Angular environmental information has a precedent in the modified embedded
atom method; see [Baskes, Phys. Rev. B 46, 2727 (1992)](https://doi.org/10.1103/PhysRevB.46.2727).
The construction below is derived here as a minimal analytic test. It is NOT
an implementation of that published MEAM parameterization, and it does not
inherit any of its material validation.

## 2. An odd, rotationally covariant local moment

Use fixed microscopic units L0, and normalized kernel

    w(r*) = exp[-kappa(r*-1)] / rho_ref.

rho_ref is the total scalar density in target FCC, fixed during every strain
and interface calculation. It is an amplitude gauge, not a specimen scale.
For each atom form

    U_ijk = sum_(j_neighbor != atom) w(r*) r*_i r*_j r*_k,
    v_i = sum_j U_ijj,
    T_ijk = U_ijk - (v_i delta_jk + v_j delta_ik + v_k delta_ij)/5.

Then sum_j T_ijj=0. T is a symmetric traceless third-rank Cartesian tensor,
carrying the seven components of an l=3 solid harmonic. It changes sign under
inversion and transforms covariantly under rotations. Its Frobenius norm is
rotationally invariant and nonnegative.

Every atom in a centrosymmetric affine FCC bulk has neighbors r and -r with
equal weights. Thus T=0 exactly for ANY affine strain, not just at a fitted
equilibrium. The same is true for an isolated atom. The angular term therefore
does not silently refit bulk cohesion or homogeneous elastic constants.

Faulted and free-surface atoms need not have this inversion symmetry. This
distinguishes local stacking/environment information from scalar coordination.
No empirical plastic strain, yield stress, or hardening variable is added.

## 3. Infinite Poisson evaluation, not a finite-neighbor potential

For d>0 and the convention exp(-i G.x), the existing plane transform is

    Q=sqrt(kappa^2+|G|^2),
    H(d,G)=2 pi kappa exp(-d Q)(1+d Q)/Q^3.

Equivalently,

    H(d,G)=2 sqrt(2 pi) kappa (d/Q)^(3/2) K_(3/2)(dQ).

Multiplication by an in-plane coordinate corresponds to i times differentiation
in its reciprocal component:

    transform[x_i ... x_k w] = i^p partial_Gi ... partial_Gk H.

With t=|G|^2 and H_j=partial_t^j H,

    partial_Gi H = 2 G_i H_1,
    partial_GiGj H = 2 delta_ij H_1 + 4 G_i G_j H_2,
    partial_GiGjGk H =
       4(delta_ij G_k+delta_ik G_j+delta_jk G_i) H_2
       +8 G_i G_j G_k H_3.

For c normal indices, multiply by d^c; normal displacement changes d and must
also differentiate this factor. The plane moment is

    U_plane(d,delta) = (1/A_atomic_cell) sum_G
                         U_hat(d,G) exp(i G.delta).

Odd Fourier coefficients are imaginary; their products with sine phases give
real moments. Throwing them away would incorrectly erase the angular signal.
The G=0 mean is retained, including its nonzero normal components.

## 4. Small exact derivative recurrence

The code represents the transform as a finite list of monomials

    exp(-dQ) d^j Q^(-k).

No symbolic expression explosion or numerical derivative is needed:

    partial_t [exp(-dQ)d^j Q^-k]
      = -1/2 exp(-dQ)
           [d^(j+1) Q^(-k-1) + k d^j Q^(-k-2)],

    partial_d [exp(-dQ)d^j Q^-k]
      = exp(-dQ)
           [j d^(j-1) Q^-k - d^j Q^(-k+1)].

Starting from H/(2 pi kappa)=exp(-dQ)(Q^-3+d Q^-2), three t derivatives and
two d derivatives produce all required coefficients. Registry derivatives
multiply each phase by i G.e_s or -(G.e_s)^2. Mixed derivatives use both rules.
Finite differences are used only as independent tests.

## 5. Per-atom interface aggregation

The perfect bulk moment is zero. For a lower-half atom at depth r,

    Delta T_r = sum_(k>=r+1)
      [T_plane(kh+delta_a,k tau+s e_s)-T_plane(kh,k tau)].

Inversion relates the upper-half atom to -Delta T_r. The energy invariant per
primitive interface cell is therefore

    I3(a,s) = 2 sum_(r>=0) ||Delta T_r||^2,
    W_research(a,s) = W_scalar(a,s) + D I3(a,s),  D>=0.

The norm is taken PER ATOM after summing its entire neighbor moment. Summing
plane norms or taking one norm of a global crystal moment would be wrong.

For q_i in {a,s},

    I3_i = 4 sum_r T_r : T_r,i,
    I3_ij = 4 sum_r [T_r,i : T_r,j + T_r : T_r,ij].

Thus the energy, forces, and symmetric Hessian remain analytic/semi-analytic.
At the pristine cut W_research(h,0)=0 and its extra force is zero. Its local
interface curvature need not be zero: a LOCAL cut breaks centrosymmetry,
unlike a homogeneous affine bulk deformation.

## 6. Parameter count and admissibility

The smallest test ties the angular radial decay to the scalar density decay
and adds only D. A second, explicitly recorded sensitivity study allows a
separate angular decay. This is another parameter and must earn its place
through independent targets and sensitivity rank, not a desirable fatigue
response. The old single-exponential scalar EAM family remains intact.

D has energy units in this specified moment normalization. Its numerical
magnitude is not directly comparable with a bond energy without the invariant
scale. Neither D nor angular decay is a mobility or a physical time scale.

At fixed radial decays the new angular column is linear in D, so it can be
profiled with the LJ and embedding coefficients using deterministic bounded
linear least squares. Bulk rows of this column are exactly zero. Interface
targets, not bulk data, identify this added degree of freedom.

## 7. Verification and limitations

The new tests check trace removal, rotation/inversion, perfect bulk, full
registry periodicity, direct-versus-reciprocal convergence, and analytic first
and second derivatives. Reciprocal truncation uses an absolute coefficient
envelope instead of cancellation of an odd shell. Depth truncation is checked
by independent direct and tolerance refinement.

The finite direct lattice is only a validation evaluator. The canonical LJ
pair, its infinite Bessel sums, scalar density sum, and per-atom F are not
modified. Production probability dynamics are not connected to this research
surface. A better static objective or a positive ISF alone is NOT adoption:
held-out GSF/opening curves, stress branches, material stability, sensitivity,
and broader atomistic validation must also be examined.

This extension is a hypothesis for the missing information, not a proof of
novelty, a published Al MEAM, or a successful aluminum calibration.

## 8. Scenario-driven first-moment companion (separate nested test)

The first I3 fits improve intrinsic-fault energy but predict a local normal
tangent around 188 GPa versus about 126 GPa in the same rigid reference.
This is an actual held-out scenario failure, not an aesthetic preference.

The next tested component is the first odd environment moment with the SAME
angular radial range, so no additional decay parameter is introduced:

    V_i = sum_neighbors w(r*) r*_i,
    I1 = 2 sum_depth ||Delta V||^2,
    W_test = W_scalar + D3 I3 + D1 I1.

Again V=0 under all centrosymmetric affine bulk deformations. At an intrinsic
fault, in-plane threefold symmetry makes the vector much less informative
than the third harmonic; at a free/open interface its normal component is
sensitive. This separates two environmental channels rather than prescribing
a plastic strain or reducing a barrier by hand.

Its plane Fourier coefficients are simply

    Vhat_x = 2i G_x H_1, Vhat_y = 2i G_y H_1, Vhat_z = d H_0.

The same exact derivative recurrence, layer differences and per-atom norm
calculation apply. I1 and I3 are individually nonnegative, but the exploratory
D1 is allowed to be signed. Therefore the ADDED energy need not be nonnegative:
total local/global stability and separation work must be tested, not assumed.
A negative D1 cannot be certified by a good least-squares fit alone.

The source local normal-interface curvature is now an explicit eleventh fit
target (10% declared discrepancy scale); it is no longer held out in this
nested study. Other opening/registry curve points and stress paths remain
validation. The convex C embedding term is also tested only as a nested
already-derived option, not a high-order interpolating potential.

All of these are research candidates with recorded parameter counts and
source mappings. None alters the original EAM classes or production registry.

## 9. Independent elasticity gate: the lowest even harmonic

The odd moments vanish in EVERY affine centrosymmetric bulk. Thus neither
I1 nor I3 can repair a wrong homogeneous elastic tensor. This is an exact
limitation, not a reason to change mobility or the weights of a fatigue fit.
After the recorded odd-family bulk/interface compromise, the next nested
test is the lowest nontrivial even solid harmonic:

    Q_ij = sum_neighbors w(r*) [r*_i r*_j - |r*|^2 delta_ij/3],
    E2,atom = D2 Q:Q.

Q=0 in isotropically scaled perfect FCC by cubic symmetry. It need not vanish
under a general affine strain. Consequently this term changes neither the
isotropic equilibrium nor cohesion but CAN supply an independent anisotropic
elastic response. It is not a phenomenological elastic or plastic strain law:
the same per-atom environment also determines its interface energy.

The reciprocal second moment follows from the already derived H:

    Uhat_ij = -(2 delta_ij H1 + 4 G_i G_j H2), i,j in plane,
    Uhat_iz = 2i d G_i H1,
    Uhat_zz = d^2 H0,
    Qhat = Uhat - trace(Uhat) I/3.

For the perfect cubic cut, Q_bulk=0 and

    I2,interface = 2 sum_depth ||Delta Q_depth||^2.

The same difference/derivative formulas apply to I2 as I3. In the [111]
homogeneous strain modes, the unchanged same-plane contribution has zero
derivative under normal strain alpha or engineering shear gamma. At alpha=0,
gamma=0 the exact convergent plane sums give

    Q_alpha = 2 sum_(l>=1) l h partial_d Q_plane(lh,l tau),
    Q_gamma = 2 sum_(l>=1) l h partial_s Q_plane(lh,l tau),
    H_alphaalpha,2 = 2 D2 ||Q_alpha||^2,
    H_gammagamma,2 = 2 D2 ||Q_gamma||^2,
    H_hydro,2 = 0.

This is a single additional coefficient with the SAME angular decay as the
odd moments, not an independently tuned range. The first sector imposes
D2>=0; any signed exploration is separately recorded and requires total
finite-wavevector and interface stability checks. A good objective alone is
still insufficient. The original scalar EAM and all production defaults remain
unchanged. This is a new angular research family, not the published Al MEAM.

The rank-2 interface implementation requires a cubic FCC background. In a
noncubic bulk Q_bulk may be nonzero: the correct difference contains
2 sum_depth[2 Q_bulk:Delta Q + ||Delta Q||^2]. This extension is not implemented;
the constructor rejects it instead of dropping the background term.

## 10. Finite-wavevector STATIC stability, not frequencies

Uniform elasticity does not exclude nonuniform displacement instabilities.
For each reference neighbor R, n=R/|R|, at fixed bulk scalar F',F'',

    phi_eff' = phi_LJ' + 2 F' f',
    phi_eff'' = phi_LJ'' + 2 F' f'',
    Phi_eff(R) = phi_eff'' n n^T + (phi_eff'/r)(I-n n^T),
    K_pair(q) = sum_R Phi_eff(R) [1-cos(q.R)],
    g(q) = sum_R f'(r) n sin(q.R),
    K_scalar(q) = F'' g(q) g(q)^T.

Let B_l(R)=partial_R T_l be the tensor-valued bond derivative. Odd moments
have even B_l, whereas the quadrupole has odd B_2. Therefore

    L_1,3(q)=sum_R B_1,3(R)[cos(q.R)-1],
    L_2(q)=sum_R B_2(R) sin(q.R),
    K_l(q)=2 D_l L_l(q)^dagger L_l(q),
    K_total=K_pair+K_scalar+sum_l K_l.

The implied contraction is over tensor components. K(0)=0. Signed coefficients
require the TOTAL matrix to be checked. No mass is introduced; these are static
force constants, not phonon frequencies or physical Hz. The independent direct
validator refines neighbor radii 5/8/12 L0 and compares K against sinusoidal
displacement energy variations. Sampling 20 points on each Gamma-X/L/K does
not prove positivity over the entire Brillouin zone or finite-amplitude stability.

Later nonlocal_v5 audit: the old cubic q-path labels used the wrong orientation
for the unchanged +ABC stack. The corrected conversion is now explicit in
`FCC111Geometry.plane_basis_in_stacked_cubic_axes()`. New radius12/20/32
Gamma-X/L/K results are saved separately under nonlocal_v5. Old vectors were
valid samples, but their path names must not be retroactively certified.
The acoustic-limit matrix and independent affine cubic tensor now agree
after correcting the normal/in-plane coupling sign; see
[NONLOCAL_INTERFACE_ELASTICITY.md](NONLOCAL_INTERFACE_ELASTICITY.md).
