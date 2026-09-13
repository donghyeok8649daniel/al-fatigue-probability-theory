# True spatial environmental-density gradient — separate analytic hypothesis

## Scope and motivation

The completed v31/v32 coefficient, radial, density-shape and gauge controls
do not pass the joint interface/finite-wavelength material discrepancy box.
This document derives one additional nonuniform-environment observable before
testing whether it helps. It does not replace LJ, the production energy,
mobility, probability equation, or any material calibration file.

Wu, Lu, Garcia-Cervera and E, *Density-gradient-corrected embedded atom method*,
Physical Review B **79**, 035124 (2009),
[DOI](https://doi.org/10.1103/PhysRevB.79.035124),
[author manuscript](https://arxiv.org/abs/0901.4594), discuss nonuniform-density
corrections to EAM with Al tests. Their leading kinetic gradient structure
motivates the invariant below. Our auxiliary lattice density is **not measured
electron density**, and this is **not their complete potential** or a calibrated
Al gradient-corrected EAM. Literature motivation is not validation.

## Two different vectors

For relative positions R_ij=r_j-r_i, a fixed normalized exponential kernel
f(r)=C exp(-k r), and x_i=sum_j f(r_ij), define

    Q_i = sum_j R_ij f(r_ij),
    g_i = grad_(r_i) x_i = -sum_j f'(r_ij) R_ij/r_ij.

The existing rank-one moment is Q, not g. For an exponential kernel g has
weights k f(r) multiplying the **unit** neighbor direction, whereas Q weights
the full distance. In a general environment these vectors need not be parallel.
No center/self atom is included. The density normalization C is chosen once
so that x_bulk=1 and is not recalibrated during strain/opening.

The minimal positive research correction is

    E_grad = D sum_i |g_i|^2/x_i,       D >= 0.

It is applied per site after summing the whole environment. It is never the
norm of a global crystal sum. In affine centrosymmetric bulk g_i=0 exactly;
therefore this term changes neither cohesion nor homogeneous elastic constants.
It can change an interface and nonuniform displacement modes. For a separated
exponential neighbor, |g|^2/x is proportional to k^2 x and tends to zero, not
infinity. The bulk normalization and the same scalar density must be used in
both numerator and denominator. D=0 restores the prior energy exactly.

Lengths in implementation are reduced by L0; g is a derivative with respect
to reduced position and D is in eV. In dimensional notation the corresponding
coefficient multiplying |grad_phys x|^2/x is D L0^2, in eV m^2. There is no
statistical correlation area, kinetic time, or empirical strength factor here.

## Exact reciprocal representation

For triangular plane area A_atomic_cell and separation d>0, the existing exact
Poisson transform is

    Q = sqrt(k^2+|G|^2),
    T(d,G) = 2 pi C k exp(-d Q) (Q^-3+d Q^-2),
    x_plane(d,delta) = A_atomic_cell^-1 sum_G T(d,G) exp(i G.delta).

Thus the plane contribution to the central-position gradient, in the ordered
components (normal,in-plane-x,in-plane-y), is

    g_plane = -A_atomic_cell^-1 sum_G
              [T_d, iG_x T, iG_y T] exp(i G.delta).

The minus sign follows from R=r_neighbor-r_center, not a convention fitted to
data. The reciprocal zero mode contributes to the normal component and must
not be discarded. Registry derivatives multiply by iG_x/iG_y; normal
derivatives differentiate exp(-dQ)(Q^-3+dQ^-2) analytically. The ten-jet of g
requires derivatives of the density through order three. Reusable monomial
recurrences in `transform_terms(0,n)`, n=0..3, evaluate these without numerical
coordinate differentiation. `gradient_plane` uses all reciprocal symmetry
vectors, with the existing ABC phase and triangular geometry.

The identity is exact. The computed reciprocal/layer sums are semi-analytic:
three successive small envelopes are **empirical stopping diagnostics**, not
a rigorous total remainder bound. The implementation reports shell/layer counts
and last envelopes. Direct sums and tighter tolerances independently test them.

## Active-interface bookkeeping

For a lower-half atom at depth j>=0, bulk inversion gives total g_bulk=0.
Only cross-interface neighbors change, so

    g_j(a,u) = sum_(l>=j+1) [
        g_plane(a+(l-1)h, l tau+u) - g_plane(lh,l tau)].

The subtracted bulk terms are fixed constants: their coordinate derivatives
must not be subtracted. Upper-half gradients have opposite sign and equal
norm/density. Hence the interface correction per atomic interface cell is

    W_grad(a,u) = 2 D sum_(j>=0) |g_j(a,u)|^2/x_j(a,u).

At perfect bulk W_grad=0 and its gradient is zero, but its **interface Hessian
is not zero**. The two pristine interface-tangent anchors must therefore include
the new column in calibration. No infinite bulk energies are subtracted.
Site density/gradient arrays may have different converged depths; the existing
sitewise jet helper pads absent changes by their correct bulk values.

For any state derivative alpha and beta, writing I=g.g and v=1/x,

    I_alpha = 2 g.g_alpha,
    I_alpha,beta = 2(g_alpha.g_beta+g.g_alpha,beta),
    v_alpha = -x_alpha/x^2,
    v_alpha,beta = 2 x_alpha x_beta/x^3-x_alpha,beta/x^2,
    (Iv)_alpha,beta = I_alpha,beta v + I_alpha v_beta
                      + I_beta v_alpha + I v_alpha,beta.

The code reuses the verified `screened_site_jet(..., power=-1)` chain rule.

## Independent finite-wavelength consequence

For a cubic bulk Bloch displacement u_i=u exp(i q.R_i), inversion symmetry
gives the linear response

    delta g_i = R_h(q) u exp(i q.R_i),
    R_h(q) = sum_R [1-cos(q.R)] Hess_R f(R).

Consequently, at x_bulk=1,

    H_grad(q) = 2 D R_h(q)^T R_h(q).

This is positive semidefinite and O(|q|^4). It is **not** the existing scalar
embedding F'' term, whose rank-one factor instead contains sum sin(q.R) grad f.
The Voronoi exponential-tail bound on Hess f is propagated as

    ||delta H_grad|| <= 2D [2 ||R_h|| delta_Rh + delta_Rh^2].

The canonical interface uses reciprocal sums. Direct atom sums for this Bloch
matrix are independent controlled validation, with a reported radius/tail.
No atom mass or Hz is involved in this static Hessian.

## Actual checks and fit protocol

`test_spatial_density_gradient.py`: three tests passed (1.96 s) before fitting:

1. A direct triangular sum independently reconstructs g, its Jacobian and
   third radial tensor; it also checks that g is not silently Q.
2. Interface energy/gradient/Hessian finite differences, perfect reference,
   and full registry periodicity.
3. Actual sinusoidal per-site energies versus the Bloch matrix, amplitude
   extrapolation, and the q^4 acoustic limit.

`run_radial_channels_v32 --spatial-gradient` tests **one** added D>=0 at the
prior shape. It retains the same seven exact bulk/interface tangent anchors,
115 inspected interface curvature rows, and four fit-q matrices. Excluded q
matrices remain outside selection. Full primal/dual LP and conditional SVD
are saved. This is a necessary fixed-shape compatibility test, not a complete
spectral certificate, material calibration, or validation of physical time.

Actual output belongs in `results/radial_channels_v32/spatial_density_gradient`.
Actual run completed in 50.292 s. The normalized minimax error changes from
12.50738498 to **12.50709992**, while the criterion is <=1. The fitted new
D=1.20897994e-5 eV does not resolve the conflict. Excluded-q relative error is
210.08196%, versus the prior 210.08016%. Equality residual is 3.91e-14 and
duality gap -1.07e-14; conditional rank rises from 3 to 4 (four singular values
1.014855, .742716, .237352, .022372), but added identifiability is not adequacy.

The active minimax conflict is explicit (H in eV/L0^2):

| Observable | Reference | Prediction | Absolute dual weight |
|---|---:|---:|---:|
| saddle_Hxx | -.18642463 | .04673852 | .120748 |
| v19_new_2_Hxx | .96429798 | -.28641201 | .183301 |
| v20_new_11_Haa | 2.86542642 | -.71839103 | .014912 |
| v32_q1_H22 | 37.17324859 | 69.26565258 | .681039 |

The additional positive gradient energy raises some nonuniform stiffnesses;
that does not automatically repair the required, partly softer, source
restoring response. Wrong signs at these particular source coordinates and
large finite-q errors are material mismatch; they do not by themselves locate
the candidate's own saddles. This fixed-shape result does not rule out all
gradient-based analytic families. The tested one-amplitude extension is **not
adopted** and is not a solution of tasks1–3.

The new D bound has zero dual multiplier and D is strictly positive at this
LP solution. Its certificate therefore remains feasible if only the D sign
constraint is removed: allowing a signed D alone cannot lower this particular
fixed-shape optimum. This statement is not extrapolated to the existing
quartic coefficient, whose lower-bound multiplier is .01497856, or to other
radial shapes. Negative-energy extensions are not silently enabled.
