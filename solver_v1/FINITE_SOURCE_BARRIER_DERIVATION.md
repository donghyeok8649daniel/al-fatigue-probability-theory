# Finite pinned-line barrier from the same outer elastic energy

This is a scoped static continuation of `finite_source_reference.py`, not a
new empirical plasticity law, a calibrated Al source, or a production PDE.
The line coefficient follows the existing infinite LJ/Bessel bulk elasticity
and two-half-crystal Schur kernel. Nothing is fitted to yield or fatigue.
The missing core/obstacle physics is not replaced by a multiplicative area.

## Geometry and energy

For a planar line with two fixed endpoints separated by L, use its tangent
angle theta, resolved Burgers magnitude b, and pressure p=tau b [J/m²].
The existing long-wave energy is

    G[curve] = integral gamma(theta) d ell - p Area,
    gamma(theta)=k(theta) log(R/r_core) [J/m],
    k(theta)=b_vector^T Re K0(n_perp) b_vector/(2 pi).

L and R/r_core are explicit source/core geometry inputs. No default physical
Al source length is known. R=L, when used in scenarios, is a declared outer
asymptotic convention fixed during variation. This is NOT A_c. The expression
omits core energy, nonlocal finite parts and actual pinning/obstacle mechanics.

With gamma pi-periodic, mirror-symmetric, and T=gamma+gamma''>0, variation
gives the anisotropic curvature balance T kappa=p. Define

    Q(theta)=gamma sin(theta)+gamma' cos(theta), Q'=T cos(theta),
    r(theta)=-gamma cos(theta)+gamma' sin(theta), r'=T sin(theta).

An equilibrium bowed curve has

    x(theta)=Q(theta)/p,
    y(theta)=[r(theta_e)-r(theta)]/p,
    d ell=T/p dtheta,
    Q(theta_e)=p L/2.

The positive-T minor branch has theta_0<pi/2; a second overhanging branch
has theta_1=pi-theta_0. Their fold is p_c=2 gamma(pi/2)/L. An overhanging
saddle cannot be represented by the single-valued graph solver y(x); the
parametric curve is essential. A failed root is not declared a spinodal.

## Actual stability within the stated line energy

For a normal variation eta with fixed endpoints, second variation is

    delta²G = integral T[(deta/dell)²-kappa² eta²] d ell
            = p integral[-theta_e,theta_e] [(deta/dtheta)²-eta²] dtheta.

For eta=alpha sin[j pi(theta+theta_e)/(2theta_e)],

    delta²G / alpha² = p theta_e [(j pi/(2theta_e))²-1].

The minor arc has index0, the major arc index1, and the fold one zero mode.
This statement is exact within the positive-T planar local-line model, NOT
atomistic/nonlocal/core stability. Independent perturbed polygonal-curve
energies verify both the positive and negative modes under grid refinement.

## A finite positive barrier, without an inserted activation volume

Subtract the major and minor stationary curves with the SAME L,p,gamma.
Using the primitives above gives cancellation-safe integrals

    DeltaG = (2/p) integral[theta0,pi-theta0]
                       [Q(theta)-Q(theta0)] T(theta) sin(theta) dtheta,
    DeltaArea = (2/p²) integral[theta0,pi-theta0]
                              Q(theta) T(theta) sin(theta) dtheta.

The first integrand is nonnegative by Q monotonicity/symmetry. This is an
analytic consequence, not clipping a computed negative barrier. Quadrature
is refined and unresolved near-fold values are refused. The exact envelope
identity is

    -d DeltaG/d tau = b DeltaArea [m³].

This stress derivative follows from the finite geometry; no characteristic
volume is introduced into W. In particular, multiplying the old per-cell
barrier by a guessed number of atoms is unnecessary and unjustified.

For constant line tension Gamma, put theta=arcsin(tau/tau_c), R=Gamma/p:

    DeltaG = Gamma²/p [pi-2theta-sin(2theta)],
    DeltaArea = R² [pi-2theta+sin(2theta)].

Near the fold, DeltaG scales as (1-tau/tau_c)^(3/2), while the derivative
scales as its square root. At fixed outer logarithm and fixed tau/tau_c,
DeltaG scales with L and b DeltaArea with L². These are derived geometry
laws, not fitted stress reductions or empirical hardening.

## What the executed scenarios do and do not establish

`run_finite_source_barrier.py` uses shared0K C11/C12/C44=114/62/32GPa and
the corresponding FCC Burgers translation. It scans explicit HYPOTHETICAL
spans0.2/1/5 micrometres, core cutoffs1/2/4b, and six subcritical fractions.
All54 actual cases, angular/quadrature refinement, derivative checks and
minor/major curve coordinates are saved. No span is labeled measured or
selected because it reproduces a desired MPa value.

The public wire activation data concern the **continuous** relaxation signal
after separating abrupt bursts. They are not automatically activation data
for an entire pinned-source emission event. Apparent rate sensitivity also
depends on prefactors, internal stress and source populations. Fitting L to
those values without matching the operative event would be circular.

The result explains how a same-energy spatial mechanism can introduce a
low-stress scale distinct from homogeneous ideal registry strength. It does
not yet calibrate actual yield: core energy, source/pin geometry and populations,
finite-part elasticity, temperature and condition-matched dissipation remain.
No Arrhenius rate, probability, normal mobility, t0 or physical PDE Hz is emitted.
The existing production energy registry, Smoluchowski dynamics and UI are unchanged.
