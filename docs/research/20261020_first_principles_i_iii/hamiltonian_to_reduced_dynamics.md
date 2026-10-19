# From the atomic Hamiltonian to reduced dynamics (PART II)

Date: 20 October 2026. This is a reconstruction, not a material calibration.
The lattice result in [PART I](lattice_energy_foundations.md) remains intact.
Here $q=(\delta a,s)$; a fixed shift of the opening coordinate to $a$
changes none of the arguments. $A(q)$, also denoted $F_{\rm eff}(q)$,
means a specified marginal free energy, not an area.

## 1. A finite microscopic system and an explicit coordinate

Start with $N$ atoms, positions $R\in\mathbb R^{3N}$, momenta $\Pi$,
positive diagonal mass matrix $m$, and

$$
H(R,\Pi;\lambda)=\tfrac12\Pi^Tm^{-1}\Pi+V(R;\lambda).
\tag{H1}
$$

The potential includes all adopted pair/embedding terms, boundaries and
constraints, not merely the bonds crossing one ideal cut. Fix a finite
volume or another confining reference for which the canonical partition
function exists. Remove translation zero modes consistently. Canonical
preparation is an assumption about an ensemble of initial conditions;
Hamilton's equations remain deterministic for each such initial condition.

For a concrete linear collective variable take the difference of the mass
centres of two disjoint finite groups, projected on a fixed unit normal
and tangent. Thus $q=CR+q_c$, where $C$ has rank two. For group masses

$$
M_+=\sum_{i\in+}m_i,\quad M_-=\sum_{i\in-}m_i,
\quad G=Cm^{-1}C^T,
\quad M_q=G^{-1}.
\tag{H2}
$$

For those orthogonal centre-of-mass difference projections,

$$
G=(M_+^{-1}+M_-^{-1})I_2,\qquad
M_q=\mu I_2,\quad \mu=\frac{M_+M_-}{M_++M_-}.
$$

This statement does **not** apply to an arbitrary average of interatomic
gaps: that observable has its own matrix $C$. Define

$$
B=m^{-1}C^TM_q,\qquad CB=I_2,
\qquad R=R_0+Bq+Ey,\quad CE=0,
$$

where the columns of $E$ span the remaining directions. The coordinate
transformation is affine and invertible on the adopted constraint space,
and $B^TmE=0$. The canonical reduced momentum is

$$
p=B^T\Pi=M_qCm^{-1}\Pi=M_q\dot q.
$$

The kinetic energy separates into $p^TM_q^{-1}p/2$ and a bath kinetic
energy with constant mass matrix. The eliminated variables include internal
vibrations, common centre-of-mass motion, all unretained relative motion,
and any local rearrangements permitted by $V$. They are not automatically
fast, equilibrated or irrelevant to fatigue.

For an infinite coherent half-crystal these group masses diverge. A finite
patch of $N_c$ identical primitive cells has energy $N_cW_\Gamma$ and
mass proportional to $N_c$. A per-cell energy cannot be combined with the
mass of an unspecified whole half-crystal. A statistical correlation area
used to count initiation opportunities is not, by definition, this coherent
patch area. These choices must precede a physical clock.

## 2. Integrating out the atoms: the first arrow

Fix the reference temperature (T>0), $\beta=(k_BT)^{-1}$, and reference
measure. Momentum integration in (H1) gives a constant for flat Cartesian
configuration space. The equilibrium density relative to flat $dq$ is

$$
Z(q)=\int dR\,\delta(\Xi(R)-q)e^{-\beta V(R)},\qquad
A(q)=-k_BT\log Z(q)+C_A.
\tag{H3}
$$

The units inside the logarithm are fixed by a reference measure; additive
constants cancel from differences. This is a constrained statistical
integral, not evaluation of $V$ at a single lattice configuration.

For the affine coordinate above the Jacobian is constant. Write, exactly,

$$
V(q,y)=U_{\rm Bessel}(q)+V_\perp(q,y),\qquad
Z_\perp(q)=\int dy\,e^{-\beta V_\perp(q,y)}.
$$

After a common reference subtraction,

$$
\boxed{A(q)-A(q_0)=U_{\rm Bessel}(q)-U_{\rm Bessel}(q_0)
             -k_BT\log\frac{Z_\perp(q)}{Z_\perp(q_0)}.}
\tag{H4}
$$

Here $U_{\rm Bessel}=N_cW_\Gamma^{\rm Bessel}$ only when that is the
declared rigid potential and patch normalization. Other contributions to
the atomic potential are then included in $V_\perp$. The equality

$$
A(q)=U_{\rm Bessel}(q)+\text{constant}
\iff Z_\perp(q)\text{ is independent of }q
\tag{H5}
$$

holds on the common accessible domain for that decomposition and measure.
The Bessel identity establishes no such independence.

For a nonlinear $\Xi$, the coarea formula reads

$$
Z(q)=\int_{\Xi^{-1}(q)}
  \frac{e^{-\beta V}}{J_\Xi}\,d\sigma_q,
\quad J_\Xi=\sqrt{\det(D\Xi D\Xi^T)}.
\tag{H6}
$$

For smooth vector fields $b_i$ satisfying
$b_i\cdot\nabla\Xi_j=\delta_{ij}$, integration by parts (with vanishing
boundary terms) gives the mean force

$$
\partial_i A=\left\langle b_i\cdot\nabla V
                              -k_BT\,\nabla\cdot b_i\right\rangle_q.
\tag{H7}
$$

For affine coordinates, $b_i=B_i$ is constant, but the conditional average
still contains relaxation and entropy. For nonlinear coordinates the mass
metric $D\Xi m^{-1}D\Xi^T$ can vary over one level set. Conditional
velocities can be a mixture of Gaussians; a covariance-derived effective
mass alone is not an exact canonical inertial reduction. Do not append a
"Fixman term" without specifying whether the measure is flat marginal,
surface constrained, or another constrained ensemble. The coarea/mean-force
framework is consistent with [Legoll--Lelièvre (2010)](https://arxiv.org/abs/0906.4865);
that paper starts from an overdamped microscopic reference and does not by
itself justify a Hamiltonian-to-Brownian limit.

### Relaxation and a controlled low-temperature expansion

Let $y_*(q)$ be an isolated, interior, nondegenerate constrained minimum;
write $U_{\rm rel}(q)=V(q,y_*)$, $H_\perp=\partial_y^2V(q,y_*)>0$,
and use a smooth measure density $w(q,y)\,dy$. With fixed bath dimension
$d$, Laplace expansion gives

$$
A(q)=U_{\rm rel}(q)+\frac{k_BT}{2}\log\det H_\perp(q)
       -k_BT\log w(q,y_*(q))
       -\frac{d k_BT}{2}\log(2\pi k_BT)+C_A+O(T^2).
\tag{H8}
$$

Use ratios of determinants/measure densities for dimensionless differences.
For the affine constant measure, the $w$ term cancels. Zero modes,
multiple minima, boundaries, changing bath dimension, instabilities and
nonuniform limits require separate treatment; then (H8) is not a blanket
formula. Neither relaxed zero-temperature energy nor rigid energy is
automatically the finite-temperature PMF.

### Exact counterexample using the same nonlinear Bessel energy

Keep the full counted Bessel $U(q)$, and add one transverse coordinate:

$$
V(q,y)=U(q)+\tfrac12\kappa(q)y^2,\qquad\kappa(q)>0.
$$

The minimum is $y_*=0$, so $U_{\rm rel}=U_{\rm rigid}=U$, yet Gaussian
integration gives **exactly**

$$
\Delta A=\Delta U+\tfrac12 k_BT
                 \log[\kappa(q)/\kappa(q_0)].
\tag{H9}
$$

This is not a quadratic replacement of the Bessel energy. In the isolated
dimensionless test $k_BT=0.2$,
$\kappa=2\exp[0.8\delta a+0.4(1-\cos 2\pi s)]$, at
$q=(0.1,0.17)$, $\Delta U=0.5170967738516374$ and
$\Delta A=0.5458266268875693$. Analytic Gaussian integration and numerical
quadrature agree at the recorded precision; quadrature error estimates are
also retained in the JSON. This proves non-identifiability of the PMF from
the rigid energy alone, not the magnitude of entropy in aluminium.

## 3. Exact probability transport before stochastic closure

For an initial density $\rho_0$, Liouville's equation is
$\partial_t\rho_t=-\{\rho_t,H\}$. Differentiating its coordinate marginal
and integrating by parts gives the exact identity

$$
\partial_t P(q,t)=-\nabla_q\cdot J(q,t),\qquad
\boxed{J=P(q,t)\,\mathbb E_{\rho_t}[\dot q\mid q].}
\tag{H10}
$$

The conditional velocity depends on the evolving hidden distribution. This
does not close as $J=-M(P\nabla A+k_BT\nabla P)$. Marginalization creates
an ensemble description; it does not produce white noise, a Markov property
or a mobility. Finite smooth Hamiltonian position paths have zero quadratic
variation on finite time intervals; a diffusion with nonzero noise has
positive quadratic variation. A singular/coarse-time limit is needed to
identify the latter with the former.

## 4. Exact conditional projection and memory

The following identities require an autonomous equilibrium reference with
invariant canonical measure, suitable domains for the Liouville operator,
and boundary terms that vanish. For $Z=(q,p)$, let

$$
Lf=\{f,H\},\quad\mathcal P f=\mathbb E_{\rm eq}[f\mid Z],
\quad\mathcal Q=I-\mathcal P.
$$

The Dyson identity yields, formally on that common domain,

$$
\dot Z(t)=e^{tL}\mathcal P LZ
 +e^{t\mathcal QL}\mathcal QLZ
 +\int_0^t e^{(t-u)L}\mathcal P L
                    e^{u\mathcal QL}\mathcal QLZ\,du.
\tag{H11}
$$

For the affine coordinate with constant $M_q$,
$\mathcal QLq=0$, $\mathcal PLq=M_q^{-1}p$, and
$\mathcal PLp=-\nabla A(q)$. Define

$$
F(u)=e^{u\mathcal QL}\mathcal QLp,\qquad
C_{ij}(z,u)=\mathbb E_{\rm eq}[F_i(u)F_j(0)\mid Z=z].
$$

Orthogonality gives $\mathcal PF(u)=0$; it does **not** assert independent
increments or Gaussianity. One can derive the memory drift explicitly.
For a smooth test function $g(Z)$, anti-self-adjointness of $L$ gives

$$
\langle g,LF_i(u)\rangle
 =-\langle Lg,F_i(u)\rangle
 =-\int dz\,\rho_Z(z)\,\partial_{p_j}g\,C_{ij}(z,u).
$$

The resolved part of $Lg$ drops out by conditional orthogonality. A second
integration by parts therefore gives (sum over (j))

$$
\boxed{\mathcal PLF_i(u)(z)
 =\rho_Z^{-1}\partial_{p_j}(\rho_Z C_{ij})
 =\partial_{p_j}C_{ij}-\beta C_{ij}(M_q^{-1}p)_j.}
\tag{H12}
$$

We used $\rho_Z\propto\exp[-\beta(A(q)+p^TM_q^{-1}p/2)]$.
In (H11) this expression is evaluated at $Z(t-u)$. It is generally
state- and momentum-dependent and includes the derivative of $C$.
Only in special cases with $C$ independent of $z$ does the memory term
reduce to a linear velocity convolution with $\Gamma(u)=\beta C(u)$.
An invariant reference and mathematical domain assumptions are part of this
identity; no existence theorem for every infinite LJ lattice is asserted.

This distinguishes conditional (Zwanzig) projection from a linear Mori
projection, whose resolved drift and correlation kernel refer to a chosen
linear observable span. One cannot combine a nonlinear PMF from one
projection with a constant convolution kernel from another without extra
conditions. See [Mori (1965)](https://doi.org/10.1143/PTP.33.423) and the
explicit nonlinear system/bath construction in
[Zwanzig (1973)](https://doi.org/10.1007/BF01008729).
The projection distinctions and general nonlinear-memory caveat are also
discussed by [Ayaz--Dalton--Netz (2022)](https://publications.imp.fu-berlin.de/2648/2/non-linear_GLEfinal.pdf);
their hybrid projection is a different construction, not the definition
of the conditional projection used in (H11)--(H12).

## 5. A solvable bath: an exact bridge with explicit restrictions

For any smooth nonlinear $U(q)$, including the full Bessel cut energy,
consider the finite Hamiltonian

$$
H=\tfrac12p^TM_q^{-1}p+U(q)+
\sum_j\left[\frac{P_j^2}{2m_j}+
 \frac{m_j\omega_j^2}{2}
 \left(x_j-\frac{c_j^Tq}{m_j\omega_j^2}\right)^2\right].
\tag{H13}
$$

The counterterm is inside the square and must be retained. At fixed $q$,
the translated Gaussian bath integral is independent of $q$, hence
$A=U+C_A$ **exactly in this model**. Eliminating each driven oscillator
by its elementary sine/cosine solution gives

$$
M_q\ddot q=-\nabla U(q)-\int_0^t\Gamma(t-u)\dot q(u)\,du
                         +\eta(t)+F_{\rm ext}(t),
\tag{H14}
$$
$$
\Gamma(t)=\sum_j\frac{c_jc_j^T}{m_j\omega_j^2}\cos\omega_jt,
$$
$$
\eta(t)=\sum_j c_j\left[
 \left(x_j(0)-\frac{c_j^Tq(0)}{m_j\omega_j^2}\right)\cos\omega_jt
 +\frac{P_j(0)}{m_j\omega_j}\sin\omega_jt\right].
\tag{H15}
$$

Prepare the bath conditionally at the specified $q(0)$ with the shifted
canonical Gaussian. Then its mean is zero and

$$
\mathbb E[\eta(t)\eta(u)^T]=k_BT\,\Gamma(t-u).
\tag{H16}
$$

Thus Gaussian random force and an equilibrium fluctuation relation follow
from this **specified bath and preparation**. Finitely many oscillators
give colored, recurrent noise and persistent memory, not Brownian noise.
Unshifted preparation instead produces the deterministic initial-slip term
$-\Gamma(t)q(0)$. The off-diagonal entries need not vanish; only a symmetry
of the full model and its coordinate transformation can enforce a vanishing
cross term at the corresponding symmetric state.

Changing $m_j,\omega_j,c_j$ leaves this $A=U+C_A$ unchanged while changing
the entire memory spectrum. Even after the first arrow is exact, the second
arrow is not identifiable from $A$ alone. This is a constructive result,
not a generic claim that every real bath is harmonic.

For externally imposed $F_{\rm ext}(t)$ acting only on $q$, (H14)--(H16)
remain exact for this shifted linear bath even under cyclic driving. For
a general driven anharmonic crystal, the reference distribution and hidden
state can evolve. A two-time, preparation-dependent memory/noise law is then
possible; an equilibrium stationary FDT cannot simply be assumed. Likewise
one equilibrium $A(q)$ cannot encode dislocation structures, residual
stress, or slow hidden history unless their conditional distribution is
actually slaved to $q$ on the required time scale.

## 6. What has and has not been established

The independent deterministic integration keeps the nonlinear Bessel normal
force, a three-mode Hamiltonian bath and a purely normal sinusoidal load.
After refinement, the largest force identity residual is
$1.23\times10^{-15}$ and the energy-minus-external-work residual is
$1.04\times10^{-14}$ on the recorded finite trajectory. Analytic Gaussian
covariance verifies (H16), and the saved recurrence/initial-slip checks
distinguish it from a diffusion. These are mathematical demonstrations,
not finite-temperature Al trajectory validation. The actual material's
coordinate, bath spectrum, preparation and elimination error remain inputs
to be established. [PART III](stochastic_limit_and_fokker_planck.md) states
the additional assumptions needed for a Fokker--Planck current.
