# Stochastic limit, probability current and the physical clock (PART III)

Date: 20 October 2026. $A(q)$ denotes the flat-coordinate PMF defined in
[PART II](hamiltonian_to_reduced_dynamics.md). All limits below are conditional
on their stated microscopic or mesoscopic model. A numerical solution of a
chosen PDE does not demonstrate those hypotheses for aluminium.

## 1. Memory becoming white is a separate limit

The exact finite harmonic bath gives a Gaussian force with colored covariance
and recurrences. Gaussianity concerns finite-dimensional distributions;
Markov behavior concerns conditional time dependence. One implies neither
white covariance nor independent increments of the integrated force.

A controlled mesoscopic family can be constructed by choosing symmetric
positive definite friction $\zeta$ and

$$
\Gamma_\varepsilon(t)=\varepsilon^{-1}\zeta e^{-t/\varepsilon},
\quad t\geq0,\qquad
\int_0^\infty\Gamma_\varepsilon(t)dt=\zeta.
\tag{S1}
$$

For continuous velocities away from the initial layer, its one-sided
convolution approaches $\zeta\dot q(t)$. If

$$
\mathbb E[\eta_\varepsilon(t)\eta_\varepsilon(u)^T]
 =k_BT\varepsilon^{-1}\zeta e^{-|t-u|/\varepsilon},
$$

the distributional covariance limit is $2k_BT\zeta\delta(t-u)$.
The factor two comes from the **two-sided** covariance integral; the
causal friction integral in (S1) has mass $\zeta$. A continuous spectral
bath or a stationary OU auxiliary force can realize this chosen family;
a fixed finite oscillator list cannot have this decay for all time.

For a non-Gaussian atomic bath, a white-noise limit needs mixing and an
appropriate functional central limit/scaling result, as well as control of
the coupling, initial layer and feedback. A short observed correlation at
one state or a nearly Gaussian histogram alone is insufficient. The bath
must relax faster than changes in the resolved state and applied load:

$$
\tau_{\rm bath}\ll\tau_q,\ \tau_{\rm load}.
\tag{S2}
$$

Slow phonons, defects and bath recurrences can violate this separation.
Neither the lattice-sum identity nor the static surface energy supplies
the mixing theorem, friction tensor or these time scales.

Under the chosen white-bath limit, constant $T$, constant symmetric
$\zeta>0$ and a constant mass matrix, the candidate inertial diffusion is

$$
dq=M_q^{-1}p\,dt,\qquad
dp=[-\nabla A(q)+F_{\rm ext}(q,t)-\zeta M_q^{-1}p]dt
                            +\sqrt{2k_BT\zeta}\,dB_t.
\tag{S3}
$$

Its Kramers equation for phase-space density $f(q,p,t)$ is

$$
\partial_t f=-\nabla_q\!\cdot(M_q^{-1}pf)
 -\nabla_p\!\cdot[(-\nabla A+F_{\rm ext}-\zeta M_q^{-1}p)f]
 +k_BT\nabla_p\!\cdot(\zeta\nabla_p f).
\tag{S4}
$$

Integrating over $p$ recovers $J=\int M_q^{-1}pf\,dp$, but does not yet
give a position-only Fick law. Initial/boundary conditions and integrability
are part of this statement. An absorbing position boundary can have an
inertial boundary layer and is not justified by formally dropping $p$
at the boundary.

## 2. Overdamping is another approximation

For a scalar mass/friction the momentum time is $\tau_p=\mu/\zeta$.
For matrices use the decay rates of $\zeta M_q^{-1}$, with consistent
units and norms. A position-only approximation additionally needs

$$
\tau_p\ll\tau_q,\ \tau_{\rm load}.
\tag{S5}
$$

Rapid bath decorrelation and rapid momentum relaxation are different
requirements. Taking colored-noise and small-mass limits together without
specifying their ratio can change the effective drift; one should not
silently exchange these limits.

For spatially varying friction, smooth symmetric positive mobility
$M(q)=\zeta(q)^{-1}$, a constant temperature, constant microscopic mass,
and a fluctuation relation $\sigma\sigma^T=2k_BT\zeta$, the controlled
small-mass limit of the corresponding inertial SDE is the Itô equation

$$
\boxed{dq_i=\left[-M_{ij}\partial_jA+M_{ij}F_{{\rm ext},j}
                  +k_BT\partial_jM_{ij}\right]dt+B_{ik}dB_k,
\qquad BB^T=2k_BT M.}
\tag{S6}
$$

Regularity, uniformly positive friction on the visited domain, appropriate
growth/nonexplosion conditions or localization, compatible initial data,
and boundary control are required. This is an instance of the
[Hottovy--McDaniel--Volpe--Wehr (2015) limit](https://arxiv.org/abs/1404.2330).
Their theorem starts from an inertial stochastic equation; it is not proof
that an unspecified atomic bath already satisfies (S3). State-dependent
mass, nonsymmetric friction or non-isothermal noise needs its own formula.

## 3. Deriving the current, including the drift correction

The Itô forward equation is

$$
\partial_tP=-\partial_i(b_iP)+\tfrac12\partial_i\partial_j
                          [(BB^T)_{ij}P].
$$

Writing it as a continuity equation defines

$$
J_i=b_iP-k_BT\partial_j(M_{ij}P).
$$

Substitution of (S6), rather than omission of the mobility derivative,
gives the explicit current

$$
\boxed{J=-M(q)\big[(\nabla A-F_{\rm ext})P+k_BT\nabla P\big].}
\tag{S7}
$$

For conservative external force $F_{\rm ext}=-\nabla V_{\rm ext}$, the
zero-current stationary density is proportional to
$e^{-\beta(A+V_{\rm ext})}$ **if it is normalizable and boundary conditions
allow it**. Under suitable no-flux/decay boundaries at fixed load, direct
integration by parts yields

$$
\frac{d}{dt}\int P(A+V_{\rm ext}+k_BT\log P)dq
 =-\int P\,[\nabla(A+V_{\rm ext})+k_BT\nabla\log P]^T
 M[\nabla(A+V_{\rm ext})+k_BT\nabla\log P]dq\leq0.
\tag{S8}
$$

Boundary removal, time-dependent loading and nonsmooth densities require
their own boundary/work terms; (S8) is not automatically a dissipation
identity for the manuscript's intact-only absorbed population.

For constant $M$, (S6) reduces to additive noise and (S7) to Rev.11's
current after consistent nondimensionalization **only if** its input energy
has already been identified with $A$ and all earlier reduction assumptions
hold. Existing Rev.11 already calls its diffusion an assumed closure;
the audit does not label this conditional algebra an error.

For a Stratonovich representation the drift is

$$
b_i^{\rm Strat}=b_i^{\rm Ito}
 -\tfrac12\sum_{jk}B_{jk}\partial_jB_{ik}.
\tag{S9}
$$

In multiple dimensions the correction depends on the chosen noise frame
$B$; it is not universally half of $k_BT\nabla\cdot M$. In one scalar
dimension $B=\sqrt{2k_BT M}$ gives that familiar half correction. The
Itô equation (S6) with its full divergence term avoids this ambiguity.

An independent one-dimensional test, on a periodic registry quotient at
fixed opening, uses $M=1.2+0.25\sin 2\pi s$. Correct Itô drift gives a
maximum Gibbs flux residual $3.56\times10^{-15}$; omitting $k_BT M'$
gives (2.465414995254872). This checks the mathematical drift convention
on that finite quotient, not a global Gibbs law on the opening domain.

## 4. Domain, loading and ensemble compatibility

A rigid cleavage energy tends to a finite separation energy as opening
tends to infinity. It is periodic in unwrapped $s$. Consequently the
integral of $e^{-\beta U_\Gamma}$ over unrestricted opening and unwrapped
registry diverges even with zero load. A tensile tilt worsens the opening
tail. Do not invoke a normalized global Gibbs measure for that domain.

An equilibrium projection must instead specify a confining finite reference,
a metastable restricted ensemble, or another normalizable construction.
A periodic registry quotient has finite phase volume but loses cumulative
slip unless a winding/cell index is retained. A reflecting numerical box is
one artificial reference; agreement with its Gibbs distribution does not
establish physical confinement or its suitability near opening. A
quasi-stationary absorbed ensemble also differs from equilibrium.

For a uniform imposed traction conjugate to the actual patch opening,
$V_{\rm ext}=-A_\Gamma T_n(t)\delta a$, hence
$F_{{\rm ext},a}=A_\Gamma T_n(t)$ and $F_{{\rm ext},s}=0$.
The latter is the purely normal loading case used in our deterministic
bath check. In Rev.11, $-f[(a-a_0)+\chi s]$ instead generates a tangential
generalized component $\chi f$. That is work-consistent only for the
declared projected loading geometry; it cannot stand in for a pure normal
traction unless $\chi=0$. A scalar experimental stress is not by itself
the two generalized forces.

An equilibrium PMF at each fixed load can be a useful reference, but its
instantaneous substitution during a cycle requires hidden-state relaxation
faster than loading. In the special shifted harmonic bath, arbitrary forcing
of $q$ leaves its explicitly derived bath covariance unchanged. In a real
anharmonic driven solid, this observation is not a license to impose a
stationary equilibrium fluctuation relation without testing it.

## 5. Nondimensionalization and identifiable time

Let $q=\ell q^*$, $A=E_0A^*$, $M=M_0M^*$, and use the same length
scale for both coordinates. With

$$
t_0=\frac{\ell^2}{M_0E_0},\qquad \theta=\frac{k_BT}{E_0},
\tag{S10}
$$

(S7) has dimensionless diffusion tensor $D^*=\theta M^*$. Unequal
coordinate scales require the corresponding tensor transformations, not a
single unexplained scalar clock. If a coherent patch has $N_c$ repeated
cells and $A_{\rm tot}=N_c A_{\rm cell}$, dividing the energy by $N_c$
requires $M_{\rm cell}=N_c M_{\rm tot}$ and
$(k_BT)_{\rm cell}=k_BT/N_c$ in the divided equation. Changing only the
energy normalization changes the noise model.

A static $U(q)$ cannot fix $M_0$. Multiplying an overdamped mobility
by (c>0) preserves the same equilibrium density while multiplying all
relaxation rates by $c$. A fully specified Hamiltonian with physical
masses, lengths and energies *does* have a Newtonian time scale; what is
missing from a static energy alone is its eliminated bath dynamics and
the justification of a particular overdamped time map.

Existing repository v31 research includes actual Al99 MD and finite-band
plane-coordinate kinetic estimates. It must not be described as "no dynamics
data." Its own report states that the low-frequency drag uncertainty permits
zero, leaving the upper mobility interval unbounded. The equivalence of
plane-coordinate PMF, temperature/energy normalization and the production
local coordinate is unresolved. Therefore those results do not currently
identify the physical $M_0$ of this manuscript. See
`solver_v1/PHYSICAL_TIME_AND_MOBILITY.md` and
`solver_v1/WEAK_KINETICS_AND_MATERIAL_V31.md` at the provenance commit.

Without that mapping, model-time first-passage results cannot be reported
as measured seconds, hertz or cycles-to-initiation. A cycle count
$n=\int f_{\rm load}(t)dt$ is meaningful only with a compatible physical
time/load-frequency map.

## 6. Executed limit illustrations and the decision

Only for a **local quadratic control** about the independently evaluated
Bessel minimum do we replace $A$ by its Hessian. This is labelled as a
local approximation, not the nonlinear lattice model. An exponential-memory
OU embedding then has an analytic Gibbs covariance; the Lyapunov residual
is $8.00\times10^{-15}$. Refined correlation times reduce the recorded
colored-to-white frequency-response error from $3.63\times10^{-3}$ to
$1.51\times10^{-4}$; refined masses reduce inertial-to-overdamped error
from $2.13\times10^{-3}$ to $8.49\times10^{-5}$. These finite parameter
comparisons illustrate two separate limits, not a theorem inferred from
three data points. The same Hessian and temperature with twice the mobility
have exactly twice the linear relaxation rates.

The scientific default is therefore: retain the Bessel energy as an
explicit rigid energetic input; retain Rev.11's Fokker--Planck equation as
a conditional model; elevate neither to material-derived stochastic fatigue
until the PMF, memory, normalization, time and event-observation checks are
passed. The more faithful intermediate candidate is the explicitly defined
PMF plus a tested inertial/memory reduction. Which approximation is adequate
is an empirical and mathematical question still to be decided for the same
atomic system and coordinate.
