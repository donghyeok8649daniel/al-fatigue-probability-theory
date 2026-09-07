# Configurational registry transport and model plastic slip

## Status and scope

This document formalizes the plastic-slip content already implicit in the
deterministic probability density $P(a,s,t)$. It adds no empirical plasticity
law. The Bessel/Poisson-summed generalized-LJ landscape, the relaxed mapping
$f^*=\kappa_{\mathrm{axial}}\sigma/E$, and the normal-opening first-passage
surface are unchanged.

The present two-row LJ landscape is a minimal mechanistic proof of concept. It
is not a quantitatively calibrated FCC-aluminum slip surface.

## One state, four distinct responses

The coordinates have the following canonical meanings.

- $a$ is the normal inter-row/inter-layer separation.
- $s$ is the unwrapped total configurational registry/slip coordinate.
- Intrabasin motion in $a$ gives normal reversible or finite-rate response.
- Intrawell motion in $s$ gives configurational/anelastic memory.
- Interwell transport in $s$ changes accumulated registry and gives model
  plastic flow when numerically resolved.
- Crossing $a=a^\dagger(s,f)$ removes intact probability and gives crack
  initiation.

Thus the slip interfaces $s=(n+1/2)b$ transfer intact probability, whereas the
opening interface $a=a^\dagger(s,f)$ absorbs it. Their fluxes are not
interchangeable and are not combined into a generic damage variable.

## Registry decomposition and well populations

Use the exact half-open decomposition

$$
s=bn+\xi,
\qquad n=\left\lfloor\frac{s+b/2}{b}\right\rfloor\in\mathbb Z,
\qquad -\frac b2\leq\xi<\frac b2.
$$

The configurational well and its unnormalized intact mass are

$$
\Omega_n=\{(a,s):(n-1/2)b\leq s<(n+1/2)b\},
$$

$$
P_n(t)=\int_{\Omega_n}P(a,s,t)\,da\,ds.
$$

On a domain containing every populated well,

$$
\sum_nP_n(t)=M_{\mathrm{intact}}(t).
$$

The implementation enumerates every represented well; it does not assume
that only $n=-1,0,+1$ exist.

## Strain and survivor conditioning

The solver reports constitutive strain for the still-intact ensemble. For any
observable $X(a,s)$ define both the unnormalized survivor moment and the
conditional expectation:

$$
M_X(t)=\int X(a,s)P(a,s,t)\,da\,ds,
\qquad
\langle X\rangle_S=\frac{M_X}{S},
$$

where $S=\int P\,da\,ds$ is the direct intact mass. The current strain fields
use $\langle\cdot\rangle_S$, not unnormalized moments. Consequently

$$
\epsilon_{\mathrm{total}}
=\epsilon_a+\epsilon_\xi+\epsilon_p,
$$

$$
\epsilon_a=\left\langle\frac{a-a_0}{a_0}\right\rangle_S,
\qquad
\epsilon_\xi=\frac{\chi}{a_0}\langle\xi\rangle_S,
\qquad
\epsilon_p=\frac{\chi b}{a_0}\langle n\rangle_S.
$$

It follows identically that

$$
\frac{\chi}{a_0}\langle s\rangle_S
=\epsilon_\xi+\epsilon_p.
$$

A nonzero floating-point $\epsilon_p$ is only a signed registry-well moment;
it is not by itself evidence of experimentally resolved aluminum plasticity.

## Interwell flux and exact well balance

Write the existing Smoluchowski equation in conservative form,

$$
\partial_tP=-\partial_aJ_a-\partial_sJ_s,
$$

$$
J_s=-M_s\left(P\,\partial_sG+k_BT\,\partial_sP\right).
$$

At $s_{n+1/2}=(n+1/2)b$, define positive flux as transport toward increasing
$s$:

$$
\mathcal J_{n+1/2}(t)
=\int J_s(a,s_{n+1/2},t)\,da.
$$

Let $\dot A_n\geq0$ be opening absorption from well $n$. With reflecting outer
$s$ boundaries, direct integration gives

$$
\frac{dP_n}{dt}
=\mathcal J_{n-1/2}-\mathcal J_{n+1/2}-\dot A_n.
$$

An open outer $s$ boundary would add its explicitly signed boundary flux; the
present solver uses reflecting outer $s$ boundaries and exposes boundary mass
as a truncation diagnostic.

The Scharfetter--Gummel implementation stores both the signed interface flux
and the two-way activity. If the one-way nonnegative rates at an interface are
$j^+$ and $j^-$, then

$$
\mathcal J=j^+-j^- ,
\qquad
\mathcal J_{\mathrm{gross}}=j^++j^-.
$$

Therefore $|\mathcal J|$ is generally not gross hopping. Symmetric thermal
hopping can have large gross activity and zero net registry change.

## Registry moment, plastic flow, and crack-loss contamination

Define the unnormalized registry moment

$$
M_n(t)=\sum_n nP_n(t).
$$

Multiplying the exact well balance by $n$ and telescoping internal interfaces
gives

$$
\boxed{
\frac{dM_n}{dt}
=\sum_n\mathcal J_{n+1/2}-\sum_n n\dot A_n
}.
$$

The first term is actual interwell registry transport. The second is selective
removal by crack opening and is not plastic flow. The unnormalized interwell
contribution to axial plastic-strain rate is

$$
\dot\epsilon_{p,\mathrm{flow}}^{\mathrm{unnorm}}
=\frac{\chi b}{a_0}\sum_n\mathcal J_{n+1/2}.
$$

For the survivor-conditioned mean $\bar n=M_n/S$, with
$\dot S=-\sum_n\dot A_n$,

$$
\frac{d\bar n}{dt}
=\frac{1}{S}\sum_n\mathcal J_{n+1/2}
+\frac{-\sum_n n\dot A_n+\bar n\sum_n\dot A_n}{S}.
$$

The corresponding reported strain-rate decomposition is

$$
\frac{d\epsilon_p}{dt}
=\frac{\chi b}{a_0S}\sum_n\mathcal J_{n+1/2}
+\frac{\chi b}{a_0S}
\left(-\sum_n n\dot A_n+\bar n\sum_n\dot A_n\right).
$$

The code exposes these as `net_plastic_flow_rate` and
`selective_opening_plastic_rate`. Only the former is model plastic flow. It
also exposes the exact cumulative identity

$$
M_n(t)-M_n(0)-\sum_{n}\int_0^t\mathcal J_{n+1/2}\,dt'
+\sum_n nA_n(t)=0
$$

through `registry_moment_balance_residual`.

## Evidence hierarchy

The model uses four levels of interpretation.

1. A merely nonzero numerical value is not plasticity evidence.
2. Converged forward-plus-backward interface activity is configurational slip
   activity.
3. Converged directional transfer, well-population balance, and a signal above
   the convergence floor constitute model plastic flow.
4. A net registry change that persists after unload and zero-stress relaxation,
   and is stable under $s$-domain enlargement, is a residual-plastic candidate.

The previous apparent $\epsilon_p\sim10^{-11}$ at constant $f^*=3$ collapsed
to approximately $4\times10^{-16}$ under aligned-grid refinement. It remains
classified as numerical leakage.

## Why periodic energy can retain accumulated slip

The local landscape obeys

$$
W(a,s+b)=W(a,s),
$$

so neighboring minima are locally energetically equivalent. Nevertheless the
unwrapped translations differ by $\Delta s=b$. A transition $n\to n+1$ can
therefore return the local registry to an equivalent minimum while changing
the accumulated relative translation. This is the minimal mechanism for slip
without an empirical evolution law. The model contains no physical hardening
law, and none is inferred from periodicity.

## Candidate crystallographic grounding of $\chi$

For one slip system let $\mathbf m$ be the unit slip direction,
$\mathbf n_{\mathrm{slip}}$ the unit slip-plane normal, and $\mathbf e$ the
unit uniaxial loading direction. With shear $\gamma$,

$$
\boldsymbol\epsilon_p
=\frac{\gamma}{2}
\left(\mathbf m\otimes\mathbf n_{\mathrm{slip}}
+\mathbf n_{\mathrm{slip}}\otimes\mathbf m\right).
$$

Its axial projection is

$$
\epsilon_{p,\mathrm{axial}}
=\mathbf e\cdot\boldsymbol\epsilon_p\mathbf e
=\gamma(\mathbf e\cdot\mathbf m)
(\mathbf e\cdot\mathbf n_{\mathrm{slip}}).
$$

If $s$ is a slip displacement across plane spacing $h$, then $\gamma=s/h$.
Comparison with $\epsilon_s=\chi s/a_0$ gives the candidate

$$
\boxed{
\chi=\frac{a_0}{h}
(\mathbf e\cdot\mathbf m)(\mathbf e\cdot\mathbf n_{\mathrm{slip}})
}.
$$

`crystallographic_slip_projection` evaluates this geometry for supplied unit
vectors but never changes default solver physics. FCC aluminum conventionally
slips in the $\{111\}\langle110\rangle$ family. The current $s$ is only one
effective slip/registry coordinate: no orientation, plane spacing, or mapping
between the two-row coordinate and a particular FCC plane has yet been
calibrated.

## Current bound-basin configurational landscape

The relevant finite-temperature diagnostic is the bound-basin potential of
mean force, using the actual opening saddle as the normal integration limit.
It is not the fixed-$a$-box free energy, and truncated Gibbs is not the
absorbing fast-process quasi-stationary distribution.

| $f^*$ | principal $s$ minimum | forward $s$ saddle | $\Delta\mathcal F_s^\dagger$ | opening barrier at $s$ saddle |
|---:|---:|---:|---:|---:|
| 0.5 | 0.00216 | 0.47630 | 0.59551 | 0.92166 |
| 2.0 | 0.01173 | 0.35177 | 0.26614 | 0.32321 |
| 3.5 | 0.03519 | 0.17179 | 0.03474 | 0.10720 |
| 4.0 | 0.06761 | 0.09497 | 0.000363 | 0.09492 |
| 4.020 | -- | -- | $3.33\times10^{-5}$ | positive |
| 4.025 | 0.07934 | 0.08113 | $1.03\times10^{-7}$ | positive |
| 4.030 | no distinct minimum/saddle | no distinct minimum/saddle | 0 | -- |

The minimum curvature is positive and the saddle curvature negative up to the
merger. The configurational spinodal is therefore bracketed by approximately
$4.025<f^*_{s,\mathrm{sp}}<4.030$. Barrier reduction, exact spinodal, and
finite-time probability transfer are three different events.

The earlier fixed-domain estimate near $f^*\simeq2$ was an artifact of the
tensile-tilted large-$a$ tail. The bound-basin audit supersedes it. The
principal $s=0$ normal-opening spinodal is about $f^*=5.252$, but the opening
threshold varies with registry and the minimum over one period is about
$3.753$. Hence opening loss can compete before the principal configurational
branch reaches its spinodal.

## Loading-time and mobility sensitivity

For the canonical $M_a=1$, $M_s=0.05$, $k_BT=0.02$, the pristine slow
relaxation time is approximately $\tau_s=0.397$ model time. A one-cycle pulse
from zero to $f^*=4.0$ and back gave:

| period | $2\pi\tau_s/T$ | absorbed opening mass | max outside-center mass | max $|\epsilon_p|$ |
|---:|---:|---:|---:|---:|
| 0.04 | 62.36 | $6.57\times10^{-9}$ | $5.56\times10^{-16}$ | $1.44\times10^{-16}$ |
| 0.4 | 6.24 | $2.74\times10^{-5}$ | $2.80\times10^{-10}$ | $7.25\times10^{-11}$ |
| 2.5 | 1.00 | $2.55\times10^{-3}$ | $1.57\times10^{-6}$ | $4.08\times10^{-7}$ |

Thus a nearly vanished static barrier does not guarantee transfer: the
distribution needs model time to reach the interface, while opening absorption
acts concurrently.

At period 2.5 and peak $f^*=4.0$, the sensitivity values
$M_s=0.01,0.05,0.25$ produced max $|\epsilon_p|$ of approximately
$6.6\times10^{-10}$, $4.0\times10^{-7}$, and $1.7\times10^{-4}$, respectively.
This is a sensitivity result, not an aluminum mobility calibration.

## Slow candidate, convergence, and order of events

A slower period-10 pulse to $f^*=4.1$ produced a candidate positive-well tail.
The observed refinement sequence was:

| $(n_a,n_s)$ | max $dt$ | absorbed mass | max outside-center mass | max $|\epsilon_p|$ |
|---:|---:|---:|---:|---:|
| 19, 30 | 0.0500 | 0.03425 | $6.40\times10^{-5}$ | $1.72\times10^{-5}$ |
| 23, 36 | 0.0333 | 0.03543 | $6.15\times10^{-5}$ | $1.65\times10^{-5}$ |
| 27, 42 | 0.0250 | 0.03746 | $4.50\times10^{-5}$ | $1.21\times10^{-5}$ |
| 31, 48 | 0.0200 | 0.03647 | $3.79\times10^{-5}$ | $1.02\times10^{-5}$ |

The sign and existence of a small interwell tail persist, but its magnitude has
not reached a convincing grid-independent plateau and an explicit reference
at this long duration is computationally impractical because of the known LJ
repulsive-region CFL stiffness. The conservative classification is therefore
**unresolved interwell transfer coupled to already resolved opening loss**.

With a 5-model-time zero-stress hold, 3-, 5-, and 7-well domains at equal
cells-per-well produced the same coarse $\epsilon_p=3.63\times10^{-5}$ and
outside-center mass $1.34\times10^{-4}$; outer-boundary mass fell from
$1.5\times10^{-13}$ to $7.7\times10^{-37}$. This rules out the immediate outer
boundary as the source at that resolution, but it does not overcome the grid
convergence failure above and therefore does not establish residual
plasticity.

Resolved opening absorption is already percent-level in this slow candidate.
The current evidence therefore supports order-of-events classification **B:
opening occurs before any independently resolved interwell plastic
transition**. The model has not demonstrated usable, converged plastic slip
before crack opening.

## Interpretation and future aluminum landscape

The correct current claims are:

- normal response, intrawell memory, interwell transfer, and crack first
  passage are mathematically distinct consequences of one $P(a,s,t)$;
- the old $10^{-11}$ plastic signal is numerical leakage;
- a slow, high-force pulse produces a candidate interwell tail only after
  opening loss is already resolved;
- no calibrated residual plasticity or experimental aluminum yield stress is
  established.

The strongest next physical validation is to replace

$$
W_{\mathrm{LJ}}(a,s)\longrightarrow W_{\mathrm{Al}}(a,s)
$$

with an EAM/MEAM or generalized stacking-fault/registry surface computed for a
declared FCC orientation and slip system. That comparison must ground the
$s$-corrugation, normal-opening surface, plane spacing $h$, and the geometric
projection $\chi$. It must not be achieved by lowering the LJ barrier,
increasing $M_s$, or adding empirical hardening merely to manufacture slip.
