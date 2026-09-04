# Candidate activated registry rare-event route

Status: **CANDIDATE / PARTIAL SUCCESS — not an active governing law.**

> **2026-09-04 supersession note.** The coherent-patch assumption $\Delta G_{s,N}=N\Delta G_s$ and the associated $N$-sensitivity table below are retained only as a historical/sensitivity diagnostic. After spatially resolving registry as $s_j$, the preferred physical candidate is kink-pair nucleation from `CANDIDATE_SPATIAL_REGISTRY_KINK_FEEDBACK.md`. A local extended system need not translate all $N$ repeats coherently, so the physical activation barrier must be computed from the spatial kink-pair saddle $\Delta G_{\mathrm{kp}}$, not assumed to equal $N\Delta G_s$. No $N$ value or coherent-patch rate is adopted as a governing law.

## 1. Question

The previously tested conservative clocks fail at laboratory fatigue frequency:

- one local normal oscillator is adiabatic at laboratory Hz;
- a finite normal elastic segment is also adiabatic unless its length is macroscopic;
- the physical registry inertia gives a THz-scale natural frequency, so direct 20-Hz registry resonance is rejected.

The remaining question is different:

> Can the existing registry barrier itself generate a slow *rare-event* clock even though all conservative vibration modes are fast?

This document tests that question without fitting a damping coefficient or a named lifetime distribution.

## 2. Existing energy barrier

Retain the historical multilayer registry energy only as an extension diagnostic,

$$
U_0(a,s)=\sum_{k\ge1}\sum_{p\in\mathbb Z}
v_{m,n}\left(\sqrt{k^2a^2+(pb+s)^2}\right),
$$

with stable registry $s_0/b=1/2$. The barrier to the symmetry saddle is

$$
\Delta G_s(a)=U_0(a,0)-U_0(a,b/2).
$$

For the current $m=12.19$, $n=6$ diagnostic double sum, the registry equilibrium is

$$
a_{0,r}/b\approx0.9910707144.
$$

Matching the dimensionless normal curvature of this registry surface to the retained physical normal stiffness gives a diagnostic energy conversion. With the current Al normal calibration,

$$
E=69\ {\rm GPa},\qquad
A_0=6.0338\times10^{-20}\ {\rm m^2},\qquad
a_0=2.8627442948\times10^{-10}\ {\rm m},
$$

the current direct-sum audit gives approximately

$$
U_{aa}=108.3424,\qquad U_{ss}=26.1868,
$$

$$
r_K=U_{ss}/U_{aa}\approx0.2417043,
$$

and an energy scale of about

$$
E_r\approx0.0699044\ {\rm eV}
$$

per reference repeat. The equilibrium registry barrier is therefore

$$
\Delta G_s(a_{0,r})\approx0.10635\ {\rm eV/repeat}.
$$

This is a *diagnostic mapping*, not yet an independently calibrated Al slip barrier.

## 3. Coherent-patch barrier — retained sensitivity toy only

If a coherent patch contains $N$ repeats and all repeats are **forced by assumption** to cross the same ideal registry barrier coherently, both its energy and inertia are extensive and the trial barrier is

$$
\Delta G_{s,N}(a)=N\,\Delta G_s(a).
$$

This construction remains useful only to show how an exponential barrier can separate THz attempt frequencies from many-cycle event probabilities. It is no longer the preferred spatial transition path.

The integer $N$ is deliberately **not calibrated here**. It is characteristic-area / coherent-event information and belongs to a later spatial calibration step.

## 4. Conditional activated-rate diagnostic

To convert a barrier into a rate, extra physics is unavoidable. The following is therefore a **conditional finite-temperature hypothesis**, not a consequence of the conservative Hamiltonian alone.

Assume:

1. a finite-temperature bath exists;
2. intrawell relaxation is fast compared with interwell escape;
3. crossings are rare enough for a local activated-rate description;
4. for this historical diagnostic only, the coherent patch follows one dominant registry saddle;
5. the harmonic registry frequency supplies an order-one attempt frequency.

Then the historical sensitivity rate is

$$
k_N(a,T)=\nu_s(a)\exp\left[-\frac{N\Delta G_s(a)}{k_BT}\right].
$$

At the reference state, the retained curvature mapping gives

$$
\nu_s\approx1.41\times10^{12}\ {\rm Hz}.
$$

No mobility or damping coefficient is fitted in this diagnostic. Nevertheless the thermal-bath and activated-rate assumptions are real assumptions and must be validated before promotion. The spatially resolved candidate must replace $N\Delta G_s$ by the computed kink-pair saddle barrier $\Delta G_{\mathrm{kp}}$.

## 5. Quasistatic cyclic normal loading

At laboratory Hz, the normal coordinate is adiabatic relative to atomic and finite-segment elastic modes. Therefore a consistent low-frequency trial uses the stable normal force root rather than inertial resonance:

$$
\phi'(\lambda_{\rm eq}(t))=q(t)=\frac{\sigma_n(t)}{E}.
$$

Map that normal opening onto the registry surface as

$$
a_r(t)=a_{0,r}\lambda_{\rm eq}(t).
$$

For any physically justified activated rate $k(a_r,T)$, the accumulated hazard has the general form

$$
H(t)=\int_0^t k(a_r(\tau),T)\,d\tau,
$$

with local intact survival

$$
S(t)=\exp[-H(t)].
$$

The old coherent-patch model is one conditional example, not the adopted rate.

## 6. P0-to-surviving-subprobability map

Let $A(a_0,t;\sigma)$ denote the quasistatic normal-spacing map for a prepared local structural spacing label $a_0$. If a physically defined activated transition is treated as an absorbing transition out of an intact state, then

$$
H(a_0,t)=\int_0^t k(A(a_0,\tau),T)\,d\tau,
$$

and the intact subprobability measure is

$$
\rho(a,t)=\int P_0(a_0)
\exp[-H(a_0,t)]
\delta[a-A(a_0,t)]\,da_0.
$$

Its total mass is

$$
S(t)=\int\rho(a,t)\,da,
$$

and the spacing density conditioned on remaining intact is

$$
P_{\rm int}(a,t)=\frac{\rho(a,t)}{S(t)}.
$$

Thus, conditional on a justified activated rate, the reduced map

$$
P_0+\sigma(0:t)\longrightarrow \rho(a,t),\ S(t),\ P_{\rm int}(a,t)
$$

is closed without resolving ordered atom trajectories.

This is not yet the same as a closed evolution law for the full material spacing marginal after plastic transitions.

## 7. Historical numerical sensitivity — not calibration

For a 300 K diagnostic, 20 Hz sinusoidal tensile loading with mean stress 100 MPa and amplitude 100 MPa, the old coherent-patch rate gives approximately:

| coherent repeats $N$ | hazard per cycle | local survival per cycle | median cycles from this rate |
|---:|---:|---:|---:|
| 7 | $4.96\times10^{-2}$ | 0.95165 | $1.40\times10^1$ |
| 8 | $9.26\times10^{-4}$ | 0.999074 | $7.48\times10^2$ |
| 9 | $1.74\times10^{-5}$ | 0.9999826 | $3.99\times10^4$ |
| 10 | $3.28\times10^{-7}$ | 0.99999967 | $2.12\times10^6$ |
| 11 | $6.19\times10^{-9}$ | 0.999999994 | $1.12\times10^8$ |

There is **no adopted $N$** and these values are **not** current lifetime predictions. Their only retained use is to demonstrate exponential separation between a THz attempt scale and rare many-cycle events.

## 8. Two critical limitations

### 8.1 A single uniform registry crossing does not by itself evolve normal P

The ideal registry energy is exactly periodic,

$$
U_0(a,s+b)=U_0(a,s).
$$

Therefore a uniform transfer from well index $z$ to $z+1$ lands in an energetically equivalent ideal state. Spatially nonuniform registry is required to create a non-equivalent compatibility/core state; see `CANDIDATE_SPATIAL_REGISTRY_KINK_FEEDBACK.md`.

### 8.2 Zero-load thermal crossings exist in the ideal periodic model

Any finite-temperature activated model can give nonzero crossing probability at zero applied stress. Equivalent registry-well hopping is not automatically damage. A completed fatigue/plastic event must be tied to a non-equivalent residual structural state, such as a spatial kink/core configuration, before its rate can be called a damage rate.

## 9. Updated verdict

The robust part that survives is

$$
\text{an activation barrier can provide a many-cycle rare-event clock without a slow natural frequency.}
$$

The coherent $N\Delta G_s$ implementation does **not** survive as the preferred physical barrier. The next candidate is the spatial registry field, where a local transition creates kink/antikink incompatibility and a non-equivalent core. Its activation barrier must be computed from the spatial energy landscape.

Therefore this file is retained as a conditional historical stepping stone, while `CANDIDATE_SPATIAL_REGISTRY_KINK_FEEDBACK.md` is the current non-mainline candidate for the missing post-transition feedback.
