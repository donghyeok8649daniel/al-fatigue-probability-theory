# Candidate activated registry rare-event route

Status: **CANDIDATE / PARTIAL SUCCESS — not an active governing law.**

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
\boxed{\Delta G_s(a_{0,r})\approx0.10635\ {\rm eV/repeat}.}
$$

This is a *diagnostic mapping*, not yet an independently calibrated Al slip barrier.

## 3. Coherent-patch barrier

If a coherent patch contains $N$ repeats and all repeats cross the same ideal registry barrier coherently, both its energy and inertia are extensive. The natural frequency therefore does not slow with $N$, but the barrier does grow:

$$
\boxed{\Delta G_{s,N}(a)=N\,\Delta G_s(a).}
$$

This distinction is important. Patch size cannot produce a slow conservative oscillator, but it can make a barrier-crossing event exponentially rare.

The integer $N$ is deliberately **not calibrated here**. It is characteristic-area / coherent-event information and belongs to a later spatial calibration step.

## 4. Conditional activated-rate diagnostic

To convert a barrier into a rate, extra physics is unavoidable. The following is therefore a **conditional finite-temperature hypothesis**, not a consequence of the conservative Hamiltonian alone.

Assume:

1. a finite-temperature bath exists;
2. intrawell relaxation is fast compared with interwell escape;
3. crossings are rare enough for a local activated-rate description;
4. the coherent patch follows one dominant registry saddle;
5. the harmonic registry frequency supplies an order-one attempt frequency.

Then use the diagnostic rate

$$
\boxed{
k_N(a,T)=\nu_s(a)\exp\left[-\frac{N\Delta G_s(a)}{k_BT}\right].
}
$$

At the reference state, the retained curvature mapping gives

$$
\nu_s\approx1.41\times10^{12}\ {\rm Hz}.
$$

No mobility or damping coefficient is fitted in this diagnostic. Nevertheless the thermal-bath and activated-rate assumptions are real assumptions and must be validated before promotion.

## 5. Quasistatic cyclic normal loading

At laboratory Hz, the normal coordinate is adiabatic relative to atomic and finite-segment elastic modes. Therefore a consistent low-frequency trial uses the stable normal force root rather than inertial resonance:

$$
\phi'(\lambda_{\rm eq}(t))=q(t)=\frac{\sigma_n(t)}{E}.
$$

Map that normal opening onto the registry surface as

$$
a_r(t)=a_{0,r}\lambda_{\rm eq}(t).
$$

Then the activated hazard accumulated over one loading history is

$$
H_N(t)=\int_0^t k_N(a_r(\tau),T)\,d\tau,
$$

with local intact survival

$$
\boxed{S_N(t)=\exp[-H_N(t)].}
$$

This supplies a slow clock through rare barrier crossing even though the conservative coordinates themselves respond quasistatically.

## 6. P0-to-surviving-subprobability map

Let $A(a_0,t;\sigma)$ denote the quasistatic normal-spacing map for a prepared local structural spacing label $a_0$. If the activated transition is treated as an absorbing transition out of an intact state, then

$$
H_N(a_0,t)=\int_0^t k_N(A(a_0,\tau),T)\,d\tau,
$$

and the intact subprobability measure is

$$
\boxed{
\rho(a,t)=\int P_0(a_0)
\exp[-H_N(a_0,t)]
\delta[a-A(a_0,t)]\,da_0.
}
$$

Its total mass is

$$
S(t)=\int\rho(a,t)\,da,
$$

and the spacing density conditioned on remaining intact is

$$
P_{\rm int}(a,t)=\frac{\rho(a,t)}{S(t)}.
$$

Thus, **conditional on $T$, $N$, the barrier surface, and the activated-rate hypothesis**, the reduced map

$$
P_0+\sigma(0:t)\longrightarrow \rho(a,t),\ S(t),\ P_{\rm int}(a,t)
$$

is closed without resolving ordered atom trajectories.

This is not yet the same as a closed evolution law for the full material spacing marginal after plastic transitions.

## 7. Numerical sensitivity — not calibration

For a 300 K diagnostic, 20 Hz sinusoidal tensile loading with mean stress 100 MPa and amplitude 100 MPa, the direct-sum barrier plus harmonic attempt-frequency audit gives approximately:

| coherent repeats $N$ | hazard per cycle | local survival per cycle | median cycles from this rate |
|---:|---:|---:|---:|
| 7 | $4.96\times10^{-2}$ | 0.95165 | $1.40\times10^1$ |
| 8 | $9.26\times10^{-4}$ | 0.999074 | $7.48\times10^2$ |
| 9 | $1.74\times10^{-5}$ | 0.9999826 | $3.99\times10^4$ |
| 10 | $3.28\times10^{-7}$ | 0.99999967 | $2.12\times10^6$ |
| 11 | $6.19\times10^{-9}$ | 0.999999994 | $1.12\times10^8$ |

The important result is not any preferred $N$. There is **no adopted $N$** here. The result only shows that an activated coherent barrier naturally separates the THz attempt scale from a many-cycle rare-event scale through the exponential barrier factor.

It also matches the expected qualitative requirement that a single local domain can have survival extremely close to one per cycle.

## 8. Two critical limitations

### 8.1 Equivalent registry wells do not by themselves evolve the normal P

The ideal registry energy is exactly periodic,

$$
U_0(a,s+b)=U_0(a,s).
$$

Therefore transfer from well index $z$ to $z+1$ changes the unwrapped registry/plastic label but lands in an energetically equivalent ideal well. In the present ideal potential this transition **does not automatically change the normal interaction branch**.

Hence activated registry crossing can create a slow transition/first-passage clock, but by itself it does not solve the stronger requirement

$$
P_a(a,0)+\sigma(0:t)\to P_a(a,t)
$$

for the *full* spacing marginal with progressive structural evolution.

A physical post-slip defect state, boundary incompatibility, hardening/storage variable, or other non-equivalent structural state would be required for that feedback. None is inserted here.

### 8.2 Zero-load thermal crossings exist in the ideal periodic model

The same activated formula gives nonzero crossing probability even at zero applied stress. That is not automatically damage: equivalent registry-well hopping in the ideal periodic surface is not the same thing as creation of an irreversible defect.

Therefore the rate cannot be labeled a fatigue-damage rate until the physical meaning of a completed transition is established.

## 9. Verdict

The audit gives a mixed result:

$$
\boxed{
\text{activated registry barrier can provide a many-cycle rare-event clock}
}
$$

without lowering any natural frequency to laboratory Hz, but

$$
\boxed{
\text{the current periodic registry transition does not by itself drive progressive }P_a\text{ evolution.}
}
$$

Therefore this route is retained as a **candidate rare-event / plastic first-passage mechanism**, not promoted as the missing $P_a$ evolution law.

The next non-arbitrary question is whether a physically justified **non-equivalent post-transition structural state** can be derived from the existing mechanics. If not, the active theory should keep the normal first-passage model and leave irreversible G3 / plastic feedback open rather than inventing a state variable.
