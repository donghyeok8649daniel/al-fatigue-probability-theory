# Final reduced P0-to-survival closure / 최종 P0-생존 축약법칙

Status: **ACTIVE REDUCED LAB-SCALE CLOSURE for the declared 1D normal-instability hypothesis.**  
The exact finite LJ chain and its $P$-$u$-$\Theta$ projection remain the microscopic reference model. This document supplies a separate controlled laboratory-time-scale reduction. It is not claimed to be a complete dislocation/slip theory of real single-crystal aluminum.

## 1. Problem that this closure must solve

The target is

$$
P_0(\lambda)+\sigma(0:t)
\longrightarrow
P_b(\lambda,t),\ S(t),\ F_{\mathrm{ci}}(t)
$$

without reconstructing every microscopic trajectory and without prescribing a named PDF or empirical fatigue-damage variable.

Previous audits established three facts:

1. the exact finite conservative chain is not autonomous in $P$ alone;
2. laboratory-frequency normal elastic motion is quasistatic on the atomic time scale and does not by itself create progressive cycle memory;
3. permanent drift of the normalized spacing PDF is not mathematically required if survivor mass is removed by first passage.

The remaining task is therefore to close the survivor escape flux, not to force permanent $P$-shape drift.

## 2. Structural P0, not instantaneous thermal displacement

$P_0(\lambda)$ is the reference-phase **structural/prestress spacing density** on the stable branch

$$
0<\lambda<\lambda_c.
$$

It is not the instantaneous finite-temperature displacement distribution. Fast thermal coordinates are eliminated from the structural PDF and enter only through the rare-event crossing law below.

Let the external reduced normal traction at the reference phase be

$$
q_{\mathrm{ref}}=\frac{\sigma_{\mathrm{ref}}}{E}.
$$

For a starting value $\lambda_0$ define the local residual conjugate bias

$$
\boxed{
q_r(\lambda_0)=\phi'(\lambda_0)-q_{\mathrm{ref}}.
}
$$

This is the minimal $P_0$-only embedding closure: $\lambda_0$ is, by construction, the stable local equilibrium at the reference load. It does not claim that the full finite-chain neighbour configuration is reconstructible from $P_0$.

## 3. Quasistatic stable-branch map

The local reduced effective energy is

$$
\psi(\lambda;q_{\mathrm{tot}})
=\phi(\lambda)-q_{\mathrm{tot}}\lambda,
$$

with

$$
q_{\mathrm{tot}}(\lambda_0,t)
=q_r(\lambda_0)+q(t),
\qquad
q(t)=\frac{\sigma(t)}{E}.
$$

In the laboratory-frequency quasistatic regime, the intact structural spacing follows the stable root

$$
\boxed{
\phi'[\Lambda(\lambda_0,t)]
=\phi'(\lambda_0)+q(t)-q_{\mathrm{ref}},
\qquad
\phi''[\Lambda(\lambda_0,t)]>0.
}
$$

At the reference phase this gives exactly

$$
\Lambda(\lambda_0,t_0)=\lambda_0.
$$

Differentiating the stable-root equation gives the closed characteristic velocity

$$
\boxed{
\dot\Lambda
=\frac{\dot q(t)}{\phi''(\Lambda)}.
}
$$

Therefore the nonabsorbing structural density obeys the closed transport equation

$$
\boxed{
\partial_tP
+\partial_\lambda
\left[
\frac{\dot q(t)}{\phi''(\lambda)}P
\right]
=0,
\qquad 0<\lambda<\lambda_c.
}
$$

This is the requested $P_0+\sigma(0:t)\to P(t)$ law under the declared quasistatic/local-prestress reduction.

## 4. Exact push-forward form of the mechanical part

On the stable branch the map is monotone. Differentiating with respect to $\lambda_0$ gives

$$
\phi''(\Lambda)\frac{\partial\Lambda}{\partial\lambda_0}
=\phi''(\lambda_0),
$$

so

$$
\boxed{
\frac{\partial\Lambda}{\partial\lambda_0}
=\frac{\phi''(\lambda_0)}{\phi''(\Lambda)}.
}
$$

Hence

$$
P(\lambda,t)
=\int P_0(\lambda_0)
\delta[\lambda-\Lambda(\lambda_0,t)]\,d\lambda_0,
$$

or, when the inverse map exists,

$$
\boxed{
P(\lambda,t)
=P_0(\lambda_0)
\frac{\phi''(\lambda)}{\phi''(\lambda_0)}.
}
$$

A closed load cycle returns this nonabsorbing structural $P$ to the same reference-phase shape if the stable branch is never lost. Fatigue accumulation therefore does not come from forcing a permanent drift into this mechanical PDF.

## 5. Mechanical instability boundary

The operational initiation boundary remains

$$
\phi''(\lambda_c)=0,
$$

with

$$
\boxed{
\lambda_c
=\left(\frac{m+1}{n+1}\right)^{1/(m-n)}.
}
$$

The corresponding maximum stable reduced traction is

$$
\boxed{
q_c=\phi'(\lambda_c).
}
$$

If

$$
q_{\mathrm{tot}}\ge q_c,
$$

no stable intact root exists and the local domain reaches the operational instability deterministically.

## 6. Characteristic cohesive area remains symbolic

Let $A_c$ be the characteristic area participating coherently in one local normal-instability event. It is **not calibrated here** and is not a FEM element area or a specimen independent-cell area.

Let $A_0$ be the existing mechanical atomic/reference area. The coherent effective mass is taken consistently with the same reduced normal coordinate as

$$
m_c=\frac{A_c}{A_0}m_a.
$$

The characteristic-domain energy scale is

$$
E_c=EA_ca_0.
$$

Because stiffness and mass both scale with $A_c$, the local small-oscillation frequency is independent of the unknown characteristic area:

$$
\omega_s^2
=\frac{EA_0}{m_aa_0}\phi''(\lambda_s)
=\frac{1}{t_0^2}\phi''(\lambda_s).
$$

Thus

$$
\boxed{
\nu_s(\lambda_s)
=\frac{\sqrt{\phi''(\lambda_s)}}{2\pi t_0}.
}
$$

## 7. Rare thermal first passage to lambda_c

For the current stable spacing $\lambda_s<\lambda_c$, define the energy climb to the operational absorbing boundary

$$
\Delta\psi_c(\lambda_s)
=
\left[
\phi(\lambda_c)-\phi'(\lambda_s)\lambda_c
\right]
-
\left[
\phi(\lambda_s)-\phi'(\lambda_s)\lambda_s
\right].
$$

The characteristic-domain barrier is

$$
\boxed{
\Delta G_c(\lambda_s)
=EA_ca_0\,\Delta\psi_c(\lambda_s).
}
$$

In the rare-event regime

$$
\Delta G_c\gg k_BT,
$$

and assuming fast thermal re-equilibration inside the intact well compared with both the fatigue period and the escape time, the positive-flux transition-state approximation gives

$$
\boxed{
k_c(\lambda,T;A_c)
=\frac{\sqrt{\phi''(\lambda)}}{2\pi t_0}
\exp\left[
-\frac{EA_ca_0\Delta\psi_c(\lambda)}{k_BT}
\right].
}
$$

This is not an arbitrary Kramers kernel or fitted fatigue rate. The barrier and prefactor are both inherited from the already declared generalized-LJ normal mechanics. The additional physical assumptions are finite temperature, local re-equilibration, and rare first passage.

## 8. Closed survivor equation

Let $P_b(\lambda,t)$ be the intact survivor subdensity. The final reduced equation is

$$
\boxed{
\partial_tP_b
+\partial_\lambda
\left[
\frac{\dot q(t)}{\phi''(\lambda)}P_b
\right]
=-k_c(\lambda,T;A_c)P_b.
}
$$

Initial condition:

$$
\boxed{
P_b(\lambda,t_0)=P_0(\lambda)
}
$$

when the reference population is initially intact.

Along a characteristic,

$$
W(\lambda_0,t)
=\exp\left[
-\int_{t_0}^{t}
k_c(\Lambda(\lambda_0,s),T;A_c)\,ds
\right].
$$

Hence the exact characteristic solution of the reduced equation is

$$
\boxed{
P_b(\lambda,t)
=\int P_0(\lambda_0)
W(\lambda_0,t)
\delta[\lambda-\Lambda(\lambda_0,t)]\,d\lambda_0.
}
$$

If the inverse map exists,

$$
\boxed{
P_b(\lambda,t)
=P_0(\lambda_0)
W(\lambda_0,t)
\frac{\phi''(\lambda)}{\phi''(\lambda_0)}.
}
$$

The local survival and cumulative first-passage probability are therefore

$$
\boxed{
S(t)=\int_0^{\lambda_c}P_b(\lambda,t)\,d\lambda
=\int P_0(\lambda_0)W(\lambda_0,t)\,d\lambda_0,
}
$$

$$
\boxed{
F_{\mathrm{ci}}(t)=1-S(t).
}
$$

This completes the requested mapping

$$
\boxed{
P_0+\sigma(0:t)
\longrightarrow
P_b(\lambda,t),\ S(t),\ F_{\mathrm{ci}}(t).
}
$$

## 9. Periodic loading

For a periodic stress waveform of period $T_f$, each reference label has one-cycle integrated hazard

$$
\mathcal H_c(\lambda_0)
=\int_0^{T_f}
k_c[\Lambda(\lambda_0,t),T;A_c]dt.
$$

For one delta reference state,

$$
S_N=\exp[-N\mathcal H_c].
$$

For a general structural $P_0$,

$$
\boxed{
S_N
=\int P_0(\lambda_0)
\exp[-N\mathcal H_c(\lambda_0)]\,d\lambda_0.
}
$$

Thus repeated cycles progressively remove the more weakly protected part of $P_0$ even when the reversible mechanical map itself returns to the same phase each cycle. The normalized survivor distribution can change by **selection**, not by an invented permanent deformation of every surviving spacing.

In the strict quasistatic phase-controlled limit,

$$
\mathcal H_c\propto\frac1f.
$$

Therefore the model predicts approximately frequency-independent life in physical time and a cycle count proportional to $f$ until quasistatic/local-equilibrium assumptions fail. This is a falsifiable prediction, not a hidden calibration rule.

## 10. Why this route survives the earlier no-go tests

### Deterministic finite-chain mixing

At laboratory frequency the pure normal elastic dynamics is quasistatic compared with atomic mechanical time scales. A label-preserving reversible cycle cannot generate new first-passage events indefinitely. This route remains a reference/no-go for high-cycle accumulation.

### Local conservative oscillator

The local oscillator supplied $P_0\to P(t)$ but became quasistatic at laboratory Hz and therefore had no cycle accumulation. The present reduction keeps its valid quasistatic mechanical limit and adds only the finite-temperature first-passage physics required for renewed rare opportunities.

### Spatial registry kink

The registry-kink audit found a rare formation scale but only very shallow post-formation trapping under the ideal pure-normal reduced surface. It is therefore not promoted as the primary long-lived memory state. It remains a transient/plasticity-extension diagnostic.

### Arbitrary diffusion

No diffusion coefficient, damping constant, Gaussian PDF, Weibull life law, or fitted fatigue kernel appears in the final equation. Thermal physics enters only through $k_BT$ and the rare-event transition-state flux derived from the same normal potential.

## 11. What is calibrated later

The following are deliberately left for experiment or higher-fidelity calculation:

- characteristic cohesive area $A_c$;
- the physical structural/prestress $P_0$;
- single-crystal loading-axis $E$ and the mechanical reference $A_0,a_0$ if the present calibration is revised;
- validity of the operational $\lambda_c$ initiation boundary;
- transmission/recrossing corrections if the rare-event TST approximation is not accurate enough;
- specimen-scale correlation volume/area and the mapping from local survival to specimen probability.

$A_c$ is especially important because it multiplies the activation barrier exponentially. It must not be chosen now merely to force an S-N curve.

## 12. Real-aluminum scope boundary

Published fatigue experiments on pure Al single crystals show strong involvement of slip bands, secondary slip, dislocation structures, lattice rotation, and subgrain evolution in actual crack initiation. Therefore the present closed law must be described as a **1D normal-instability hypothesis/submodel**, not as a proven complete microscopic theory of aluminum fatigue.

This distinction is scientifically useful: the reduced law is now closed and testable. If its temperature, frequency, mean-stress, and $P_0$ predictions fail against dedicated normal-loading experiments, the missing state is likely plastic/dislocation structure rather than another arbitrary probability closure.

## 13. Final verdict

The exact finite-chain probability projection alone cannot satisfy the requested $P_0$-only autonomous evolution without lost correlation information. However, after the laboratory-frequency quasistatic reduction and an explicit finite-temperature rare-first-passage assumption, the requested reduced law **can** be closed:

$$
\boxed{
\partial_tP_b
+\partial_\lambda
\left[
\frac{\dot\sigma(t)/E}{\phi''(\lambda)}P_b
\right]
=-
\frac{\sqrt{\phi''(\lambda)}}{2\pi t_0}
\exp\left[
-\frac{EA_ca_0\Delta\psi_c(\lambda)}{k_BT}
\right]P_b.
}
$$

Together with $P_b(\lambda,t_0)=P_0(\lambda)$ and the absorbing deterministic condition at $q_{\mathrm{tot}}\ge q_c$, this is the final reduced local initiation equation of the current 1D normal hypothesis.
