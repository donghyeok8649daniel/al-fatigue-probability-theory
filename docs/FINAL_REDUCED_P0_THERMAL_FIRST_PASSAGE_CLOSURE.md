# Final reduced P0-to-survival closure / 최종 P0-생존 축약법칙

Status: **ACTIVE REDUCED LAB-SCALE CLOSURE for the declared 1D normal-instability hypothesis.**

The exact finite LJ chain and its $P$-$u$-$\Theta$ projection remain the microscopic reference model. This document supplies a separate controlled laboratory-time-scale reduction. It is not claimed to be a complete dislocation/slip theory of real single-crystal aluminum.

## 1. Target

The reduced laboratory model must predict

$$
\boxed{
P_0(\lambda)+\sigma(0:t)
\longrightarrow
P_b(\lambda,t),\ S(t),\ F_{\mathrm{ci}}(t)
}
$$

without reconstructing every microscopic trajectory and without prescribing a named PDF or empirical fatigue-damage variable.

Previous audits established that the exact finite conservative chain is not autonomous in one-point $P_0$, laboratory-frequency normal elastic motion is quasistatic on the atomic time scale, and permanent drift of a normalized spacing PDF is not mathematically required if survivor mass is removed by first passage.

## 2. Structural P0 embedding

$P_0(\lambda)$ is the reference-phase **structural/prestress spacing density** on the intact stable branch, not an instantaneous thermal-displacement distribution. Let

$$
q_{\mathrm{ref}}=\frac{\sigma_{\mathrm{ref}}}{E},
\qquad
q(t)=\frac{\sigma(t)}{E}.
$$

For starting label $\lambda_0$, define

$$
\boxed{
q_r(\lambda_0)=\phi'(\lambda_0)-q_{\mathrm{ref}}.
}
$$

This is the explicit $P_0$-only embedding hypothesis: $\lambda_0$ is a local stable equilibrium at the reference phase. It does not reconstruct the missing finite-chain neighbour ordering.

## 3. Quasistatic stable-branch map

The intact reduced characteristic satisfies

$$
\boxed{
\phi'[\Lambda(\lambda_0,t)]
=\phi'(\lambda_0)+q(t)-q_{\mathrm{ref}},
\qquad
\phi''[\Lambda(\lambda_0,t)]>0.
}
$$

At the declared initial/reference phase,

$$
\boxed{
\Lambda(\lambda_0,0)=\lambda_0.
}
$$

Differentiating gives

$$
\boxed{
\dot\Lambda=\frac{\dot q(t)}{\phi''(\Lambda)}.
}
$$

Therefore the nonabsorbing structural density obeys

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

The mechanical push-forward is

$$
P(\lambda,t)
=\int P_0(\lambda_0)
\delta[\lambda-\Lambda(\lambda_0,t)]d\lambda_0.
$$

Where the inverse map exists,

$$
\boxed{
P(\lambda,t)
=P_0(\lambda_0)
\frac{\phi''(\lambda)}{\phi''(\lambda_0)}.
}
$$

A closed load cycle returns this nonabsorbing structural map to the same reference-phase shape if the stable branch is never lost.

## 4. Operational instability boundary

The retained generalized-LJ tangent stiffness vanishes at

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

The maximum stable reduced traction is

$$
\boxed{q_c=\phi'(\lambda_c).}
$$

If the local total reduced traction reaches or exceeds $q_c$, the stable intact root disappears and that characteristic is absorbed deterministically. The use of $\lambda_c$ is operational and must ultimately be checked against an experimental initiation definition.

## 5. Characteristic cohesive area

Let $A_c$ be the characteristic area participating coherently in one local normal-instability event. It is not calibrated here and is not a FEM element area or specimen independent-cell area. The event-energy scaling

$$
\Delta G_c\propto A_c
$$

is an explicit coherent-area hypothesis rather than a hidden fitted fatigue parameter.

With the same reduced coordinate, take

$$
m_c=\frac{A_c}{A_0}m_a.
$$

Because stiffness and mass both scale with $A_c$, the small-oscillation attempt frequency is

$$
\boxed{
\nu_s(\lambda)=\frac{\sqrt{\phi''(\lambda)}}{2\pi t_0},
\qquad
t_0=\sqrt{\frac{m_a a_0}{EA_0}}.
}
$$

## 6. Rare thermal first passage

For stable $\lambda<\lambda_c$, define the effective-potential climb to the operational dividing surface at the same reduced traction,

$$
\Delta\psi_c(\lambda)
=
\left[
\phi(\lambda_c)-\phi'(\lambda)\lambda_c
\right]
-
\left[
\phi(\lambda)-\phi'(\lambda)\lambda
\right].
$$

Then

$$
\boxed{
\Delta G_c(\lambda)=EA_ca_0\Delta\psi_c(\lambda).
}
$$

Under the declared high-barrier, harmonic-intrawell, fast-re-equilibration assumptions, and **away from the immediate spinodal neighbourhood**, the positive-flux transition-state approximation is

$$
\boxed{
k_c(\lambda,T;A_c)
=\frac{\sqrt{\phi''(\lambda)}}{2\pi t_0}
\exp\left[
-\frac{EA_ca_0\Delta\psi_c(\lambda)}{k_BT}
\right].
}
$$

This is not a fitted fatigue kernel and not an exact stochastic reduction. A transmission/recrossing factor may be required by higher-fidelity dynamics. The harmonic prefactor is not extrapolated through $\phi''\to0$; deterministic absorption is used when the stable characteristic reaches $\lambda_c$.

## 7. Closed survivor equation

The intact survivor subdensity satisfies

$$
\boxed{
\partial_tP_b
+\partial_\lambda
\left[
\frac{\dot q(t)}{\phi''(\lambda)}P_b
\right]
=-k_c(\lambda,T;A_c)P_b,
\qquad \lambda<\lambda_c.
}
$$

with

$$
\boxed{P_b(\lambda,0)=P_0(\lambda).}
$$

Along a stable characteristic,

$$
W(\lambda_0,t)
=\exp\left[
-\int_0^t k_c(\Lambda(\lambda_0,s),T;A_c)ds
\right].
$$

Hence, before deterministic absorption of that characteristic,

$$
\boxed{
P_b(\lambda,t)
=\int P_0(\lambda_0)
W(\lambda_0,t)
\delta[\lambda-\Lambda(\lambda_0,t)]d\lambda_0.
}
$$

The local survival and cumulative initiation probability are

$$
\boxed{
S(t)=\int_{0}^{\lambda_c}P_b(\lambda,t)d\lambda,
\qquad
F_{\mathrm{ci}}(t)=1-S(t).
}
$$

This completes the requested reduced mapping from $P_0$ and the stress history.

## 8. Periodic loading

For a periodic stress waveform of period $T_f=1/f$, each structural label that remains mechanically stable has

$$
\mathcal H_c(\lambda_0)
=\int_0^{T_f}k_c[\Lambda(\lambda_0,t),T;A_c]dt.
$$

After $N$ identical cycles,

$$
\boxed{
S_N
=\int_{\mathrm{stable}} P_0(\lambda_0)
\exp[-N\mathcal H_c(\lambda_0)]d\lambda_0.
}
$$

Thus the reversible mechanical map may return after every cycle while survivor mass decreases. A heterogeneous survivor population changes by selection rather than by an imposed permanent distortion of every surviving spacing.

In the strict quasistatic phase-controlled limit,

$$
\mathcal H_c\propto\frac1f.
$$

The model therefore predicts approximately frequency-independent local life in physical time and a cycle count proportional to $f$ until the quasistatic or fast-equilibration assumptions fail. This is a falsification signature.

## 9. Why earlier routes are not the mainline

- **Exact finite-chain mixing:** reference truth, but not autonomous in one-point $P_0$.
- **Local conservative oscillator:** gives reversible $P_0\to P(t)$ but no laboratory-Hz cycle accumulation.
- **Collective normal elastic mode:** too fast in a local characteristic region.
- **Permanent normalized-P drift:** not mathematically required for cumulative first passage.
- **Registry-kink storage:** rare transient states exist, but long-lived trapping was not established in the ideal pure-normal surface.
- **Arbitrary diffusion/Kramers fitting:** rejected; no fitted diffusion coefficient, damping constant, named PDF, or empirical life kernel is used in the active equation.

## 10. What remains for calibration and falsification

The following remain open by design:

1. characteristic cohesive area $A_c$;
2. physical structural/prestress $P_0$;
3. validation or revision of the operational $\lambda_c$ boundary;
4. transmission/recrossing correction to the positive-flux TST approximation;
5. loading-axis calibration $E,a_0,A_0,m,n$ where required;
6. specimen-scale correlation area/volume and local-to-specimen survival mapping;
7. comparison with actual pure-Al single-crystal fatigue, where slip-band and dislocation evolution may dominate.

The peak-hazard asymptotic and temperature-slope relation provide a future route to identify $A_c$ from local-hazard data rather than by forcing an S-N fit.

## 11. Final verdict

The exact finite-chain projection alone cannot satisfy an autonomous $P_0$-only high-cycle law because projection discards ordering/correlation information. After the laboratory-frequency quasistatic reduction and the explicitly declared finite-temperature transition-state hypothesis, however, the local survivor problem is mathematically closed:

$$
\boxed{
P_0+\sigma(0:t)
\to P_b(\lambda,t)
\to S(t)
\to F_{\mathrm{ci}}(t).
}
$$

This is the successful endpoint of the present closure search. The remaining work is parameter identification, transmission testing, specimen-scale scaling, and experimental falsification—not invention of another arbitrary probability closure.
