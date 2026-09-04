# Final research verdict — reduced 1D fatigue-initiation theory

Status: **ACTIVE SYNTHESIS / current final mainline.**

This document records the endpoint of the 2026-09-04 closure search. It does not claim that the 1D normal-instability hypothesis is already experimentally validated for real single-crystal Al. It records which mathematical routes failed, which reduction survived, and what must be calibrated later.

## 1. Required mapping

The laboratory-scale theory must start from the structural spacing density and the applied normal-stress history,

$$
P_0(\lambda)+\sigma(0:t),
$$

and predict local intact survival without reconstructing all atomic labels,

$$
\boxed{
P_0+\sigma(0:t)
\longrightarrow
P_b(\lambda,t),\ S(t),\ F_{\mathrm{ci}}(t).
}
$$

Here

$$
S(t)=\int_0^{\lambda_c}P_b(\lambda,t)d\lambda,
\qquad
F_{\mathrm{ci}}(t)=1-S(t).
$$

Permanent cycle-to-cycle drift of a normalized PDF is not required. What must accumulate is first-passage loss from the intact population.

## 2. Results that are now ruled out as the main fatigue clock

### Exact finite conservative chain as a P0-only solver

The finite generalized-LJ chain is the microscopic reference truth, but the one-point initial marginal $P_0$ loses ordering and correlation information. Different microscopic arrangements can share the same $P_0$ and even the same one-point phase-space density while producing different future marginals. The exact projected hierarchy therefore does not close autonomously in $P_0$.

### Laboratory-frequency conservative normal oscillation

The atomic and elastic normal time scales are far faster than ordinary fatigue loading. In the laboratory-frequency limit the pure normal conservative response approaches a label-preserving quasistatic cycle. Repeating that cycle does not create new first-passage labels indefinitely.

### Collective normal elastic modes

Reducing the lowest longitudinal elastic mode to tens of hertz would require macroscopic lengths far beyond the local characteristic region. A small characteristic domain raises, rather than lowers, the elastic-mode frequency. Collective normal elasticity is therefore not the missing slow clock.

### Permanent P-shape drift as a mathematical requirement

It is not required. A periodic survivor-conditioned shape can coexist with decreasing survivor mass. Conversely, a periodic snapshot marginal alone does not determine cumulative first passage because it does not record which labels have already crossed.

### Registry-kink storage as the primary long-lived memory

Spatial registry resolution produced a real rare formation-energy scale, but the retained kink-pair branch is extremely shallow in the separation/migration direction. It is not established as a long-lived room-temperature memory state under the present ideal pure-normal model.

A further coupled check removes an earlier over-interpretation: the approximately 14 percent relaxed normal opening near the registry saddle must not be compared directly with the normal-only chain value $\lambda_c$. In the same $U_0(a,s)$ surface, the registry-dependent relaxed spacing remains below that same surface's normal curvature-loss point. Thus the current kink candidate is not promoted as the crack-initiation trigger.

## 3. Surviving mechanical reduction

The structural $P_0$ is interpreted as a reference-phase prestress distribution, not an instantaneous thermal-displacement distribution. For each initial label $\lambda_0$, define

$$
q_r(\lambda_0)=\phi'(\lambda_0)-q_{\mathrm{ref}},
\qquad
q(t)=\frac{\sigma(t)}{E}.
$$

At laboratory frequency, intact normal mechanics follows the stable quasistatic branch

$$
\boxed{
\phi'[\Lambda(\lambda_0,t)]
=\phi'(\lambda_0)+q(t)-q_{\mathrm{ref}},
\qquad
\phi''(\Lambda)>0.
}
$$

Differentiation gives

$$
\boxed{
\dot\Lambda
=\frac{\dot q(t)}{\phi''(\Lambda)}.
}
$$

Therefore the reversible nonabsorbing structural density is closed directly from $P_0$:

$$
\boxed{
\partial_tP
+\partial_\lambda
\left[
\frac{\dot q(t)}{\phi''(\lambda)}P
\right]=0.
}
$$

This solves the requested $P_0+\sigma(0:t)\to P(t)$ mechanical mapping under the declared local-prestress/quasistatic reduction.

## 4. Minimal additional physics required for high-cycle first passage

A perfectly label-preserving quasistatic cycle cannot generate indefinitely renewed first-passage opportunities. Some non-deterministic fast degrees of freedom must therefore remain after the structural coordinate is reduced.

The current minimum closure is **finite-temperature transition-state first passage**, not an arbitrary diffusion law. Fast thermal degrees of freedom are assumed to re-equilibrate inside the intact well on a time scale much shorter than a fatigue period and the rare escape time. This separation is physically plausible for Al: measured and ab-initio studies report phonon and electron-phonon relaxation on microscopic time scales many orders of magnitude below ordinary mechanical-fatigue periods, but the effective reduced-coordinate rate still requires validation.

The operational normal-instability boundary remains

$$
\phi''(\lambda_c)=0,
\qquad
\lambda_c=\left(\frac{m+1}{n+1}\right)^{1/(m-n)}.
$$

For a stable current spacing $\lambda$, define the effective-potential climb

$$
\Delta\psi_c(\lambda)
=
\left[\phi(\lambda_c)-\phi'(\lambda)\lambda_c\right]
-
\left[\phi(\lambda)-\phi'(\lambda)\lambda\right].
$$

Let $A_c$ be the characteristic cohesive area of one local event. It remains symbolic and must be calibrated later. The event barrier is

$$
\boxed{
\Delta G_c(\lambda)=EA_ca_0\Delta\psi_c(\lambda).
}
$$

If the coherent effective mass scales with the same area, the harmonic attempt frequency is

$$
\boxed{
\nu_s(\lambda)=\frac{\sqrt{\phi''(\lambda)}}{2\pi t_0}.
}
$$

The rare positive-flux transition-state rate is therefore

$$
\boxed{
k_c(\lambda,T;A_c)
=\nu_s(\lambda)
\exp\left[-\frac{EA_ca_0\Delta\psi_c(\lambda)}{k_BT}\right].
}
$$

This is a controlled approximation. A transmission/recrossing coefficient may be required by higher-fidelity dynamics; it is not silently set by an S-N fit.

## 5. Final reduced survivor equation

The active lab-scale closure is

$$
\boxed{
\partial_tP_b
+\partial_\lambda
\left[
\frac{\dot\sigma(t)/E}{\phi''(\lambda)}P_b
\right]
=-k_c(\lambda,T;A_c)P_b,
}
$$

with

$$
\boxed{P_b(\lambda,t_0)=P_0(\lambda)}
$$

for an initially intact reference population.

Along a mechanical characteristic,

$$
W(\lambda_0,t)
=\exp\left[-\int_{t_0}^{t}k_c(\Lambda(\lambda_0,s),T;A_c)ds\right],
$$

and therefore

$$
\boxed{
S(t)=\int P_0(\lambda_0)W(\lambda_0,t)d\lambda_0,
\qquad
F_{\mathrm{ci}}(t)=1-S(t).
}
$$

The theory is thus closed from the requested initial PDF and load history once $A_c$, material parameters, temperature and the operational threshold are supplied.

## 6. Periodic loading and high local survival

For a stress period $T_f$,

$$
\mathcal H_c(\lambda_0)
=\int_0^{T_f}k_c[\Lambda(\lambda_0,t),T;A_c]dt.
$$

Hence

$$
\boxed{
S_N
=\int P_0(\lambda_0)
\exp[-N\mathcal H_c(\lambda_0)]d\lambda_0.
}
$$

For a narrow initial state this reduces to

$$
S_N\simeq e^{-N\mathcal H_c}.
$$

Therefore a local survival probability extremely close to unity per cycle is not a problem; it is exactly the expected high-cycle regime. Spatial multiplication over many characteristic regions is a later specimen-scale calibration and is not inserted here.

## 7. Sensitivity audit, not calibration

Using the current historical Al bridge $E=69$ GPa, $a_0=2.8627442948\times10^{-10}$ m, $A_0=6.0338\times10^{-20}$ m$^2$, $T=300$ K and a $100\pm100$ MPa, 20 Hz sinusoid, the atomic-reference barrier to $\lambda_c$ is only about $0.0216$ eV at zero stress and $0.0193$ eV at 200 MPa. Therefore an atomic-area event would not be a rare high-cycle event. This is not a failure of the survival equation; it proves that the characteristic cohesive area cannot be identified with one atomic reference area without validation.

Writing $A_c=R A_0$, the same diagnostic gives an extreme exponential sensitivity. Representative results are:

| $R=A_c/A_0$ | one-cycle integrated hazard | one-cycle escape probability | local median cycles |
|---:|---:|---:|---:|
| 30 | $9.89$ | $0.99995$ | $7.0\times10^{-2}$ |
| 40 | $4.67\times10^{-3}$ | $4.66\times10^{-3}$ | $1.48\times10^2$ |
| 50 | $2.30\times10^{-6}$ | $2.30\times10^{-6}$ | $3.02\times10^5$ |
| 60 | $1.16\times10^{-9}$ | $1.16\times10^{-9}$ | $5.97\times10^8$ |

These values are a **sensitivity sweep only**. No value of $R$ is adopted. The sweep quantifies why $A_c$ must be determined independently rather than hidden in a fitted fatigue law.

## 8. What is successful, and what remains unproven

### Successful mathematical closure

The requested autonomous reduced mapping is now explicit:

$$
\boxed{
P_0+\sigma(0:t)
\to P_b(\lambda,t)
\to S(t)
\to F_{\mathrm{ci}}(t).
}
$$

It requires no full microscopic ordering, no named spacing PDF, no Weibull life law, no scalar damage variable, no permanent normalized-P drift, and no arbitrary diffusion coefficient.

### Remaining physical validation

The following are not solved by algebra and must remain open until independent data or higher-fidelity calculations exist:

1. $A_c$, the characteristic cohesive area of one local rare event;
2. physical construction/measurement of structural $P_0$;
3. transmission/recrossing correction to transition-state flux;
4. experimental validity of $\lambda_c$ as the local crack-initiation dividing surface;
5. specimen correlation area/volume and local-to-specimen survival scaling;
6. whether real pure-Al single-crystal fatigue under the target loading is dominated by this normal-instability route or by dislocation/slip-band evolution.

## 9. Final scientific conclusion

The conservative finite-chain projection is indispensable as a microscopic reference but cannot, by itself, provide a $P_0$-only high-cycle fatigue law. The laboratory-frequency mechanical part is quasistatic and reversible. High-cycle accumulation must therefore be placed in **survivor first passage**, not in an invented irreversible drift of the normalized spacing PDF.

Under the explicitly declared finite-temperature rare-event assumptions, the current generalized-LJ mechanics supplies both the stable-branch transport and the activation barrier. The resulting survivor equation is a closed, falsifiable 1D normal-instability fatigue-initiation submodel.

This is the current successful endpoint of the reduced-theory search. The next stage is no longer closure invention; it is independent calibration and falsification of $A_c$, $P_0$, the transition-state transmission factor, and the operational initiation boundary.