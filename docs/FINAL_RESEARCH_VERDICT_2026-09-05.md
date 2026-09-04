# Final research verdict — 2026-09-05

Status: **ACTIVE FINAL SYNTHESIS.**

This document supersedes the 2026-09-04 endpoint as the current interpretation. It preserves the successful reduced normal-instability closure but adds the final cycle-frequency consistency result.

## 1. What is solved

Within the declared one-dimensional normal-instability hypothesis, the laboratory reduced mapping is closed from the structural initial spacing density and prescribed normal-stress history:

$$
P_0(\lambda)+\sigma(0:t)
\longrightarrow
P_b(\lambda,t),\ S(t),\ F_{\mathrm{ci}}(t).
$$

The reversible structural characteristic is

$$
\phi'[\Lambda(\lambda_0,t)]
=\phi'(\lambda_0)+q(t)-q_{\mathrm{ref}},
\qquad
\phi''(\Lambda)>0,
$$

and the survivor equation is

$$
\partial_tP_b
+\partial_\lambda
\left[
\frac{\dot q(t)}{\phi''(\lambda)}P_b
\right]
=-k_c(\lambda,T;A_c)P_b.
$$

The thermal sink is an explicit transition-state approximation derived from the retained normal potential:

$$
k_c
=\frac{\sqrt{\phi''(\lambda)}}{2\pi t_0}
\exp\left[-\frac{EA_ca_0\Delta\psi_c(\lambda)}{k_BT}\right].
$$

No named fatigue-life distribution, arbitrary diffusion coefficient, fitted scalar damage law, or imposed permanent normalized-$P$ drift is required.

## 2. What the exact microscopic layer proves

The finite conservative generalized-LJ chain remains the microscopic reference. Its probability projection is exact, but a one-point initial marginal $P_0$ does not retain neighbour ordering and higher correlations. Therefore it is not an autonomous exact $P_0$-only solver of the microscopic future.

This result justifies separating the exact microscopic reference from the reduced laboratory model.

## 3. Final pure-normal cycle-accumulation limitation

The strict pure-normal reduced state has a hard limitation.

If thermal renewal is absent and every intact structural label is periodic under an identical closed stress cycle, deterministic first passage cannot continue creating new labels indefinitely. The deterministic crossing set saturates after the first completed cycle.

If instead renewed crossing opportunity comes only from a fast stationary thermal bath, a fixed phase-shaped waveform gives

$$
\mathcal H_f
=\frac1f\int_0^1 k_c[\Lambda_*(\theta),T;A_c]d\theta,
$$

so

$$
\mathcal H_f\propto\frac1f.
$$

Thus median cycle count scales as $N_{50}\propto f$, while median elapsed time is approximately frequency independent in the strict quasistatic/fast-equilibration limit.

Therefore

$$
\boxed{
\text{pure normal}
+\text{quasistatic reversible mechanics}
+\text{no slow internal state}
}
$$

cannot by itself supply a genuinely cycle-controlled high-cycle fatigue mechanism.

## 4. External consistency check

A historical high-purity-aluminum study reported room-temperature fatigue strength to be insensitive to frequency over 25--1440 cycles/min. Those endpoints are approximately 0.4167 and 24 Hz. The strict fast-equilibrium normal model predicts a one-cycle hazard ratio of

$$
\frac{24}{0.4167}\approx57.6
$$

between the low- and high-frequency cases at the same phase-shaped stress waveform. This is a strong falsification warning, although the historical test is not an exact single-crystal match.

Room-temperature aluminum single-crystal studies further report secondary slip, persistent-slip-band-related extrusions/intrusions, irreversible slip, subgrain development, crystal-orientation dependence, and dislocation-density-related crack-initiation signatures. These observations identify cycle-evolving microstructure that the strict pure-normal state does not contain.

## 5. Scientific scope that survives

The current result is therefore not “the complete pure-Al fatigue theory is solved.” The defensible result is:

$$
\boxed{
\text{closed 1D normal-instability first-passage submodel}
+\text{no-go theorem for strict pure-normal cycle accumulation}.
}
$$

The normal submodel remains useful because it gives a mechanically defined local spacing transport, operational instability surface, characteristic energy climb, local survivor probability, and explicit falsification predictions.

## 6. What is still calibration versus what requires new physics

The following remain calibration/validation tasks inside the normal submodel:

1. characteristic cohesive area $A_c$;
2. physical structural/prestress $P_0$;
3. validation or revision of $\lambda_c$;
4. transition-state transmission/recrossing;
5. loading-axis material calibration;
6. specimen correlation area/volume and local-to-specimen scaling;
7. temperature, frequency, mean-stress, and amplitude validation.

A complete room-temperature pure-Al fatigue model requires more than calibration. It requires a physically derived cycle-evolving internal state that changes under repeated loading and feeds back on the normal-spacing state or its initiation barrier. Published aluminum evidence points most strongly to slip/dislocation/subgrain structure.

The user-level problem may still begin from $P_0+\sigma(0:t)$ if the added internal state's initial condition is fixed by the declared specimen preparation. But the constitutive state cannot in general remain the one-point normal-spacing marginal alone.

## 7. Final decision

If the project retains the strict no-shear/no-slip restriction, stop extending the closure and present the theory as the normal-instability submodel plus the no-go result.

If the objective remains a complete room-temperature pure-aluminum fatigue theory, the next research problem is to derive a cycle-evolving slip/dislocation state from mechanics and couple it back to the normal-spacing probability. It must not be introduced as an empirical scalar damage variable.
