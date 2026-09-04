# Final reduced closure verdict

This report accompanies `simulations/audit_final_reduced_closure_sensitivity.py`.

## Result

The reduced laboratory model is closed by combining:

1. quasistatic stable-branch transport of the structural spacing distribution;
2. a finite-temperature positive-flux transition-state approximation to the declared normal-instability boundary;
3. survivor mass loss as the cumulative first-passage memory.

On the intact domain $\lambda<\lambda_c$,

$$
\partial_tP_b
+\partial_\lambda
\left[
\frac{\dot\sigma/E}{\phi''(\lambda)}P_b
\right]
=-k_c(\lambda,T;A_c)P_b,
$$

with

$$
k_c
=\frac{\sqrt{\phi''(\lambda)}}{2\pi t_0}
\exp\left[-\frac{EA_ca_0\Delta\psi_c(\lambda)}{k_BT}\right].
$$

The rate expression is used only under the declared high-barrier, harmonic-well, fast-intrawell-equilibration approximation and away from the immediate spinodal neighbourhood. If the quasistatic stable branch itself reaches $\lambda_c$, the local state is absorbed deterministically rather than extrapolating the harmonic prefactor through $\phi''\to0$. A dynamical transmission/recrossing correction may be required by higher-fidelity calculations.

No value of $A_c$ is selected here.

## Sensitivity

For the historical Al bridge and a 300 K, 20 Hz, $100\pm100$ MPa protocol, the effective-potential climb associated with one reference atomic area is only about 0.02 eV. The event is therefore not rare at that scale. Increasing the coherent area multiplies the energy climb and changes the local hazard exponentially.

The sensitivity sweep gives approximately:

| $A_c/A_0$ | escape per cycle | local median cycles |
|---:|---:|---:|
| 30 | 0.99995 | 0.070 |
| 40 | $4.66\times10^{-3}$ | $1.48\times10^2$ |
| 50 | $2.30\times10^{-6}$ | $3.02\times10^5$ |
| 60 | $1.16\times10^{-9}$ | $5.97\times10^8$ |

These values are not a fit and do not identify the physical characteristic area. Rows that do not satisfy the rare-event assumption are **failure-of-assumption diagnostics**, not quantitative fatigue-life predictions. The rare-event rows demonstrate the desired high-cycle structure: one local event can have survival extremely close to one per cycle while cumulative first passage remains finite over many cycles.

## Rejected alternatives summarized

- Exact finite-chain projection: reference truth, but not autonomous in $P_0$.
- Conservative local oscillator: closes $P_0\to P$, but becomes quasistatic and non-accumulating at laboratory frequency.
- Collective normal modes: too fast for a small characteristic domain.
- Permanent normalized-$P$ drift: not mathematically required.
- Spatial registry kink: rare transient candidate, but not established as a long-lived memory or as the same-potential normal-instability trigger.
- Arbitrary diffusion/damping: not required in the final reduced equation.

## Interpretation

The successful mathematical endpoint is not a complete experimentally validated theory of real Al fatigue. It is a closed and falsifiable **1D normal-instability initiation submodel under a declared transition-state approximation**. The next unknowns are calibration/validation quantities: $A_c$, physical $P_0$, transmission/recrossing correction, the operational initiation boundary, and specimen-scale correlation structure.
