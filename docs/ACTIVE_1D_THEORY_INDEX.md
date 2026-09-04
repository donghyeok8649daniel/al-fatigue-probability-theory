# Active 1D theory index

The authoritative mainline is split into two layers: an exact microscopic reference and a reduced laboratory-time-scale initiation model.

1. `../README_EQUATION_INDEX.md` — equation/symbol navigation;
2. `EQUATION_SUMMARY_1D_P_U_THETA.md` — exact microscopic governing equations;
3. `VARIABLE_INDEX_1D_P_U_THETA.md` — bilingual primary symbols;
4. `AUXILIARY_SYMBOL_INDEX_1D.md` — derivation/reduced-closure symbols;
5. `MASTER_1D_P_U_THETA_FORMULATION.md` — exact differential derivation;
6. `MILESTONE25_EXACT_INTEGRAL_REPRESENTATION.md` — exact push-forward/integral representations;
7. `CRACK_INITIATION_DEFINITION.md` — first-passage initiation definition;
8. `FINAL_REDUCED_P0_THERMAL_FIRST_PASSAGE_CLOSURE.md` — active reduced lab-scale closure;
9. `PERIODIC_P_FIRST_PASSAGE_SEPARATION.md` — periodic PDF versus cumulative first passage;
10. `PEAK_HAZARD_ASYMPTOTICS_AND_IDENTIFIABILITY.md` — tensile-peak hazard asymptotic and future $A_c$ identification route;
11. `EXPERIMENTAL_FALSIFICATION_PLAN_REDUCED_CLOSURE.md` — ordered validation/rejection tests.

## Layer A — exact microscopic reference

$$
\text{finite LJ chain}
\to \Phi^q
\to F(\lambda,c,\tau)
\to \{P,u,\Theta\}
$$

The finite microscopic chain is closed. Its projected one-point hierarchy is exact but not autonomous because neighbour ordering and higher correlations are discarded. Given the full initial-state measure $\mu_0$, the projected probability is obtained exactly by push-forward.

The active mean-rate symbol is

$$
u(\lambda,\tau)=\mathrm E[c\mid\lambda,\tau].
$$

The general variance balance is

$$
D_\tau\Theta
+2\Theta\partial_\lambda u
+\frac1P\partial_\lambda(PC_3)
=2\Psi,
$$

with

$$
\Psi(\lambda,\tau)=\mathrm{Cov}(c,\ddot\lambda\mid\lambda,\tau).
$$

No zero closure is adopted without additional justification.

## Layer B — active laboratory closure

$P_0(\lambda)$ is a structural/prestress spacing density at the declared reference phase. For each starting label,

$$
q_r(\lambda_0)=\phi'(\lambda_0)-q_{\mathrm{ref}},
$$

and the intact quasistatic characteristic satisfies

$$
\phi'[\Lambda(\lambda_0,t)]
=\phi'(\lambda_0)+q(t)-q_{\mathrm{ref}},
\qquad
\phi''(\Lambda)>0,
$$

with

$$
\Lambda(\lambda_0,0)=\lambda_0.
$$

Therefore

$$
\dot\Lambda=\frac{\dot q(t)}{\phi''(\Lambda)},
$$

and the reversible structural density obeys

$$
\partial_tP
+\partial_\lambda
\left[
\frac{\dot q(t)}{\phi''(\lambda)}P
\right]
=0.
$$

The operational instability boundary is

$$
\phi''(\lambda_c)=0,
\qquad
\lambda_c=\left(\frac{m+1}{n+1}\right)^{1/(m-n)}.
$$

Let $A_c$ be the coherent local event area. For stable $\lambda<\lambda_c$,

$$
\Delta G_c(\lambda)
=EA_ca_0\Delta\psi_c(\lambda).
$$

Away from the immediate spinodal neighbourhood, the declared high-barrier/harmonic-well/fast-equilibration transition-state approximation is

$$
k_c(\lambda,T;A_c)
=\frac{\sqrt{\phi''(\lambda)}}{2\pi t_0}
\exp\left[-\frac{\Delta G_c(\lambda)}{k_BT}\right].
$$

This rate is not exact. A transmission/recrossing correction may be required. If the stable characteristic reaches $\lambda_c$, that characteristic is absorbed deterministically rather than extrapolating the harmonic prefactor through $\phi''\to0$.

The active survivor equation is

$$
\partial_tP_b
+\partial_\lambda
\left[
\frac{\dot q(t)}{\phi''(\lambda)}P_b
\right]
=-k_c(\lambda,T;A_c)P_b,
\qquad \lambda<\lambda_c,
$$

with

$$
P_b(\lambda,0)=P_0(\lambda).
$$

Hence

$$
S(t)=\int_0^{\lambda_c}P_b(\lambda,t)d\lambda,
\qquad
F_{\mathrm{ci}}(t)=1-S(t),
$$

and the active target is

$$
P_0+\sigma(0:t)
\to P_b(\lambda,t),S(t),F_{\mathrm{ci}}(t).
$$

For periodic loading,

$$
\mathcal H_c(\lambda_0)
=\int_0^{T_f}k_c[\Lambda(\lambda_0,t),T;A_c]dt,
$$

and

$$
S_N
=\int_{\mathrm{stable}}P_0(\lambda_0)
\exp[-N\mathcal H_c(\lambda_0)]d\lambda_0.
$$

Thus permanent drift of the normalized structural PDF is not required; cumulative history is carried by survivor loss and selection.

## Superseded or non-mainline routes

The following are not active governing assumptions unless separately re-justified:

- named Gaussian/Weibull life or spacing distributions;
- Taylor/harmonic full-distribution closures;
- arbitrary diffusion or mobility closures;
- fitted Kramers/fatigue kernels;
- registry/slip $s$ as a required coordinate for the normal-only mainline;
- FCC reconstruction;
- independent statistical-cell/FEM-element probability products.

Registry/kink calculations remain a plasticity/defect extension; long-lived trapping and direct crack triggering were not established in the ideal pure-normal reduced surface.

## Remaining calibration and falsification

1. Determine $A_c$ independently.
2. Construct or measure structural $P_0$.
3. Validate or revise $\lambda_c$ as the operational dividing surface.
4. Quantify transition-state transmission/recrossing.
5. Validate temperature, frequency, mean-stress, amplitude, and initial-state predictions.
6. Determine specimen correlation area/volume and local-to-specimen survival scaling later.
7. Test whether real pure-Al single-crystal fatigue is governed by this normal route or by slip-band/dislocation evolution.

The closure problem is therefore solved only **within the declared 1D normal-instability hypothesis and transition-state approximation**. The next stage is calibration and falsification, not another arbitrary probability closure.
