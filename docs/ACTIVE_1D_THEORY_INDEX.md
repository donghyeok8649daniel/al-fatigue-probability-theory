# Active 1D theory index

The authoritative mathematical mainline is now split into two layers: an exact microscopic reference layer and a reduced laboratory-time-scale initiation layer.

1. `../README_EQUATION_INDEX.md` — root equation/symbol navigation and notation policy;
2. `EQUATION_SUMMARY_1D_P_U_THETA.md` — compact exact microscopic governing-equation sheet;
3. `VARIABLE_INDEX_1D_P_U_THETA.md` — authoritative bilingual primary symbol/term dictionary;
4. `AUXILIARY_SYMBOL_INDEX_1D.md` — derivation-only and reduced-closure auxiliary symbols;
5. `MASTER_1D_P_U_THETA_FORMULATION.md` — full exact differential derivation;
6. `MILESTONE25_EXACT_INTEGRAL_REPRESENTATION.md` — exact full-flow, Volterra, characteristic, and first-passage integral representations;
7. `MILESTONE24_THETA_HISTORY_STATE.md` — history-bearing $P$-$u$-$\Theta$ interpretation;
8. `CRACK_INITIATION_DEFINITION.md` — kinetic first-passage initiation definition;
9. `FINAL_REDUCED_P0_THERMAL_FIRST_PASSAGE_CLOSURE.md` — active reduced lab-scale $P_0\to P_b,S,F_{\mathrm{ci}}$ closure;
10. `PERIODIC_P_FIRST_PASSAGE_SEPARATION.md` — exact separation between periodic PDF shape and cumulative first passage;
11. `theory/normal_lj_chain.py` — closed finite microscopic chain;
12. `theory/normal_lj_distribution_shape.py` — exact smooth density-shape reconstruction;
13. `theory/normal_lj_moment_hierarchy.py` — corrected exact $\Theta$ balance and spacing-coordinate kinetic metric.

## Mandatory notation rule

A new mathematical symbol is not part of the active theory unless `VARIABLE_INDEX_1D_P_U_THETA.md` or `AUXILIARY_SYMBOL_INDEX_1D.md` is updated simultaneously with:

$$
\boxed{
\text{equation definition}
+\text{English term}
+\text{Korean term}
+\text{mathematical meaning}
+\text{physical meaning}
+\text{unit/scaling}
+\text{status}
+\text{dependencies}
}
$$

The same glyph must not silently change meaning between files. If a mathematical defining relation exists, prose alone is insufficient.

GitHub Markdown math uses `$$ ... $$` for display math and `$...$` for inline math. Active files do not use `\[` / `\]` delimiters.

## Layer A — exact microscopic reference

$$
\boxed{
\text{finite LJ chain}
\to \Phi^q
\to F(\lambda,c,\tau)
\to\{P,u,\Theta\}
\to\{\bar a,\bar U,\text{first passage}\}
}
$$

The finite microscopic chain is closed. The reduced differential description obtained from it is exact but hierarchical. Given the finite-chain flow $\Phi^q$ and an initial full-state measure $\mu_0$, projected fields are obtained by exact push-forward integrals.

This layer remains the reference truth for testing any reduced closure. It is not, by itself, an autonomous $P_0$-only laboratory fatigue law.

## Layer B — active reduced laboratory closure

Laboratory-frequency audits show that the pure normal elastic dynamics is quasistatic relative to the atomic mechanical time scale. The active reduced model therefore uses the stable branch

$$
\phi'[\Lambda(\lambda_0,t)]
=\phi'(\lambda_0)+q(t)-q_{\mathrm{ref}}.
$$

Its characteristic velocity is

$$
\dot\Lambda=\frac{\dot q(t)}{\phi''(\Lambda)}.
$$

The nonabsorbing structural density therefore satisfies

$$
\boxed{
\partial_tP
+\partial_\lambda
\left[
\frac{\dot q(t)}{\phi''(\lambda)}P
\right]
=0.
}
$$

Finite-temperature rare first passage to the existing operational boundary $\lambda_c$ supplies the survivor sink. With characteristic cohesive area $A_c$,

$$
\Delta G_c(\lambda)
=EA_ca_0\Delta\psi_c(\lambda),
$$

$$
k_c(\lambda,T;A_c)
=\frac{\sqrt{\phi''(\lambda)}}{2\pi t_0}
\exp\left[-\frac{\Delta G_c(\lambda)}{k_BT}\right],
$$

and the active reduced survivor equation is

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

With $P_b(\lambda,t_0)=P_0(\lambda)$,

$$
S(t)=\int P_b(\lambda,t)d\lambda,
\qquad
F_{\mathrm{ci}}(t)=1-S(t).
$$

Thus the currently active laboratory-scale target is

$$
\boxed{
P_0+\sigma(0:t)
\to P_b(\lambda,t),S(t),F_{\mathrm{ci}}(t).
}
$$

This layer is a controlled reduced hypothesis. It assumes a structural/prestress $P_0$, quasistatic normal mechanics, fast thermal re-equilibration inside the intact well, and a rare-event regime. It does not claim to reconstruct the finite-chain neighbour ordering.

## Superseded or non-mainline routes

The following material may remain in the repository for historical, verification, or future-extension purposes but is not part of the active 1D governing model unless explicitly re-justified:

- named Gaussian/Weibull life or spacing distributions;
- harmonic/Taylor full-distribution closures;
- an unjustified Smoluchowski/Fokker--Planck mobility closure or an arbitrary diffusion coefficient;
- arbitrary Kramers rates or fitted fatigue kernels;
- registry/slip $s$ as a required coordinate for the normal-only mainline;
- FCC lattice geometry as a required basis of the present reduced chain;
- independent statistical-cell or FEM-element probability products.

The active thermal sink is not the rejected arbitrary-rate route: its energy barrier and attempt frequency are derived from the retained normal potential and inertial scale, while $A_c$ remains explicit for later calibration.

## Corrections that must be preserved

### 1. Active mean-rate symbol is $u$

Greek $\nu$ is reserved and is **not** the exact-chain conditional mean-rate field. The active microscopic definition is

$$
\boxed{
u(\lambda,\tau)=\mathbb E[c\mid\lambda,\tau]}
$$

The reduced thermal-closure attempt frequency is written $\nu_s$ and must not be confused with $u$.

### 2. General $\Theta$ equation

For the actual spatial LJ chain,

$$
\boxed{
D_\tau\Theta
+2\Theta\partial_\lambda u
+\frac1P\partial_\lambda(PC_3)
=2\Psi
}
$$

where

$$
\boxed{
\Psi(\lambda,\tau)=\operatorname{Cov}(c,\ddot\lambda\mid\lambda,\tau).
}
$$

The zero-right-hand-side equation is conditional on $\Psi=0$ and is not the general chain identity.

### 3. $\Theta$ is not the complete kinetic energy

$$
\boxed{
\mathbb E[c^2\mid\lambda,\tau]=u^2+\Theta
}
$$

is exact, but the chain kinetic energy in spacing coordinates is

$$
\boxed{
T^*=\frac12\boldsymbol c^T\mathbf G_\lambda\boldsymbol c,
\qquad
\mathbf G_\lambda=\mathbf L^T\mathbf L.
}
$$

Therefore exact total kinetic energy also requires cross-spacing rate correlations.

## Integral solution level

For deterministic finite-chain trajectories generated from $\Gamma_0$,

$$
P(\lambda,\tau)
=\frac1M\sum_i\int
\delta[\lambda-\Lambda_i(\tau;\Gamma_0)]\,\mu_0(d\Gamma_0).
$$

Analogous exact projection integrals define $u$, $\Theta$, $C_3$, and $\Psi$. Consequently,

$$
\boxed{
\text{non-closure of the exact reduced hierarchy}
\neq
\text{absence of an exact microscopic solution representation}.
}
$$

The new laboratory closure does not replace this identity. It is an asymptotic model that deliberately trades microscopic ordering information for the structural $P_0$ prestress embedding and a rare thermal escape law.

## Remaining open physics and calibration

1. The characteristic cohesive area $A_c$ is not yet calibrated.
2. The physical structural/prestress $P_0$ must be measured or constructed independently.
3. The operational $\lambda_c$ boundary and transition-state approximation require experimental or higher-fidelity validation.
4. Specimen-scale correlation area/volume and local-to-specimen survival scaling remain later calibration tasks.
5. Real pure-Al single-crystal fatigue is experimentally associated with slip-band/dislocation evolution, so the present normal-only closure is a testable normal-instability hypothesis rather than a complete dislocation theory.
6. Experimental validation of temperature, frequency, mean-stress, and initial-state predictions remains required.

These remaining tasks are calibration and model-scope questions, not an unresolved mathematical $P_0\to P_b,S,F_{\mathrm{ci}}$ closure problem.