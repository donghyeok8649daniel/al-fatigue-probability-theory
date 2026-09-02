# Active 1D theory index

The authoritative mathematical mainline is now:

1. `../README_EQUATION_INDEX.md` — root equation/symbol navigation and notation policy;
2. `EQUATION_SUMMARY_1D_P_U_THETA.md` — compact active governing-equation sheet;
3. `VARIABLE_INDEX_1D_P_U_THETA.md` — authoritative bilingual primary symbol/term dictionary with equation, mathematical, physical, unit, status, and dependency definitions;
4. `AUXILIARY_SYMBOL_INDEX_1D.md` — derivation-only auxiliary symbols with the same bilingual definition contract;
5. `MASTER_1D_P_U_THETA_FORMULATION.md` — full active differential derivation;
6. `MILESTONE25_EXACT_INTEGRAL_REPRESENTATION.md` — exact full-flow, Volterra, characteristic, and first-passage integral representations;
7. `MILESTONE24_THETA_HISTORY_STATE.md` — history-bearing $P$–$u$–$\Theta$ interpretation;
8. `CRACK_INITIATION_DEFINITION.md` — kinetic first-passage initiation definition;
9. `theory/normal_lj_chain.py` — closed finite microscopic chain;
10. `theory/normal_lj_distribution_shape.py` — exact smooth density-shape reconstruction;
11. `theory/normal_lj_moment_hierarchy.py` — corrected exact $\Theta$ balance and spacing-coordinate kinetic metric.

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

## Active chain

$$
\boxed{
\text{LJ chain mechanics}
\to \Phi^q
\to F(\lambda,c,\tau)
\to\{P,u,\Theta\}
\to\{\bar a,\bar U,\text{first passage}\}
}
$$

The finite microscopic chain is closed. The reduced differential description is exact but hierarchical. This does not prevent an exact integral solution representation: given the finite-chain flow $\Phi^q$ and an initial full-state measure $\mu_0$, all projected fields are obtained by push-forward integrals.

## Superseded or non-mainline routes

The following material may remain in the repository for historical, verification, or future-extension purposes but is not part of the active 1D paper-level governing model unless explicitly re-justified:

- Boltzmann/Gibbs initial distributions;
- named Gaussian/Weibull life or spacing distributions;
- harmonic/Taylor full-distribution closures;
- Smoluchowski/Fokker--Planck mobility closure and its $J^2/(MP)$ dissipation formula;
- registry/slip $s$ as a required coordinate for the normal-only mainline;
- FCC lattice geometry as a required basis of the present reduced chain;
- independent statistical-cell or FEM-element probability products.

## Corrections that must be preserved

### 1. Active mean-rate symbol is $u$

Greek $\nu$ is reserved and is **not** the active mean-rate field. The active definition is

$$
\boxed{
 u(\lambda,\tau)=\mathbb E[c\mid\lambda,\tau]
}
$$

Any legacy occurrence of LaTeX `\nu` where the conditional mean spacing rate is intended is a typography error and must be corrected to $u$.

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
\Psi(\lambda,\tau)=\operatorname{Cov}(c,\ddot\lambda\mid\lambda,\tau)
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
\mathbf G_\lambda=\mathbf L^T\mathbf L
}
$$

Therefore exact total kinetic energy also requires cross-spacing rate correlations.

## Integral solution level

For deterministic finite-chain trajectories generated from $\Gamma_0$,

$$
\boxed{
P(\lambda,\tau)
=\frac1M\sum_i\int
\delta[\lambda-\Lambda_i(\tau;\Gamma_0)]\,\mu_0(d\Gamma_0)
}
$$

and analogous exact projection integrals define $u$, $\Theta$, $C_3$, and $\Psi$. Consequently,

$$
\boxed{
\text{non-closure of a three-field PDE}
\neq
\text{absence of an exact integral solution representation}
}
$$

The specimen survival formula exists once an initial realization measure is declared:

$$
\boxed{
S_{\rm spec}(\tau)
=\int
\mathbf 1\left[
\max_i\sup_{s\in[\tau_0,\tau]}\Lambda_i(s;\Gamma_0)<\lambda_c
\right]\mu_0(d\Gamma_0)
}
$$

Thus the specimen-scale open issue is primarily the physical construction and validation of $\mu_0$ and the represented correlation scale, not the absence of an integral survival equation.

## Remaining open physics

1. G3 irreversibility is not generated by the conservative baseline.
2. Microscopic same-force history dependence has not yet been bridged to laboratory Hz-scale fatigue memory.
3. The physical initial/specimen measure $\mu_0$ and correlation scale have not yet been derived/calibrated.
4. Experimental validation remains required.

These are physical open problems, not reasons to replace the active $P$–$u$–$\Theta$ mathematics with an assumed PDF or empirical damage variable.
