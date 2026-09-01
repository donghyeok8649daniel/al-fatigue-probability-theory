# Two-coordinate $(a,s)$ density identity — numerical verification

This report verifies the new theory-core identity

$$
\nabla\ln P
=
\boldsymbol\Theta^{-1}
\left[
\boldsymbol{\mathcal A}
-D_t\mathbf u
-\nabla\cdot\boldsymbol\Theta
\right]
$$

on a rectangular $(a,s)$ grid with a deliberately nonzero cross covariance $\Theta_{as}$.

This is a **synthetic algebra/numerics verification**, not a physical aluminum fatigue result. The target density and the energy-like surface used for the G2 quadrature check are analytic test fields. No Boltzmann distribution, Fokker--Planck dynamics, Gaussian/Weibull physical PDF, Markov bath, or state-independence assumption is used.

## Verification case

- grid: $181\times221$ points;
- non-diagonal conditional velocity covariance;
- conditional velocity correlation $\Theta_{as}/\sqrt{\Theta_{aa}\Theta_{ss}}$: 0.1616 to 0.1701;
- nonzero prescribed $D_tu_a$ and $D_tu_s$;
- conditional acceleration is constructed from the rearranged exact identity so that the target $P$ is an exact compatible solution before numerical differentiation/quadrature error.

## Results

| diagnostic | result |
|---|---:|
| max $|\partial_a g_s-\partial_s g_a|$ | $1.50\times10^{-10}$ |
| RMS compatibility curl | $4.93\times10^{-12}$ |
| rectangular-path log-density mismatch | $7.11\times10^{-15}$ |
| reconstructed-density $L^1$ error | $4.73\times10^{-6}$ |
| G2 test mean energy, target | 0.0820867796 |
| G2 test mean energy, reconstructed | 0.0820863366 |
| relative G2 mean-energy error | $5.40\times10^{-6}$ |

The nonzero $\Theta_{as}$ case therefore reconstructs the prescribed density to numerical quadrature accuracy and satisfies the two-dimensional integrability condition.

## What this proves

The tensor generalization of the old one-dimensional $\Theta$ shape identity is algebraically and numerically self-consistent:

$$
\boldsymbol\Theta\nabla\ln P
=
\boldsymbol{\mathcal A}-D_t\mathbf u-\nabla\cdot\boldsymbol\Theta.
$$

The off-diagonal covariance cannot in general be discarded. The numerical implementation handles it explicitly and recovers the correct density when $\Theta_{as}\ne0$.

The same reconstructed density can be inserted into G2,

$$
\bar U=\iint\Delta U_0P\,da\,ds,
$$

and the mean-energy quadrature is consistent to the same numerical accuracy in this verification case.

## What this does **not** prove

This does not yet show that a real single-crystal fatigue state has the synthetic $P$, $\Theta$, or $\mathcal A$ used here. A physical test requires a mechanics-generated ensemble of coupled opening/registry states $(a_\alpha,s_\alpha)$ from which

$$
\mathbf u(a,s,t),\qquad
\boldsymbol\Theta(a,s,t),\qquad
\boldsymbol{\mathcal A}(a,s,t)
$$

are measured without imposing a named distribution.

That physical step is currently blocked by one missing piece: the project has an exact intrinsic energy $U_0(a,s)$ and separate normal/registry numerical models, but it does not yet have a single deterministic spatially coupled $(a,s)$ mechanics model whose represented patches generate the required two-coordinate ensemble. Introducing such a model requires explicit choices for inertia, neighboring-patch constraints/couplings, boundary loading, and any irreversible mechanism. Those choices must be justified rather than inserted only to broaden the PDF.

## Next implementation target

Build the minimal deterministic coupled $(a,s)$ lattice/patch mechanics using the existing exact $U_0(a,s)$, then extract

1. empirical $P_M(a,s,t)$;
2. $u_a,u_s$;
3. $\Theta_{aa},\Theta_{as},\Theta_{ss}$;
4. $\mathcal A_a,\mathcal A_s$;
5. compatibility curl;
6. reconstructed $P$ vs empirical $P_M$;
7. G1 mean spacing and G2 mean intrinsic energy;
8. registry-well population for plasticity.

G3 remains separate until the irreversible mechanism is physically fixed.
