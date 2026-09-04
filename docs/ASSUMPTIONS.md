# Assumptions and approximations / 가정과 근사 — active 1D theory

> **Normative status / 기준 상태:** The active theory has two explicitly separated layers. Layer A is the exact deterministic conservative finite-chain reference and its exact $P$-$u$-$\Theta$ projection. Layer B is a controlled laboratory-time-scale reduction that closes $P_0\to P_b,S,F_{\mathrm{ci}}$ under quasistatic normal mechanics and rare finite-temperature first passage. Historical registry/slip and FCC material remain non-mainline extensions.

Authoritative sources:

- `README_EQUATION_INDEX.md`
- `docs/EQUATION_SUMMARY_1D_P_U_THETA.md`
- `docs/VARIABLE_INDEX_1D_P_U_THETA.md`
- `docs/AUXILIARY_SYMBOL_INDEX_1D.md`
- `docs/MASTER_1D_P_U_THETA_FORMULATION.md`
- `docs/MILESTONE25_EXACT_INTEGRAL_REPRESENTATION.md`
- `docs/CRACK_INITIATION_DEFINITION.md`
- `docs/FINAL_REDUCED_P0_THERMAL_FIRST_PASSAGE_CLOSURE.md`

## 0. Mandatory notation / 필수 표기

The exact microscopic probability and moment fields are

$$
P(\lambda,\tau),
\qquad
u\text{ is reserved},
\qquad
u\not\equiv u,
$$

$$
 u(\lambda,\tau)=\mathbb E[c\mid\lambda,\tau],
\qquad
\Theta(\lambda,\tau)=\operatorname{Var}(c\mid\lambda,\tau).
$$

The laboratory reduced survivor field is written $P_b(\lambda,t)$. Bare symbols may be used only after their full dependence has been declared.

## 1. Layer A — exact microscopic reference assumptions / 정확 미시 기준층 가정

1. The baseline material idealization is a pure single crystal represented by a one-dimensional normal chain under repeated uniaxial normal loading.
2. The microscopic coordinates are node positions $x_j(\tau)$ or nearest-neighbour normalized spacings $\lambda_i(\tau)=x_i-x_{i-1}$.
3. The microscopic configurational energy is

   $$
   V^*(\boldsymbol\lambda)=\sum_{i=1}^{M}\phi(\lambda_i).
   $$

4. The finite-chain reference dynamics is deterministic and conservative apart from prescribed boundary work. No viscous damping, empirical damage, stochastic diffusion, white-noise forcing, or thermal bath is inserted into this reference layer.
5. A single deterministic ideal initial state may have $\lambda_i(0)=1$ and $\dot\lambda_i(0)=0$. A broader full-state measure $\mu_0$ must be declared physically.
6. The microscopic one-point probability is mechanically generated as an empirical spatial counting measure or ensemble push-forward. No Gaussian, Weibull, Gibbs, Boltzmann, or other named family is imposed on $P(\lambda,\tau)$.
7. Neighbour independence is not assumed. The exact moment hierarchy retains neighbour joint information.
8. The exact general variance equation is

   $$
   D_\tau\Theta
   +2\Theta\partial_\lambda u
   +\frac1P\partial_\lambda(PC_3)
   =2\Psi,
   $$

   with

   $$
   \Psi(\lambda,\tau)=\operatorname{Cov}(c,\ddot\lambda\mid\lambda,\tau).
   $$

9. $\Theta$ is a conditional spacing-rate variance and is neither a fatigue-damage scalar nor the complete chain kinetic-energy density.
10. Same-load non-retracing of $(P,u,\Theta)$ establishes dynamic history dependence, not irreversible thermodynamic dissipation.
11. The conservative reference layer has $\dot D_{\mathrm{irr}}=0$ unless a separately derived irreversible microscopic force is introduced.
12. The operational local initiation threshold is the loss of positive tangent stiffness,

   $$
   \phi''(\lambda_c)=0,
   $$

   followed by first passage through $\lambda_c$.
13. Local first-passage fraction and specimen-to-specimen initiation probability remain distinct. Specimen probability requires a separate correlation/realization model.
14. The exact finite-chain projection does not close autonomously from $P_0$ alone because microscopic ordering/correlation information is lost. This is an established limitation of Layer A, not a numerical defect.

## 2. Layer B — reduced laboratory-time-scale assumptions / 축약 실험실 시간척도 가정

The final reduced closure is not claimed to be an exact projection of Layer A. It is an explicit asymptotic/coarse model with the following assumptions.

1. **Structural initial density.** $P_0(\lambda)$ is a structural/prestress spacing density defined at a declared reference load phase. It is not an instantaneous thermal-vibration PDF.
2. **Local prestress embedding.** Each reference label $\lambda_0$ is interpreted as a stable local equilibrium through

   $$
   q_r(\lambda_0)=\phi'(\lambda_0)-q_{\mathrm{ref}}.
   $$

   This is the explicit $P_0$-only closure assumption. It does not reconstruct full neighbour ordering.
3. **Laboratory-frequency quasistatics.** Normal elastic relaxation is much faster than the laboratory fatigue period, so the intact structural spacing follows

   $$
   \phi'[\Lambda(\lambda_0,t)]
   =\phi'(\lambda_0)+q(t)-q_{\mathrm{ref}},
   \qquad
   \phi''(\Lambda)>0.
   $$

4. **Finite-temperature eliminated fast modes.** Atomic thermal degrees of freedom are not represented inside $P_0$. They are treated as a fast local thermal reservoir for the rare crossing calculation.
5. **Fast intrawell re-equilibration.** Thermal equilibration within a stable intact well is assumed fast compared with both the fatigue period and the mean escape time.
6. **Rare-event regime.** The transition-state sink is used only when

   $$
   \Delta G_c\gg k_BT.
   $$

   Near deterministic loss of stability, direct stable-branch failure replaces this asymptotic rate.
7. **Mechanically derived transition-state rate.** The reduced sink is

   $$
   k_c(\lambda,T;A_c)
   =\frac{\sqrt{\phi''(\lambda)}}{2\pi t_0}
   \exp\left[-\frac{EA_ca_0\Delta\psi_c(\lambda)}{k_BT}\right].
   $$

   The barrier and attempt frequency are inherited from the retained generalized-LJ normal mechanics. No free diffusion coefficient, fitted damping constant, arbitrary Kramers prefactor, or empirical S-N kernel is inserted.
8. **Characteristic cohesive area remains symbolic.** $A_c$ is a physical characteristic area controlling the activation energy. It is not calibrated at the present theory stage, is not a FEM element area, and is not automatically an independent-cell area.
9. **Absorbing local event.** Once a local domain reaches the operational instability boundary it leaves the intact survivor population. Post-crack propagation is outside the present model.
10. **No required permanent PDF drift.** Reversible mechanical $P$ may be periodic. Irreversibility in the reduced initiation model is the loss of survivor mass through first passage.
11. **Selection is allowed.** Even if the mechanical map is cycle-periodic, the normalized survivor distribution may change because vulnerable parts of the initial structural population are preferentially removed.

Under these assumptions the active reduced equation is

$$
\partial_tP_b
+\partial_\lambda
\left[
\frac{\dot q(t)}{\phi''(\lambda)}P_b
\right]
=-k_c(\lambda,T;A_c)P_b,
$$

with

$$
P_b(\lambda,t_0)=P_0(\lambda),
\qquad
S(t)=\int P_b(\lambda,t)d\lambda,
\qquad
F_{\mathrm{ci}}(t)=1-S(t).
$$

Thus the laboratory reduced model is closed as

$$
P_0+\sigma(0:t)
\longrightarrow
P_b(\lambda,t),S(t),F_{\mathrm{ci}}(t).
$$

## 3. Assumptions that remain NOT active / 계속 비활성인 가정

The present 1D mainline does **not** assume:

- Gaussian/Weibull spacing or fatigue-life distributions;
- independent neighbouring spacings in the exact finite chain;
- independent statistical cells or independent FEM element failure;
- an arbitrary Smoluchowski mobility or arbitrary stochastic diffusion coefficient;
- fitted viscous damping introduced to create a fatigue clock;
- arbitrary Kramers/Arrhenius rates unrelated to the retained mechanics;
- empirical scalar fatigue damage;
- FCC reconstruction;
- registry/slip $s$ or well index $z$ as a required state of the normal-only closure;
- permanent hysteresis-energy storage in the same conservative pair potential.

## 4. Registry extension status / registry 확장 상태

The spatial registry/kink calculations are retained as a future plasticity diagnostic. The ideal pure-normal reduced surface showed a finite rare formation scale but only shallow post-formation trapping after correction of earlier bounded-minimization artifacts. Registry is therefore not the active long-lived memory variable.

For the ideal symmetric baseline, a uniform registry state remains invariant until an independently justified symmetry-breaking mechanism is introduced.

## 5. Remaining calibration and validation / 남은 보정 및 검증

The mathematical $P_0\to P_b,S,F_{\mathrm{ci}}$ closure is no longer the open problem. The remaining issues are:

- determine the physical characteristic cohesive area $A_c$ independently;
- measure or construct the structural/prestress $P_0$;
- validate the operational $\lambda_c$ initiation boundary;
- test transition-state recrossing/prefactor corrections where required;
- determine specimen correlation area/volume and local-to-specimen scaling later;
- compare predicted temperature, frequency, mean-stress, and initial-state trends with experiment;
- determine whether actual pure-Al single-crystal fatigue is dominated by the present normal-instability channel or by slip/dislocation structure.

The last point is a scope test. The reduced law is a closed and falsifiable **1D normal-instability fatigue-initiation hypothesis**, not a claim that all real aluminum fatigue physics has been solved.
