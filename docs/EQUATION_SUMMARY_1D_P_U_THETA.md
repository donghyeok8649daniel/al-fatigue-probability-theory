# Active 1D equation summary — $P$–$u$–$\Theta$

This file is the compact authoritative equation sheet for the active normal-only theory.  
이 파일은 현재 **1D normal-only 이론의 기준 수식 정리본**이다.

All symbols are defined in `VARIABLE_INDEX_1D_P_U_THETA.md` with bilingual mathematical and physical definitions.  
모든 기호는 `VARIABLE_INDEX_1D_P_U_THETA.md`에서 영·한 수학적/물리적 정의를 갖는다.

**Status labels:** MODEL, DEFINITION, EXACT, CONDITIONAL, OPEN.

## E01. Physical–nondimensional bridge / 물리–무차원 연결

$$
\boxed{
 t_0=\sqrt{\frac{m_a a_0}{EA_0}},
 \qquad
 \tau=\frac{t}{t_0}
}
$$

$$
\boxed{
U_{\rm ref}=EA_0a_0,
\qquad
F_{\rm ref}=EA_0,
\qquad
q(\tau)=\frac{F_{\rm ext}}{EA_0}=\frac{\sigma_n(t)}{E}
}
$$

For sinusoidal loading,

$$
\boxed{
\sigma_n(t)=\sigma_m+\sigma_a\sin(2\pi f t),
\qquad
\omega^*=2\pi f t_0
}
$$

**Status:** DEFINITION under the current calibration bridge.

## E02. Generalized-LJ normal energy / generalized-LJ 수직 에너지

$$
\boxed{
\phi(\lambda)
=
\frac{\lambda^{-m}}{m(m-n)}
-
\frac{\lambda^{-n}}{n(m-n)},
\qquad m>n>1
}
$$

$$
\boxed{
\phi'(\lambda)
=\frac{\lambda^{-n-1}-\lambda^{-m-1}}{m-n}
}
$$

$$
\boxed{
\phi''(\lambda)
=\frac{(m+1)\lambda^{-m-2}-(n+1)\lambda^{-n-2}}{m-n}
}
$$

$$
\boxed{
\phi'(1)=0,
\qquad
\phi''(1)=1
}
$$

$$
\boxed{
V^*(\boldsymbol\lambda)=\sum_{i=1}^{M}\phi(\lambda_i)
}
$$

**Status:** $\phi$ is MODEL; derivatives and consequences are EXACT under that model.

## E03. Microscopic chain equations / 미시 사슬 운동식

$$
\boxed{
\lambda_i=x_i-x_{i-1},
\qquad
a_i=a_0\lambda_i,
\qquad
c_i=\dot\lambda_i
}
$$

Interior nodes:

$$
\boxed{
\ddot x_j
=\phi'(\lambda_{j+1})-\phi'(\lambda_j),
\qquad j=1,\ldots,M-1
}
$$

Loaded end:

$$
\boxed{
\ddot x_M=-\phi'(\lambda_M)+q(\tau)
}
$$

Bulk spacings:

$$
\boxed{
\ddot\lambda_i
=\phi'(\lambda_{i+1})
-2\phi'(\lambda_i)
+\phi'(\lambda_{i-1}),
\qquad i=2,\ldots,M-1
}
$$

Boundary spacings:

$$
\boxed{
\ddot\lambda_1=\phi'(\lambda_2)-\phi'(\lambda_1)
}
$$

$$
\boxed{
\ddot\lambda_M
=q(\tau)+\phi'(\lambda_{M-1})-2\phi'(\lambda_M)
}
$$

**Status:** EXACT under the active finite-chain model and boundary law.

## E04. Mechanical energy and spacing-coordinate mass metric / 기계에너지와 간격좌표 질량메트릭

$$
\boxed{
E_{\rm mech}^*=T^*+V^*,
\qquad
T^*=\frac12\sum_{j=1}^{M}\dot x_j^2
}
$$

$$
\boxed{
\frac{dE_{\rm mech}^*}{d\tau}=q(\tau)\dot x_M
}
$$

With $\boldsymbol x=\mathbf L\boldsymbol\lambda$,

$$
\boxed{
T^*=\frac12\boldsymbol c^T\mathbf G_\lambda\boldsymbol c,
\qquad
\mathbf G_\lambda=\mathbf L^T\mathbf L
}
$$

$$
\boxed{
(G_\lambda)_{k\ell}=M-\max(k,\ell)+1
}
$$

**Status:** EXACT. Therefore one-point $\Theta$ alone is not the complete chain kinetic energy.

## E05. Empirical phase-space state / 경험적 위상공간 상태

$$
\boxed{
F_M(\lambda,c,\tau)
=\frac1M\sum_{i=1}^{M}
\delta[\lambda-\lambda_i(\tau)]
\delta[c-c_i(\tau)]
}
$$

$$
\boxed{
P_M(\lambda,\tau)
=\int F_M\,dc
=\frac1M\sum_i\delta[\lambda-\lambda_i(\tau)]
}
$$

No Gaussian, Weibull, Boltzmann, or other named PDF is assumed.

**Status:** DEFINITION.

## E06. Exact projected phase-space transport / 정확한 투영 위상공간 수송

$$
\boxed{
A(\lambda,c,\tau)
=\mathbb E[\ddot\lambda_i\mid\lambda_i=\lambda,c_i=c]
}
$$

$$
\boxed{
\partial_\tau F+\partial_\lambda(cF)+\partial_c(AF)=0
}
$$

**Status:** EXACT projected identity when $A$ is the true mechanics-generated conditional acceleration; not an autonomous one-point closure.

## E07. Exact raw moment hierarchy / 정확한 원시 모멘트 계층

$$
\boxed{
R_r(\lambda,\tau)=\int c^rF(\lambda,c,\tau)\,dc
}
$$

$$
\boxed{
B_r(\lambda,\tau)=\int c^{r-1}AF\,dc,
\qquad r\ge1
}
$$

$$
\boxed{
\partial_\tau R_r+\partial_\lambda R_{r+1}=rB_r,
\qquad r=0,1,2,\ldots
}
$$

assuming finite moments and vanishing velocity-boundary terms.

**Status:** EXACT.

## E08. Reduced fields / 축약장

$$
\boxed{
P(\lambda,\tau)=\int F(\lambda,c,\tau)\,dc
}
$$

The active mean-rate symbol is $u$; Greek $\nu$ is reserved and is not used for this state field.

$$
\boxed{
 u(\lambda,\tau)=\mathbb E[c\mid\lambda,\tau]
}
$$

$$
\boxed{
\Theta(\lambda,\tau)
=\operatorname{Var}(c\mid\lambda,\tau)
=\mathbb E[(c-u)^2\mid\lambda,\tau]
}
$$

$$
\boxed{
C_3(\lambda,\tau)=\mathbb E[(c-u)^3\mid\lambda,\tau]
}
$$

$$
\boxed{
\mathcal A(\lambda,\tau)=\mathbb E[\ddot\lambda_i\mid\lambda_i=\lambda,\tau]
}
$$

$$
\boxed{
\Psi(\lambda,\tau)
=\operatorname{Cov}(c,\ddot\lambda\mid\lambda,\tau)
=\mathbb E[(c-u)\ddot\lambda\mid\lambda,\tau]
}
$$

$$
\boxed{
D_\tau=\partial_\tau+u\partial_\lambda
}
$$

## E09. Exact $P$–$u$–$\Theta$ equations / 정확한 $P$–$u$–$\Theta$ 식

$$
\boxed{
\partial_\tau P+\partial_\lambda(Pu)=0
}
$$

$$
\boxed{
D_\tau u
=\mathcal A-\frac1P\partial_\lambda(P\Theta)
}
$$

$$
\boxed{
D_\tau\Theta
+2\Theta\partial_\lambda u
+\frac1P\partial_\lambda(PC_3)
=2\Psi
}
$$

The shorter zero-source $\Theta$ equation is valid only under the additional condition $\Psi=0$.

**Status:** EXACT general equations; the $\Psi=0$ version is CONDITIONAL.

## E10. Exact neighbour-statistics form / 정확한 이웃상관 형태

$$
\boxed{
m_+
=\frac1P\int\phi'(\lambda')P_2^+(\lambda,\lambda',\tau)\,d\lambda'
}
$$

$$
\boxed{
m_-
=\frac1P\int\phi'(\lambda')P_2^-(\lambda,\lambda',\tau)\,d\lambda'
}
$$

$$
\boxed{
\mathcal A_{\rm bulk}=m_++m_- -2\phi'(\lambda)
}
$$

$$
\boxed{
\Psi_{\rm bulk}
=\frac1P\iint
(c-u)\phi'(\lambda')[F_2^++F_2^-]\,dc\,d\lambda'
}
$$

**Status:** EXACT for bulk spacings; no neighbour-independence assumption.

## E11. Instantaneous density-shape relation / 순간 확률밀도 형상식

$$
\boxed{
\Theta\partial_\lambda\ln P
=\mathcal A-D_\tau u-\partial_\lambda\Theta
}
$$

For smooth $P>0$ and $\Theta>0$,

$$
\boxed{
\partial_\lambda\ln P
=\frac{\mathcal A-D_\tau u}{\Theta}
-\partial_\lambda\ln\Theta
}
$$

and

$$
\boxed{
P(\lambda,\tau)
=\frac{\mathcal N_P(\tau)}{\Theta(\lambda,\tau)}
\exp\left[
\int_{\lambda_*}^{\lambda}
\frac{\mathcal A(\eta,\tau)-D_\tau u(\eta,\tau)}{\Theta(\eta,\tau)}\,d\eta
\right]
}
$$

where $\mathcal N_P$ is fixed by $\int P\,d\lambda=1$.

At $\Theta=0$ the divided form is invalid; use the undivided moment and transport equations.

**Status:** EXACT instantaneous reconstruction under smoothness/positivity conditions.

## E12. Exact Volterra time-integral form / 정확한 볼테라 시간적분형

$$
\boxed{
P(\lambda,\tau)
=P_0(\lambda)
-\partial_\lambda\int_{\tau_0}^{\tau}P(\lambda,s)u(\lambda,s)\,ds
}
$$

$$
\boxed{
Pu
=P_0u_0
-\partial_\lambda\int_{\tau_0}^{\tau}P(u^2+\Theta)(\lambda,s)\,ds
+\int_{\tau_0}^{\tau}P\mathcal A(\lambda,s)\,ds
}
$$

$$
\boxed{
\begin{aligned}
P(u^2+\Theta)(\lambda,\tau)
={}&P_0(u_0^2+\Theta_0)(\lambda)\\
&-\partial_\lambda\int_{\tau_0}^{\tau}P(u^3+3u\Theta+C_3)(\lambda,s)\,ds\\
&+2\int_{\tau_0}^{\tau}P(u\mathcal A+\Psi)(\lambda,s)\,ds.
\end{aligned}
}
$$

**Status:** EXACT hierarchical integral equations.

## E13. Characteristic integral form / 특성곡선 적분형

Define

$$
\boxed{
\frac{dX}{ds}=u(X(s),s),
\qquad X(\tau_0)=\alpha
}
$$

and

$$
\boxed{
\mathcal I_u(s;\alpha)
=\int_{\tau_0}^{s}\partial_\lambda u(X(r),r)\,dr
}
$$

Then

$$
\boxed{
P(X(\tau),\tau)=P_0(\alpha)e^{-\mathcal I_u(\tau;\alpha)}
}
$$

$$
\boxed{
 u(X(\tau),\tau)
=u_0(\alpha)
+\int_{\tau_0}^{\tau}
\left[\mathcal A-\frac1P\partial_\lambda(P\Theta)\right]_{(X(s),s)}\,ds
}
$$

Define

$$
\boxed{
S_\Theta=2\Psi-\frac1P\partial_\lambda(PC_3)
}
$$

Then

$$
\boxed{
\Theta(X(\tau),\tau)
=e^{-2\mathcal I_u(\tau;\alpha)}
\left[
\Theta_0(\alpha)
+\int_{\tau_0}^{\tau}e^{2\mathcal I_u(s;\alpha)}S_\Theta(X(s),s)\,ds
\right]
}
$$

**Status:** EXACT characteristic representation when a smooth characteristic map exists.

## E14. Exact full-flow push-forward / 전체 미시흐름의 정확한 push-forward

$$
\boxed{
\Gamma=(x_1,\ldots,x_M,\dot x_1,\ldots,\dot x_M)
}
$$

$$
\boxed{
\Gamma(\tau)=\Phi_{\tau,\tau_0}^{q}(\Gamma_0)
}
$$

$$
\boxed{
\int\mu_0(d\Gamma_0)=1
}
$$

$$
\boxed{
\Lambda_i(\tau;\Gamma_0)=x_i(\tau;\Gamma_0)-x_{i-1}(\tau;\Gamma_0)
}
$$

$$
\boxed{
C_i=\frac{d\Lambda_i}{d\tau},
\qquad
A_i=\frac{d^2\Lambda_i}{d\tau^2}
}
$$

$$
\boxed{
F(\lambda,c,\tau)
=\frac1M\sum_i\int
\delta[\lambda-\Lambda_i]
\delta[c-C_i]\,\mu_0(d\Gamma_0)
}
$$

$$
\boxed{
P(\lambda,\tau)
=\frac1M\sum_i\int\delta[\lambda-\Lambda_i]\,\mu_0(d\Gamma_0)
}
$$

$$
\boxed{
Pu
=\frac1M\sum_i\int C_i\delta[\lambda-\Lambda_i]\,\mu_0(d\Gamma_0)
}
$$

$$
\boxed{
P(u^2+\Theta)
=\frac1M\sum_i\int C_i^2\delta[\lambda-\Lambda_i]\,\mu_0(d\Gamma_0)
}
$$

Also,

$$
\boxed{
P\mathcal A
=\frac1M\sum_i\int A_i\delta[\lambda-\Lambda_i]\,\mu_0(d\Gamma_0)
}
$$

$$
\boxed{
PC_3
=\frac1M\sum_i\int(C_i-u)^3\delta[\lambda-\Lambda_i]\,\mu_0(d\Gamma_0)
}
$$

$$
\boxed{
P\Psi
=\frac1M\sum_i\int(C_i-u)A_i\delta[\lambda-\Lambda_i]\,\mu_0(d\Gamma_0)
}
$$

**Status:** EXACT under the closed finite-chain model and declared $\mu_0$.

## E15. G1 — mean spacing / 평균 간격

$$
\boxed{
\bar\lambda(\tau)=\int_0^\infty\lambda P(\lambda,\tau)\,d\lambda
}
$$

$$
\boxed{
\bar a(t)=a_0\bar\lambda(t/t_0)
}
$$

$$
\boxed{
\bar a(\tau)=\frac{a_0}{M}\sum_i\int\Lambda_i(\tau;\Gamma_0)\,\mu_0(d\Gamma_0)
}
$$

**Status:** G1 DEFINITION plus exact equivalent projection form.

## E16. G2 — mean intrinsic configurational energy / 평균 고유 배치에너지

$$
\boxed{
\Delta\phi(\lambda)=\phi(\lambda)-\phi(1)
}
$$

$$
\boxed{
\bar U(\tau)
=U_{\rm ref}\int_0^\infty\Delta\phi(\lambda)P(\lambda,\tau)\,d\lambda
}
$$

$$
\boxed{
\bar U(\tau)
=\frac{U_{\rm ref}}{M}\sum_i\int[\phi(\Lambda_i)-\phi(1)]\,\mu_0(d\Gamma_0)
}
$$

**Status:** EXACT under the same nearest-neighbour energy used in the microscopic equations.

## E17. G3 — irreversible hysteresis energy / 비가역 히스테리시스 에너지

$$
\boxed{
E_{\rm hyst}(t)=\int_0^t\dot D_{\rm irr}(t')\,dt',
\qquad
\dot D_{\rm irr}\ge0
}
$$

Present conservative baseline:

$$
\boxed{
\dot D_{\rm irr}=0,
\qquad
E_{\rm hyst}=0
}
$$

If a future physical irreversible node force is derived,

$$
\boxed{
\dot D_{\rm irr}^*=-\sum_jr_j^{\rm irr}\dot x_j\ge0
}
$$

$$
\boxed{
\frac{dE_{\rm mech}^*}{d\tau}=q\dot x_M-\dot D_{\rm irr}^*
}
$$

**Status:** G3 observable DEFINED; physical irreversible law OPEN.

## E18. G4 — mechanical first passage / 기계적 최초통과 균열개시

$$
\boxed{
\phi''(\lambda_c)=0
}
$$

$$
\boxed{
\lambda_c=\left(\frac{m+1}{n+1}\right)^{1/(m-n)}
}
$$

$$
\boxed{
\tau_i^c=\inf\{\tau\ge\tau_0:\lambda_i(\tau)\ge\lambda_c\}
}
$$

For the survivor phase-space subdensity $F_b$,

$$
\boxed{
F_b(\lambda_c,c,\tau)=0\qquad(c<0)
}
$$

$$
\boxed{
j_{\rm esc}=\int_0^\infty cF_b(\lambda_c^-,c,\tau)\,dc
}
$$

$$
\boxed{
\dot S=-j_{\rm esc},
\qquad
F_{\rm ci}^{\rm local}=1-S,
\qquad
h=\frac{j_{\rm esc}}S=-\frac{d}{d\tau}\ln S
}
$$

Full-flow local survival:

$$
\boxed{
S_{\rm local}(\tau)
=\frac1M\sum_i\int
\mathbf1\left[
\sup_{s\in[\tau_0,\tau]}\Lambda_i(s;\Gamma_0)<\lambda_c
\right]\mu_0(d\Gamma_0)
}
$$

Specimen survival:

$$
\boxed{
S_{\rm spec}(\tau)
=\int
\mathbf1\left[
\max_i\sup_{s\in[\tau_0,\tau]}\Lambda_i(s;\Gamma_0)<\lambda_c
\right]\mu_0(d\Gamma_0)
}
$$

$$
\boxed{
F_{\rm ci}^{\rm spec}=1-S_{\rm spec}
}
$$

**Status:** exact first-passage mathematics once $\mu_0$ is declared; physical specimen-scale $\mu_0$ and correlation calibration remain OPEN.

## E19. Final mathematical structure / 최종 수학 구조

$$
\boxed{
\text{closed finite LJ ODE}
\Longrightarrow
\Phi^q
\Longrightarrow
F
\Longrightarrow
(P,u,\Theta)
}
$$

$$
\boxed{
(P,u,\Theta)\text{ satisfy exact hierarchical reduced equations}
}
$$

$$
\boxed{
\text{no autonomous three-field closure}
\neq
\text{no exact integral solution representation}
}
$$

Remaining open problems are physical: the irreversible G3 mechanism, physical specimen measure $\mu_0$/correlation scale, microscopic-to-laboratory fatigue time-scale bridge, and experimental validation.
