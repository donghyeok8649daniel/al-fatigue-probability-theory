# Active 1D Equation Summary / 활성 1D 수식 정리

이 문서는 현재 normal-only $P$-$u$-$\Theta$ 이론의 기준 수식 정리본이다.

모든 display 수식은 단순 `$$ ... $$` 형식으로 작성하며, 렌더러에서 자주 깨지는 고급 매크로와 정렬 environment는 사용하지 않는다.

## E01. Physical and nondimensional variables / 물리·무차원 변수

$$
t_0=\sqrt{\frac{m_a a_0}{EA_0}}
$$

$$
\tau=\frac{t}{t_0}
$$

$$
\lambda=\frac{a}{a_0}
$$

$$
F_{\mathrm{ref}}=EA_0
$$

$$
U_{\mathrm{ref}}=EA_0a_0
$$

$$
q(\tau)=\frac{F_{\mathrm{ext}}(t)}{EA_0}=\frac{\sigma_n(t)}{E}
$$

사인 하중은

$$
\sigma_n(t)=\sigma_m+\sigma_a\sin(2\pi f t)
$$

$$
\omega^*=2\pi f t_0
$$

으로 쓴다.

## E02. Generalized-LJ normal energy / generalized-LJ 수직 에너지

$$
\phi(\lambda)=\frac{\lambda^{-m}}{m(m-n)}-\frac{\lambda^{-n}}{n(m-n)}
$$

$$
\phi'(\lambda)=\frac{\lambda^{-n-1}-\lambda^{-m-1}}{m-n}
$$

$$
\phi''(\lambda)=\frac{(m+1)\lambda^{-m-2}-(n+1)\lambda^{-n-2}}{m-n}
$$

$$
\phi'(1)=0
$$

$$
\phi''(1)=1
$$

$$
V^*(\lambda_1,\ldots,\lambda_M)=\sum_{i=1}^{M}\phi(\lambda_i)
$$

Status: $\phi$는 MODEL이고, 그 미분식은 해당 모델 아래 EXACT이다.

## E03. Microscopic chain equations / 미시 사슬 운동식

$$
\lambda_i=x_i-x_{i-1}
$$

$$
a_i=a_0\lambda_i
$$

$$
c_i=\dot\lambda_i
$$

내부 node는

$$
\ddot x_j=\phi'(\lambda_{j+1})-\phi'(\lambda_j)
$$

$$
j=1,\ldots,M-1
$$

을 만족한다.

loaded end는

$$
\ddot x_M=-\phi'(\lambda_M)+q(\tau)
$$

이다.

bulk spacing은

$$
\ddot\lambda_i=\phi'(\lambda_{i+1})-2\phi'(\lambda_i)+\phi'(\lambda_{i-1})
$$

$$
i=2,\ldots,M-1
$$

이다.

경계 spacing은

$$
\ddot\lambda_1=\phi'(\lambda_2)-\phi'(\lambda_1)
$$

$$
\ddot\lambda_M=q(\tau)+\phi'(\lambda_{M-1})-2\phi'(\lambda_M)
$$

이다.

## E04. Mechanical energy / 기계에너지

$$
T^*=\frac{1}{2}\sum_{j=1}^{M}\dot x_j^2
$$

$$
E_{\mathrm{mech}}^*=T^*+V^*
$$

$$
\frac{dE_{\mathrm{mech}}^*}{d\tau}=q(\tau)\dot x_M
$$

spacing coordinate를 쓰면

$$
x_j=\sum_{k=1}^{j}\lambda_k
$$

이고, $L_{jk}=1$ for $k\le j$, otherwise $0$로 두면

$$
x=L\lambda
$$

$$
G_\lambda=L^TL
$$

$$
T^*=\frac{1}{2}c^T G_\lambda c
$$

$$
(G_\lambda)_{k\ell}=M-\max(k,\ell)+1
$$

이다. 따라서 $\frac12(u^2+\Theta)$는 전체 chain kinetic energy가 아니다.

## E05. Empirical phase-space state / 경험적 위상공간 상태

$$
F_M(\lambda,c,\tau)=\frac{1}{M}\sum_{i=1}^{M}\delta(\lambda-\lambda_i(\tau))\delta(c-c_i(\tau))
$$

$$
P_M(\lambda,\tau)=\int F_M(\lambda,c,\tau)\,dc
$$

$$
P_M(\lambda,\tau)=\frac{1}{M}\sum_i\delta(\lambda-\lambda_i(\tau))
$$

Gaussian, Weibull, Boltzmann 등 named PDF를 가정하지 않는다.

## E06. Exact projected phase-space transport / 정확한 투영 수송

$$
A(\lambda,c,\tau)=\mathrm{E}[\ddot\lambda_i\mid\lambda_i=\lambda,c_i=c,\tau]
$$

$$
\partial_\tau F+\partial_\lambda(cF)+\partial_c(AF)=0
$$

이 식은 true conditional acceleration $A$를 사용할 때 EXACT projected identity이며 autonomous closure는 아니다.

## E07. Raw moment hierarchy / 원시 모멘트 계층

$$
R_r(\lambda,\tau)=\int c^rF(\lambda,c,\tau)\,dc
$$

$$
B_r(\lambda,\tau)=\int c^{r-1}A(\lambda,c,\tau)F(\lambda,c,\tau)\,dc
$$

$$
\partial_\tau R_r+\partial_\lambda R_{r+1}=rB_r
$$

velocity-boundary term이 사라지고 필요한 모멘트가 존재할 때 EXACT이다.

## E08. Reduced fields / 축약장

$$
P(\lambda,\tau)=\int F(\lambda,c,\tau)\,dc
$$

$$
u(\lambda,\tau)=\mathrm{E}[c\mid\lambda,\tau]
$$

$$
\Theta(\lambda,\tau)=\mathrm{Var}(c\mid\lambda,\tau)
$$

$$
\Theta(\lambda,\tau)=\mathrm{E}[(c-u)^2\mid\lambda,\tau]
$$

$$
C_3(\lambda,\tau)=\mathrm{E}[(c-u)^3\mid\lambda,\tau]
$$

$$
\mathcal A(\lambda,\tau)=\mathrm{E}[\ddot\lambda_i\mid\lambda_i=\lambda,\tau]
$$

$$
\Psi(\lambda,\tau)=\mathrm{Cov}(c,\ddot\lambda\mid\lambda,\tau)
$$

$$
\Psi(\lambda,\tau)=\mathrm{E}[(c-u)\ddot\lambda\mid\lambda,\tau]
$$

$$
D_\tau=\partial_\tau+u\partial_\lambda
$$

## E09. Exact reduced equations / 정확한 축약식

Continuity:

$$
\partial_\tau P+\partial_\lambda(Pu)=0
$$

Mean spacing-rate balance:

$$
D_\tau u=\mathcal A-\frac{1}{P}\partial_\lambda(P\Theta)
$$

Variance balance:

$$
D_\tau\Theta+2\Theta\partial_\lambda u+\frac{1}{P}\partial_\lambda(PC_3)=2\Psi
$$

$\Psi=0$은 자동 조건이 아니다.

## E10. Neighbour-conditioned mechanics / 이웃 조건부 역학

$$
m_+(\lambda,\tau)=\frac{1}{P}\int\phi'(\lambda')P_2^+(\lambda,\lambda',\tau)\,d\lambda'
$$

$$
m_-(\lambda,\tau)=\frac{1}{P}\int\phi'(\lambda')P_2^-(\lambda,\lambda',\tau)\,d\lambda'
$$

$$
\mathcal A_{\mathrm{bulk}}=m_++m_- -2\phi'(\lambda)
$$

$$
\Psi_{\mathrm{bulk}}=\frac{1}{P}\iint(c-u)\phi'(\lambda')[F_2^++F_2^-]\,dc\,d\lambda'
$$

neighbour independence는 사용하지 않는다.

## E11. Density-shape identity / 확률밀도 형상식

$$
\Theta\partial_\lambda\ln P=\mathcal A-D_\tau u-\partial_\lambda\Theta
$$

$P>0$ 및 $\Theta>0$인 smooth region에서는

$$
\partial_\lambda\ln P=\frac{\mathcal A-D_\tau u}{\Theta}-\partial_\lambda\ln\Theta
$$

이고

$$
P(\lambda,\tau)=\frac{\mathcal N_P(\tau)}{\Theta(\lambda,\tau)}\exp\left(\int_{\lambda_*}^{\lambda}\frac{\mathcal A(\eta,\tau)-D_\tau u(\eta,\tau)}{\Theta(\eta,\tau)}\,d\eta\right)
$$

이다.

정규화는

$$
\int_0^\infty P(\lambda,\tau)\,d\lambda=1
$$

로 정한다.

## E12. Exact Volterra forms / 정확한 볼테라 적분형

$$
P(\lambda,\tau)=P_0(\lambda)-\partial_\lambda\int_{\tau_0}^{\tau}P(\lambda,s)u(\lambda,s)\,ds
$$

$$
P(\lambda,\tau)u(\lambda,\tau)=P_0(\lambda)u_0(\lambda)-\partial_\lambda\int_{\tau_0}^{\tau}P(u^2+\Theta)(\lambda,s)\,ds+\int_{\tau_0}^{\tau}P(\lambda,s)\mathcal A(\lambda,s)\,ds
$$

$$
P(u^2+\Theta)(\lambda,\tau)=P_0(u_0^2+\Theta_0)(\lambda)-\partial_\lambda\int_{\tau_0}^{\tau}P(u^3+3u\Theta+C_3)(\lambda,s)\,ds+2\int_{\tau_0}^{\tau}P(u\mathcal A+\Psi)(\lambda,s)\,ds
$$

## E13. Characteristic forms / 특성곡선 적분형

$$
\frac{dX}{ds}=u(X(s),s)
$$

$$
X(\tau_0)=\alpha
$$

$$
\mathcal I_u(s;\alpha)=\int_{\tau_0}^{s}\partial_\lambda u(X(r),r)\,dr
$$

$$
P(X(\tau),\tau)=P_0(\alpha)e^{-\mathcal I_u(\tau;\alpha)}
$$

$$
S_\Theta=2\Psi-\frac{1}{P}\partial_\lambda(PC_3)
$$

$$
\Theta(X(\tau),\tau)=e^{-2\mathcal I_u(\tau;\alpha)}\left(\Theta_0(\alpha)+\int_{\tau_0}^{\tau}e^{2\mathcal I_u(s;\alpha)}S_\Theta(X(s),s)\,ds\right)
$$

## E14. Exact full-flow push-forward / 전체 미시흐름 적분표현

$$
\Gamma=(x_1,\ldots,x_M,\dot x_1,\ldots,\dot x_M)
$$

$$
\Gamma(\tau)=\Phi^q_{\tau,\tau_0}(\Gamma_0)
$$

$$
\int\mu_0(d\Gamma_0)=1
$$

$$
\Lambda_i(\tau;\Gamma_0)=x_i(\tau;\Gamma_0)-x_{i-1}(\tau;\Gamma_0)
$$

$$
C_i=\frac{d\Lambda_i}{d\tau}
$$

$$
A_i=\frac{d^2\Lambda_i}{d\tau^2}
$$

$$
F(\lambda,c,\tau)=\frac{1}{M}\sum_i\int\delta(\lambda-\Lambda_i)\delta(c-C_i)\,\mu_0(d\Gamma_0)
$$

$$
P(\lambda,\tau)=\frac{1}{M}\sum_i\int\delta(\lambda-\Lambda_i)\,\mu_0(d\Gamma_0)
$$

$$
P u=\frac{1}{M}\sum_i\int C_i\delta(\lambda-\Lambda_i)\,\mu_0(d\Gamma_0)
$$

$$
P(u^2+\Theta)=\frac{1}{M}\sum_i\int C_i^2\delta(\lambda-\Lambda_i)\,\mu_0(d\Gamma_0)
$$

## E15. G1 / 평균 간격

$$
\bar\lambda(\tau)=\int_0^\infty\lambda P(\lambda,\tau)\,d\lambda
$$

$$
\bar a(t)=a_0\bar\lambda(t/t_0)
$$

## E16. G2 / 평균 고유 배치에너지

$$
\Delta\phi(\lambda)=\phi(\lambda)-\phi(1)
$$

$$
\bar U(\tau)=U_{\mathrm{ref}}\int_0^\infty\Delta\phi(\lambda)P(\lambda,\tau)\,d\lambda
$$

## E17. G3 / 비가역 히스테리시스 에너지

$$
E_{\mathrm{hyst}}(t)=\int_0^t\dot D_{\mathrm{irr}}(t')\,dt'
$$

$$
\dot D_{\mathrm{irr}}\ge0
$$

현재 conservative baseline에서는

$$
\dot D_{\mathrm{irr}}=0
$$

$$
E_{\mathrm{hyst}}=0
$$

이다.

## E18. G4 / 기계적 최초통과

$$
\phi''(\lambda_c)=0
$$

$$
\lambda_c=\left(\frac{m+1}{n+1}\right)^{1/(m-n)}
$$

$$
\tau_i^c=\inf\{\tau\ge\tau_0:\lambda_i(\tau)\ge\lambda_c\}
$$

생존 phase-space subdensity $F_b$에 대해 incoming $c<0$에서는

$$
F_b(\lambda_c,c,\tau)=0
$$

을 적용한다.

$$
j_{\mathrm{esc}}=\int_0^\infty cF_b(\lambda_c^-,c,\tau)\,dc
$$

$$
\frac{dS}{d\tau}=-j_{\mathrm{esc}}
$$

$$
F_{\mathrm{ci}}^{\mathrm{local}}=1-S
$$

$$
h=\frac{j_{\mathrm{esc}}}{S}=-\frac{d}{d\tau}\ln S
$$

Specimen survival은

$$
S_{\mathrm{spec}}(\tau)=\int I\left[\max_i\sup_{s\in[\tau_0,\tau]}\Lambda_i(s;\Gamma_0)<\lambda_c\right]\,\mu_0(d\Gamma_0)
$$

이다.

## E19. Closure status / 폐쇄 상태

finite microscopic ODE는 closed deterministic system이다.

$$
P,u,\Theta\to C_3,\Psi,P_2^\pm,F_2^\pm,\ldots
$$

따라서 3-field PDE는 exact하지만 autonomous하지 않다. 그러나 closed microscopic flow를 사용하면 exact push-forward integral representation은 존재한다.
