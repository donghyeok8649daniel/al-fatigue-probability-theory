# Master 1D $P$-$u$-$\Theta$ Formulation / 기준 1D 유도

이 문서는 활성 normal-only 이론의 기준 미분형 유도문서다. 수식 렌더링 안정성을 위해 단순한 LaTeX 문법만 사용한다.

## 1. Physical scaling / 물리 스케일

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

## 2. Finite generalized-LJ chain / 유한 generalized-LJ 사슬

$M+1$개의 node $x_0,\ldots,x_M$을 두고 $x_0=0$으로 둔다.

$$
\lambda_i=x_i-x_{i-1}
$$

$$
a_i=a_0\lambda_i
$$

활성 normalized interaction energy는

$$
\phi(\lambda)=\frac{\lambda^{-m}}{m(m-n)}-\frac{\lambda^{-n}}{n(m-n)}
$$

이다.

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

전체 configurational energy는

$$
V^*=\sum_{i=1}^{M}\phi(\lambda_i)
$$

이다.

내부 node equation은

$$
\ddot x_j=\phi'(\lambda_{j+1})-\phi'(\lambda_j)
$$

$$
j=1,\ldots,M-1
$$

이고 loaded end는

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

을 만족한다.

경계 spacing은

$$
\ddot\lambda_1=\phi'(\lambda_2)-\phi'(\lambda_1)
$$

$$
\ddot\lambda_M=q(\tau)+\phi'(\lambda_{M-1})-2\phi'(\lambda_M)
$$

이다.

## 3. Exact mechanical energy / 정확한 기계에너지

$$
T^*=\frac{1}{2}\sum_{j=1}^{M}\dot x_j^2
$$

$$
E_{\mathrm{mech}}^*=T^*+V^*
$$

운동방정식을 대입하면

$$
\frac{dE_{\mathrm{mech}}^*}{d\tau}=q(\tau)\dot x_M
$$

을 얻는다. 현재 baseline은 prescribed external work를 제외하면 conservative다.

spacing coordinate로는

$$
x_j=\sum_{k=1}^{j}\lambda_k
$$

이다. $L_{jk}=1$ for $k\le j$, otherwise $0$로 두면

$$
x=L\lambda
$$

$$
G_\lambda=L^TL
$$

$$
T^*=\frac{1}{2}c^TG_\lambda c
$$

$$
(G_\lambda)_{k\ell}=M-\max(k,\ell)+1
$$

이다. 따라서 one-point variance $\Theta$만으로 전체 kinetic energy를 표현할 수 없다.

## 4. Mechanics-generated empirical probability / 역학에서 생성되는 경험적 확률

$$
c_i=\dot\lambda_i
$$

finite deterministic chain의 empirical phase-space measure는

$$
F_M(\lambda,c,\tau)=\frac{1}{M}\sum_{i=1}^{M}\delta(\lambda-\lambda_i(\tau))\delta(c-c_i(\tau))
$$

이다.

spacing marginal은

$$
P_M(\lambda,\tau)=\int F_M(\lambda,c,\tau)\,dc
$$

$$
P_M(\lambda,\tau)=\frac{1}{M}\sum_i\delta(\lambda-\lambda_i(\tau))
$$

이다. named PDF를 가정하지 않는다.

## 5. Exact phase-space transport / 정확한 위상공간 수송

empirical acceleration flux를

$$
\mathcal G_M(\lambda,c,\tau)=\frac{1}{M}\sum_i\ddot\lambda_i\delta(\lambda-\lambda_i)\delta(c-c_i)
$$

로 정의한다.

분포론적으로

$$
\partial_\tau F_M+\partial_\lambda(cF_M)+\partial_c\mathcal G_M=0
$$

이다.

smooth representation에서는

$$
A(\lambda,c,\tau)=\mathrm{E}[\ddot\lambda_i\mid\lambda_i=\lambda,c_i=c,\tau]
$$

로 두고

$$
\mathcal G=AF
$$

$$
\partial_\tau F+\partial_\lambda(cF)+\partial_c(AF)=0
$$

을 얻는다. 이 식은 true conditional acceleration을 사용하면 exact projected identity다.

## 6. Raw moment hierarchy / 원시 모멘트 계층

$$
R_r(\lambda,\tau)=\int c^rF(\lambda,c,\tau)\,dc
$$

$$
B_r(\lambda,\tau)=\int c^{r-1}A(\lambda,c,\tau)F(\lambda,c,\tau)\,dc
$$

velocity-space boundary term이 사라진다고 하면

$$
\partial_\tau R_r+\partial_\lambda R_{r+1}=rB_r
$$

이다.

## 7. Reduced fields / 축약장

$$
P(\lambda,\tau)=R_0
$$

$$
u(\lambda,\tau)=\mathrm{E}[c\mid\lambda,\tau]=\frac{R_1}{P}
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

raw moments는

$$
R_0=P
$$

$$
R_1=Pu
$$

$$
R_2=P(u^2+\Theta)
$$

$$
R_3=P(u^3+3u\Theta+C_3)
$$

이다.

## 8. Exact continuity and mean-flow equations / 정확한 연속·평균류 식

$r=0$에서

$$
\partial_\tau P+\partial_\lambda(Pu)=0
$$

을 얻는다.

확률류는

$$
J=Pu
$$

이다.

$r=1$에서

$$
\partial_\tau(Pu)+\partial_\lambda[P(u^2+\Theta)]=P\mathcal A
$$

을 얻는다.

$$
D_\tau=\partial_\tau+u\partial_\lambda
$$

를 사용하면

$$
D_\tau u=\mathcal A-\frac{1}{P}\partial_\lambda(P\Theta)
$$

이다.

## 9. Correct exact variance equation / 정확한 분산식

$r=2$에서

$$
\partial_\tau[P(u^2+\Theta)]+\partial_\lambda[P(u^3+3u\Theta+C_3)]=2P\mathrm{E}[c\ddot\lambda\mid\lambda,\tau]
$$

을 얻는다.

$$
\mathrm{E}[c\ddot\lambda\mid\lambda,\tau]=u\mathcal A+\Psi
$$

이므로

$$
D_\tau\Theta+2\Theta\partial_\lambda u+\frac{1}{P}\partial_\lambda(PC_3)=2\Psi
$$

이다.

$\Psi=0$은 additional condition이며 일반적인 spatial LJ chain에서 자동으로 성립하지 않는다.

## 10. Neighbour-conditioned acceleration / 이웃 조건부 가속도

bulk spacing에서

$$
\ddot\lambda_i=\phi'(\lambda_{i+1})-2\phi'(\lambda_i)+\phi'(\lambda_{i-1})
$$

이다.

ordered joint densities $P_2^+$와 $P_2^-$를 사용해

$$
m_+=\frac{1}{P}\int\phi'(\lambda')P_2^+(\lambda,\lambda',\tau)\,d\lambda'
$$

$$
m_-=\frac{1}{P}\int\phi'(\lambda')P_2^-(\lambda,\lambda',\tau)\,d\lambda'
$$

로 두면

$$
\mathcal A_{\mathrm{bulk}}=m_++m_- -2\phi'(\lambda)
$$

이다.

central rate를 포함한 $F_2^+$와 $F_2^-$를 쓰면

$$
\Psi_{\mathrm{bulk}}=\frac{1}{P}\iint(c-u)\phi'(\lambda')[F_2^++F_2^-]\,dc\,d\lambda'
$$

이다. neighbour independence는 필요하지 않다.

## 11. Density-shape relation / 확률밀도 형상식

mean-flow equation을 전개하면

$$
\frac{1}{P}\partial_\lambda(P\Theta)=\partial_\lambda\Theta+\Theta\partial_\lambda\ln P
$$

이다.

따라서

$$
\Theta\partial_\lambda\ln P=\mathcal A-D_\tau u-\partial_\lambda\Theta
$$

이다.

$P>0$ 및 $\Theta>0$인 smooth interval에서는

$$
\partial_\lambda\ln P=\frac{\mathcal A-D_\tau u}{\Theta}-\partial_\lambda\ln\Theta
$$

이고

$$
P(\lambda,\tau)=\frac{\mathcal N_P(\tau)}{\Theta(\lambda,\tau)}\exp\left(\int_{\lambda_*}^{\lambda}\frac{\mathcal A(\eta,\tau)-D_\tau u(\eta,\tau)}{\Theta(\eta,\tau)}\,d\eta\right)
$$

이다.

$$
\int_0^\infty P(\lambda,\tau)\,d\lambda=1
$$

로 $\mathcal N_P$를 정한다. $\Theta=0$에서는 divided form을 사용하지 않는다.

## 12. Exact time-integral forms / 정확한 시간적분형

$$
P(\lambda,\tau)=P_0(\lambda)-\partial_\lambda\int_{\tau_0}^{\tau}P(\lambda,s)u(\lambda,s)\,ds
$$

$$
P(\lambda,\tau)u(\lambda,\tau)=P_0(\lambda)u_0(\lambda)-\partial_\lambda\int_{\tau_0}^{\tau}P(u^2+\Theta)(\lambda,s)\,ds+\int_{\tau_0}^{\tau}P(\lambda,s)\mathcal A(\lambda,s)\,ds
$$

$$
P(u^2+\Theta)(\lambda,\tau)=P_0(u_0^2+\Theta_0)(\lambda)-\partial_\lambda\int_{\tau_0}^{\tau}P(u^3+3u\Theta+C_3)(\lambda,s)\,ds+2\int_{\tau_0}^{\tau}P(u\mathcal A+\Psi)(\lambda,s)\,ds
$$

## 13. Characteristic forms / 특성곡선형

$$
\frac{dX}{ds}=u(X(s),s)
$$

$$
X(\tau_0)=\alpha
$$

$$
\mathcal I_u(s;\alpha)=\int_{\tau_0}^{s}\partial_\lambda u(X(r),r)\,dr
$$

continuity equation으로부터

$$
P(X(\tau),\tau)=P_0(\alpha)e^{-\mathcal I_u(\tau;\alpha)}
$$

이다.

$$
S_\Theta=2\Psi-\frac{1}{P}\partial_\lambda(PC_3)
$$

로 두면

$$
\Theta(X(\tau),\tau)=e^{-2\mathcal I_u(\tau;\alpha)}\left(\Theta_0(\alpha)+\int_{\tau_0}^{\tau}e^{2\mathcal I_u(s;\alpha)}S_\Theta(X(s),s)\,ds\right)
$$

이다.

## 14. Same-load history dependence / 동일하중 이력의존성

$$
q(\tau_L)=q(\tau_U)=q^*
$$

$$
\dot q(\tau_L)>0
$$

$$
\dot q(\tau_U)<0
$$

로 두고

$$
\mathcal R_2(\tau)=\{P(\lambda,\tau),u(\lambda,\tau),\Theta(\lambda,\tau)\}
$$

를 비교한다.

$$
\mathcal R_2(\tau_L)\ne\mathcal R_2(\tau_U)
$$

이면 instantaneous load $q$만의 memoryless map으로 reduced state를 나타낼 수 없다. 이 사실은 dynamic history dependence이며 irreversible dissipation 자체를 증명하지 않는다.

## 15. G1 mean spacing / 평균 간격

$$
\bar\lambda(\tau)=\int_0^\infty\lambda P(\lambda,\tau)\,d\lambda
$$

$$
\bar a(t)=a_0\bar\lambda(t/t_0)
$$

boundary flux가 사라지면

$$
\frac{d\bar\lambda}{d\tau}=\int_0^\infty P(\lambda,\tau)u(\lambda,\tau)\,d\lambda
$$

이다.

## 16. G2 mean intrinsic configurational energy / 평균 고유 배치에너지

$$
\Delta\phi(\lambda)=\phi(\lambda)-\phi(1)
$$

$$
\bar U(\tau)=U_{\mathrm{ref}}\int_0^\infty\Delta\phi(\lambda)P(\lambda,\tau)\,d\lambda
$$

boundary flux가 사라지면

$$
\frac{d\bar U}{d\tau}=U_{\mathrm{ref}}\int_0^\infty\phi'(\lambda)P(\lambda,\tau)u(\lambda,\tau)\,d\lambda
$$

이다.

## 17. G3 irreversible history / 비가역 이력

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

future irreversible node force가 유도된다면

$$
\dot D_{\mathrm{irr}}^*=-\sum_jr_j^{\mathrm{irr}}\dot x_j
$$

이고 dissipative interpretation을 위해

$$
\dot D_{\mathrm{irr}}^*\ge0
$$

이어야 한다.

## 18. G4 first passage / 최초통과 균열개시

local tangent-stiffness threshold를

$$
\phi''(\lambda_c)=0
$$

로 정의한다.

$$
\lambda_c=\left(\frac{m+1}{n+1}\right)^{1/(m-n)}
$$

각 spacing의 first-passage time은

$$
\tau_i^c=\inf\{\tau\ge\tau_0:\lambda_i(\tau)\ge\lambda_c\}
$$

이다.

instantaneous tail은

$$
Q_c(\tau)=\int_{\lambda_c}^{\infty}P(\lambda,\tau)\,d\lambda
$$

이며 cumulative first passage가 아니다.

survivor subdensity $F_b$의 right boundary에서 incoming $c<0$에 대해

$$
F_b(\lambda_c,c,\tau)=0
$$

을 적용한다.

$$
j_{\mathrm{esc}}=\int_0^\infty cF_b(\lambda_c^-,c,\tau)\,dc
$$

$$
S(\tau)=\int_0^{\lambda_c}\int_{-\infty}^{\infty}F_b(\lambda,c,\tau)\,dc\,d\lambda
$$

$$
\frac{dS}{d\tau}=-j_{\mathrm{esc}}
$$

$$
F_{\mathrm{ci}}^{\mathrm{local}}=1-S
$$

$$
h_\tau=\frac{j_{\mathrm{esc}}}{S}=-\frac{d}{d\tau}\ln S
$$

## 19. Local and specimen probability / 국소·시편 확률

one realization에서

$$
\tau_{\mathrm{spec}}^c=\min_i\tau_i^c
$$

이다.

ensemble measure $\mu_0$가 선언되면 specimen survival은

$$
S_{\mathrm{spec}}(\tau)=\int I\left[\max_i\sup_{s\in[\tau_0,\tau]}\Lambda_i(s;\Gamma_0)<\lambda_c\right]\,\mu_0(d\Gamma_0)
$$

이다.

수학적 survival formula는 존재하지만 실제 재료의 $\mu_0$와 spatial correlation scale은 OPEN이다.

## 20. Closure status / 폐쇄 상태

full microscopic state는

$$
\Gamma=(x_1,\ldots,x_M,\dot x_1,\ldots,\dot x_M)
$$

이며 finite ODE 아래 closed deterministic state다.

projected fields는

$$
P,u,\Theta\to C_3,\Psi,P_2^\pm,F_2^\pm,\ldots
$$

의 hierarchy를 가진다.

따라서 3-field PDE는 exact하지만 autonomous하지 않다. 그러나 full deterministic flow를 투영하면 exact integral representation이 존재한다.
