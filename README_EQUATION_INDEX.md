# Equation and Symbol Index / 수식·기호 Index

이 파일은 `theory-core`의 활성 수학 문서 진입점이다.

## 1. Authoritative documents / 기준 문서

1. [`docs/EQUATION_SUMMARY_1D_P_U_THETA.md`](docs/EQUATION_SUMMARY_1D_P_U_THETA.md) — 핵심 수식 정리본
2. [`docs/VARIABLE_INDEX_1D_P_U_THETA.md`](docs/VARIABLE_INDEX_1D_P_U_THETA.md) — 영·한 기호 사전
3. [`docs/AUXILIARY_SYMBOL_INDEX_1D.md`](docs/AUXILIARY_SYMBOL_INDEX_1D.md) — 보조기호 사전
4. [`docs/MASTER_1D_P_U_THETA_FORMULATION.md`](docs/MASTER_1D_P_U_THETA_FORMULATION.md) — 전체 미분형 유도
5. [`docs/MILESTONE25_EXACT_INTEGRAL_REPRESENTATION.md`](docs/MILESTONE25_EXACT_INTEGRAL_REPRESENTATION.md) — 정확한 적분 표현
6. [`docs/CRACK_INITIATION_DEFINITION.md`](docs/CRACK_INITIATION_DEFINITION.md) — 최초통과 균열개시 정의

## 2. Mandatory symbol-definition contract / 기호 정의 규칙

새 기호를 도입할 때는 해당 기호의 다음 항목을 Index에 함께 추가한다.

- Symbol / 기호
- Equation definition / 수식적 정의
- English term / 영문 명칭
- Korean term / 한국어 명칭
- Mathematical definition / 수학적 정의
- Physical definition / 물리적 정의
- Unit or scaling / 단위 또는 스케일
- Status / 상태
- Dependencies / 선행 정의

수식으로 정의할 수 있는 기호는 문장 설명만으로 끝내지 않는다. 같은 기호를 서로 다른 수학적 객체에 중복 사용하지 않는다.

## 3. Markdown math safety rule / Markdown 수식 안전 규칙

이 문서군은 렌더러 호환성을 위해 보수적인 LaTeX 부분집합만 사용한다.

- Inline math: `$ ... $`
- Display math: `$$ ... $$`
- 사용 금지: `\operatorname`, `\boxed`, `\text`, `\begin`, `\end`, `\mathbb`, `\mathbf`, `\boldsymbol`, `\mathsf`
- 기대값: `\mathrm{E}`
- 분산: `\mathrm{Var}`
- 공분산: `\mathrm{Cov}`
- 지시함수: `I[condition]`
- display 수식은 Markdown table 안에 넣지 않는다.
- 여러 줄 수식은 정렬 environment 대신 여러 개의 독립 `$$ ... $$` 블록으로 나눈다.
- 수식 안에 긴 영어/한국어 문장을 넣지 않는다.

## 4. Active mathematical chain / 활성 수학 체계

$$
\mathrm{1D\ LJ\ chain}
\to
\Phi^q_{\tau,\tau_0}
\to
F(\lambda,c,\tau)
\to
P(\lambda,\tau),u(\lambda,\tau),\Theta(\lambda,\tau)
$$

$$
P,u,\Theta
\to
\bar a,\bar U,S,F_{ci}
$$

물리 좌표와 무차원 좌표의 관계는

$$
\lambda=\frac{a}{a_0}
$$

$$
\tau=\frac{t}{t_0}
$$

이다. 따라서 물리적 분포와 무차원 분포는 서로 다른 밀도이며 Jacobian을 포함해 연결한다.

$$
P_a(a,t)\,da=P_\lambda(\lambda,t)\,d\lambda
$$

$$
P_\lambda(\lambda,t)=a_0 P_a(a_0\lambda,t)
$$

## 5. Governing-equation index / 지배방정식 Index

### E01. Generalized-LJ energy

$$
\phi(\lambda)=\frac{\lambda^{-m}}{m(m-n)}-\frac{\lambda^{-n}}{n(m-n)}
$$

### E02. Bulk spacing dynamics

$$
\ddot\lambda_i=\phi'(\lambda_{i+1})-2\phi'(\lambda_i)+\phi'(\lambda_{i-1})
$$

### E03. Empirical phase-space measure

$$
F_M(\lambda,c,\tau)=\frac{1}{M}\sum_i\delta(\lambda-\lambda_i)\delta(c-c_i)
$$

### E04. Projected phase-space transport

$$
\partial_\tau F+\partial_\lambda(cF)+\partial_c(AF)=0
$$

### E05. Continuity

$$
\partial_\tau P+\partial_\lambda(Pu)=0
$$

### E06. Conditional mean spacing-rate balance

$$
D_\tau u=\mathcal A-\frac{1}{P}\partial_\lambda(P\Theta)
$$

### E07. Exact variance balance

$$
D_\tau\Theta+2\Theta\partial_\lambda u+\frac{1}{P}\partial_\lambda(PC_3)=2\Psi
$$

### E08. Density-shape identity

$$
\Theta\partial_\lambda\ln P=\mathcal A-D_\tau u-\partial_\lambda\Theta
$$

### E09. Instantaneous integral form of the density

$$
P(\lambda,\tau)=\frac{\mathcal N_P(\tau)}{\Theta(\lambda,\tau)}\exp\left(\int_{\lambda_*}^{\lambda}\frac{\mathcal A(\eta,\tau)-D_\tau u(\eta,\tau)}{\Theta(\eta,\tau)}\,d\eta\right)
$$

이 식은 smooth region에서 $P>0$ 및 $\Theta>0$일 때만 나눗셈 형태로 사용한다.

### E10. Exact full-flow push-forward

$$
F(\lambda,c,\tau)=\frac{1}{M}\sum_i\int\delta(\lambda-\Lambda_i(\tau;\Gamma_0))\delta(c-C_i(\tau;\Gamma_0))\,\mu_0(d\Gamma_0)
$$

### E11. Characteristic integral for the variance

$$
\mathcal I_u(s;\alpha)=\int_{\tau_0}^{s}\partial_\lambda u(X(r),r)\,dr
$$

$$
S_\Theta=2\Psi-\frac{1}{P}\partial_\lambda(PC_3)
$$

$$
\Theta(X(\tau),\tau)=e^{-2\mathcal I_u(\tau;\alpha)}\left(\Theta_0(\alpha)+\int_{\tau_0}^{\tau}e^{2\mathcal I_u(s;\alpha)}S_\Theta(X(s),s)\,ds\right)
$$

### E12. G1 mean spacing

$$
\bar a=a_0\int\lambda P(\lambda,\tau)\,d\lambda
$$

### E13. G2 mean intrinsic configurational energy

$$
\bar U=U_{\mathrm{ref}}\int[\phi(\lambda)-\phi(1)]P(\lambda,\tau)\,d\lambda
$$

### E14. G3 irreversible history

$$
E_{\mathrm{hyst}}(t)=\int_0^t\dot D_{\mathrm{irr}}(t')\,dt'
$$

$$
\dot D_{\mathrm{irr}}\ge0
$$

현재 보존계 baseline에서는

$$
\dot D_{\mathrm{irr}}=0
$$

이다.

### E15. Local mechanical threshold

$$
\phi''(\lambda_c)=0
$$

$$
\lambda_c=\left(\frac{m+1}{n+1}\right)^{1/(m-n)}
$$

### E16. Specimen survival

$$
S_{\mathrm{spec}}(\tau)=\int I\left[\max_i\sup_{s\in[\tau_0,\tau]}\Lambda_i(s;\Gamma_0)<\lambda_c\right]\,\mu_0(d\Gamma_0)
$$

## 6. Core definitions / 핵심 정의

### Normalized spacing / 무차원 간격

$$
\lambda=\frac{a}{a_0}
$$

### Spacing rate / 간격 변화율

$$
c=\frac{d\lambda}{d\tau}
$$

### Marginal density / 주변밀도

$$
P(\lambda,\tau)=\int F(\lambda,c,\tau)\,dc
$$

### Conditional mean / 조건부 평균

$$
u(\lambda,\tau)=\mathrm{E}[c\mid\lambda,\tau]
$$

### Conditional variance / 조건부 분산

$$
\Theta(\lambda,\tau)=\mathrm{Var}(c\mid\lambda,\tau)=\mathrm{E}[(c-u)^2\mid\lambda,\tau]
$$

### Conditional acceleration / 조건부 가속도

$$
\mathcal A(\lambda,\tau)=\mathrm{E}[\ddot\lambda\mid\lambda,\tau]
$$

### Rate-acceleration covariance / 속도-가속도 공분산

$$
\Psi(\lambda,\tau)=\mathrm{Cov}(c,\ddot\lambda\mid\lambda,\tau)
$$

### Probability current / 확률류

$$
J(\lambda,\tau)=P(\lambda,\tau)u(\lambda,\tau)
$$

## 7. Scope / 범위

현재 수학적으로는 finite deterministic chain, exact projected transport, moment hierarchy, instantaneous density-shape relation, exact push-forward integral, characteristic integral, 그리고 first-passage survival 표현까지 정의돼 있다.

아직 물리적으로 열려 있는 항목은 다음과 같다.

- 비가역 G3 microscopic mechanism
- 실제 시편의 초기 측도 $\mu_0$와 spatial correlation scale
- 원자 시간척도와 실험실 피로 시간척도의 연결
- 실험 검증 및 calibration
