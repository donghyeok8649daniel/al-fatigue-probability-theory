# Variable Definitions — 1D Layer-LJ Closure Falsification

## Classification labels

- **EXACT / IDENTITY** — exact under the explicitly stated equations and boundary conditions.
- **DEFINITION** — mathematical definition.
- **CONTROLLED APPROXIMATION** — reduction requiring numerical or theoretical validation.
- **NUMERICAL DIAGNOSTIC** — quantity used to compare a finite calculation with the closure.

## Variables

| Symbol | Definition | Meaning | Unit | Classification |
|---|---|---|---|---|
| $M$ | number of represented layer spacings in a deterministic snapshot | finite spatial sample size | dimensionless | DEFINITION |
| $\lambda_i(t)$ | $a_i(t)/a_0$ | local dimensionless layer separation | dimensionless | DEFINITION |
| $\mu_{\rm sim}(t)$ | $M^{-1}\sum_i\lambda_i(t)$ | empirical mean stretch | dimensionless | DEFINITION |
| $\mathcal E_{\rm sim}(t)$ | $M^{-1}\sum_i\psi(\lambda_i(t))$ | empirical mean shifted configurational energy | dimensionless | DEFINITION |
| $m_3(t)$ | $\int(\lambda-\mu)^3p_\lambda\,d\lambda$ | third central spacing moment | dimensionless | DEFINITION |
| $\gamma_1(t)$ | $m_3/\operatorname{Var}(\lambda)^{3/2}$ | spacing skewness | dimensionless | DEFINITION |
| $D_{\rm KS}$ | supremum distance between empirical and closure CDFs | full-distribution discrepancy diagnostic | dimensionless | NUMERICAL DIAGNOSTIC |
| $F_{\rm emp}(\lambda)$ | empirical cumulative spacing distribution | finite deterministic CDF | dimensionless | DEFINITION |
| $F_{\rm closure}(\lambda)$ | CDF of the solved two-moment closure | predicted CDF | dimensionless | CONTROLLED APPROXIMATION |
| $v_\lambda(\lambda,t)$ | conditional spacing-space velocity | transport velocity in normalized spacing space | s$^{-1}$ if $t$ is physical | DEFINITION |
| $T$ | loading period | period of imposed sinusoidal forcing | s or normalized time | DEFINITION |

## Closure used in the comparison

$$
p_\lambda(\lambda,t)=Z^{-1}\exp[-\alpha(t)\lambda-\beta(t)\psi(\lambda)].
$$

The multipliers are determined by

$$
\int\lambda p_\lambda\,d\lambda=\mu_{\rm sim}(t)
$$

and

$$
\int\psi(\lambda)p_\lambda\,d\lambda=\mathcal E_{\rm sim}(t).
$$

They are not histogram-fit parameters.

## Third central moment

$$
m_3(t)=\int[\lambda-\mu(t)]^3p_\lambda(\lambda,t)\,d\lambda.
$$

Under

$$
\partial_t p_\lambda+\partial_\lambda(p_\lambda v_\lambda)=0
$$

with vanishing boundary flux,

$$
\dot m_3=3\operatorname{Cov}\left((\lambda-\mu)^2,v_\lambda\right).
$$

This is an **EXACT / IDENTITY** under the stated kinematic conditions.

---

# 한국어 번역 — 1D Layer-LJ Closure 반증시험 변수정의

## 분류 라벨

- **EXACT / IDENTITY** — 명시된 equation과 boundary condition 아래 정확히 성립.
- **DEFINITION** — 수학적 정의.
- **CONTROLLED APPROXIMATION** — 이론적 또는 수치적 검증이 필요한 축약.
- **NUMERICAL DIAGNOSTIC** — finite calculation과 closure를 비교하기 위한 진단량.

## 변수

| 기호 | 정의 | 의미 | 단위 | 분류 |
|---|---|---|---|---|
| $M$ | deterministic snapshot에 포함된 represented layer spacing 수 | 유한 spatial sample size | 무차원 | DEFINITION |
| $\lambda_i(t)$ | $a_i(t)/a_0$ | 국부 무차원 layer separation | 무차원 | DEFINITION |
| $\mu_{\rm sim}(t)$ | $M^{-1}\sum_i\lambda_i(t)$ | empirical mean stretch | 무차원 | DEFINITION |
| $\mathcal E_{\rm sim}(t)$ | $M^{-1}\sum_i\psi(\lambda_i(t))$ | empirical mean shifted configurational energy | 무차원 | DEFINITION |
| $m_3(t)$ | $\int(\lambda-\mu)^3p_\lambda\,d\lambda$ | spacing의 third central moment | 무차원 | DEFINITION |
| $\gamma_1(t)$ | $m_3/\operatorname{Var}(\lambda)^{3/2}$ | spacing skewness | 무차원 | DEFINITION |
| $D_{\rm KS}$ | empirical CDF와 closure CDF의 supremum distance | full-distribution discrepancy diagnostic | 무차원 | NUMERICAL DIAGNOSTIC |
| $F_{\rm emp}(\lambda)$ | empirical cumulative spacing distribution | finite deterministic CDF | 무차원 | DEFINITION |
| $F_{\rm closure}(\lambda)$ | solved two-moment closure의 CDF | predicted CDF | 무차원 | CONTROLLED APPROXIMATION |
| $v_\lambda(\lambda,t)$ | conditional spacing-space velocity | normalized spacing space의 transport velocity | physical $t$이면 s$^{-1}$ | DEFINITION |
| $T$ | loading period | imposed sinusoidal forcing의 주기 | s 또는 normalized time | DEFINITION |

## 비교에 사용하는 closure

$$
p_\lambda(\lambda,t)=Z^{-1}\exp[-\alpha(t)\lambda-\beta(t)\psi(\lambda)]
$$

이다.

multiplier는

$$
\int\lambda p_\lambda\,d\lambda=\mu_{\rm sim}(t)
$$

및

$$
\int\psi(\lambda)p_\lambda\,d\lambda=\mathcal E_{\rm sim}(t)
$$

로 결정한다.

histogram fitting parameter가 아니다.

## Third central moment

$$
m_3(t)=\int[\lambda-\mu(t)]^3p_\lambda(\lambda,t)\,d\lambda
$$

이다.

$$
\partial_t p_\lambda+\partial_\lambda(p_\lambda v_\lambda)=0
$$

및 vanishing boundary flux 아래에서

$$
\dot m_3=3\operatorname{Cov}\left((\lambda-\mu)^2,v_\lambda\right)
$$

가 성립한다.

이 식은 stated kinematic condition 아래 **EXACT / IDENTITY**다.
