# Variable Definitions — Continuous-Time 1D LJ Energy Feasibility

This file defines the symbols introduced by the continuous-time energy-feasibility formulation.

| Symbol | Definition | Physical meaning | Unit / type | Classification |
|---|---|---|---|---|
| $t$ | Physical time | Fundamental evolution coordinate | s | DEFINITION |
| $a$ | Local normal spacing | Normal interatomic separation coordinate | m | DEFINITION |
| $a_0$ | Reference equilibrium spacing | Scale used to nondimensionalize spacing | m | EMPIRICAL INPUT or calibration output |
| $\lambda$ | $a/a_0$ | Dimensionless local normal stretch | dimensionless | DEFINITION |
| $P(\lambda,t)$ | Normalized spacing density | Continuous-time probability density in normal-stretch space | dimensionless density in $\lambda$ | DEFINITION |
| $\mu(t)$ | $\int\lambda P(\lambda,t)d\lambda$ | Mean normal stretch | dimensionless | DEFINITION |
| $\phi(\lambda)$ | normalized generalized-LJ pair energy | Fixed 1D microscopic normal energy | dimensionless | DEFINITION within model |
| $\psi(\lambda)$ | $\phi(\lambda)-\phi(1)$ | LJ energy shifted so equilibrium has zero energy | dimensionless | DEFINITION |
| $\mathcal E(t)$ | $\int\psi(\lambda)P(\lambda,t)d\lambda$ | Mean configurational LJ energy per represented spacing | dimensionless | DEFINITION |
| $D_\psi(\lambda\mid\mu)$ | $\psi(\lambda)-\psi(\mu)-\psi'(\mu)(\lambda-\mu)$ | Convexity/Bregman remainder measuring energy associated with spread away from the mean | dimensionless | DEFINITION |
| $\lambda_c$ | first tensile solution of $\phi''(\lambda_c)=0$ | Idealized 1D LJ tangent-stability limit | dimensionless | mechanically derived model output |
| $\lambda_L(t)$ | hard lower support bound | Maximum admissible reverse-compression extent | dimensionless | PHYSICAL CONSTRAINT; must be derived or measured, not fitted |
| $\mathcal A_{\rm safe}(t)$ | admissible probability measures supported on $[\lambda_L(t),\lambda_c]$ with stated normalization and mean | Set of crack-free distributions allowed by the current constraints | set of measures | DEFINITION |
| $\mathcal E_{\rm safe}^{\min}(t)$ | $\psi(\mu(t))$ | Minimum energy of a safe distribution at the stated mean | dimensionless | EXACT under convexity/support assumptions |
| $\mathcal E_{\rm safe}^{\max}(t)$ | endpoint-chord expression | Maximum energy of a safe distribution at the stated mean | dimensionless | EXACT under convexity/support assumptions |
| $w_L(t)$ | $(\lambda_c-\mu)/(\lambda_c-\lambda_L)$ | Endpoint probability weight at $\lambda_L$ for the maximum-energy safe measure | dimensionless | EXACT derived quantity |
| $w_c(t)$ | $(\mu-\lambda_L)/(\lambda_c-\lambda_L)$ | Endpoint probability weight at $\lambda_c$ for the maximum-energy safe measure | dimensionless | EXACT derived quantity |
| $M_E(t)$ | $\mathcal E_{\rm safe}^{\max}(t)-\mathcal E(t)$ | Remaining safe energy margin | dimensionless | DEFINITION |
| $\tau_E$ | $\inf\{t:\mathcal E(t)>\mathcal E_{\rm safe}^{\max}(t)\}$ | First time at which no crack-free distribution satisfies the stated constraints | s | DEFINITION |
| $Q_c(t)$ | $\int_{\lambda_c}^{\infty}P(\lambda,t)d\lambda$ | Instantaneous unstable tensile-tail mass | dimensionless | DEFINITION |
| $N_E$ | $f\tau_E$ for constant $f$ only | Optional experimental conversion of a time to an equivalent cycle count | dimensionless | DERIVED LABEL, not a state variable |

## Important distinction

The symbol $N$ may still appear in a finite empirical density such as

$$
P_N(a,t)=\frac1N\sum_i\delta(a-a_i(t)),
$$

where it means **finite sample/system size**. It must not be confused with fatigue cycle count. The active evolution variable is $t$.

---

# 한국어 번역 — 연속시간 1D LJ 에너지 실현가능성 변수정의

이 문서는 연속시간 energy-feasibility formulation에서 새로 도입된 기호를 정의한다.

| 기호 | 정의 | 물리적 의미 | 단위 / 종류 | 분류 |
|---|---|---|---|---|
| $t$ | 물리적 시간 | 근본 evolution coordinate | s | DEFINITION |
| $a$ | 국부 수직 spacing | 수직 원자간 거리 좌표 | m | DEFINITION |
| $a_0$ | 기준 평형 spacing | spacing을 무차원화하는 기준길이 | m | EMPIRICAL INPUT 또는 calibration output |
| $\lambda$ | $a/a_0$ | 무차원 국부 수직 stretch | dimensionless | DEFINITION |
| $P(\lambda,t)$ | 정규화된 spacing density | normal-stretch space의 연속시간 확률밀도 | $\lambda$에 대한 dimensionless density | DEFINITION |
| $\mu(t)$ | $\int\lambda P(\lambda,t)d\lambda$ | 평균 수직 stretch | dimensionless | DEFINITION |
| $\phi(\lambda)$ | normalized generalized-LJ pair energy | 고정된 1D microscopic normal energy | dimensionless | model 내부 DEFINITION |
| $\psi(\lambda)$ | $\phi(\lambda)-\phi(1)$ | equilibrium energy를 0으로 이동한 LJ energy | dimensionless | DEFINITION |
| $\mathcal E(t)$ | $\int\psi(\lambda)P(\lambda,t)d\lambda$ | 대표 spacing 하나당 평균 configurational LJ energy | dimensionless | DEFINITION |
| $D_\psi(\lambda\mid\mu)$ | $\psi(\lambda)-\psi(\mu)-\psi'(\mu)(\lambda-\mu)$ | 평균에서 분포가 퍼지며 생기는 에너지를 나타내는 convexity/Bregman remainder | dimensionless | DEFINITION |
| $\lambda_c$ | $\phi''(\lambda_c)=0$의 첫 tensile solution | 이상화된 1D LJ tangent-stability limit | dimensionless | mechanics-derived model output |
| $\lambda_L(t)$ | hard lower support bound | 허용 가능한 최대 reverse-compression 범위 | dimensionless | PHYSICAL CONSTRAINT; fitting이 아니라 유도 또는 측정 필요 |
| $\mathcal A_{\rm safe}(t)$ | $[\lambda_L(t),\lambda_c]$에 support를 가지며 정규화와 평균조건을 만족하는 probability measure 집합 | 현재 제약 아래 허용되는 crack-free distribution 집합 | measure set | DEFINITION |
| $\mathcal E_{\rm safe}^{\min}(t)$ | $\psi(\mu(t))$ | 주어진 평균에서 safe distribution이 가질 수 있는 최소에너지 | dimensionless | convexity/support 조건 아래 EXACT |
| $\mathcal E_{\rm safe}^{\max}(t)$ | endpoint chord 식 | 주어진 평균에서 safe distribution이 가질 수 있는 최대에너지 | dimensionless | convexity/support 조건 아래 EXACT |
| $w_L(t)$ | $(\lambda_c-\mu)/(\lambda_c-\lambda_L)$ | maximum-energy safe measure의 $\lambda_L$ endpoint weight | dimensionless | EXACT derived quantity |
| $w_c(t)$ | $(\mu-\lambda_L)/(\lambda_c-\lambda_L)$ | maximum-energy safe measure의 $\lambda_c$ endpoint weight | dimensionless | EXACT derived quantity |
| $M_E(t)$ | $\mathcal E_{\rm safe}^{\max}(t)-\mathcal E(t)$ | 남아 있는 safe energy margin | dimensionless | DEFINITION |
| $\tau_E$ | $\inf\{t:\mathcal E(t)>\mathcal E_{\rm safe}^{\max}(t)\}$ | 주어진 제약을 만족하는 crack-free distribution이 더 이상 존재할 수 없는 최초시간 | s | DEFINITION |
| $Q_c(t)$ | $\int_{\lambda_c}^{\infty}P(\lambda,t)d\lambda$ | 순간 unstable tensile-tail mass | dimensionless | DEFINITION |
| $N_E$ | 일정한 $f$에서만 $f\tau_E$ | 시간을 equivalent cycle count로 환산한 선택적 실험표시 | dimensionless | DERIVED LABEL, state variable 아님 |

## 중요한 구분

다음과 같은 finite empirical density에서는

$$
P_N(a,t)=\frac1N\sum_i\delta(a-a_i(t))
$$

$N$이 계속 등장할 수 있다. 여기서 $N$은 **finite sample/system size**를 뜻하며 fatigue cycle count와 혼동하면 안 된다. 활성 evolution variable은 $t$이다.
