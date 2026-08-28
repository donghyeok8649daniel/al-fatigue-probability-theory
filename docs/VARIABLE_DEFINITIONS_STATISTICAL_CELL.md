# Variable Definitions — 1D Statistical Cell and Dependence

## Classification labels

- **EXACT / IDENTITY** — exact under the stated second-order stationary 1D assumptions.
- **DEFINITION** — a mathematical definition.
- **ESTIMATOR / DIAGNOSTIC** — exact formula evaluated using empirical finite-snapshot correlations rather than the true process correlation.
- **PHYSICAL INPUT / UNRESOLVED** — must be fixed by physics outside the current 1D axial correlation derivation.

## Variables

| Symbol | Definition | Meaning | Classification |
|---|---|---|---|
| $\rho_k$ | $\operatorname{Cov}(\lambda_i,\lambda_{i+k})/\sigma_\lambda^2$ | lag-$k$ spacing correlation | DEFINITION |
| $\tau_M$ | $1+2\sum_{k=1}^{M-1}(1-k/M)\rho_k$ | finite correlation factor for the variance of the sample mean | EXACT / IDENTITY |
| $M_{\rm eff}$ | $M/\tau_M$ | variance-equivalent number of independent spacings | DEFINITION |
| $a_0$ | equilibrium represented-layer spacing | axial microscopic length scale | PHYSICAL INPUT |
| $\ell_{\rm stat}^{(2)}$ | $a_0\tau_M=Ma_0/M_{\rm eff}$ | second-moment / variance-equivalent axial statistical length | DEFINITION |
| $P_1(\lambda)$ | one-point spacing density | marginal spacing state | DEFINITION |
| $P_2(\lambda,\lambda';k)$ | lag-$k$ joint spacing density | pair dependence state | DEFINITION |
| $A_0$ | area entering $E_0=EA_0a_0$ | mechanical energy-calibration patch area | PHYSICAL INPUT / UNRESOLVED |
| $A_{\rm stat}$ | not yet derived in 1D | transverse statistical correlation/independence area | PHYSICAL INPUT / UNRESOLVED |
| $V_{\rm stat}$ | $A_{\rm stat}\ell_{\rm stat}$ when a transverse theory exists | future statistical correlation volume | OUTSIDE CURRENT 1D SCOPE |

## Exact dependence criteria

Complete identical dependence:

$$
\boxed{
\mathbb E[(X-Y)^2]=0
\iff
X=Y\text{ almost surely}.
}
$$

Full independence:

$$
\boxed{
P_2(\lambda,\lambda';k)
=P_1(\lambda)P_1(\lambda').
}
$$

Zero covariance alone is not an independence criterion.

## Exact variance identity

For a second-order stationary sequence,

$$
\boxed{
\operatorname{Var}(\bar\lambda_M)
=
\frac{\sigma_\lambda^2}{M}\tau_M.
}
$$

When empirical $\rho_k$ are substituted, $\tau_M$, $M_{\rm eff}$, and $\ell_{\rm stat}^{(2)}$ are estimators/diagnostics rather than exact material constants.

---

# 한국어 번역 — 1D 통계셀 및 종속성 변수정의

## 분류 라벨

- **EXACT / IDENTITY** — 명시한 second-order stationary 1D 조건에서 정확.
- **DEFINITION** — 수학적 정의.
- **ESTIMATOR / DIAGNOSTIC** — true process correlation 대신 finite snapshot의 empirical correlation을 exact 식에 넣은 결과.
- **PHYSICAL INPUT / UNRESOLVED** — 현재 1D axial correlation 유도 바깥의 물리에서 정해야 하는 값.

## 변수

| 기호 | 정의 | 의미 | 분류 |
|---|---|---|---|
| $\rho_k$ | $\operatorname{Cov}(\lambda_i,\lambda_{i+k})/\sigma_\lambda^2$ | lag-$k$ spacing correlation | DEFINITION |
| $\tau_M$ | $1+2\sum_{k=1}^{M-1}(1-k/M)\rho_k$ | sample mean 분산에 들어가는 finite correlation factor | EXACT / IDENTITY |
| $M_{\rm eff}$ | $M/\tau_M$ | 분산 기준으로 동등한 독립 spacing 개수 | DEFINITION |
| $a_0$ | equilibrium represented-layer spacing | 축방향 microscopic length scale | PHYSICAL INPUT |
| $\ell_{\rm stat}^{(2)}$ | $a_0\tau_M=Ma_0/M_{\rm eff}$ | second-moment / variance-equivalent axial statistical length | DEFINITION |
| $P_1(\lambda)$ | one-point spacing density | marginal spacing state | DEFINITION |
| $P_2(\lambda,\lambda';k)$ | lag-$k$ joint spacing density | pair dependence state | DEFINITION |
| $A_0$ | $E_0=EA_0a_0$에 들어가는 면적 | mechanical energy-calibration patch area | PHYSICAL INPUT / UNRESOLVED |
| $A_{\rm stat}$ | 1D에서는 아직 유도하지 않음 | 횡방향 statistical correlation/independence area | PHYSICAL INPUT / UNRESOLVED |
| $V_{\rm stat}$ | 횡방향 이론이 생긴 뒤 $A_{\rm stat}\ell_{\rm stat}$ | 미래의 statistical correlation volume | OUTSIDE CURRENT 1D SCOPE |

## 정확한 종속성 기준

완전히 동일한 종속:

$$
\boxed{
\mathbb E[(X-Y)^2]=0
\iff
X=Y\text{ almost surely}
}
$$

완전 독립:

$$
\boxed{
P_2(\lambda,\lambda';k)
=P_1(\lambda)P_1(\lambda')
}
$$

covariance가 0인 것만으로 독립이라고 할 수 없다.

## 정확한 평균분산 식

second-order stationary sequence에 대해

$$
\boxed{
\operatorname{Var}(\bar\lambda_M)
=
\frac{\sigma_\lambda^2}{M}\tau_M
}
$$

이다.

empirical $\rho_k$를 넣으면 $\tau_M$, $M_{\rm eff}$, $\ell_{\rm stat}^{(2)}$는 exact material constant가 아니라 estimator/diagnostic이다.
