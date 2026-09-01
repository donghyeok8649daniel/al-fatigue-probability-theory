# Variable Definitions — 1D Layer-LJ Spatial Correlation

## Classification labels

- **EXACT / IDENTITY** — exact under the stated finite 1D definition.
- **DEFINITION** — chosen mathematical definition.
- **CONTROLLED NUMERICAL PROTOCOL** — numerical comparison protocol.
- **NUMERICAL RESULT / DIAGNOSTIC** — computed result, not a material theorem.

## Variables

| Symbol | Definition | Meaning | Unit | Classification |
|---|---|---|---|---|
| $M$ | number of represented spacings | finite 1D system size | count | DEFINITION |
| $k$ | integer lag | separation in layer-spacing index | count | DEFINITION |
| $\eta$ | $k/M$ | scaled spatial lag | dimensionless | DEFINITION |
| $\mu$ | $M^{-1}\sum_i\lambda_i$ | mean normalized layer spacing | dimensionless | DEFINITION |
| $C_k$ | lag-$k$ covariance | spatial ordering measure | dimensionless$^2$ | DEFINITION |
| $C_0$ | variance of $\lambda_i$ | zero-lag covariance | dimensionless$^2$ | EXACT / IDENTITY |
| $\rho_k$ | $C_k/C_0$ | normalized spatial correlation | dimensionless | DEFINITION |
| $k_0$ | first interpolated zero crossing of $\rho_k$ | correlation sign-change scale | spacing count | NUMERICAL DIAGNOSTIC |
| $k_0/M$ | scaled first zero crossing | system-size-normalized correlation scale | dimensionless | NUMERICAL DIAGNOSTIC |
| $P_2(\lambda,\lambda',t)$ | adjacent-spacing joint density | candidate pair state carrying ordering | normalized joint-density units | DEFINITION |

The dynamically matched sweep uses

$$
\boxed{\omega M=0.62}
$$

and samples at $t_s=2T$.

## Random-permutation null

For any centered finite spacing sample and any nonzero lag under a uniform random permutation,

$$
\boxed{\mathbb E_{\rm perm}[\rho_k]=-\frac{1}{M-1}}.
$$

This is **EXACT / IDENTITY** for the stated permutation ensemble.

---

# 한국어 번역 — 1D Layer-LJ 공간상관 변수정의

## 분류 라벨

- **EXACT / IDENTITY** — 명시된 finite 1D 정의 아래 정확.
- **DEFINITION** — 선택한 수학적 정의.
- **CONTROLLED NUMERICAL PROTOCOL** — 수치 비교 protocol.
- **NUMERICAL RESULT / DIAGNOSTIC** — 계산결과이며 material theorem이 아님.

## 변수

| 기호 | 정의 | 의미 | 단위 | 분류 |
|---|---|---|---|---|
| $M$ | represented spacing 수 | finite 1D system size | count | DEFINITION |
| $k$ | 정수 lag | layer-spacing index separation | count | DEFINITION |
| $\eta$ | $k/M$ | scaled spatial lag | 무차원 | DEFINITION |
| $\mu$ | $M^{-1}\sum_i\lambda_i$ | 평균 normalized layer spacing | 무차원 | DEFINITION |
| $C_k$ | lag-$k$ covariance | spatial ordering measure | 무차원$^2$ | DEFINITION |
| $C_0$ | $\lambda_i$ variance | zero-lag covariance | 무차원$^2$ | EXACT / IDENTITY |
| $\rho_k$ | $C_k/C_0$ | normalized spatial correlation | 무차원 | DEFINITION |
| $k_0$ | $\rho_k$의 첫 interpolated zero crossing | correlation sign-change scale | spacing count | NUMERICAL DIAGNOSTIC |
| $k_0/M$ | scaled first zero crossing | system-size-normalized correlation scale | 무차원 | NUMERICAL DIAGNOSTIC |
| $P_2(\lambda,\lambda',t)$ | adjacent-spacing joint density | ordering을 담는 candidate pair state | normalized joint-density unit | DEFINITION |

Dynamically matched sweep는

$$
\boxed{\omega M=0.62}
$$

를 사용하고 $t_s=2T$에서 sample한다.

## Random-permutation null

centered finite spacing sample을 uniform random permutation했을 때 모든 nonzero lag에서

$$
\boxed{\mathbb E_{\rm perm}[\rho_k]=-\frac{1}{M-1}}
$$

가 성립한다.

이는 stated permutation ensemble에 대해 **EXACT / IDENTITY**다.
