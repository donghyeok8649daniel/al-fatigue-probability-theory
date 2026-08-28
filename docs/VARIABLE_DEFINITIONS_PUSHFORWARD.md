# Variable Definitions — 1D Layer-LJ Push-Forward Distribution

## Classification labels

- **EXACT / IDENTITY** — exact under the stated finite/continuum definitions.
- **DEFINITION** — chosen mathematical definition.
- **LINEARIZED MECHANICS** — exact result of the linearized 1D layer-LJ equation.
- **NONLINEAR EXPANSION** — Taylor expansion of the calibrated LJ force.
- **NUMERICAL DIAGNOSTIC** — comparison with deterministic snapshots.

## Variables

| Symbol | Definition | Meaning | Unit | Classification |
|---|---|---|---|---|
| $\xi$ | continuum layer label in $[0,1]$ | normalized spatial label | dimensionless | DEFINITION |
| $\Lambda(\xi,t)$ | continuum normalized spacing field | deterministic layer-spacing waveform | dimensionless | DEFINITION |
| $p_\lambda(\lambda,t)$ | $\int_0^1\delta[\lambda-\Lambda(\xi,t)]d\xi$ | one-point spacing density as a spatial push-forward | inverse stretch | EXACT / IDENTITY |
| $u_i$ | $\lambda_i-1$ | small spacing perturbation about equilibrium | dimensionless | DEFINITION |
| $q$ | spatial wavenumber | discrete normal-mode wavenumber | radians per layer index | DEFINITION |
| $\omega_q$ | $2|\sin(q/2)|$ | linearized mode angular frequency | dimensionless inverse time | LINEARIZED MECHANICS |
| $A$ | first-harmonic amplitude | amplitude of $\cos\vartheta$ spacing component | dimensionless | DEFINITION |
| $B$ | second-harmonic amplitude | amplitude of $\cos2\vartheta$ spacing component | dimensionless | DEFINITION |
| $c_1$ | $\phi''(1)$ | linear LJ force coefficient | dimensionless | NONLINEAR EXPANSION |
| $c_2$ | $\phi'''(1)/2$ | quadratic LJ force coefficient | dimensionless | NONLINEAR EXPANSION |
| $c_3$ | $\phi''''(1)/6$ | cubic LJ force coefficient | dimensionless | NONLINEAR EXPANSION |

For the active calibration $m=12.19$, $n=6$,

$$
\boxed{c_1=1,\qquad c_2=-10.595,\qquad c_3=62.97935.}
$$

## Push-forward formula

For a piecewise monotone spacing field,

$$
\boxed{
p_\lambda(\lambda,t)
=
\sum_{\xi_j:\Lambda(\xi_j,t)=\lambda}
\frac{1}{|\partial_\xi\Lambda(\xi_j,t)|}.
}
$$

This formula is not a probability-family assumption. It follows from the delta-function representation of the spatial push-forward.

## Single-mode reference

For

$$
\Lambda=\mu+A\cos\vartheta,
$$

uniform spatial phase gives

$$
\boxed{
p(\lambda)
=
\frac{\mathbf 1_{|\lambda-\mu|<|A|}}
{\pi\sqrt{A^2-(\lambda-\mu)^2}}.
}
$$

The resulting variance is $A^2/2$ and the skewness is zero.

## Two-harmonic reference

For

$$
\Lambda=\mu+A\cos\vartheta+B\cos2\vartheta,
$$

$$
\operatorname{Var}(\lambda)=\frac{A^2+B^2}{2},
$$

$$
\mu_3=\frac34A^2B,
$$

and

$$
\boxed{|\gamma_1|\le\sqrt{2/3}.}
$$

---

# 한국어 번역 — 1D Layer-LJ Push-Forward Distribution 변수정의

## 분류 라벨

- **EXACT / IDENTITY** — 명시된 finite/continuum 정의 아래 정확.
- **DEFINITION** — 선택한 수학적 정의.
- **LINEARIZED MECHANICS** — 선형화된 1D layer-LJ 방정식의 정확한 결과.
- **NONLINEAR EXPANSION** — calibration된 LJ force의 Taylor expansion.
- **NUMERICAL DIAGNOSTIC** — deterministic snapshot과의 비교.

## 변수

| 기호 | 정의 | 의미 | 단위 | 분류 |
|---|---|---|---|---|
| $\xi$ | $[0,1]$의 continuum layer label | normalized spatial label | 무차원 | DEFINITION |
| $\Lambda(\xi,t)$ | continuum normalized spacing field | deterministic layer-spacing waveform | 무차원 | DEFINITION |
| $p_\lambda(\lambda,t)$ | $\int_0^1\delta[\lambda-\Lambda(\xi,t)]d\xi$ | spatial push-forward로 표현한 one-point spacing density | inverse stretch | EXACT / IDENTITY |
| $u_i$ | $\lambda_i-1$ | equilibrium 주변의 작은 spacing perturbation | 무차원 | DEFINITION |
| $q$ | spatial wavenumber | discrete normal-mode wavenumber | layer index당 radian | DEFINITION |
| $\omega_q$ | $2|\sin(q/2)|$ | 선형화된 mode angular frequency | 무차원 inverse time | LINEARIZED MECHANICS |
| $A$ | first-harmonic amplitude | $\cos\vartheta$ spacing component의 amplitude | 무차원 | DEFINITION |
| $B$ | second-harmonic amplitude | $\cos2\vartheta$ spacing component의 amplitude | 무차원 | DEFINITION |
| $c_1$ | $\phi''(1)$ | linear LJ force coefficient | 무차원 | NONLINEAR EXPANSION |
| $c_2$ | $\phi'''(1)/2$ | quadratic LJ force coefficient | 무차원 | NONLINEAR EXPANSION |
| $c_3$ | $\phi''''(1)/6$ | cubic LJ force coefficient | 무차원 | NONLINEAR EXPANSION |

현재 calibration $m=12.19$, $n=6$에서

$$
\boxed{c_1=1,\qquad c_2=-10.595,\qquad c_3=62.97935}
$$

이다.

## Push-forward 식

piecewise monotone spacing field에서는

$$
\boxed{
p_\lambda(\lambda,t)
=
\sum_{\xi_j:\Lambda(\xi_j,t)=\lambda}
\frac{1}{|\partial_\xi\Lambda(\xi_j,t)|}
}
$$

이다.

이 식은 probability-family assumption이 아니다. spatial push-forward의 delta-function 표현에서 직접 나온다.

## Single-mode 기준

$$
\Lambda=\mu+A\cos\vartheta
$$

이고 spatial phase가 uniform하면

$$
\boxed{
p(\lambda)
=
\frac{\mathbf 1_{|\lambda-\mu|<|A|}}
{\pi\sqrt{A^2-(\lambda-\mu)^2}}
}
$$

이다.

이때 variance는 $A^2/2$이고 skewness는 0이다.

## Two-harmonic 기준

$$
\Lambda=\mu+A\cos\vartheta+B\cos2\vartheta
$$

이면

$$
\operatorname{Var}(\lambda)=\frac{A^2+B^2}{2},
$$

$$
\mu_3=\frac34A^2B,
$$

그리고

$$
\boxed{|\gamma_1|\le\sqrt{2/3}}
$$

이다.
