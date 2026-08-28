# Variable Definitions — 1D Layer-LJ Distribution Closure

## Classification

- **EXACT / IDENTITY** — exact within the explicitly stated mathematical model.
- **DEFINITION** — notation or mathematical definition.
- **ASSUMPTION** — additional physical/statistical assumption.
- **CONTROLLED APPROXIMATION** — approximation that must be validated or falsified.
- **EMPIRICAL INPUT / CALIBRATION** — externally supplied or previously calibrated material input.

## Variables

| Symbol | Definition | Meaning | Unit | Classification |
|---|---|---|---|---|
| $a_i(t)$ | neighboring represented-layer separation | local 1D normal layer spacing | m | DEFINITION |
| $a_0$ | reference layer spacing | equilibrium/reference spacing | m | CALIBRATION |
| $\lambda_i$ | $a_i/a_0$ | normalized local layer spacing | dimensionless | DEFINITION |
| $P_a(a,t)$ | physical spacing density | probability density in physical spacing | m$^{-1}$ | DEFINITION |
| $p_\lambda(\lambda,t)$ | $a_0P_a(a_0\lambda,t)$ | density in normalized spacing | dimensionless density | DEFINITION |
| $M$ | number of represented layer gaps in the finite ensemble | finite system size, not fatigue cycles | dimensionless | DEFINITION |
| $L(t)$ | $\sum_{i=1}^M\lambda_i=M\mu(t)$ | total normalized length constraint | dimensionless | DEFINITION |
| $E_c(t)$ | $\sum_{i=1}^M\psi(\lambda_i)=M\mathcal E(t)$ | total shifted configurational energy in normalized units | dimensionless | DEFINITION |
| $\phi(\lambda)$ | normalized generalized-LJ energy | calibrated layer-LJ shape used by active code | dimensionless | DEFINITION / CALIBRATED MODEL |
| $\psi(\lambda)$ | $\phi(\lambda)-\phi(1)$ | shifted layer-LJ energy | dimensionless | DEFINITION |
| $\mu(t)$ | $\int\lambda p_\lambda d\lambda$ | mean normalized layer spacing | dimensionless | DEFINITION |
| $\mathcal E(t)$ | $\int\psi p_\lambda d\lambda$ | mean shifted configurational energy | dimensionless | DEFINITION |
| $\Omega_M(L,E)$ | constrained configuration-space measure | density of reduced layer-spacing states at fixed $L,E$ | convention-dependent | DEFINITION |
| $\alpha(t)$ | length-conjugate saddle-point multiplier | enforces $\mu(t)$ | dimensionless in normalized form | CONTROLLED APPROXIMATION variable |
| $\beta(t)$ | energy-conjugate saddle-point multiplier | enforces $\mathcal E(t)$ | dimensionless in normalized form | CONTROLLED APPROXIMATION variable |
| $Z(\alpha,\beta)$ | $\int_0^\infty e^{-\alpha\lambda-\beta\psi(\lambda)}d\lambda$ | closure partition integral | dimensionless | DEFINITION |
| $\lambda_c$ | $\phi''(\lambda_c)=0$ | current idealized layer-opening tangent-instability stretch | dimensionless | EXACT under stated LJ shape |
| $Q_c(t)$ | $\int_{\lambda_c}^\infty p_\lambda d\lambda$ | instantaneous closure probability above $\lambda_c$ | dimensionless | DEFINITION |

## Distribution closure

The finite-$M$ constrained density of states is

$$
\Omega_M(L,E)
=
\int
\delta\left(\sum_i\lambda_i-L\right)
\delta\left(\sum_i\psi(\lambda_i)-E\right)
\,d^M\lambda.
$$

**ASSUMPTION:** at a stated physical time, reduced configurations compatible with the instantaneous $(L,E)$ constraints receive equal base measure.

Under that assumption the finite-$M$ one-spacing marginal is

$$
p_M(\lambda\mid L,E)
=
\frac{
\Omega_{M-1}(L-\lambda,E-\psi(\lambda))
}{
\Omega_M(L,E)
}.
$$

The large-$M$ saddle-point approximation gives

$$
\boxed{
p_\lambda(\lambda,t)
=
Z^{-1}
\exp[-\alpha(t)\lambda-\beta(t)\psi(\lambda)].
}
$$

The constraints are

$$
\mu(t)
=-\partial_\alpha\ln Z,
$$

$$
\mathcal E(t)
=-\partial_\beta\ln Z.
$$

The symbol $\beta$ is **not defined as thermodynamic inverse temperature**. Such an identification would require a separate derivation.

---

# 한국어 번역 — 1D Layer-LJ Distribution Closure 변수정의

## 분류

- **EXACT / IDENTITY** — 명시한 수학모델 내부에서 정확.
- **DEFINITION** — 기호 또는 수학적 정의.
- **ASSUMPTION** — 추가 물리/통계 가정.
- **CONTROLLED APPROXIMATION** — 검증 또는 반증해야 하는 근사.
- **EMPIRICAL INPUT / CALIBRATION** — 외부에서 공급되거나 기존에 calibration된 재료입력.

## 변수

| 기호 | 정의 | 의미 | 단위 | 분류 |
|---|---|---|---|---|
| $a_i(t)$ | 인접 represented layer 사이 간격 | 국부 1D 수직 layer spacing | m | DEFINITION |
| $a_0$ | reference layer spacing | 평형/기준 spacing | m | CALIBRATION |
| $\lambda_i$ | $a_i/a_0$ | normalized local layer spacing | 무차원 | DEFINITION |
| $P_a(a,t)$ | physical spacing density | physical spacing 좌표의 probability density | m$^{-1}$ | DEFINITION |
| $p_\lambda(\lambda,t)$ | $a_0P_a(a_0\lambda,t)$ | normalized spacing density | 무차원 density | DEFINITION |
| $M$ | finite ensemble의 represented layer gap 수 | finite system size이며 fatigue cycle이 아님 | 무차원 | DEFINITION |
| $L(t)$ | $\sum_{i=1}^M\lambda_i=M\mu(t)$ | total normalized length constraint | 무차원 | DEFINITION |
| $E_c(t)$ | $\sum_{i=1}^M\psi(\lambda_i)=M\mathcal E(t)$ | normalized unit의 total shifted configurational energy | 무차원 | DEFINITION |
| $\phi(\lambda)$ | normalized generalized-LJ energy | active code에서 사용하는 calibrated layer-LJ shape | 무차원 | DEFINITION / CALIBRATED MODEL |
| $\psi(\lambda)$ | $\phi(\lambda)-\phi(1)$ | shifted layer-LJ energy | 무차원 | DEFINITION |
| $\mu(t)$ | $\int\lambda p_\lambda d\lambda$ | 평균 normalized layer spacing | 무차원 | DEFINITION |
| $\mathcal E(t)$ | $\int\psi p_\lambda d\lambda$ | 평균 shifted configurational energy | 무차원 | DEFINITION |
| $\Omega_M(L,E)$ | constrained configuration-space measure | fixed $L,E$에서 reduced layer-spacing state density | convention-dependent | DEFINITION |
| $\alpha(t)$ | length-conjugate saddle-point multiplier | $\mu(t)$를 강제 | normalized form에서 무차원 | CONTROLLED APPROXIMATION variable |
| $\beta(t)$ | energy-conjugate saddle-point multiplier | $\mathcal E(t)$를 강제 | normalized form에서 무차원 | CONTROLLED APPROXIMATION variable |
| $Z(\alpha,\beta)$ | $\int_0^\infty e^{-\alpha\lambda-\beta\psi(\lambda)}d\lambda$ | closure partition integral | 무차원 | DEFINITION |
| $\lambda_c$ | $\phi''(\lambda_c)=0$ | 현재 이상화된 layer-opening tangent-instability stretch | 무차원 | stated LJ shape 아래 EXACT |
| $Q_c(t)$ | $\int_{\lambda_c}^\infty p_\lambda d\lambda$ | $\lambda_c$ 위 instantaneous closure probability | 무차원 | DEFINITION |

## Distribution closure

finite-$M$ constrained density of states는

$$
\Omega_M(L,E)
=
\int
\delta\left(\sum_i\lambda_i-L\right)
\delta\left(\sum_i\psi(\lambda_i)-E\right)
\,d^M\lambda
$$

이다.

**ASSUMPTION:** 주어진 물리시간에서 instantaneous $(L,E)$ constraint와 양립하는 reduced configuration에 동일한 base measure를 부여한다.

이 가정 아래 finite-$M$ one-spacing marginal은

$$
p_M(\lambda\mid L,E)
=
\frac{
\Omega_{M-1}(L-\lambda,E-\psi(\lambda))
}{
\Omega_M(L,E)
}
$$

이다.

large-$M$ saddle-point approximation은

$$
\boxed{
p_\lambda(\lambda,t)
=
Z^{-1}
\exp[-\alpha(t)\lambda-\beta(t)\psi(\lambda)]
}
$$

를 준다.

constraint는

$$
\mu(t)
=-\partial_\alpha\ln Z,
$$

$$
\mathcal E(t)
=-\partial_\beta\ln Z
$$

이다.

$\beta$는 **thermodynamic inverse temperature로 정의되지 않는다.** 그 해석을 하려면 별도의 유도가 필요하다.
