# Variable Definitions — Quasistatic Protocol and Ensemble Distinction

## Classification labels

- **EXACT STATIC RESULT**: follows directly from the stated homogeneous 1D force-controlled potential and stable-branch conditions.
- **DEFINITION**: introduced notation or diagnostic quantity.
- **CONTROLLED NUMERICAL DIAGNOSTIC**: protocol-dependent numerical quantity used to test a hypothesis; not a material constant.
- **FUTURE PHYSICAL ENSEMBLE OBJECT**: formally defined target whose physical ensemble still requires justification.

## Static force-control variables

| Symbol | Meaning | Classification |
|---|---|---|
| $M$ | number of represented 1D layer spacings | DEFINITION |
| $\lambda_i$ | normalized spacing $a_i/a_0$ | DEFINITION |
| $f$ | dimensionless tensile end force; under the current stress mapping $f=\sigma/E$ | DEFINITION / calibrated mapping |
| $f_c$ | maximum stable dimensionless tensile force where $\phi''(\lambda_c)=0$ | EXACT model quantity |
| $\lambda_c$ | tangent-instability stretch | EXACT model quantity |
| $\lambda_s(f)$ | unique stable root of $\phi'(\lambda)=f$ for $0\le f\le f_c$ | EXACT STATIC RESULT |
| $\Pi$ | force-controlled total potential $\sum_i[\phi(\lambda_i)-f\lambda_i]$ | DEFINITION |

For $0\le f<f_c$, stability implies

$$
\phi''(\lambda_s)>0,
$$

and the exact homogeneous static state is

$$
\lambda_i=\lambda_s(f)\quad\forall i.
$$

## Protocol variables

| Symbol | Meaning | Classification |
|---|---|---|
| $\omega$ | reduced angular drive frequency used by the deterministic chain integrator | DEFINITION |
| $N$ | integer cycle index used only as bookkeeping for the imposed periodic drive | DEFINITION |
| $\alpha$ | $\omega M$, a drive-rate/system-transit protocol diagnostic | CONTROLLED NUMERICAL DIAGNOSTIC |
| $C_0$ | empirical spatial spacing variance in one deterministic snapshot | DEFINITION |
| $\sqrt{C_0}$ | RMS nonuniformity amplitude of that snapshot | DEFINITION |
| $\rho_k$ | normalized open-chain lag-$k$ spacing correlation $C_k/C_0$ | DEFINITION |
| $\widehat\tau_M^{(+)}$ | first-positive-lobe finite-snapshot correlation-factor estimator | ESTIMATOR / DIAGNOSTIC |
| $\widehat M_{\rm eff}^{(+)}$ | $M/\widehat\tau_M^{(+)}$ | ESTIMATOR / DIAGNOSTIC |

A central Milestone 15 warning is

$$
C_0\to0
$$

does not force the normalized shape $\rho_k$ to vanish. Therefore a finite or system-scale $\widehat M_{\rm eff}^{(+)}$ is not physically meaningful as a material statistical-cell count when the underlying fluctuation amplitude is collapsing.

## Two different meanings of P

### Single-trajectory spatial empirical measure

$$
P_{M,\mathrm{spatial}}^{\mathrm{traj}}(\lambda,t)
=
\frac1M\sum_{i=1}^{M}\delta[\lambda-\lambda_i(t)].
$$

This is an **EXACT DEFINITION** for one deterministic trajectory. It contains spatial heterogeneity present in that trajectory, but it does not itself establish physical randomness.

### Ensemble-averaged probability target

$$
P_{\mathrm{ens}}(\lambda,t)
=
\mathbb E_{\Gamma_0}
\left[
\frac1M\sum_i\delta(\lambda-\lambda_i(t;\Gamma_0))
\right].
$$

This is a **FUTURE PHYSICAL ENSEMBLE OBJECT**. The initial phase-space ensemble $\Gamma_0$ must be physically specified before numerical values are interpreted as aluminum fatigue probabilities.

Candidate ensemble sources are not interchangeable and must be separately justified. Examples include a fixed-length canonical ensemble, a controlled metastable intact-basin ensemble, or independently supported material heterogeneity.

## Current exact cycle-boundary reference

For the Milestone 13/15 protocol with zero mean force and sinusoidal loading,

$$
f(t)=A(t)f_a\sin(\omega t).
$$

At an exact integer cycle,

$$
\omega t=2\pi N,
$$

so

$$
f=0.
$$

Therefore

$$
\lambda_s=1,
\qquad
C_0^{\rm qs}=0.
$$

Any nonzero $C_0$ at that snapshot is a dynamical residual relative to the exact static reference.

---

# 변수 정의 — 준정적 프로토콜과 ensemble 구분

## 분류 라벨

- **정확한 정적 결과(EXACT STATIC RESULT)**: 명시한 균질 1D force-controlled potential과 안정 branch 조건에서 직접 따라오는 결과.
- **정의(DEFINITION)**: 새로 도입한 기호 또는 진단량.
- **통제된 수치 진단(CONTROLLED NUMERICAL DIAGNOSTIC)**: 특정 protocol 가설을 검사하기 위한 수치량이며 물질상수가 아님.
- **미래 물리 ensemble 객체(FUTURE PHYSICAL ENSEMBLE OBJECT)**: 수학적으로는 정의하지만 어떤 물리 ensemble을 쓸지 아직 검증해야 하는 대상.

## 정적 force-control 변수

| 기호 | 의미 | 분류 |
|---|---|---|
| $M$ | represented 1D layer spacing 수 | 정의 |
| $\lambda_i$ | normalized spacing $a_i/a_0$ | 정의 |
| $f$ | dimensionless tensile end force; 현재 stress mapping에서 $f=\sigma/E$ | 정의 / calibration mapping |
| $f_c$ | $\phi''(\lambda_c)=0$에서의 최대 안정 tensile force | 정확한 모델량 |
| $\lambda_c$ | tangent-instability stretch | 정확한 모델량 |
| $\lambda_s(f)$ | $0\le f\le f_c$에서 $\phi'(\lambda)=f$의 유일한 안정 root | 정확한 정적 결과 |
| $\Pi$ | force-controlled potential $\sum_i[\phi(\lambda_i)-f\lambda_i]$ | 정의 |

$0\le f<f_c$에서

$$
\phi''(\lambda_s)>0
$$

이고 정확한 균질 정적 상태는

$$
\lambda_i=\lambda_s(f)\quad\forall i
$$

이다.

## 프로토콜 변수

| 기호 | 의미 | 분류 |
|---|---|---|
| $\omega$ | deterministic chain integrator의 reduced angular drive frequency | 정의 |
| $N$ | imposed periodic drive의 정수 cycle bookkeeping index | 정의 |
| $\alpha$ | $\omega M$, drive-rate/system-transit protocol 진단량 | 통제된 수치 진단 |
| $C_0$ | deterministic snapshot 하나의 empirical spatial spacing variance | 정의 |
| $\sqrt{C_0}$ | 그 snapshot의 RMS nonuniformity amplitude | 정의 |
| $\rho_k$ | normalized open-chain lag-$k$ correlation $C_k/C_0$ | 정의 |
| $\widehat\tau_M^{(+)}$ | first-positive-lobe finite-snapshot correlation-factor estimator | estimator / 진단 |
| $\widehat M_{\rm eff}^{(+)}$ | $M/\widehat\tau_M^{(+)}$ | estimator / 진단 |

Milestone 15의 핵심 경고는

$$
C_0\to0
$$

이어도 normalized shape $\rho_k$가 사라질 필요는 없다는 것이다. 따라서 fluctuation amplitude가 무너지는 상황에서 finite/system-scale $\widehat M_{\rm eff}^{(+)}$를 물질 고유 statistical-cell count로 해석하면 안 된다.

## 서로 다른 두 P

### single-trajectory spatial empirical measure

$$
P_{M,\mathrm{spatial}}^{\mathrm{traj}}(\lambda,t)
=
\frac1M\sum_{i=1}^{M}\delta[\lambda-\lambda_i(t)]
$$

이다. 한 deterministic trajectory에 대한 **정확한 정의**다. 그 trajectory에 존재하는 spatial heterogeneity는 담지만, 그 자체가 물리적 randomness를 입증하지는 않는다.

### ensemble-averaged probability target

$$
P_{\mathrm{ens}}(\lambda,t)
=
\mathbb E_{\Gamma_0}
\left[
\frac1M\sum_i\delta(\lambda-\lambda_i(t;\Gamma_0))
\right]
$$

이다. 이는 **미래 물리 ensemble 객체**다. initial phase-space ensemble $\Gamma_0$를 물리적으로 정의한 뒤에야 수치값을 Al fatigue probability로 해석할 수 있다.

fixed-length canonical ensemble, controlled metastable intact-basin ensemble, 독립적으로 뒷받침된 material heterogeneity 등은 서로 같은 것이 아니며 각각 별도로 정당화해야 한다.

## 현재 cycle-boundary의 정확한 기준

Milestone 13/15의 zero-mean sinusoidal loading에서는

$$
f(t)=A(t)f_a\sin(\omega t)
$$

이고 정수 cycle에서

$$
\omega t=2\pi N
$$

이므로

$$
f=0.
$$

따라서

$$
\lambda_s=1,
\qquad
C_0^{\rm qs}=0
$$

이다. 이 snapshot에서 nonzero $C_0$는 정확한 정적 기준에 대한 동적 잔류량이다.
