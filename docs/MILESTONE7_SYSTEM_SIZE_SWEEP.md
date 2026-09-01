# Milestone 7 — Dynamically Matched System-Size Sweep of the 1D Layer-LJ Closure

## Scope

The active model remains strictly one-dimensional, normal-only, and based on the calibrated generalized Lennard-Jones interaction between represented material layers.

Milestone 6 showed that the two-moment closure

$$
p_\lambda(\lambda,t)
=
Z^{-1}
\exp[-\alpha(t)\lambda-\beta(t)\psi(\lambda)]
$$

does not reproduce the full deterministic spacing distribution at $M=31$ represented spacings, even when the empirical mean stretch and configurational energy are matched exactly.

The present question is narrower:

$$
\boxed{
\text{Does the mismatch vanish as }M\text{ increases?}
}
$$

## 1. Why a naive fixed-frequency sweep is not sufficient

The chain is driven from a boundary. Therefore changing the chain length while keeping the same dimensionless angular frequency changes the ratio between loading period and acoustic transit time.

For a sample taken at two periods,

$$
t_s=2T=\frac{4\pi}{\omega}.
$$

The relevant normalized transit ratio scales as

$$
\boxed{
\frac{t_s}{M}
=
\frac{4\pi}{\omega M}.
}
$$

Thus a fixed-$\omega$ sweep changes the number of chain transits and reflections before the snapshot. Such a sweep mixes finite-$M$ statistics with a different dynamical state.

This observation is an **IDENTITY / SCALING ARGUMENT** under the normalized 1D chain geometry.

## 2. Controlled dynamic-similarity rule

To isolate represented-system-size effects more cleanly, the sweep keeps

$$
\boxed{
\omega M=\chi
}
$$

constant.

The reference state uses

$$
M_0=31,
\qquad
\omega_0=0.02,
$$

so

$$
\boxed{
\chi=\omega_0M_0=0.62.
}
$$

For every represented spacing count $M$,

$$
\boxed{
\omega(M)=\frac{0.62}{M}.
}
$$

The normal force amplitude remains

$$
F_a^*=0.03,
$$

the loading ramp remains two periods, and the snapshot is taken at

$$
t_s=2T.
$$

Because $\omega M$ is fixed,

$$
\frac{t_s}{M}
=
\frac{4\pi}{0.62}
\approx20.26834
$$

is also fixed.

**Classification:** this is a **CONTROLLED NUMERICAL SCALING PROTOCOL**, not a material law.

## 3. Numerical-resolution correction used in this sweep

Large $M$ produces narrow closure states with large $\alpha$ and $\beta$. A fixed global quadrature can miss the narrow probability peak.

The closure solver was therefore corrected before interpreting the sweep. For sharply concentrated states, the numerical integration is centered on the mode $\lambda_*$ determined by

$$
\alpha+\beta\psi'(\lambda_*)=0,
$$

with a width obtained from

$$
s_*=
[\beta\psi''(\lambda_*)]^{-1/2}.
$$

This correction changes only the numerical resolution of the same closure. It adds no physical parameter.

## 4. Dynamically matched results

The tested chains contain 32, 64, 128, and 256 represented layer nodes, corresponding to

$$
M=31,\ 63,\ 127,\ 255
$$

spacings.

| $M$ | $\omega$ | variance error | empirical skewness | closure skewness | $D_{\rm KS}$ |
|---:|---:|---:|---:|---:|---:|
| 31 | $2.0000\times10^{-2}$ | $1.65\%$ | $0.55592$ | $0.10321$ | $0.12933$ |
| 63 | $9.84127\times10^{-3}$ | $1.79\%$ | $0.69280$ | $0.08552$ | $0.14510$ |
| 127 | $4.88189\times10^{-3}$ | $1.82\%$ | $0.76759$ | $0.07692$ | $0.15824$ |
| 255 | $2.43137\times10^{-3}$ | $1.82\%$ | $0.80674$ | $0.07269$ | $0.16531$ |

The exact machine-readable values are stored in

`results/data/normal_lj_closure_system_size.json`.

## 5. Main result

The Kolmogorov distance does not decrease over the tested dynamically matched sequence.

Instead,

$$
D_{\rm KS}
:
0.1293
\rightarrow
0.1451
\rightarrow
0.1582
\rightarrow
0.1653.
$$

At the same time, the deterministic skewness increases,

$$
\gamma_{1,\rm sim}
:
0.5559
\rightarrow
0.8067,
$$

while the two-moment closure remains much less skewed,

$$
\gamma_{1,\rm closure}
:
0.1032
\rightarrow
0.0727.
$$

Therefore the tested data do not support the hypothesis

$$
\boxed{
D_{\rm KS}(M)\to0
\text{ merely because }M\to\infty.
}
$$

This is a **NUMERICAL FALSIFICATION RESULT on the tested dynamically matched sequence**. It is not a proof for every possible large-$M$ driving protocol.

## 6. Exploratory inverse-$M$ extrapolation

A linear fit in $1/M$ over only these four points gives

$$
D_{\rm KS}(M\to\infty)
\approx0.1682,
$$

$$
\gamma_{1,\rm sim}(M\to\infty)
\approx0.8377,
$$

and

$$
\gamma_{1,\rm closure}(M\to\infty)
\approx0.06845.
$$

These values are explicitly classified as

**EXPLORATORY NUMERICAL EXTRAPOLATION — NOT A THEOREM.**

The important result is not the intercept itself. The important result is that no trend toward zero mismatch is observed.

## 7. Scientific interpretation

Finite-$M$ saddle-point error alone is no longer a convincing explanation of the closure failure.

The deterministic chain carries information that is not specified by

$$
\mu(t)
\quad\text{and}\quad
\mathcal E(t).
$$

Two leading possibilities remain:

1. a one-point shape variable such as
   $$
   m_3(t)
   =
   \int(\lambda-\mu)^3p_\lambda\,d\lambda;
   $$
2. persistent spatial information between neighboring layer spacings.

The second possibility is especially natural because the exact 1D mechanics contains

$$
\ddot\lambda_i
\propto
\phi'(\lambda_{i+1})
-2\phi'(\lambda_i)
+\phi'(\lambda_{i-1}),
$$

so neighboring spacings are dynamically coupled even though the one-point closure does not retain their ordering.

## 8. Next decisive observable

Before adding a third-moment multiplier to the probability law, measure the nearest-neighbor spacing covariance

$$
\boxed{
C_1(t)
=
\frac{1}{M-1}
\sum_{i=1}^{M-1}
[\lambda_i(t)-\mu(t)]
[\lambda_{i+1}(t)-\mu(t)].
}
$$

More generally,

$$
\boxed{
C_k(t)
=
\frac{1}{M-k}
\sum_{i=1}^{M-k}
[\lambda_i(t)-\mu(t)]
[\lambda_{i+k}(t)-\mu(t)].
}
$$

These are direct observables of the existing deterministic 1D layer-LJ mechanics. No fitted constitutive law is introduced.

The next question is

$$
\boxed{
M\uparrow
\quad\Longrightarrow\quad
C_k(t)\ ?
}
$$

If a nonzero correlation structure survives the dynamically matched large-$M$ sequence, the equal-base-measure one-point closure is missing spatial state information.

---

# 한국어 번역 — 1D Layer-LJ Closure의 동적 유사 System-Size Sweep

## 범위

활성 모델은 계속 엄격하게 1차원, 수직변형 전용이며 represented material layer 사이의 calibration된 generalized Lennard-Jones 상호작용을 사용한다.

Milestone 6에서는 two-moment closure

$$
p_\lambda(\lambda,t)
=
Z^{-1}
\exp[-\alpha(t)\lambda-\beta(t)\psi(\lambda)]
$$

가 empirical mean stretch와 configurational energy를 정확히 맞춰도 $M=31$ deterministic spacing distribution 전체를 재현하지 못한다는 것을 확인했다.

이번 질문은 더 좁다.

$$
\boxed{
M\text{을 증가시키면 mismatch가 사라지는가?}
}
$$

## 1. 단순 fixed-frequency sweep가 충분하지 않은 이유

chain은 boundary에서 구동된다. 따라서 dimensionless angular frequency를 고정한 채 chain length만 바꾸면 loading period와 acoustic transit time의 비가 달라진다.

두 period 뒤에 sample한다면

$$
t_s=2T=\frac{4\pi}{\omega}
$$

이고 normalized transit ratio는

$$
\boxed{
\frac{t_s}{M}
=
\frac{4\pi}{\omega M}
}
$$

로 scaling한다.

따라서 fixed-$\omega$ sweep는 snapshot 이전 chain 왕복 및 reflection 횟수를 바꾼다. 이런 sweep는 finite-$M$ statistics와 서로 다른 dynamical state를 섞는다.

이 관찰은 normalized 1D chain geometry 아래 **IDENTITY / SCALING ARGUMENT**다.

## 2. Controlled dynamic-similarity rule

represented-system-size effect를 더 깨끗하게 분리하기 위해

$$
\boxed{
\omega M=\chi
}
$$

를 일정하게 유지한다.

reference state는

$$
M_0=31,
\qquad
\omega_0=0.02
$$

이므로

$$
\boxed{
\chi=\omega_0M_0=0.62
}
$$

다.

각 represented spacing 수 $M$에서

$$
\boxed{
\omega(M)=\frac{0.62}{M}
}
$$

로 둔다.

normal force amplitude는

$$
F_a^*=0.03
$$

으로 유지하고, loading ramp는 두 period, snapshot은

$$
t_s=2T
$$

에서 취한다.

$\omega M$이 일정하므로

$$
\frac{t_s}{M}
=
\frac{4\pi}{0.62}
\approx20.26834
$$

도 일정하다.

**분류:** 이는 material law가 아니라 **CONTROLLED NUMERICAL SCALING PROTOCOL**이다.

## 3. 이번 sweep에 사용한 numerical-resolution correction

큰 $M$에서는 $\alpha,\beta$가 큰 매우 좁은 closure state가 나타난다. fixed global quadrature는 이 좁은 probability peak를 놓칠 수 있다.

따라서 sweep를 해석하기 전에 closure solver의 numerical resolution을 보정했다. sharply concentrated state에서는

$$
\alpha+\beta\psi'(\lambda_*)=0
$$

으로 정해지는 mode $\lambda_*$를 중심으로 적분하고,

$$
s_*=
[\beta\psi''(\lambda_*)]^{-1/2}
$$

에서 numerical width를 얻는다.

이 보정은 동일한 closure의 numerical resolution만 바꾸며 새로운 physical parameter를 추가하지 않는다.

## 4. Dynamically matched 결과

시험한 chain은 32, 64, 128, 256개의 represented layer node를 가지며 spacing 수는

$$
M=31,\ 63,\ 127,\ 255
$$

이다.

| $M$ | $\omega$ | variance error | empirical skewness | closure skewness | $D_{\rm KS}$ |
|---:|---:|---:|---:|---:|---:|
| 31 | $2.0000\times10^{-2}$ | $1.65\%$ | $0.55592$ | $0.10321$ | $0.12933$ |
| 63 | $9.84127\times10^{-3}$ | $1.79\%$ | $0.69280$ | $0.08552$ | $0.14510$ |
| 127 | $4.88189\times10^{-3}$ | $1.82\%$ | $0.76759$ | $0.07692$ | $0.15824$ |
| 255 | $2.43137\times10^{-3}$ | $1.82\%$ | $0.80674$ | $0.07269$ | $0.16531$ |

정확한 machine-readable 값은

`results/data/normal_lj_closure_system_size.json`

에 저장한다.

## 5. 핵심 결과

시험한 dynamically matched sequence에서 Kolmogorov distance는 감소하지 않았다.

오히려

$$
D_{\rm KS}
:
0.1293
\rightarrow
0.1451
\rightarrow
0.1582
\rightarrow
0.1653
$$

로 움직인다.

동시에 deterministic skewness는

$$
\gamma_{1,\rm sim}
:
0.5559
\rightarrow
0.8067
$$

로 증가하지만 two-moment closure는

$$
\gamma_{1,\rm closure}
:
0.1032
\rightarrow
0.0727
$$

수준에 머문다.

따라서 시험한 데이터는

$$
\boxed{
M\to\infty
\text{라는 이유만으로 }
D_{\rm KS}(M)\to0
}
$$

이라는 가설을 지지하지 않는다.

이는 시험한 dynamically matched sequence에 대한 **NUMERICAL FALSIFICATION RESULT**다. 모든 possible large-$M$ driving protocol에 대한 증명은 아니다.

## 6. Exploratory inverse-$M$ extrapolation

네 점에 대해 단순히 $1/M$ linear fit을 하면

$$
D_{\rm KS}(M\to\infty)
\approx0.1682,
$$

$$
\gamma_{1,\rm sim}(M\to\infty)
\approx0.8377,
$$

그리고

$$
\gamma_{1,\rm closure}(M\to\infty)
\approx0.06845
$$

가 나온다.

이 값은 명시적으로

**EXPLORATORY NUMERICAL EXTRAPOLATION — NOT A THEOREM**

으로 분류한다.

중요한 것은 intercept 값 자체가 아니라 mismatch가 0으로 향하는 경향이 관찰되지 않는다는 점이다.

## 7. 과학적 해석

finite-$M$ saddle-point error만으로 closure failure를 설명하기는 이제 설득력이 약하다.

deterministic chain에는

$$
\mu(t)
\quad\text{와}\quad
\mathcal E(t)
$$

가 지정하지 않는 정보가 남아 있다.

주요 후보는 두 가지다.

1. one-point shape variable
   $$
   m_3(t)
   =
   \int(\lambda-\mu)^3p_\lambda\,d\lambda;
   $$
2. neighboring layer spacing 사이의 지속적인 spatial information.

두 번째 후보가 특히 자연스러운 이유는 exact 1D mechanics가

$$
\ddot\lambda_i
\propto
\phi'(\lambda_{i+1})
-2\phi'(\lambda_i)
+\phi'(\lambda_{i-1})
$$

처럼 neighboring spacing을 직접 coupling하기 때문이다. 반면 one-point closure는 spacing의 spatial ordering을 보존하지 않는다.

## 8. 다음 결정적 observable

probability law에 third-moment multiplier를 바로 추가하기 전에 nearest-neighbor spacing covariance를 측정한다.

$$
\boxed{
C_1(t)
=
\frac{1}{M-1}
\sum_{i=1}^{M-1}
[\lambda_i(t)-\mu(t)]
[\lambda_{i+1}(t)-\mu(t)]
}
$$

를 정의한다.

더 일반적으로

$$
\boxed{
C_k(t)
=
\frac{1}{M-k}
\sum_{i=1}^{M-k}
[\lambda_i(t)-\mu(t)]
[\lambda_{i+k}(t)-\mu(t)]
}
$$

이다.

이 값들은 기존 deterministic 1D layer-LJ mechanics에서 직접 측정되는 observable이다. fitted constitutive law는 추가하지 않는다.

다음 질문은

$$
\boxed{
M\uparrow
\quad\Longrightarrow\quad
C_k(t)\ ?
}
$$

이다.

dynamically matched large-$M$ sequence에서 nonzero correlation structure가 유지된다면 equal-base-measure one-point closure가 spatial state information을 놓치고 있다는 뜻이다.
