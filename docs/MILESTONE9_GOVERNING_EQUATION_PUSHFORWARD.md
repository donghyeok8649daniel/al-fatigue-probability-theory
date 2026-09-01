# Milestone 9 — Governing-Equation Push-Forward Clue for the 1D Spacing Distribution

> **SUPERSEDED AS ACTIVE GLOBAL-DISTRIBUTION ROUTE:** Sections below that use linear modes, finite harmonics, or a Taylor expansion about $\lambda=1$ are retained only as historical/local diagnostics. They are not used for the active full-support distribution theory. The current exact nonlinear route is `MILESTONE10_EXACT_DISTRIBUTION_TRANSPORT.md`.

<!-- M9_SUPERSEDED_EN -->

## Scope

The active model remains strictly one-dimensional, normal-only, and layer based. This milestone does not choose a new fitted probability family. Instead it asks what can be said about the form of the one-point spacing density directly from the existing deterministic layer-spacing field and the 1D layer-LJ governing equation.

The normalized layer spacing is

$$
\lambda_i(t)=a_i(t)/a_0.
$$

For interior spacings the dimensionless nearest-neighbor governing equation is

$$
\boxed{
\ddot\lambda_i
=
\phi'(\lambda_{i+1})
-2\phi'(\lambda_i)
+\phi'(\lambda_{i-1}).
}
$$

## 1. Exact push-forward representation

Introduce a continuum layer label $\xi\in[0,1]$ and let $\Lambda(\xi,t)$ denote the spacing field obtained in the represented large-system limit. Uniform weighting of the layer label gives the one-point spacing density as the push-forward measure

$$
\boxed{
p_\lambda(\lambda,t)
=
\int_0^1
\delta[\lambda-\Lambda(\xi,t)]\,d\xi.
}
$$

This is a **DEFINITION / EXACT PUSH-FORWARD REPRESENTATION** of the one-point density once the deterministic spacing field is given. It is not a statistical-equilibrium assumption.

If $\Lambda(\xi,t)$ is piecewise monotone, the delta-function change-of-variables identity gives

$$
\boxed{
p_\lambda(\lambda,t)
=
\sum_{\xi_j:\,\Lambda(\xi_j,t)=\lambda}
\frac{1}{|\partial_\xi\Lambda(\xi_j,t)|}.
}
$$

Therefore the shape of $p_\lambda$ is controlled by the spatial waveform and its slope. In particular, extrema of the spacing field, where $\partial_\xi\Lambda=0$, naturally create enhanced probability density near the corresponding spacing values.

This is the first new structural clue: **the functional form of the one-point distribution is a kinematic image of the spatial mode structure, not an independently chosen probability family.**

## 2. Linearized governing equation

Write

$$
\lambda_i=1+u_i.
$$

Because the calibrated generalized LJ model satisfies

$$
\phi'(1)=0,
\qquad
\phi''(1)=1,
$$

linearization gives

$$
\boxed{
\ddot u_i
=u_{i+1}-2u_i+u_{i-1}.
}
$$

A normal mode

$$
u_i(t)=A\cos(qi-\omega_q t+\theta)
$$

has the exact linear dispersion relation

$$
\boxed{
\omega_q^2=4\sin^2(q/2).
}
$$

## 3. Single-mode distribution implied by the governing equation

If the represented spacing field at a fixed time is dominated by one coherent linear mode and the spatial phase is sampled uniformly, then

$$
\Lambda=\mu+A\cos\vartheta,
\qquad
\vartheta\sim\mathrm{Uniform}(0,2\pi).
$$

The push-forward formula gives

$$
\boxed{
p_{\rm 1mode}(\lambda)
=
\frac{
\mathbf 1_{|\lambda-\mu|<|A|}
}{
\pi\sqrt{A^2-(\lambda-\mu)^2}
}.
}
$$

This arcsine form is **not introduced as a fitted distribution**. It is the exact spatial push-forward of a single linear normal mode.

Its moments are

$$
\operatorname{Var}(\lambda)=A^2/2,
$$

$$
\mu_3=0,
$$

and

$$
\frac{\mu_4}{\operatorname{Var}(\lambda)^2}=\frac32.
$$

Thus a single linear mode can explain strong non-Gaussian endpoint accumulation but cannot explain nonzero skewness.

## 4. The calibrated LJ equation predicts nonlinear harmonic distortion

Expand the force around equilibrium:

$$
\phi'(1+u)
=
c_1u+c_2u^2+c_3u^3+O(u^4).
$$

For the normalized generalized LJ form,

$$
c_1=1,
$$

$$
\boxed{
c_2=\frac12\phi'''(1)
=-\frac{m+n+3}{2},
}
$$

and

$$
\boxed{
c_3=\frac16\phi''''(1).}
$$

With the active calibration $m=12.19$, $n=6$,

$$
\boxed{
c_2=-10.595,}
$$

$$
\boxed{
c_3=62.97935.}
$$

Hence the governing equation itself contains strong quadratic and cubic waveform distortion even before any empirical fatigue law is introduced.

A minimal distorted spatial waveform is

$$
\Lambda(\vartheta)
=
\mu+A\cos\vartheta+B\cos2\vartheta.
$$

Its push-forward density is exactly

$$
\boxed{
p(\lambda)
=
\frac{1}{2\pi}
\sum_{\vartheta_j:\,\Lambda(\vartheta_j)=\lambda}
\frac{1}{| -A\sin\vartheta_j-2B\sin2\vartheta_j|}.
}
$$

The first low moments are

$$
\boxed{
\operatorname{Var}(\lambda)=\frac{A^2+B^2}{2},
}
$$

and

$$
\boxed{
\mu_3=\frac{3}{4}A^2B.
}
$$

Therefore skewness can arise directly from nonlinear harmonic distortion of the mechanically generated spacing field.

The two-harmonic family also has the exact bound

$$
\boxed{
|\gamma_1|\le\sqrt{\frac23}\approx0.816497.
}
$$

This bound is useful as a falsification diagnostic: a deterministic snapshot with larger $|\gamma_1|$ cannot be represented by only the first two spatial harmonics.

## 5. Comparison with the existing deterministic snapshots

Using no histogram fit, a single-mode arcsine reference can be fixed by the measured mean and variance through

$$
A=\sqrt{2\operatorname{Var}(\lambda)}.
$$

For the previously studied 32-node snapshots, the Kolmogorov distances are approximately

| deterministic snapshot | exponential mean-energy closure | single-mode push-forward |
|---|---:|---:|
| slow, $t\approx10T$ | $0.159$ | $0.260$ |
| dynamic, $t=2T$ | $0.146$ | $0.197$ |

Thus the pure single-mode arcsine law is **not** a better full-distribution closure for those snapshots. This is a useful negative result, not a reason to discard the push-forward structure.

The slow snapshot has empirical skewness

$$
\gamma_{1,\rm sim}\approx1.062,
$$

which even exceeds the two-harmonic upper bound $\sqrt{2/3}$. Therefore that state necessarily contains richer spatial structure than a first-plus-second-harmonic waveform.

## 6. What is learned about the form of $P$

The governing equations support the following hierarchy of statements.

1. **EXACT:** the one-point density is the spatial push-forward of the deterministic spacing field.
2. **EXACT under piecewise monotonicity:** the density is a sum of inverse spatial slopes over all preimages.
3. **LINEARIZED MECHANICS:** coherent normal modes generate arcsine-type endpoint structure.
4. **NONLINEAR LJ MECHANICS:** the large quadratic and cubic coefficients generate harmonic distortion and hence skewness/multimodality without inserting a named statistical distribution.
5. **NUMERICAL FALSIFICATION:** one mode is insufficient for the existing driven snapshots, and the slow snapshot is too skewed even for a two-harmonic representation.

Therefore the next probability-form search should not be

$$
\text{choose another }p(\lambda)\text{ family}.
$$

It should be

$$
\boxed{
\text{derive/resolve the mechanically excited spatial mode content}
\rightarrow
\Lambda(\xi,t)
\rightarrow
p_\lambda(\lambda,t)\text{ by push-forward}.
}
$$

This also explains why the previously observed $O(M)$ spatial correlation is directly relevant to the distribution form.

## 7. Next target

The next exact derivation should combine the push-forward representation with the pair-state hierarchy. In particular, derive how the pair field or mode amplitudes evolve from

$$
\ddot\lambda_i
=
\phi'(\lambda_{i+1})-2\phi'(\lambda_i)+\phi'(\lambda_{i-1}),
$$

and determine whether a small number of mechanically selected modes can reproduce the observed $P_1$ and $C_k$ simultaneously.

No empirical relaxation constant, damage variable, Gaussian assumption, Weibull assumption, or fitted correlation length is introduced here.

---

# 한국어 번역 — 1D Spacing Distribution 형식에 대한 지배방정식 Push-Forward 단서

> **활성 전역분포 경로에서 제외됨:** 아래의 linear mode, finite harmonic, $\lambda=1$ 주변 Taylor expansion 부분은 historical/local diagnostic으로만 보존한다. 전체 support의 활성 distribution theory에는 사용하지 않는다. 현재 exact nonlinear 경로는 `MILESTONE10_EXACT_DISTRIBUTION_TRANSPORT.md`다.

<!-- M9_SUPERSEDED_KO -->

## 범위

활성 모델은 계속 엄격하게 1차원, normal-only, layer 기반이다. 이번 milestone에서는 새로운 fitted probability family를 선택하지 않는다. 대신 현재 가지고 있는 deterministic layer-spacing field와 1D layer-LJ 지배방정식만으로 one-point spacing density의 형식에 대해 무엇을 말할 수 있는지 본다.

normalized layer spacing은

$$
\lambda_i(t)=a_i(t)/a_0
$$

이고 interior spacing의 dimensionless nearest-neighbor 지배방정식은

$$
\boxed{
\ddot\lambda_i
=
\phi'(\lambda_{i+1})
-2\phi'(\lambda_i)
+\phi'(\lambda_{i-1})
}
$$

이다.

## 1. 정확한 push-forward 표현

continuum layer label $\xi\in[0,1]$을 도입하고 represented large-system limit에서 spacing field를 $\Lambda(\xi,t)$라고 하자. layer label에 uniform weight를 주면 one-point spacing density는

$$
\boxed{
p_\lambda(\lambda,t)
=
\int_0^1
\delta[\lambda-\Lambda(\xi,t)]\,d\xi
}
$$

라는 push-forward measure가 된다.

이는 deterministic spacing field가 주어졌을 때 one-point density의 **DEFINITION / EXACT PUSH-FORWARD REPRESENTATION**이다. statistical-equilibrium assumption이 아니다.

$\Lambda(\xi,t)$가 piecewise monotone이면 delta-function change-of-variables identity로

$$
\boxed{
p_\lambda(\lambda,t)
=
\sum_{\xi_j:\,\Lambda(\xi_j,t)=\lambda}
\frac{1}{|\partial_\xi\Lambda(\xi_j,t)|}
}
$$

를 얻는다.

따라서 $p_\lambda$의 shape은 spatial waveform과 그 기울기에 의해 정해진다. 특히 $\partial_\xi\Lambda=0$인 spacing-field extremum에서는 해당 spacing 값 근처의 probability density가 자연스럽게 증가한다.

이게 새로운 첫 구조적 단서다. **one-point distribution의 함수형식은 독립적으로 고른 probability family가 아니라 spatial mode structure의 kinematic image다.**

## 2. 선형화된 지배방정식

$$
\lambda_i=1+u_i
$$

로 두자.

calibrated generalized LJ model은

$$
\phi'(1)=0,
\qquad
\phi''(1)=1
$$

을 만족하므로 선형화하면

$$
\boxed{
\ddot u_i
=u_{i+1}-2u_i+u_{i-1}
}
$$

를 얻는다.

normal mode

$$
u_i(t)=A\cos(qi-\omega_q t+\theta)
$$

의 exact linear dispersion은

$$
\boxed{
\omega_q^2=4\sin^2(q/2)
}
$$

이다.

## 3. 지배방정식이 주는 single-mode distribution

고정 시간에서 represented spacing field가 하나의 coherent linear mode에 의해 지배되고 spatial phase를 uniform하게 sample한다면

$$
\Lambda=\mu+A\cos\vartheta,
\qquad
\vartheta\sim\mathrm{Uniform}(0,2\pi)
$$

이다.

push-forward formula로

$$
\boxed{
p_{\rm 1mode}(\lambda)
=
\frac{
\mathbf 1_{|\lambda-\mu|<|A|}
}{
\pi\sqrt{A^2-(\lambda-\mu)^2}
}
}
$$

를 얻는다.

이 arcsine form은 fitted distribution으로 도입한 것이 아니다. single linear normal mode의 exact spatial push-forward다.

moment는

$$
\operatorname{Var}(\lambda)=A^2/2,
$$

$$
\mu_3=0,
$$

그리고

$$
\frac{\mu_4}{\operatorname{Var}(\lambda)^2}=\frac32
$$

이다.

따라서 single linear mode는 강한 non-Gaussian endpoint accumulation은 설명할 수 있지만 nonzero skewness는 설명할 수 없다.

## 4. Calibrated LJ 지배방정식은 nonlinear harmonic distortion을 예측한다

평형점 주변 force를

$$
\phi'(1+u)
=
c_1u+c_2u^2+c_3u^3+O(u^4)
$$

로 전개한다.

normalized generalized LJ에서는

$$
c_1=1,
$$

$$
\boxed{
c_2=\frac12\phi'''(1)
=-\frac{m+n+3}{2}
}
$$

이고

$$
\boxed{
c_3=\frac16\phi''''(1)}
$$

이다.

현재 calibration $m=12.19$, $n=6$를 넣으면

$$
\boxed{c_2=-10.595}
$$

및

$$
\boxed{c_3=62.97935}
$$

를 얻는다.

따라서 경험적 fatigue law를 하나도 넣지 않아도 지배방정식 자체에 강한 quadratic/cubic waveform distortion 항이 존재한다.

가장 단순한 distorted spatial waveform은

$$
\Lambda(\vartheta)
=
\mu+A\cos\vartheta+B\cos2\vartheta
$$

이다.

그 push-forward density는 정확히

$$
\boxed{
p(\lambda)
=
\frac{1}{2\pi}
\sum_{\vartheta_j:\,\Lambda(\vartheta_j)=\lambda}
\frac{1}{| -A\sin\vartheta_j-2B\sin2\vartheta_j|}
}
$$

이고 low moment는

$$
\boxed{
\operatorname{Var}(\lambda)=\frac{A^2+B^2}{2}
}
$$

및

$$
\boxed{
\mu_3=\frac{3}{4}A^2B
}
$$

이다.

따라서 skewness는 mechanically generated spacing field의 nonlinear harmonic distortion만으로 직접 발생할 수 있다.

또한 two-harmonic family는 정확히

$$
\boxed{
|\gamma_1|\le\sqrt{\frac23}\approx0.816497
}
$$

라는 bound를 가진다.

이 bound는 좋은 falsification diagnostic이다. deterministic snapshot의 $|\gamma_1|$가 이보다 크면 first+second harmonic만으로는 그 상태를 표현할 수 없다.

## 5. 기존 deterministic snapshot과의 비교

histogram fitting 없이 measured mean과 variance만 사용하면 single-mode amplitude는

$$
A=\sqrt{2\operatorname{Var}(\lambda)}
$$

로 정해진다.

이전에 사용한 32-node snapshot에서 Kolmogorov distance는 대략 다음과 같다.

| deterministic snapshot | exponential mean-energy closure | single-mode push-forward |
|---|---:|---:|
| slow, $t\approx10T$ | $0.159$ | $0.260$ |
| dynamic, $t=2T$ | $0.146$ | $0.197$ |

따라서 pure single-mode arcsine law는 이 snapshot들의 full-distribution closure로 더 좋은 결과를 주지 않는다. 이는 유용한 negative result이며 push-forward structure 자체를 버릴 이유는 아니다.

slow snapshot의 empirical skewness는

$$
\gamma_{1,\rm sim}\approx1.062
$$

로 two-harmonic upper bound $\sqrt{2/3}$보다도 크다. 따라서 이 상태에는 first+second harmonic보다 더 복잡한 spatial structure가 반드시 존재한다.

## 6. $P$의 형식에 대해 새로 얻은 것

지배방정식으로부터 다음 hierarchy를 얻는다.

1. **EXACT:** one-point density는 deterministic spacing field의 spatial push-forward다.
2. **piecewise monotonicity 아래 EXACT:** density는 같은 spacing 값을 만드는 모든 preimage에서 inverse spatial slope를 합친 형태다.
3. **LINEARIZED MECHANICS:** coherent normal mode는 arcsine-type endpoint structure를 만든다.
4. **NONLINEAR LJ MECHANICS:** 큰 quadratic/cubic coefficient가 harmonic distortion과 이에 따른 skewness/multimodality를 만든다. named statistical distribution을 넣을 필요가 없다.
5. **NUMERICAL FALSIFICATION:** 기존 driven snapshot에는 one mode가 부족하고, slow snapshot은 two-harmonic representation으로도 skewness가 너무 크다.

따라서 다음 probability-form 탐색은

$$
\text{다른 }p(\lambda)\text{ family를 고르는 것}
$$

이 아니라

$$
\boxed{
\text{mechanically excited spatial mode content를 유도/해석}
\rightarrow
\Lambda(\xi,t)
\rightarrow
p_\lambda(\lambda,t)\text{를 push-forward로 계산}
}
$$

하는 방향이어야 한다.

이 결과는 앞서 관찰한 $O(M)$ spatial correlation이 왜 distribution 형식과 직접 연결되는지도 설명한다.

## 7. 다음 목표

다음 exact derivation은 push-forward representation과 pair-state hierarchy를 결합해야 한다. 특히

$$
\ddot\lambda_i
=
\phi'(\lambda_{i+1})-2\phi'(\lambda_i)+\phi'(\lambda_{i-1})
$$

에서 pair field 또는 mode amplitude가 어떻게 진화하는지 유도하고, mechanically selected mode 몇 개만으로 관찰된 $P_1$과 $C_k$를 동시에 재현할 수 있는지 확인한다.

여기에는 empirical relaxation constant, damage variable, Gaussian assumption, Weibull assumption, fitted correlation length를 도입하지 않는다.
