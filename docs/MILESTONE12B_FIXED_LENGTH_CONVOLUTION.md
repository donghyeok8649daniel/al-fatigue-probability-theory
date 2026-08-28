# Milestone 12B — Exact Fixed-Length Canonical $P$ as a Convolution Recursion

## Purpose

Milestone 12 derived the exact finite-$M$ canonical one-spacing marginal at fixed total normalized length $L$,

$$
P_M(\lambda\mid L,\chi)
=
\frac{e^{-\chi\phi(\lambda)}Z_{M-1}(L-\lambda,\chi)}{Z_M(L,\chi)}.
$$

The remaining partition function is not a fitted unknown function. Because the active 1D model has an additive spacing energy

$$
V(\lambda_1,\ldots,\lambda_M)
=E_0\sum_{i=1}^M\phi(\lambda_i),
$$

it obeys an exact one-dimensional convolution recursion.

## Exact recursion

Define the one-spacing Boltzmann weight

$$
z_\chi(\lambda)
=
\exp[-\chi\phi(\lambda)]\mathbf 1_{\lambda>0}.
$$

Then

$$
Z_1(L,\chi)=z_\chi(L),
$$

and for $M\ge2$,

$$
\boxed{
Z_M(L,\chi)
=
\int_0^L
z_\chi(\lambda)
Z_{M-1}(L-\lambda,\chi)
\,d\lambda.
}
$$

Equivalently,

$$
\boxed{Z_M=z_\chi^{*M},}
$$

where $*$ is convolution on the positive-spacing axis.

Therefore the exact equilibrium function form is computationally explicit:

$$
\boxed{
P_M(\lambda\mid L,\chi)
=
\frac{
z_\chi(\lambda)
Z_{M-1}(L-\lambda,\chi)
}{Z_M(L,\chi)}.
}
$$

No Gaussian, Weibull, polynomial potential, harmonic mode, or saddle-point approximation is required.

## Symmetry and exact mean

All $M$ spacings are exchangeable in this equilibrium ensemble. Since

$$
\sum_{i=1}^M\lambda_i=L
$$

holds exactly,

$$
\boxed{
\mathbb E[\lambda_i]=\frac{L}{M}
}
$$

for every spacing. This is an exact diagnostic for any numerical convolution implementation.

For $M=2$,

$$
\boxed{
P_2(\lambda\mid L,\chi)
\propto
\exp\{-\chi[\phi(\lambda)+\phi(L-\lambda)]\},
\qquad0<\lambda<L,
}
$$

which is symmetric under $\lambda\mapsto L-\lambda$.

## Numerical implementation

`theory/normal_lj_physical_distribution.py` evaluates the recursion on a uniform $\lambda$ grid using FFT convolution. The potential is shifted numerically as

$$
\psi(\lambda)=\phi(\lambda)-\phi(1)
$$

before exponentiation. For fixed $M$, this multiplies $Z_M$ by a global constant only and therefore cancels exactly from the normalized marginal $P_M$.

The FFT/grid calculation is a **NUMERICAL APPROXIMATION of the exact convolution integral** and must be checked by grid refinement. The unit tests verify the exact $M=2$ formula and the mean identity $\langle\lambda\rangle=L/M$ for a finite-$M$ case.

## Physical interpretation

This result provides a concrete candidate function $P$ whenever the system can physically be treated as canonical and fixed-length. It does not say that a cyclically driven conservative chain is automatically canonical. The ensemble must be tested against deterministic mechanics.

The next comparison should therefore use measured $L(t)$ and a physically justified or independently inferred $\chi$, calculate the full finite-$M$ canonical $P_M$, and compare it directly with the deterministic spacing distribution. No distribution-shape parameter should be fitted.

---

# 한국어 번역 — Exact Fixed-Length Canonical $P$의 Convolution Recursion

## 목적

Milestone 12에서는 total normalized length $L$가 고정된 finite-$M$ canonical one-spacing marginal

$$
P_M(\lambda\mid L,\chi)
=
\frac{e^{-\chi\phi(\lambda)}Z_{M-1}(L-\lambda,\chi)}{Z_M(L,\chi)}
$$

를 유도했다.

여기 남은 partition function은 fitting해야 하는 미지함수가 아니다. active 1D model의 spacing energy가

$$
V(\lambda_1,\ldots,\lambda_M)
=E_0\sum_{i=1}^M\phi(\lambda_i)
$$

처럼 additive이므로 exact 1D convolution recursion을 만족한다.

## Exact recursion

one-spacing Boltzmann weight를

$$
z_\chi(\lambda)
=
\exp[-\chi\phi(\lambda)]\mathbf 1_{\lambda>0}
$$

로 정의한다.

그러면

$$
Z_1(L,\chi)=z_\chi(L)
$$

이고 $M\ge2$에서

$$
\boxed{
Z_M(L,\chi)
=
\int_0^L
z_\chi(\lambda)
Z_{M-1}(L-\lambda,\chi)
\,d\lambda
}
$$

이다.

즉

$$
\boxed{Z_M=z_\chi^{*M}}
$$

이며 $*$는 positive-spacing axis의 convolution이다.

따라서 equilibrium에서 $P$의 함수형은 실제 계산 가능한 형태로

$$
\boxed{
P_M(\lambda\mid L,\chi)
=
\frac{
z_\chi(\lambda)
Z_{M-1}(L-\lambda,\chi)
}{Z_M(L,\chi)}
}
$$

가 된다.

Gaussian, Weibull, polynomial potential, harmonic mode, saddle-point approximation은 필요하지 않다.

## 대칭성과 exact mean

이 equilibrium ensemble에서 모든 $M$ spacing은 exchangeable하다. 또한

$$
\sum_{i=1}^M\lambda_i=L
$$

이 exact이므로 모든 spacing에 대해

$$
\boxed{
\mathbb E[\lambda_i]=\frac{L}{M}
}
$$

이다. 이는 convolution 수치구현의 exact diagnostic으로 쓸 수 있다.

$M=2$에서는

$$
\boxed{
P_2(\lambda\mid L,\chi)
\propto
\exp\{-\chi[\phi(\lambda)+\phi(L-\lambda)]\},
\qquad0<\lambda<L
}
$$

이고 $\lambda\mapsto L-\lambda$에 대해 대칭이다.

## 수치 구현

`theory/normal_lj_physical_distribution.py`에서는 uniform $\lambda$ grid에서 FFT convolution으로 recursion을 계산한다. numerical exponentiation에서는

$$
\psi(\lambda)=\phi(\lambda)-\phi(1)
$$

로 potential을 shift한다. fixed $M$에서 이 shift는 $Z_M$에 global constant만 곱하므로 normalized marginal $P_M$에서는 exact하게 소거된다.

FFT/grid 계산은 **exact convolution integral의 NUMERICAL APPROXIMATION**이므로 grid refinement가 필요하다. unit test에서는 exact $M=2$ 식과 finite-$M$의 mean identity $\langle\lambda\rangle=L/M$을 검증한다.

## 물리적 해석

이 결과는 system을 canonical fixed-length로 물리적으로 취급할 수 있을 때 구체적인 $P$ 함수 후보를 준다. cyclically driven conservative chain이 자동으로 canonical이라는 뜻은 아니다. ensemble validity를 deterministic mechanics와 비교해야 한다.

따라서 다음 비교에서는 measured $L(t)$와 physically justified 또는 independently inferred $\chi$를 사용해 full finite-$M$ canonical $P_M$을 계산하고 deterministic spacing distribution과 직접 비교해야 한다. distribution-shape parameter는 fitting하지 않는다.
