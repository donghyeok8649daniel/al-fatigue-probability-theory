# Milestone 5 — Derived Distribution Closure for the 1D Layer-LJ Theory

## Status

This milestone introduces the first explicit mathematical form for the active spacing distribution without choosing a named probability family.

The active model remains strictly **one-dimensional, normal-only, and layer based**. The reduced microscopic coordinate is the normal separation between neighboring represented layers,

$$
a_i(t)>0.
$$

The effective normal interaction between layers is represented by the already calibrated generalized Lennard-Jones energy. No shear coordinate, three-dimensional FCC kinematics, cycle-dependent damage variable, or fitted probability family is introduced here.

The main result is the large-system closure

$$
\boxed{
p_\lambda(\lambda,t)
=
\frac{1}{Z(t)}
\exp\left[-\alpha(t)\lambda-\beta(t)\psi(\lambda)\right],
}
$$

where $\lambda=a/a_0$ and $\psi$ is the shifted normalized layer-LJ energy.

This result is **not exact nonequilibrium cyclic dynamics**. It is a **CONTROLLED APPROXIMATION** obtained from a precisely stated fixed-length/fixed-configurational-energy ensemble and a large-$M$ saddle-point reduction.

## 1. Physical and normalized spacing densities

Let $P_a(a,t)$ denote the physical spacing density,

$$
\int_0^\infty P_a(a,t)\,da=1.
$$

Introduce

$$
\lambda=\frac{a}{a_0},
$$

and the corresponding dimensionless density

$$
\boxed{
p_\lambda(\lambda,t)=a_0P_a(a_0\lambda,t).
}
$$

Then

$$
\int_0^\infty p_\lambda(\lambda,t)\,d\lambda=1.
$$

The mean normalized layer spacing is

$$
\boxed{
\mu(t)=\int_0^\infty \lambda p_\lambda(\lambda,t)\,d\lambda.
}
$$

## 2. Layer-LJ energy used by the closure

The current normalized generalized-LJ form is

$$
\phi(\lambda)
=
\frac{\lambda^{-m}}{m(m-n)}
-
\frac{\lambda^{-n}}{n(m-n)},
$$

with

$$
m=12.19,
\qquad
n=6.
$$

Define the shifted energy

$$
\boxed{
\psi(\lambda)=\phi(\lambda)-\phi(1).
}
$$

Therefore

$$
\psi(1)=0,
$$

and the mean shifted configurational energy is

$$
\boxed{
\mathcal E(t)
=
\int_0^\infty
\psi(\lambda)p_\lambda(\lambda,t)\,d\lambda.
}
$$

The local tangent-instability stretch remains

$$
\boxed{
\lambda_c
=
\left(\frac{m+1}{n+1}\right)^{1/(m-n)}
\approx1.1077715386.
}
$$

## 3. Finite-$M$ density-of-states construction

Consider $M$ positive layer spacings $\lambda_1,\ldots,\lambda_M$ at a fixed physical time $t$.

Define the instantaneous total normalized length and configurational energy constraints

$$
L(t)=\sum_{i=1}^{M}\lambda_i,
$$

$$
E_c(t)=\sum_{i=1}^{M}\psi(\lambda_i).
$$

Equivalently,

$$
L(t)=M\mu(t),
\qquad
E_c(t)=M\mathcal E(t).
$$

Define the density of admissible configurations

$$
\boxed{
\Omega_M(L,E)
=
\int_{(0,\infty)^M}
\delta\left(\sum_i\lambda_i-L\right)
\delta\left(\sum_i\psi(\lambda_i)-E\right)
\,d^M\lambda.
}
$$

### Ensemble assumption

**ASSUMPTION:** conditional on the instantaneous values $L(t)$ and $E_c(t)$, configurations on this constraint manifold are assigned equal base measure in the reduced layer-spacing coordinates.

This is the additional statistical closure assumption of this milestone. It is not claimed to follow automatically from the deterministic cyclic dynamics.

Under this assumption, the exact one-spacing marginal for finite $M$ is

$$
\boxed{
p_M(\lambda\mid L,E)
=
\frac{
\Omega_{M-1}\left(L-\lambda,E-\psi(\lambda)\right)
}{
\Omega_M(L,E)
}.
}
$$

## 4. Large-$M$ saddle-point closure

For large $M$, expand

$$
\ln\Omega_{M-1}
\left(L-\lambda,E-\psi(\lambda)\right)
$$

about the macroscopic state $(L,E)$.

To first order,

$$
\ln\Omega_{M-1}
\approx
C
-
\alpha\lambda
-
\beta\psi(\lambda),
$$

where $\alpha$ and $\beta$ are the conjugate derivatives of the large-system log density of states.

Therefore

$$
\boxed{
p_\lambda(\lambda,t)
=
\frac{
\exp[-\alpha(t)\lambda-\beta(t)\psi(\lambda)]
}{
Z(\alpha(t),\beta(t))
},
}
$$

with

$$
\boxed{
Z(\alpha,\beta)
=
\int_0^\infty
\exp[-\alpha\lambda-\beta\psi(\lambda)]\,d\lambda.
}
$$

For the present unbounded domain, normalizability requires

$$
\alpha>0,
\qquad
\beta>0.
$$

The parameter $\beta$ is an energy-conjugate multiplier of this closure. It must **not** be identified with $1/(k_BT)$ unless a separate equilibrium derivation justifies that identification.

## 5. The multipliers are determined by the moments

The two multipliers are not fitted fatigue parameters.

They are fixed by

$$
\boxed{
\mu(t)
=
-\frac{\partial\ln Z}{\partial\alpha},
}
$$

and

$$
\boxed{
\mathcal E(t)
=
-\frac{\partial\ln Z}{\partial\beta}.
}
$$

Thus the closure has the computational structure

$$
\boxed{
\mu(t),\mathcal E(t)
\longrightarrow
\alpha(t),\beta(t)
\longrightarrow
p_\lambda(\lambda,t).
}
$$

## 6. Exact derivative identity inside the exponential family

For this exponential family,

$$
d\mu
=
-\operatorname{Var}(\lambda)d\alpha
-\operatorname{Cov}(\lambda,\psi)d\beta.
$$

At fixed mean $d\mu=0$,

$$
\frac{d\alpha}{d\beta}\bigg|_\mu
=
-
\frac{\operatorname{Cov}(\lambda,\psi)}
{\operatorname{Var}(\lambda)}.
$$

Therefore

$$
\boxed{
\frac{d\mathcal E}{d\beta}\bigg|_\mu
=
-
\left[
\operatorname{Var}(\psi)
-
\frac{
\operatorname{Cov}(\lambda,\psi)^2
}{
\operatorname{Var}(\lambda)
}
\right]
\le0.
}
$$

The sign follows from Cauchy-Schwarz.

Hence, within this closure,

$$
\boxed{
\mu=\text{constant},
\quad
\mathcal E\uparrow
\quad\Longrightarrow\quad
\beta\downarrow.
}
$$

This statement is exact **inside the adopted exponential-family closure**.

## 7. Normal-opening tail

Define

$$
\boxed{
Q_c(t)
=
\int_{\lambda_c}^{\infty}
p_\lambda(\lambda,t)\,d\lambda.
}
$$

Substitution of the closure gives

$$
\boxed{
Q_c(t)
=
\frac{
\displaystyle
\int_{\lambda_c}^{\infty}
\exp[-\alpha(t)\lambda-\beta(t)\psi(\lambda)]\,d\lambda
}{
\displaystyle
\int_0^{\infty}
\exp[-\alpha(t)\lambda-\beta(t)\psi(\lambda)]\,d\lambda
}.
}
$$

There is not yet a general proof that $Q_c$ must be monotone in $\mathcal E$ for every admissible potential and mean because $\alpha$ changes together with $\beta$.

For the current calibrated generalized-LJ shape at $\mu=1$, however, the numerical closure gives monotone tail growth over the tested energy range.

## 8. Reference numerical result

The solver uses Gauss-Legendre quadrature after the infinite-domain map

$$
\lambda=\frac{x}{1-x},
\qquad
0<x<1,
$$

so no finite tensile or compression cutoff is inserted.

At $\mu=1$:

| $\mathcal E$ | $\alpha$ | $\beta$ | $\operatorname{Var}(\lambda)$ | $Q_c$ |
|---:|---:|---:|---:|---:|
| $2\times10^{-4}$ | 10.5403 | 2479.69 | $4.0917\times10^{-4}$ | $4.6411\times10^{-5}$ |
| $5\times10^{-4}$ | 10.4455 | 978.517 | $1.0604\times10^{-3}$ | $5.1948\times10^{-3}$ |
| $10^{-3}$ | 10.2323 | 476.217 | $2.2480\times10^{-3}$ | $2.7631\times10^{-2}$ |
| $2\times10^{-3}$ | 9.63809 | 223.561 | $4.7677\times10^{-3}$ | $6.8314\times10^{-2}$ |
| $4\times10^{-3}$ | 8.53618 | 100.057 | $9.3832\times10^{-3}$ | $1.1464\times10^{-1}$ |

Thus, for the sampled closure states,

$$
\boxed{
\mathcal E\uparrow
\quad\Rightarrow\quad
\operatorname{Var}(\lambda)\uparrow,
\quad
Q_c\uparrow.
}
$$

The second implication is a **NUMERICAL RESULT FOR THE TESTED RANGE**, not a general theorem.

## 9. Quadrature convergence

At

$$
\mu=1,
\qquad
\mathcal E=10^{-3},
$$

the computed tail is

$$
Q_c=0.0276314374023
$$

at quadrature order 320 and

$$
Q_c=0.0276314375062
$$

at order 640.

The absolute difference is approximately

$$
1.04\times10^{-10}.
$$

The numerical tail result is therefore converged far more tightly than the physical-model uncertainty.

## 10. What has and has not been achieved

### Achieved

1. A specific mathematical distribution form has been derived from stated 1D layer-spacing constraints rather than chosen as Gaussian or Weibull.
2. The two distribution parameters are determined by $\mu(t)$ and $\mathcal E(t)$ rather than fitted to fatigue-life data.
3. Increasing stored configurational energy at fixed mean necessarily lowers $\beta$ inside the closure.
4. For the present layer-LJ shape, the tested states broaden and develop a rapidly growing normal-opening tail.

### Not achieved

1. The equal-measure fixed-$(L,E)$ ensemble has not yet been derived from the exact driven cyclic dynamics.
2. The time law $\mathcal E(t)$ has not yet been derived from external work and energy redistribution.
3. $Q_c(t)>0$ is an instantaneous unstable-tail statement, not yet a cumulative crack-initiation probability.
4. The current closure does not prove a universal 100% crack-initiation energy threshold.

## 11. Next falsification step

The highest-value next calculation is direct.

Run the deterministic 1D layer-LJ chain, measure at each selected time

$$
\mu_{\rm sim}(t),
\qquad
\mathcal E_{\rm sim}(t),
$$

solve the closure using exactly those two measured moments, and compare

$$
p_{\lambda,\rm closure}(\lambda,t)
$$

against the empirical spacing distribution from the same simulation.

If the shapes disagree systematically, the closure is rejected or enlarged. If they agree after a mechanically justified mixing regime develops, the distribution form gains direct microscopic support.

---

# 한국어 번역 — 1D Layer-LJ 이론의 유도된 분포 Closure

## 상태

이번 마일스톤에서는 특정 named probability family를 임의로 선택하지 않고 활성 spacing distribution의 첫 번째 명시적 수학형태를 도입한다.

활성 모델은 계속해서 **1차원, 수직변형 전용, layer 기반**이다. 축약된 microscopic coordinate는 인접한 represented layer 사이의 수직간격

$$
a_i(t)>0
$$

이다.

layer 사이의 유효 수직상호작용은 이미 calibration한 generalized Lennard-Jones energy로 표현한다. 여기서는 shear coordinate, 3D FCC kinematics, cycle-dependent damage variable 또는 fitted probability family를 도입하지 않는다.

핵심 결과는 large-system closure

$$
\boxed{
p_\lambda(\lambda,t)
=
\frac{1}{Z(t)}
\exp\left[-\alpha(t)\lambda-\beta(t)\psi(\lambda)\right]
}
$$

이다. 여기서 $\lambda=a/a_0$이고 $\psi$는 shift한 normalized layer-LJ energy다.

이 결과는 **정확한 nonequilibrium cyclic dynamics가 아니다.** 명확하게 선언한 fixed-length/fixed-configurational-energy ensemble과 large-$M$ saddle-point reduction으로 얻은 **CONTROLLED APPROXIMATION**이다.

## 1. 물리적 spacing density와 normalized density

$P_a(a,t)$를 physical spacing density라고 하자.

$$
\int_0^\infty P_a(a,t)\,da=1.
$$

다음 normalized spacing을 도입한다.

$$
\lambda=\frac{a}{a_0}.
$$

이에 대응하는 dimensionless density는

$$
\boxed{
p_\lambda(\lambda,t)=a_0P_a(a_0\lambda,t)
}
$$

이다.

따라서

$$
\int_0^\infty p_\lambda(\lambda,t)\,d\lambda=1
$$

이다.

평균 normalized layer spacing은

$$
\boxed{
\mu(t)=\int_0^\infty \lambda p_\lambda(\lambda,t)\,d\lambda
}
$$

이다.

## 2. Closure에 사용하는 Layer-LJ energy

현재 normalized generalized-LJ form은

$$
\phi(\lambda)
=
\frac{\lambda^{-m}}{m(m-n)}
-
\frac{\lambda^{-n}}{n(m-n)}
$$

이고

$$
m=12.19,
\qquad
n=6
$$

이다.

shifted energy를

$$
\boxed{
\psi(\lambda)=\phi(\lambda)-\phi(1)
}
$$

로 정의한다.

그러면

$$
\psi(1)=0
$$

이고 평균 shifted configurational energy는

$$
\boxed{
\mathcal E(t)
=
\int_0^\infty
\psi(\lambda)p_\lambda(\lambda,t)\,d\lambda
}
$$

이다.

local tangent-instability stretch는 계속

$$
\boxed{
\lambda_c
=
\left(\frac{m+1}{n+1}\right)^{1/(m-n)}
\approx1.1077715386
}
$$

이다.

## 3. 유한 $M$ density-of-states 구성

물리시간 $t$에서 양의 layer spacing $M$개

$$
\lambda_1,\ldots,\lambda_M
$$

를 생각한다.

순간 total normalized length와 configurational energy constraint를

$$
L(t)=\sum_{i=1}^{M}\lambda_i,
$$

$$
E_c(t)=\sum_{i=1}^{M}\psi(\lambda_i)
$$

로 둔다.

동등하게

$$
L(t)=M\mu(t),
\qquad
E_c(t)=M\mathcal E(t)
$$

이다.

admissible configuration density를

$$
\boxed{
\Omega_M(L,E)
=
\int_{(0,\infty)^M}
\delta\left(\sum_i\lambda_i-L\right)
\delta\left(\sum_i\psi(\lambda_i)-E\right)
\,d^M\lambda
}
$$

로 정의한다.

### Ensemble assumption

**ASSUMPTION:** 순간적인 $L(t)$와 $E_c(t)$가 주어졌을 때 reduced layer-spacing coordinate의 constraint manifold 위 configuration에 동일한 base measure를 부여한다.

이것이 이번 milestone에서 추가되는 statistical closure assumption이다. deterministic cyclic dynamics에서 자동으로 따라온다고 주장하지 않는다.

이 가정 아래 finite $M$의 정확한 one-spacing marginal은

$$
\boxed{
p_M(\lambda\mid L,E)
=
\frac{
\Omega_{M-1}\left(L-\lambda,E-\psi(\lambda)\right)
}{
\Omega_M(L,E)
}
}
$$

이다.

## 4. Large-$M$ saddle-point closure

큰 $M$에서

$$
\ln\Omega_{M-1}
\left(L-\lambda,E-\psi(\lambda)\right)
$$

을 macroscopic state $(L,E)$ 주변에서 전개한다.

1차까지 쓰면

$$
\ln\Omega_{M-1}
\approx
C-\alpha\lambda-\beta\psi(\lambda)
$$

이다. $\alpha$, $\beta$는 large-system log density of states의 conjugate derivative다.

따라서

$$
\boxed{
p_\lambda(\lambda,t)
=
\frac{
\exp[-\alpha(t)\lambda-\beta(t)\psi(\lambda)]
}{
Z(\alpha(t),\beta(t))
}
}
$$

이고

$$
\boxed{
Z(\alpha,\beta)
=
\int_0^\infty
\exp[-\alpha\lambda-\beta\psi(\lambda)]\,d\lambda
}
$$

이다.

현재의 unbounded domain에서는 normalizability를 위해

$$
\alpha>0,
\qquad
\beta>0
$$

가 필요하다.

$\beta$는 이 closure의 energy-conjugate multiplier다. 별도의 equilibrium derivation이 없는 한 $1/(k_BT)$라고 해석하면 안 된다.

## 5. Multiplier는 moment로 결정된다

두 multiplier는 fitted fatigue parameter가 아니다.

$$
\boxed{
\mu(t)
=
-\frac{\partial\ln Z}{\partial\alpha}
}
$$

및

$$
\boxed{
\mathcal E(t)
=
-\frac{\partial\ln Z}{\partial\beta}
}
$$

로 결정된다.

따라서 computational structure는

$$
\boxed{
\mu(t),\mathcal E(t)
\longrightarrow
\alpha(t),\beta(t)
\longrightarrow
p_\lambda(\lambda,t)
}
$$

이다.

## 6. Exponential family 내부의 정확한 derivative identity

이 exponential family에서는

$$
d\mu
=
-\operatorname{Var}(\lambda)d\alpha
-\operatorname{Cov}(\lambda,\psi)d\beta
$$

이다.

평균이 고정되어 $d\mu=0$이면

$$
\frac{d\alpha}{d\beta}\bigg|_\mu
=
-
\frac{\operatorname{Cov}(\lambda,\psi)}
{\operatorname{Var}(\lambda)}
$$

이다.

따라서

$$
\boxed{
\frac{d\mathcal E}{d\beta}\bigg|_\mu
=
-
\left[
\operatorname{Var}(\psi)
-
\frac{
\operatorname{Cov}(\lambda,\psi)^2
}{
\operatorname{Var}(\lambda)
}
\right]
\le0
}
$$

이다.

부호는 Cauchy-Schwarz에서 나온다.

따라서 이 closure 내부에서는

$$
\boxed{
\mu=\text{constant},
\quad
\mathcal E\uparrow
\quad\Longrightarrow\quad
\beta\downarrow
}
$$

이다.

이 명제는 **채택한 exponential-family closure 내부에서 정확하다.**

## 7. Normal-opening tail

다음을 정의한다.

$$
\boxed{
Q_c(t)
=
\int_{\lambda_c}^{\infty}
p_\lambda(\lambda,t)\,d\lambda
}
$$

closure를 대입하면

$$
\boxed{
Q_c(t)
=
\frac{
\displaystyle
\int_{\lambda_c}^{\infty}
\exp[-\alpha(t)\lambda-\beta(t)\psi(\lambda)]\,d\lambda
}{
\displaystyle
\int_0^{\infty}
\exp[-\alpha(t)\lambda-\beta(t)\psi(\lambda)]\,d\lambda
}
}
$$

이다.

$\beta$와 함께 $\alpha$도 변하므로 모든 admissible potential과 mean에서 $Q_c$가 $\mathcal E$에 대해 반드시 단조증가한다는 일반증명은 아직 없다.

하지만 현재 calibrated generalized-LJ shape와 $\mu=1$에서는 시험한 energy 범위에서 numerical closure가 monotone tail growth를 준다.

## 8. 기준 수치결과

solver는 infinite-domain map

$$
\lambda=\frac{x}{1-x},
\qquad
0<x<1
$$

뒤 Gauss-Legendre quadrature를 사용한다. 따라서 finite tensile/compression cutoff를 넣지 않는다.

$\mu=1$에서 다음을 얻었다.

| $\mathcal E$ | $\alpha$ | $\beta$ | $\operatorname{Var}(\lambda)$ | $Q_c$ |
|---:|---:|---:|---:|---:|
| $2\times10^{-4}$ | 10.5403 | 2479.69 | $4.0917\times10^{-4}$ | $4.6411\times10^{-5}$ |
| $5\times10^{-4}$ | 10.4455 | 978.517 | $1.0604\times10^{-3}$ | $5.1948\times10^{-3}$ |
| $10^{-3}$ | 10.2323 | 476.217 | $2.2480\times10^{-3}$ | $2.7631\times10^{-2}$ |
| $2\times10^{-3}$ | 9.63809 | 223.561 | $4.7677\times10^{-3}$ | $6.8314\times10^{-2}$ |
| $4\times10^{-3}$ | 8.53618 | 100.057 | $9.3832\times10^{-3}$ | $1.1464\times10^{-1}$ |

즉 sampled closure state에서는

$$
\boxed{
\mathcal E\uparrow
\quad\Rightarrow\quad
\operatorname{Var}(\lambda)\uparrow,
\quad
Q_c\uparrow
}
$$

가 관찰된다.

두 번째 implication은 **시험한 범위에 대한 NUMERICAL RESULT**이며 일반정리가 아니다.

## 9. Quadrature convergence

$$
\mu=1,
\qquad
\mathcal E=10^{-3}
$$

에서 quadrature order 320은

$$
Q_c=0.0276314374023
$$

을 주고 order 640은

$$
Q_c=0.0276314375062
$$

를 준다.

절대차이는 약

$$
1.04\times10^{-10}
$$

이다.

따라서 numerical tail result는 physical-model uncertainty보다 훨씬 작은 수준까지 수렴했다.

## 10. 달성한 것과 아직 달성하지 못한 것

### 달성

1. Gaussian이나 Weibull을 고른 것이 아니라 명시한 1D layer-spacing constraint에서 특정 분포형태를 유도했다.
2. 두 distribution parameter는 fatigue-life data fitting이 아니라 $\mu(t)$와 $\mathcal E(t)$로 결정된다.
3. 고정평균에서 stored configurational energy가 증가하면 closure 내부에서 $\beta$가 반드시 감소한다.
4. 현재 layer-LJ shape에서 시험한 상태는 broadening과 빠른 normal-opening tail growth를 보인다.

### 아직 미달성

1. equal-measure fixed-$(L,E)$ ensemble을 정확한 driven cyclic dynamics로부터 아직 유도하지 않았다.
2. 외부 work와 내부 energy redistribution에서 $\mathcal E(t)$의 시간법칙을 아직 유도하지 않았다.
3. $Q_c(t)>0$는 instantaneous unstable-tail 명제이며 cumulative crack-initiation probability가 아니다.
4. 현재 closure는 universal 100% crack-initiation energy threshold를 증명하지 않는다.

## 11. 다음 반증단계

다음 계산은 직접적으로 할 수 있다.

deterministic 1D layer-LJ chain을 돌리고 선택한 각 시간에서

$$
\mu_{\rm sim}(t),
\qquad
\mathcal E_{\rm sim}(t)
$$

를 측정한다.

그 두 measured moment를 그대로 closure에 넣어

$$
p_{\lambda,\rm closure}(\lambda,t)
$$

를 구하고 같은 simulation의 empirical spacing distribution과 비교한다.

shape가 체계적으로 다르면 closure를 기각하거나 확장한다. mechanics로 정당화할 수 있는 mixing regime 이후에 서로 맞는다면 이 분포형태는 직접적인 microscopic support를 얻게 된다.
