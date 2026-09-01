# Milestone 4 — Continuous-Time Energy Feasibility of the 1D Normal-LJ Distribution

## Status

This milestone replaces cycle-index evolution as the primary formulation. The active state is a continuous-time normal-spacing density

$$
P(a,t).
$$

Cycle count may be recovered later from the loading history, but it is not an independent state variable.

The active 1D chain uses the normalized spacing

$$
\lambda=\frac{a}{a_0}
$$

and the fixed generalized Lennard-Jones energy

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

For convenience define the shifted energy

$$
\boxed{
\psi(\lambda)=\phi(\lambda)-\phi(1)
}
$$

so that $\psi(1)=0$.

## 1. Continuous-time state constraints

At every time $t$, define

$$
\boxed{
\int_0^\infty P(\lambda,t)\,d\lambda=1
}
$$

and

$$
\boxed{
\mu(t)=\int_0^\infty \lambda P(\lambda,t)\,d\lambda.
}
$$

The mean LJ configurational energy per represented spacing is

$$
\boxed{
\mathcal E(t)
=
\int_0^\infty \psi(\lambda)P(\lambda,t)\,d\lambda.
}
$$

The research hypothesis is that $\mu(t)$ may remain nearly conserved while $\mathcal E(t)$ increases through redistribution of probability mass.

No named probability family is assumed.

## 2. Exact energy-broadening identity

Inside a convex interval, define the Bregman-type remainder

$$
D_\psi(\lambda\mid\mu)
=
\psi(\lambda)
-
\psi(\mu)
-
\psi'(\mu)(\lambda-\mu).
$$

Because

$$
\int(\lambda-\mu)P(\lambda,t)\,d\lambda=0,
$$

one obtains the exact identity

$$
\boxed{
\mathcal E(t)-\psi(\mu(t))
=
\int D_\psi(\lambda\mid\mu(t))P(\lambda,t)\,d\lambda.
}
$$

If $\psi''>0$ on the occupied interval, then

$$
D_\psi\ge0.
$$

Therefore, at approximately fixed mean spacing, configurational energy above $\psi(\mu)$ is mathematically associated with distributional spreading away from the mean.

**Classification: EXACT / IDENTITY under the stated definitions.**

## 3. Why normalization + mean + energy are not enough

The generalized LJ repulsion diverges as

$$
\lambda\rightarrow0^+.
$$

Suppose only that a crack-free state satisfies

$$
0<\lambda\le\lambda_c
$$

and has fixed mean $\mu<\lambda_c$.

For any $\epsilon\in(0,\mu)$, construct the two-point probability measure

$$
P_\epsilon
=
w_-\delta(\lambda-\epsilon)
+w_c\delta(\lambda-\lambda_c),
$$

with

$$
w_c=\frac{\mu-\epsilon}{\lambda_c-\epsilon},
\qquad
w_-=1-w_c.
$$

This measure is normalized, has mean $\mu$, and contains no probability beyond $\lambda_c$.

However,

$$
\mathcal E_\epsilon
=
w_-\psi(\epsilon)+w_c\psi(\lambda_c)
\rightarrow\infty
$$

as

$$
\epsilon\rightarrow0^+.
$$

Therefore

$$
\boxed{
\text{normalization + mean + energy alone cannot force a tensile crack tail.}
}
$$

Arbitrarily large energy can mathematically be hidden in sufficiently strong reverse compression.

**Classification: EXACT impossibility result for the generalized LJ potential.**

## 4. The additional condition

The missing condition must limit the reverse-compression route.

The most direct form is a mechanically justified lower support bound

$$
\boxed{
\lambda\ge\lambda_L(t)>0.
}
$$

This is not allowed to be an arbitrary fitting parameter. A later step must derive or measure $\lambda_L(t)$ from 1D normal mechanics, boundary conditions, energy accessibility, or experimental constraints.

For now it is an explicit input to the mathematical feasibility theorem.

## 5. Crack-free admissible set

The 1D LJ tangent-stability limit is

$$
\boxed{
\phi''(\lambda_c)=0,
}
$$

with

$$
\lambda_c\approx1.1077715386.
$$

Define the crack-free admissible set at time $t$ as

$$
\boxed{
\mathcal A_{\rm safe}(t)
=
\left\{
P:\
P\ge0,
\int P=1,
\int\lambda P=\mu(t),
\operatorname{supp}P\subset[\lambda_L(t),\lambda_c]
\right\}.
}
$$

For the generalized LJ potential,

$$
\phi''(\lambda)>0
$$

for

$$
0<\lambda<\lambda_c.
$$

Thus the energy is convex throughout the crack-free interval.

## 6. Exact minimum safe energy

By Jensen's inequality,

$$
\boxed{
\mathcal E_{\rm safe}^{\min}(t)
=
\psi(\mu(t)).
}
$$

The minimum is attained by

$$
P(\lambda,t)=\delta(\lambda-\mu(t)).
$$

## 7. Exact maximum safe energy

A convex graph lies below the secant chord connecting the endpoints of the interval. Therefore, for every

$$
\lambda\in[\lambda_L,\lambda_c],
$$

$$
\psi(\lambda)
\le
\frac{\lambda_c-\lambda}{\lambda_c-\lambda_L}\psi(\lambda_L)
+
\frac{\lambda-\lambda_L}{\lambda_c-\lambda_L}\psi(\lambda_c).
$$

Taking the expectation and using the fixed mean gives

$$
\boxed{
\mathcal E_{\rm safe}^{\max}(t)
=
\frac{\lambda_c-\mu(t)}{\lambda_c-\lambda_L(t)}
\psi(\lambda_L(t))
+
\frac{\mu(t)-\lambda_L(t)}{\lambda_c-\lambda_L(t)}
\psi(\lambda_c).
}
$$

This maximum is attained by the endpoint measure

$$
\boxed{
P^*(\lambda,t)
=
w_L(t)\delta(\lambda-\lambda_L(t))
+w_c(t)\delta(\lambda-\lambda_c),
}
$$

where

$$
w_c(t)
=
\frac{\mu(t)-\lambda_L(t)}{\lambda_c-\lambda_L(t)},
$$

and

$$
w_L(t)=1-w_c(t).
$$

Therefore the maximum is not merely a loose inequality; it is the exact extremum over all probability measures satisfying the stated constraints.

**Classification: EXACT theorem under the lower-support and crack-free-support conditions.**

## 8. Exact feasibility interval

A crack-free probability measure with the stated normalization, mean, and support exists if and only if

$$
\boxed{
\psi(\mu(t))
\le
\mathcal E(t)
\le
\mathcal E_{\rm safe}^{\max}(t).
}
$$

Every intermediate energy is attainable by mixing the minimum-energy delta measure with the maximum-energy endpoint measure while preserving the same mean.

Define the energy margin

$$
\boxed{
M_E(t)
=
\mathcal E_{\rm safe}^{\max}(t)-\mathcal E(t).
}
$$

Then

- $M_E(t)>0$: a crack-free distribution is still mathematically feasible;
- $M_E(t)=0$: the state is on the energy-feasibility boundary;
- $M_E(t)<0$: no probability distribution confined to $[\lambda_L(t),\lambda_c]$ can satisfy the measured mean and energy.

## 9. Continuous-time first passage

Define

$$
\boxed{
\tau_E
=
\inf\left\{
t\ge0:
\mathcal E(t)>
\mathcal E_{\rm safe}^{\max}(t)
\right\}.
}
$$

At $t\ge\tau_E$, the assumptions

$$
\lambda\ge\lambda_L(t)
$$

and

$$
\lambda\le\lambda_c
$$

cannot both remain true for the entire distribution.

If $\lambda_L(t)$ is a valid hard lower bound, the only remaining route is

$$
\boxed{
\int_{\lambda_c}^{\infty}P(\lambda,t)\,d\lambda>0.
}
$$

Thus the unstable tensile tail becomes mathematically unavoidable.

This is a first-passage time in physical time $t$, not in cycle count.

For a constant loading frequency one may later report

$$
N_E=f\tau_E
$$

only as a derived experimental label, not as the fundamental evolution coordinate.

## 10. Numerical reference

The repository code evaluates the exact bound for several illustrative values of $\lambda_L$. These values are parameter-sweep examples only and are not Al material constants.

At

$$
\mu=1,
$$

the dimensionless safe-energy ceilings are approximately

| illustrative $\lambda_L$ | $\mathcal E_{\rm safe}^{\max}$ |
| ---: | ---: |
| 0.90 | 0.00704373 |
| 0.95 | 0.00215815 |
| 0.98 | 0.000650058 |
| 0.99 | 0.000296088 |

The strong dependence on $\lambda_L$ confirms that the compression constraint is not a detail. It is the essential additional physical condition required to convert stored energy into a forced tensile-tail statement.

## 11. Next problem

The next research task is now sharply defined:

$$
\boxed{
\text{derive }\lambda_L(t)\text{ from 1D LJ normal mechanics.}
}
$$

Possible derivation routes must remain one-dimensional and normal-only. Candidates include finite-chain energy accessibility, imposed-force bounds, boundary conditions, and direct experimental constraints on maximum reverse compression.

No 3D FCC model, shear variable, fitted damping law, or empirical damage law is required for this step.

---

# 한국어 번역 — 1D Normal-LJ 분포의 연속시간 에너지 실현가능성

## 상태

이 마일스톤에서는 cycle index를 주된 evolution 변수로 사용하지 않는다. 활성 상태는 연속시간 수직 원자간격 밀도

$$
P(a,t)
$$

이다.

cycle 수는 필요하면 나중에 loading history로부터 환산할 수 있지만 독립적인 상태변수가 아니다.

활성 1D chain은 normalized spacing

$$
\lambda=\frac{a}{a_0}
$$

와 고정 generalized Lennard-Jones energy

$$
\phi(\lambda)
=
\frac{\lambda^{-m}}{m(m-n)}
-
\frac{\lambda^{-n}}{n(m-n)}
$$

를 사용하며

$$
m=12.19,
\qquad
n=6
$$

이다.

편의를 위해 equilibrium energy를 뺀

$$
\boxed{
\psi(\lambda)=\phi(\lambda)-\phi(1)
}
$$

을 정의하여 $\psi(1)=0$으로 둔다.

## 1. 연속시간 상태 제약

모든 시각 $t$에서

$$
\boxed{
\int_0^\infty P(\lambda,t)\,d\lambda=1
}
$$

이고

$$
\boxed{
\mu(t)=\int_0^\infty \lambda P(\lambda,t)\,d\lambda
}
$$

를 정의한다.

대표 spacing 하나당 평균 LJ configurational energy는

$$
\boxed{
\mathcal E(t)
=
\int_0^\infty \psi(\lambda)P(\lambda,t)\,d\lambda
}
$$

이다.

연구 가설은 $\mu(t)$가 거의 보존되는 동안 확률질량의 재분배를 통해 $\mathcal E(t)$가 증가할 수 있다는 것이다.

특정 named probability family는 가정하지 않는다.

## 2. 정확한 에너지-분포확산 identity

convex interval 안에서

$$
D_\psi(\lambda\mid\mu)
=
\psi(\lambda)
-
\psi(\mu)
-
\psi'(\mu)(\lambda-\mu)
$$

를 정의한다.

평균 정의로부터

$$
\int(\lambda-\mu)P(\lambda,t)\,d\lambda=0
$$

이므로 정확히

$$
\boxed{
\mathcal E(t)-\psi(\mu(t))
=
\int D_\psi(\lambda\mid\mu(t))P(\lambda,t)\,d\lambda
}
$$

을 얻는다.

점유구간에서 $\psi''>0$이면

$$
D_\psi\ge0
$$

이다.

따라서 평균 spacing이 거의 고정되어 있을 때 $\psi(\mu)$를 초과하는 configurational energy는 수학적으로 평균으로부터 분포가 퍼지는 것과 연결된다.

**분류: 주어진 정의 아래 EXACT / IDENTITY.**

## 3. 정규화 + 평균 + 에너지만으로 부족한 이유

generalized LJ repulsion은

$$
\lambda\rightarrow0^+
$$

에서 발산한다.

crack-free 상태에 대해 단지

$$
0<\lambda\le\lambda_c
$$

이고 평균이 $\mu<\lambda_c$라고만 하자.

임의의 $\epsilon\in(0,\mu)$에 대해

$$
P_\epsilon
=
w_-\delta(\lambda-\epsilon)
+w_c\delta(\lambda-\lambda_c)
$$

를 만들고

$$
w_c=\frac{\mu-\epsilon}{\lambda_c-\epsilon},
\qquad
w_-=1-w_c
$$

로 두면 이 분포는 정규화되고 평균이 $\mu$이며 $\lambda_c$를 넘는 확률질량이 없다.

그런데

$$
\mathcal E_\epsilon
=
w_-\psi(\epsilon)+w_c\psi(\lambda_c)
\rightarrow\infty
$$

가

$$
\epsilon\rightarrow0^+
$$

에서 성립한다.

따라서

$$
\boxed{
\text{정규화 + 평균 + 에너지만으로 tensile crack tail을 강제할 수 없다.}
}
$$

아주 강한 역압축에 임의로 큰 에너지를 수학적으로 저장할 수 있기 때문이다.

**분류: generalized LJ potential에 대한 EXACT impossibility result.**

## 4. 추가로 필요한 조건

빠져 있는 조건은 reverse-compression 경로를 제한해야 한다.

가장 직접적인 형태는 역학적으로 정당화된 lower support bound

$$
\boxed{
\lambda\ge\lambda_L(t)>0
}
$$

이다.

이 값은 임의 fitting parameter가 될 수 없다. 다음 단계에서 1D normal mechanics, boundary condition, energy accessibility 또는 experimental constraint로부터 $\lambda_L(t)$를 유도하거나 측정해야 한다.

현재 단계에서는 수학적 feasibility theorem의 명시적 input으로 둔다.

## 5. crack-free admissible set

1D LJ tangent-stability limit는

$$
\boxed{
\phi''(\lambda_c)=0
}
$$

이고

$$
\lambda_c\approx1.1077715386
$$

이다.

시각 $t$에서 crack-free admissible set을

$$
\boxed{
\mathcal A_{\rm safe}(t)
=
\left\{
P:\
P\ge0,
\int P=1,
\int\lambda P=\mu(t),
\operatorname{supp}P\subset[\lambda_L(t),\lambda_c]
\right\}
}
$$

으로 정의한다.

generalized LJ potential은

$$
0<\lambda<\lambda_c
$$

에서

$$
\phi''(\lambda)>0
$$

이므로 crack-free interval 전체에서 convex하다.

## 6. 정확한 최소 safe energy

Jensen inequality로부터

$$
\boxed{
\mathcal E_{\rm safe}^{\min}(t)
=
\psi(\mu(t))
}
$$

이다.

최솟값은

$$
P(\lambda,t)=\delta(\lambda-\mu(t))
$$

에서 달성된다.

## 7. 정확한 최대 safe energy

convex function의 graph는 interval 양 끝점을 잇는 secant chord 아래에 있으므로

$$
\lambda\in[\lambda_L,\lambda_c]
$$

에서

$$
\psi(\lambda)
\le
\frac{\lambda_c-\lambda}{\lambda_c-\lambda_L}\psi(\lambda_L)
+
\frac{\lambda-\lambda_L}{\lambda_c-\lambda_L}\psi(\lambda_c)
$$

이다.

기댓값을 취하고 평균조건을 사용하면

$$
\boxed{
\mathcal E_{\rm safe}^{\max}(t)
=
\frac{\lambda_c-\mu(t)}{\lambda_c-\lambda_L(t)}
\psi(\lambda_L(t))
+
\frac{\mu(t)-\lambda_L(t)}{\lambda_c-\lambda_L(t)}
\psi(\lambda_c)
}
$$

를 얻는다.

이 최댓값은 endpoint measure

$$
\boxed{
P^*(\lambda,t)
=
w_L(t)\delta(\lambda-\lambda_L(t))
+w_c(t)\delta(\lambda-\lambda_c)
}
$$

에서 실제로 달성되며

$$
w_c(t)
=
\frac{\mu(t)-\lambda_L(t)}{\lambda_c-\lambda_L(t)}
$$

이고

$$
w_L(t)=1-w_c(t)
$$

이다.

따라서 이 식은 단순한 loose bound가 아니라 명시된 제약을 만족하는 모든 probability measure에 대한 정확한 extremum이다.

**분류: lower-support 및 crack-free-support 조건 아래 EXACT theorem.**

## 8. 정확한 feasibility interval

정규화, 평균 및 support 조건을 만족하는 crack-free probability measure가 존재할 필요충분조건은

$$
\boxed{
\psi(\mu(t))
\le
\mathcal E(t)
\le
\mathcal E_{\rm safe}^{\max}(t)
}
$$

이다.

같은 평균을 유지한 채 최소-energy delta measure와 최대-energy endpoint measure를 혼합하면 중간의 모든 energy를 만들 수 있다.

energy margin을

$$
\boxed{
M_E(t)
=
\mathcal E_{\rm safe}^{\max}(t)-\mathcal E(t)
}
$$

으로 정의한다.

- $M_E(t)>0$: crack-free distribution이 아직 수학적으로 가능하다.
- $M_E(t)=0$: energy-feasibility boundary에 있다.
- $M_E(t)<0$: $[\lambda_L(t),\lambda_c]$ 안에 갇힌 어떤 probability distribution도 해당 평균과 에너지를 동시에 만족할 수 없다.

## 9. 연속시간 first passage

$$
\boxed{
\tau_E
=
\inf\left\{
t\ge0:
\mathcal E(t)>
\mathcal E_{\rm safe}^{\max}(t)
\right\}
}
$$

로 정의한다.

$t\ge\tau_E$에서는

$$
\lambda\ge\lambda_L(t)
$$

와

$$
\lambda\le\lambda_c
$$

를 분포 전체가 동시에 만족할 수 없다.

$\lambda_L(t)$가 실제 hard lower bound라면 남은 유일한 경로는

$$
\boxed{
\int_{\lambda_c}^{\infty}P(\lambda,t)\,d\lambda>0
}
$$

이므로 unstable tensile tail이 수학적으로 필수가 된다.

이것은 cycle count가 아니라 물리적 시간 $t$에 대한 first-passage time이다.

주파수가 일정한 실험에서만 필요하면 나중에

$$
N_E=f\tau_E
$$

를 실험 표시값으로 환산할 수 있지만 이것을 근본 evolution coordinate로 사용하지 않는다.

## 10. 수치 reference

repository code는 여러 illustrative $\lambda_L$에 대해 정확한 bound를 계산한다. 이 값들은 parameter-sweep 예시일 뿐 Al material constant가 아니다.

$$
\mu=1
$$

에서 dimensionless safe-energy ceiling은 대략 다음과 같다.

| illustrative $\lambda_L$ | $\mathcal E_{\rm safe}^{\max}$ |
| ---: | ---: |
| 0.90 | 0.00704373 |
| 0.95 | 0.00215815 |
| 0.98 | 0.000650058 |
| 0.99 | 0.000296088 |

$\lambda_L$에 대한 강한 의존성은 compression constraint가 사소한 detail이 아니라는 것을 보여준다. 저장에너지로부터 tensile tail의 필연성을 말하기 위해 필요한 핵심 추가 물리조건이다.

## 11. 다음 문제

다음 연구과제는 이제 명확하다.

$$
\boxed{
\text{1D LJ normal mechanics로부터 }\lambda_L(t)\text{를 유도한다.}
}
$$

가능한 유도경로는 1차원 및 normal-only로 유지해야 한다. finite-chain energy accessibility, imposed-force bound, boundary condition, maximum reverse compression에 대한 직접 실험제약 등이 후보가 될 수 있다.

이 단계에는 3D FCC model, shear variable, fitted damping law 또는 empirical damage law가 필요하지 않다.
