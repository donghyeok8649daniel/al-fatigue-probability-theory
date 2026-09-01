# Milestone 2 — Continuous-Time 1D Normal-Deformation State

## Status

The active theory uses physical time $t$ as the evolution coordinate.

The working chain is

$$
\boxed{
\sigma_n(t)
\rightarrow
\{a_i(t)\}
\rightarrow
P(a,t)
\rightarrow
\mu(t),\mathcal E(t)
\rightarrow
\text{normal-opening feasibility}.
}
$$

Fatigue cycle count is not a state variable. If a constant-frequency experiment later requires a cycle label, it may be derived from $N=ft$.

## 1. One-dimensional normal loading only

The imposed loading is a normal stress history

$$
\sigma_n(t).
$$

For a sinusoidal experiment one may use

$$
\sigma_n(t)=\sigma_m+\sigma_a\sin(\omega t),
$$

but the theory itself is formulated for a general time history.

No shear coordinate or three-dimensional crystal kinematics are part of the active mainline.

## 2. Microscopic normal spacing state

The microscopic coordinates are the 1D normal spacings

$$
a_i(t)=x_{i+1}(t)-x_i(t).
$$

For a finite set of $M$ spacings,

$$
P_M(a,t)
=
\frac1M\sum_{i=1}^{M}\delta(a-a_i(t)).
$$

Here $M$ is a finite system/sample count, not fatigue cycle count.

The continuum state is

$$
\boxed{P(a,t).}
$$

No Gaussian, Weibull, or other named family is assumed.

## 3. Exact kinematics

For deterministic spacing trajectories,

$$
\boxed{
\partial_tP+\partial_a(Pv)=0,
}
$$

where

$$
v(a,t)=\langle\dot a_i\mid a_i=a\rangle.
$$

This is an exact kinematic identity. It does not by itself close the dynamics.

## 4. Fixed generalized-LJ normal energy

The active 1D microscopic energy is derived from the fixed generalized Lennard-Jones interaction.

Using

$$
\lambda_i=\frac{a_i}{a_0},
$$

the reduced pair energy is

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

The parameters do not evolve with time or loading history.

The normalization gives

$$
\phi'(1)=0,
\qquad
\phi''(1)=1.
$$

The idealized tensile tangent-stability limit is

$$
\boxed{
\phi''(\lambda_c)=0,
}
$$

with

$$
\lambda_c\approx1.1077715386.
$$

## 5. Continuous-time moments used by the new theory

Define

$$
\boxed{
\mu(t)=\int\lambda P(\lambda,t)\,d\lambda
}
$$

and the shifted LJ energy

$$
\psi(\lambda)=\phi(\lambda)-\phi(1).
$$

Then

$$
\boxed{
\mathcal E(t)
=
\int\psi(\lambda)P(\lambda,t)\,d\lambda.
}
$$

The current working idea is that $\mu(t)$ may remain approximately conserved while $\mathcal E(t)$ increases through redistribution and broadening of $P$.

This does not require specifying the full evolution law for $P$ first.

## 6. Existing null test

The existing 32-atom perfect-chain reference calculation at 100 MPa normal stress amplitude remains essentially reversible and does not cross $\lambda_c$.

This is a required null/falsification result. Cyclic loading alone must not create artificial fatigue.

## 7. Stronger energy-feasibility route

The active next step is no longer to demand a discrete condition such as $P_{N+1}\neq P_N$.

Instead, at every physical time $t$, ask whether a crack-free distribution satisfying the known normalization, mean, energy, and support constraints can exist.

The exact result is developed in `docs/MILESTONE4_TIME_ENERGY_FEASIBILITY.md`.

The key additional condition is a physically justified compression bound

$$
\boxed{
\lambda\ge\lambda_L(t)>0.
}
$$

Without such a condition, LJ reverse compression can absorb arbitrarily large energy and energy alone cannot force a tensile tail.

## 8. Crack-initiation variable

The instantaneous unstable tail is

$$
Q_c(t)
=
\int_{\lambda_c}^{\infty}P(\lambda,t)\,d\lambda.
$$

The active energy-feasibility first-passage variable is

$$
\boxed{
\tau_E
=
\inf\left\{
t:\mathcal E(t)>\mathcal E_{\rm safe}^{\max}(t)
\right\}.
}
$$

If the lower compression bound is physically valid, crossing this ceiling makes a fully crack-free support mathematically impossible.

---

# 한국어 번역 — 연속시간 1D 수직변형 상태

## 상태

활성 이론은 물리적 시간 $t$를 evolution coordinate로 사용한다.

작동 연결은

$$
\boxed{
\sigma_n(t)
\rightarrow
\{a_i(t)\}
\rightarrow
P(a,t)
\rightarrow
\mu(t),\mathcal E(t)
\rightarrow
\text{normal-opening feasibility}
}
$$

이다.

fatigue cycle count는 상태변수가 아니다. 일정 주파수 실험에서 cycle label이 필요하면 나중에 $N=ft$로 환산한다.

## 1. 1차원 수직하중만 사용

외부하중은 normal stress history

$$
\sigma_n(t)
$$

이다.

sinusoidal experiment에서는

$$
\sigma_n(t)=\sigma_m+\sigma_a\sin(\omega t)
$$

를 사용할 수 있지만 이론 자체는 일반적인 time history에 대해 작성한다.

shear coordinate와 3차원 crystal kinematics는 활성 mainline에 포함하지 않는다.

## 2. 미시 수직-spacing 상태

microscopic coordinate는 1D normal spacing

$$
a_i(t)=x_{i+1}(t)-x_i(t)
$$

이다.

유한한 $M$개의 spacing에 대해

$$
P_M(a,t)
=
\frac1M\sum_{i=1}^{M}\delta(a-a_i(t))
$$

를 정의한다.

여기서 $M$은 fatigue cycle count가 아니라 finite system/sample count다.

continuum state는

$$
\boxed{P(a,t)}
$$

이다.

Gaussian, Weibull 또는 다른 named family를 가정하지 않는다.

## 3. 정확한 kinematics

deterministic spacing trajectory에서는

$$
\boxed{
\partial_tP+\partial_a(Pv)=0
}
$$

이고

$$
v(a,t)=\langle\dot a_i\mid a_i=a\rangle
$$

이다.

이것은 정확한 kinematic identity지만 dynamics를 자동으로 close하지는 않는다.

## 4. 고정 generalized-LJ 수직에너지

활성 1D microscopic energy는 고정 generalized Lennard-Jones interaction에서 유도한다.

$$
\lambda_i=\frac{a_i}{a_0}
$$

를 사용하면 reduced pair energy는

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

parameter는 시간이나 loading history에 따라 변하지 않는다.

normalization으로

$$
\phi'(1)=0,
\qquad
\phi''(1)=1
$$

이고, 이상화된 tensile tangent-stability limit는

$$
\boxed{
\phi''(\lambda_c)=0
}
$$

이며

$$
\lambda_c\approx1.1077715386
$$

이다.

## 5. 새 이론에서 사용하는 연속시간 moment

$$
\boxed{
\mu(t)=\int\lambda P(\lambda,t)\,d\lambda
}
$$

를 정의하고 shifted LJ energy를

$$
\psi(\lambda)=\phi(\lambda)-\phi(1)
$$

로 둔다.

그러면

$$
\boxed{
\mathcal E(t)
=
\int\psi(\lambda)P(\lambda,t)\,d\lambda
}
$$

이다.

현재 working idea는 $\mu(t)$가 거의 보존되는 동안 $P$의 redistribution과 broadening으로 $\mathcal E(t)$가 증가할 수 있다는 것이다.

이를 위해 먼저 $P$의 전체 evolution law를 특정할 필요는 없다.

## 6. 기존 null test

기존 32-atom perfect-chain 100 MPa normal stress amplitude reference calculation은 거의 가역적이며 $\lambda_c$를 넘지 않는다.

이는 반드시 유지해야 하는 null/falsification result다. loading이 cyclic이라는 사실만으로 artificial fatigue가 생기면 안 된다.

## 7. 더 강한 energy-feasibility 경로

활성 다음 단계는 더 이상 $P_{N+1}\neq P_N$ 같은 discrete condition을 요구하는 것이 아니다.

대신 모든 물리적 시간 $t$에서 알려진 정규화, 평균, 에너지, support 조건을 만족하는 crack-free distribution이 존재할 수 있는지를 묻는다.

정확한 결과는 `docs/MILESTONE4_TIME_ENERGY_FEASIBILITY.md`에 정리한다.

핵심 추가조건은 물리적으로 정당화된 compression bound

$$
\boxed{
\lambda\ge\lambda_L(t)>0
}
$$

이다.

이 조건이 없으면 LJ reverse compression이 임의로 큰 energy를 흡수할 수 있으므로 energy만으로 tensile tail을 강제할 수 없다.

## 8. 균열개시 변수

순간 unstable tail은

$$
Q_c(t)
=
\int_{\lambda_c}^{\infty}P(\lambda,t)\,d\lambda
$$

이다.

활성 energy-feasibility first-passage variable은

$$
\boxed{
\tau_E
=
\inf\left\{
t:\mathcal E(t)>\mathcal E_{\rm safe}^{\max}(t)
\right\}
}
$$

이다.

lower compression bound가 물리적으로 유효하다면 이 ceiling을 넘는 순간 fully crack-free support는 수학적으로 불가능해진다.
