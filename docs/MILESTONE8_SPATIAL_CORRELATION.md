# Milestone 8 — Spatial Correlation in the 1D Layer-LJ State

## Scope

The active theory remains strictly one-dimensional and normal-only. The local state variable is the normalized layer spacing

$$
\lambda_i(t)=\frac{a_i(t)}{a_0}.
$$

The current one-point closure uses only the unordered spacing density $p_\lambda(\lambda,t)$. This milestone tests whether deterministic 1D layer-LJ mechanics retains spatial ordering information that cannot be represented by that one-point density.

## Exact finite-chain definitions

For $M$ represented spacings,

$$
\mu(t)=\frac1M\sum_{i=1}^M\lambda_i(t),
$$

$$
\boxed{
C_k(t)=\frac{1}{M-k}\sum_{i=1}^{M-k}[\lambda_i(t)-\mu(t)][\lambda_{i+k}(t)-\mu(t)]
}
$$

and

$$
\boxed{\rho_k(t)=\frac{C_k(t)}{C_0(t)}}.
$$

Since

$$
C_0(t)=\frac1M\sum_i[\lambda_i(t)-\mu(t)]^2,
$$

$C_0$ is exactly the empirical spacing variance.

## Exact non-identifiability of ordering from the one-point density

The empirical one-point density

$$
p_M(\lambda,t)=\frac1M\sum_{i=1}^M\delta[\lambda-\lambda_i(t)]
$$

is invariant under any permutation $\pi$ of the layer-spacing labels:

$$
\boxed{p_M^\pi(\lambda,t)=p_M(\lambda,t)}.
$$

All one-point moments are therefore unchanged by permuting exactly the same spacing values, while $C_k$ generally changes. Hence

$$
\boxed{p_M(\lambda,t)\text{ does not determine spatial ordering}.}
$$

This is an **EXACT / IDENTITY-level structural result** and does not depend on the proposed exponential closure.

## Exact random-permutation reference

Let $d_i=\lambda_i-\mu$, so that $\sum_i d_i=0$. Under a uniformly random permutation, two distinct elements satisfy

$$
\mathbb E[d_{\pi(i)}d_{\pi(j)}]=-rac{C_0}{M-1},\qquad i\ne j.
$$

Therefore every nonzero lag obeys

$$
\boxed{\mathbb E_{\rm perm}[\rho_k]=-rac{1}{M-1}}.
$$

This is exact for the stated finite permutation ensemble.

## Dynamically matched protocol

The existing controlled system-size protocol is retained:

$$
F_a^*=0.03,
$$

$$
\boxed{\omega M=0.62},
$$

and the phase-locked snapshot is taken at $t_s=2T$. The tested represented spacing counts are

$$
M=31,\ 63,\ 127,\ 255.
$$

This is a **CONTROLLED NUMERICAL PROTOCOL**, not a material law.

## Numerical result

Nearest-neighbor correlation is

| $M$ | $\rho_1$ | random-permutation expectation |
|---:|---:|---:|
| 31 | 0.933439 | -0.033333 |
| 63 | 0.966289 | -0.016129 |
| 127 | 0.982820 | -0.007937 |
| 255 | 0.991302 | -0.003937 |

A fixed lag becomes a vanishing fraction of the chain as $M$ grows, so the more meaningful comparison uses

$$
\eta=\frac{k}{M}.
$$

At similar $\eta$, the four correlation profiles nearly collapse:

| $\eta$ | $M=31$ | $M=63$ | $M=127$ | $M=255$ |
|---:|---:|---:|---:|---:|
| 0.05 | 0.8618 | 0.8955 | 0.8929 | 0.8819 |
| 0.10 | 0.7852 | 0.7808 | 0.7560 | 0.7530 |
| 0.20 | 0.5258 | 0.4740 | 0.4916 | 0.4757 |
| 0.30 | 0.2250 | 0.1707 | 0.1668 | 0.1638 |
| 0.40 | -0.1113 | -0.1623 | -0.1884 | -0.1877 |

The first zero crossing satisfies

| $M$ | $k_0/M$ |
|---:|---:|
| 31 | 0.3560 |
| 63 | 0.3513 |
| 127 | 0.3481 |
| 255 | 0.3463 |

and the positive-correlation area divided by $M$ approaches about

$$
\boxed{0.184}.
$$

These are **NUMERICAL RESULTS** of the stated reduced-model protocol, not a theorem for aluminum fatigue.

## Interpretation and next state

The driven 1D chain retains an $O(M)$ coherent spatial structure. A one-point density forgets layer order and therefore cannot distinguish the actual deterministic state from a permutation of exactly the same spacing values.

The physically motivated next statistical object is therefore a neighboring-spacing joint density

$$
\boxed{P_2(\lambda,\lambda',t)}.
$$

It can carry, for example,

$$
C_1(t)=\iint(\lambda-\mu)(\lambda'-\mu)P_2(\lambda,\lambda',t)\,d\lambda\,d\lambda'
$$

up to the finite open-chain endpoint convention.

The next task is to derive the minimum pair-state evolution required by the 1D layer-LJ equations. No empirical correlation length or fitted relaxation law is introduced.

---

# 한국어 번역 — 1D Layer-LJ 상태의 공간상관

## 범위

활성 이론은 계속 엄격하게 1차원 수직변형만 다룬다. 국부 상태변수는 normalized layer spacing

$$
\lambda_i(t)=\frac{a_i(t)}{a_0}
$$

이다.

현재 one-point closure는 순서를 제거한 spacing density $p_\lambda(\lambda,t)$만 사용한다. 이번 마일스톤에서는 deterministic 1D layer-LJ mechanics가 이 one-point density로 표현할 수 없는 공간적 순서정보를 유지하는지 검사한다.

## 정확한 finite-chain 정의

$M$개의 represented spacing에 대해

$$
\mu(t)=\frac1M\sum_{i=1}^M\lambda_i(t),
$$

$$
\boxed{C_k(t)=\frac{1}{M-k}\sum_{i=1}^{M-k}[\lambda_i(t)-\mu(t)][\lambda_{i+k}(t)-\mu(t)]}
$$

및

$$
\boxed{\rho_k(t)=\frac{C_k(t)}{C_0(t)}}
$$

를 정의한다.

$$
C_0(t)=\frac1M\sum_i[\lambda_i(t)-\mu(t)]^2
$$

이므로 $C_0$는 정확히 empirical spacing variance다.

## one-point density로 ordering을 결정할 수 없다는 정확한 결과

empirical one-point density

$$
p_M(\lambda,t)=\frac1M\sum_{i=1}^M\delta[\lambda-\lambda_i(t)]
$$

는 layer-spacing label의 임의 permutation $\pi$에 대해

$$
\boxed{p_M^\pi(\lambda,t)=p_M(\lambda,t)}
$$

로 불변이다.

따라서 정확히 같은 spacing 값들의 순서를 바꾸면 모든 one-point moment는 그대로지만 $C_k$는 일반적으로 바뀐다. 즉

$$
\boxed{p_M(\lambda,t)\text{만으로 spatial ordering을 결정할 수 없다}.}
$$

이는 proposed exponential closure와 무관한 **EXACT / IDENTITY 수준의 구조적 결과**다.

## 정확한 random-permutation 기준

$d_i=\lambda_i-\mu$로 두면 $\sum_i d_i=0$이다. uniform random permutation에서 서로 다른 두 원소는

$$
\mathbb E[d_{\pi(i)}d_{\pi(j)}]=-rac{C_0}{M-1},\qquad i\ne j
$$

를 만족한다.

따라서 모든 nonzero lag에서

$$
\boxed{\mathbb E_{\rm perm}[\rho_k]=-rac{1}{M-1}}
$$

이다. 이는 stated finite permutation ensemble에 대해 정확하다.

## Dynamically matched protocol

기존 controlled system-size protocol을 그대로 유지한다.

$$
F_a^*=0.03,
$$

$$
\boxed{\omega M=0.62}
$$

이며 phase-locked snapshot은 $t_s=2T$에서 취한다. 시험한 represented spacing 수는

$$
M=31,\ 63,\ 127,\ 255
$$

이다.

이는 material law가 아니라 **CONTROLLED NUMERICAL PROTOCOL**이다.

## 수치결과

nearest-neighbor correlation은 다음과 같다.

| $M$ | $\rho_1$ | random-permutation expectation |
|---:|---:|---:|
| 31 | 0.933439 | -0.033333 |
| 63 | 0.966289 | -0.016129 |
| 127 | 0.982820 | -0.007937 |
| 255 | 0.991302 | -0.003937 |

$M$이 커질수록 fixed lag는 전체 chain에서 차지하는 비율이 0으로 가므로 더 의미 있는 비교는

$$
\eta=\frac{k}{M}
$$

를 사용한다.

비슷한 $\eta$에서 네 correlation profile은 거의 같은 곡선으로 collapse한다.

| $\eta$ | $M=31$ | $M=63$ | $M=127$ | $M=255$ |
|---:|---:|---:|---:|---:|
| 0.05 | 0.8618 | 0.8955 | 0.8929 | 0.8819 |
| 0.10 | 0.7852 | 0.7808 | 0.7560 | 0.7530 |
| 0.20 | 0.5258 | 0.4740 | 0.4916 | 0.4757 |
| 0.30 | 0.2250 | 0.1707 | 0.1668 | 0.1638 |
| 0.40 | -0.1113 | -0.1623 | -0.1884 | -0.1877 |

첫 zero crossing은

| $M$ | $k_0/M$ |
|---:|---:|
| 31 | 0.3560 |
| 63 | 0.3513 |
| 127 | 0.3481 |
| 255 | 0.3463 |

으로 접근하고 positive-correlation area를 $M$으로 나눈 값은 약

$$
\boxed{0.184}
$$

로 접근한다.

이는 stated reduced-model protocol의 **NUMERICAL RESULT**이며 실제 aluminum fatigue에 대한 theorem은 아니다.

## 해석과 다음 상태

Driven 1D chain은 $O(M)$ 규모의 coherent spatial structure를 유지한다. one-point density는 layer 순서를 잊기 때문에 실제 deterministic state와 정확히 같은 spacing 값을 섞은 상태를 구분할 수 없다.

따라서 물리적으로 자연스러운 다음 statistical object는 neighboring-spacing joint density

$$
\boxed{P_2(\lambda,\lambda',t)}
$$

이다.

finite open-chain endpoint convention을 제외하면

$$
C_1(t)=\iint(\lambda-\mu)(\lambda'-\mu)P_2(\lambda,\lambda',t)\,d\lambda\,d\lambda'
$$

로 표현할 수 있다.

다음 목표는 경험적 correlation length나 fitted relaxation law를 넣는 것이 아니라 1D layer-LJ equation에서 필요한 최소 pair-state evolution을 유도하는 것이다.
