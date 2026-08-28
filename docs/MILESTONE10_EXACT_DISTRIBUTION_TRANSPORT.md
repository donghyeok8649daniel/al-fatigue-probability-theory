# Milestone 10 — Exact Nonlinear Transport Clue for the 1D Spacing Distribution

## Status and correction of the previous route

The active model remains strictly one-dimensional, normal-only, and layer based.

This milestone **supersedes the harmonic/Taylor route as an active derivation of the full spacing distribution**. A linear mode, an arcsine push-forward, or a Taylor expansion of the LJ force can remain historical local diagnostics, but they are not used to determine the full $P(\lambda,t)$ because:

1. a finite harmonic ansatz inserts spatial periodicity that is not guaranteed by the driven boundary-value problem;
2. a Taylor expansion about $\lambda=1$ is local and is not controlled over the full support, especially the tensile and compression tails;
3. the target theory must retain the original nonlinear generalized-LJ force over the entire positive-spacing domain.

No Taylor expansion, harmonic ansatz, Gaussian/Weibull family, damping law, or empirical fatigue-damage variable is used below.

## 1. Microscopic governing equation

For interior normalized layer spacings,

$$
\boxed{
\ddot\lambda_i
=
\phi'(\lambda_{i+1})
-2\phi'(\lambda_i)
+\phi'(\lambda_{i-1}).
}
$$

This is the original nonlinear 1D layer-LJ equation. Boundary spacings require the actual boundary force/fixation law and are not silently replaced by the bulk equation.

## 2. Exact finite-$M$ empirical density and flux

For $M$ represented spacings define

$$
\boxed{
P_M(\lambda,t)
=
\frac1M\sum_{i=1}^{M}
\delta[\lambda-\lambda_i(t)].
}
$$

Define the spacing-space probability flux

$$
\boxed{
J_M(\lambda,t)
=
\frac1M\sum_{i=1}^{M}
\dot\lambda_i(t)
\delta[\lambda-\lambda_i(t)].
}
$$

Differentiating the empirical measure in the distributional sense gives the exact identity

$$
\boxed{
\partial_tP_M+\partial_\lambda J_M=0.
}
$$

**Classification: EXACT / IDENTITY.**

This is already a direct clue about the form of $P$: the density is not selected from a named probability family. Its shape is transported, compressed, and expanded in spacing space by the mechanically generated flux $J$.

If a smooth continuum density exists and $P>0$, define the conditional mean spacing velocity

$$
\boxed{
\bar v(\lambda,t)=\frac{J(\lambda,t)}{P(\lambda,t)}
=\mathbb E[\dot\lambda_i\mid\lambda_i=\lambda].
}
$$

Then

$$
\boxed{
\partial_tP+\partial_\lambda(P\bar v)=0.
}
$$

Where $\bar v$ is sufficiently smooth, a characteristic $\lambda(t)$ obeys

$$
\dot\lambda=\bar v(\lambda,t),
$$

and along that characteristic

$$
\boxed{
\frac{d}{dt}\ln P
=-\partial_\lambda\bar v.
}
$$

Thus local convergence of the spacing-space flow increases density, while local divergence decreases density. This is a transport statement, not a fitted distribution law.

## 3. Why $P(\lambda,t)$ alone is not the natural complete state

The microscopic mechanics is second order in time. Introduce spacing velocity

$$
v_i(t)=\dot\lambda_i(t)
$$

and the finite-$M$ empirical phase-space measure

$$
\boxed{
F_M(\lambda,v,t)
=
\frac1M\sum_i
\delta[\lambda-\lambda_i(t)]
\delta[v-v_i(t)].
}
$$

Then

$$
P_M(\lambda,t)=\int F_M(\lambda,v,t)\,dv,
$$

and

$$
J_M(\lambda,t)=\int vF_M(\lambda,v,t)\,dv.
$$

Define the acceleration flux

$$
\boxed{
G_M(\lambda,v,t)
=
\frac1M\sum_i
\ddot\lambda_i(t)
\delta[\lambda-\lambda_i(t)]
\delta[v-v_i(t)].
}
$$

Direct differentiation gives

$$
\boxed{
\partial_tF_M
+\partial_\lambda(vF_M)
+\partial_vG_M
=0.
}
$$

**Classification: EXACT / IDENTITY for the finite empirical measure.**

If a smooth phase-space density exists, one may write

$$
G=F\,\bar a(\lambda,v,t),
$$

where

$$
\bar a
=
\mathbb E[\ddot\lambda_i\mid\lambda_i=\lambda,\dot\lambda_i=v].
$$

The resulting kinetic transport equation is

$$
\boxed{
\partial_tF
+\partial_\lambda(vF)
+\partial_v(\bar aF)
=0.
}
$$

This shows a structural limitation of a one-point spacing density: two states can have the same $P(\lambda,t)$ while having different conditional velocity distributions and therefore different future evolution.

## 4. Exact one-point moment hierarchy

Define

$$
K(\lambda,t)=\int v^2F(\lambda,v,t)\,dv.
$$

Integrating the phase-space transport equation over velocity yields again

$$
\partial_tP+\partial_\lambda J=0.
$$

Multiplying by $v$ and integrating gives

$$
\boxed{
\partial_tJ+\partial_\lambda K=A_1,
}
$$

where

$$
A_1(\lambda,t)
=
\int G(\lambda,v,t)\,dv
=P(\lambda,t)\mathbb E[\ddot\lambda_i\mid\lambda_i=\lambda].
$$

Therefore even the evolution of the probability flux requires both a velocity second moment and a conditional acceleration.

For the finite empirical raw moment

$$
M_r(t)=\frac1M\sum_i\lambda_i^r,
$$

the exact kinematic identities are

$$
\boxed{
\dot M_r
=r\left\langle\lambda^{r-1}v\right\rangle,
}
$$

and

$$
\boxed{
\ddot M_r
=r(r-1)\left\langle\lambda^{r-2}v^2\right\rangle
+r\left\langle\lambda^{r-1}a\right\rangle.
}
$$

No closure assumption is involved.

## 5. The original nonlinear LJ law introduces neighbor joint states exactly

For an interior spacing with value $\lambda$,

$$
\ddot\lambda_i
=
\phi'(\lambda_{i+1})
-2\phi'(\lambda_i)
+\phi'(\lambda_{i-1}).
$$

Taking a conditional expectation at fixed central spacing gives

$$
\boxed{
\mathbb E[\ddot\lambda_i\mid\lambda_i=\lambda]
=
\mathbb E[\phi'(\lambda_{i+1})+\phi'(\lambda_{i-1})\mid\lambda_i=\lambda]
-2\phi'(\lambda).
}
$$

Let $P_2^+(\lambda,\lambda',t)$ and $P_2^-(\lambda,\lambda',t)$ denote the central/right and central/left neighboring-spacing joint densities. Then the acceleration source in the $J$ equation can be written schematically as

$$
\boxed{
A_1(\lambda,t)
=
\int_0^\infty
\phi'(\lambda')
\left[P_2^+(\lambda,\lambda',t)+P_2^-(\lambda,\lambda',t)\right]d\lambda'
-2\phi'(\lambda)P(\lambda,t),
}
$$

up to the separately handled boundary contribution.

This is the key governing-equation clue: **the exact nonlinear mechanics does not naturally close on $P(\lambda,t)$ alone. It generates a hierarchy involving velocity information and neighboring-spacing joint information.**

## 6. What this says about the form of $P$

The new conclusion is not a named analytic density. It is a structural form/evolution constraint:

$$
\boxed{
P(\lambda,t)
\text{ is the marginal of a transported spacing-velocity state,}
}
$$

with

$$
\boxed{
\partial_tP=-\partial_\lambda J,
}
$$

and $J$ itself is driven by the exact full nonlinear LJ neighbor statistics.

Therefore any valid explicit form for $P$ must be compatible with all of the following simultaneously:

1. exact spacing-space continuity;
2. the conditional velocity distribution;
3. the full nonlinear LJ force, without expansion about $\lambda=1$;
4. neighboring-spacing joint statistics already observed to remain strongly correlated;
5. boundary loading and fixation.

A global Taylor expansion or a finite harmonic representation is not required and is not used in the active route.

## 7. Next target

The next derivation should determine the minimal exact state among

$$
F_1(\lambda,v,t),
$$

$$
P_2(\lambda,\lambda',t),
$$

or a combined neighboring phase-space state, and derive where the next member of the hierarchy enters.

The objective is not to fit a closure immediately. It is first to expose the exact hierarchy generated by the original nonlinear 1D layer-LJ governing equations and identify which information can be eliminated without changing the predicted $P(\lambda,t)$ and crack-relevant tail.

---

# 한국어 번역 — 1D Spacing Distribution의 정확한 비선형 Transport 단서

## 상태 및 이전 경로의 정정

활성 모델은 계속 엄격한 1차원, normal-only, layer 기반이다.

이번 milestone부터 **harmonic/Taylor 경로를 전체 spacing distribution의 활성 유도법에서 제외한다.** 선형 mode, arcsine push-forward, LJ force의 Taylor 전개는 과거의 국소 진단으로 보존할 수는 있지만 전체 $P(\lambda,t)$를 결정하는 데 사용하지 않는다. 이유는 다음과 같다.

1. 유한 개 harmonic을 가정하면 driven boundary-value problem에서 보장되지 않은 공간 주기성을 미리 넣게 된다.
2. $\lambda=1$ 주변 Taylor 전개는 국소적이며 tensile/compression tail을 포함한 전체 support에서 controlled approximation이 아니다.
3. 목표 이론은 양의 spacing 전체 영역에서 원래 generalized-LJ 비선형 force를 그대로 유지해야 한다.

아래에서는 Taylor 전개, harmonic ansatz, Gaussian/Weibull family, damping law, empirical fatigue-damage variable을 사용하지 않는다.

## 1. 미시 지배방정식

내부 normalized layer spacing은

$$
\boxed{
\ddot\lambda_i
=
\phi'(\lambda_{i+1})
-2\phi'(\lambda_i)
+\phi'(\lambda_{i-1})
}
$$

를 따른다.

이는 원래의 비선형 1D layer-LJ 식이다. 경계 spacing은 실제 boundary force/fixation law를 따라야 하며 bulk equation으로 조용히 대체하지 않는다.

## 2. 정확한 finite-$M$ empirical density와 flux

$M$개의 represented spacing에 대해

$$
\boxed{
P_M(\lambda,t)
=
\frac1M\sum_{i=1}^{M}
\delta[\lambda-\lambda_i(t)]
}
$$

를 정의한다.

spacing-space probability flux는

$$
\boxed{
J_M(\lambda,t)
=
\frac1M\sum_{i=1}^{M}
\dot\lambda_i(t)
\delta[\lambda-\lambda_i(t)]
}
$$

이다.

empirical measure를 distribution 의미에서 시간미분하면 정확히

$$
\boxed{
\partial_tP_M+\partial_\lambda J_M=0
}
$$

을 얻는다.

**분류: EXACT / IDENTITY.**

이 식 자체가 $P$의 형식에 대한 직접적인 단서다. density는 named probability family에서 선택되는 것이 아니라 mechanics가 만드는 flux $J$에 의해 spacing space에서 이동·압축·팽창한다.

매끄러운 continuum density가 존재하고 $P>0$이면 conditional mean spacing velocity를

$$
\boxed{
\bar v(\lambda,t)=\frac{J(\lambda,t)}{P(\lambda,t)}
=\mathbb E[\dot\lambda_i\mid\lambda_i=\lambda]
}
$$

로 정의할 수 있다.

그러면

$$
\boxed{
\partial_tP+\partial_\lambda(P\bar v)=0
}
$$

이다.

$\bar v$가 충분히 매끄러운 영역에서 characteristic은

$$
\dot\lambda=\bar v(\lambda,t)
$$

를 따르고 그 characteristic을 따라

$$
\boxed{
\frac{d}{dt}\ln P
=-\partial_\lambda\bar v
}
$$

이다.

즉 spacing-space flow가 국소적으로 모이면 density가 증가하고 퍼지면 density가 감소한다. 이는 transport statement이며 fitted distribution law가 아니다.

## 3. 왜 $P(\lambda,t)$만으로는 자연스러운 완전상태가 아닌가

미시 mechanics는 시간에 대해 2차식이다. spacing velocity

$$
v_i(t)=\dot\lambda_i(t)
$$

를 도입하고 finite-$M$ empirical phase-space measure를

$$
\boxed{
F_M(\lambda,v,t)
=
\frac1M\sum_i
\delta[\lambda-\lambda_i(t)]
\delta[v-v_i(t)]
}
$$

로 정의한다.

그러면

$$
P_M(\lambda,t)=\int F_M(\lambda,v,t)\,dv
$$

이고

$$
J_M(\lambda,t)=\int vF_M(\lambda,v,t)\,dv
$$

이다.

acceleration flux는

$$
\boxed{
G_M(\lambda,v,t)
=
\frac1M\sum_i
\ddot\lambda_i(t)
\delta[\lambda-\lambda_i(t)]
\delta[v-v_i(t)]
}
$$

로 둔다.

직접 미분하면

$$
\boxed{
\partial_tF_M
+\partial_\lambda(vF_M)
+\partial_vG_M
=0
}
$$

을 얻는다.

**분류: finite empirical measure에 대한 EXACT / IDENTITY.**

매끄러운 phase-space density가 존재하면

$$
G=F\,\bar a(\lambda,v,t)
$$

로 쓸 수 있고

$$
\bar a
=
\mathbb E[\ddot\lambda_i\mid\lambda_i=\lambda,\dot\lambda_i=v]
$$

이다.

따라서 kinetic transport equation은

$$
\boxed{
\partial_tF
+\partial_\lambda(vF)
+\partial_v(\bar aF)
=0
}
$$

이다.

이는 one-point spacing density의 구조적 한계를 보여준다. 같은 $P(\lambda,t)$를 가지더라도 conditional velocity distribution이 다르면 이후 evolution이 달라질 수 있다.

## 4. 정확한 one-point moment hierarchy

$$
K(\lambda,t)=\int v^2F(\lambda,v,t)\,dv
$$

를 정의한다.

phase-space transport를 velocity에 대해 적분하면 다시

$$
\partial_tP+\partial_\lambda J=0
$$

을 얻는다.

$v$를 곱해 적분하면

$$
\boxed{
\partial_tJ+\partial_\lambda K=A_1
}
$$

이고

$$
A_1(\lambda,t)
=
\int G(\lambda,v,t)\,dv
=P(\lambda,t)\mathbb E[\ddot\lambda_i\mid\lambda_i=\lambda]
$$

이다.

따라서 probability flux의 evolution조차 velocity second moment와 conditional acceleration을 필요로 한다.

finite empirical raw moment

$$
M_r(t)=\frac1M\sum_i\lambda_i^r
$$

에 대해서는 정확히

$$
\boxed{
\dot M_r
=r\left\langle\lambda^{r-1}v\right\rangle
}
$$

및

$$
\boxed{
\ddot M_r
=r(r-1)\left\langle\lambda^{r-2}v^2\right\rangle
+r\left\langle\lambda^{r-1}a\right\rangle
}
$$

가 성립한다.

closure assumption은 없다.

## 5. 원래 비선형 LJ 식은 정확히 neighbor joint state를 요구한다

내부 spacing에 대해

$$
\ddot\lambda_i
=
\phi'(\lambda_{i+1})
-2\phi'(\lambda_i)
+\phi'(\lambda_{i-1})
$$

이다.

central spacing이 $\lambda$라는 조건에서 conditional expectation을 취하면

$$
\boxed{
\mathbb E[\ddot\lambda_i\mid\lambda_i=\lambda]
=
\mathbb E[\phi'(\lambda_{i+1})+\phi'(\lambda_{i-1})\mid\lambda_i=\lambda]
-2\phi'(\lambda)
}
$$

이다.

$P_2^+(\lambda,\lambda',t)$와 $P_2^-(\lambda,\lambda',t)$를 각각 central/right 및 central/left neighboring-spacing joint density라고 하면 $J$ equation의 acceleration source는 경계항을 별도로 처리할 때 개략적으로

$$
\boxed{
A_1(\lambda,t)
=
\int_0^\infty
\phi'(\lambda')
\left[P_2^+(\lambda,\lambda',t)+P_2^-(\lambda,\lambda',t)\right]d\lambda'
-2\phi'(\lambda)P(\lambda,t)
}
$$

로 쓸 수 있다.

이게 핵심 지배방정식 단서다. **정확한 비선형 mechanics는 자연스럽게 $P(\lambda,t)$ 하나에서 닫히지 않고 velocity 정보와 neighboring-spacing joint 정보를 포함하는 hierarchy를 만든다.**

## 6. 이것이 $P$의 형식에 대해 말해주는 것

새 결론은 named analytic density가 아니다. 대신 구조적인 형식/evolution constraint다.

$$
\boxed{
P(\lambda,t)
\text{는 transported spacing-velocity state의 marginal이다}
}
$$

이고

$$
\boxed{
\partial_tP=-\partial_\lambda J
}
$$

이며 $J$는 원래의 full nonlinear LJ neighbor statistics에 의해 움직인다.

따라서 유효한 explicit $P$ 형식은 최소한 다음을 동시에 만족해야 한다.

1. 정확한 spacing-space continuity;
2. conditional velocity distribution;
3. $\lambda=1$ 주변 전개를 하지 않은 full nonlinear LJ force;
4. 이미 강하게 남는다고 확인된 neighboring-spacing joint statistics;
5. boundary loading과 fixation.

전역 Taylor 전개나 유한 harmonic 표현은 필요하지 않으며 active route에서 사용하지 않는다.

## 7. 다음 목표

다음 derivation에서는

$$
F_1(\lambda,v,t),
$$

$$
P_2(\lambda,\lambda',t),
$$

또는 neighboring phase-space state 중 어떤 것이 최소 exact state인지 확인하고 hierarchy의 다음 member가 정확히 어디서 등장하는지 유도해야 한다.

목표는 곧바로 fitted closure를 만드는 것이 아니다. 먼저 원래 nonlinear 1D layer-LJ 지배방정식이 만드는 exact hierarchy를 드러내고, 어떤 정보를 제거해도 predicted $P(\lambda,t)$와 crack-relevant tail이 변하지 않는지 확인하는 것이다.
