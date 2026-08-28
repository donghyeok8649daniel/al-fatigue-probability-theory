# Milestone 11 — Exact Instantaneous Shape Constraint for $P(\lambda,t)$

## Scope

The priority of this milestone is the functional form of the one-point normalized layer-spacing density $P(\lambda,t)$. The derivation uses only the exact one-point moment equations obtained from the full nonlinear 1D layer-LJ dynamics. It introduces no Taylor expansion of the LJ force, finite harmonic ansatz, Gaussian/Weibull family, empirical damping, or fatigue-damage variable.

The finite-$M$ empirical density is still exactly

$$
P_M(\lambda,t)=\frac1M\sum_i\delta[\lambda-\lambda_i(t)].
$$

A smooth $P(\lambda,t)$ below is a continuum representation of this empirical measure, not a claim that a finite represented chain literally has a smooth density.

## 1. Exact first two phase-space moments

Let $F(\lambda,v,t)$ be the spacing-velocity phase-space density. Define

$$
P(\lambda,t)=\int F(\lambda,v,t)\,dv,
$$

$$
u(\lambda,t)=\mathbb E[v\mid\lambda],
$$

and

$$
\Theta(\lambda,t)=\operatorname{Var}(v\mid\lambda).
$$

Then

$$
J=P u,
$$

and the second velocity-moment density is

$$
K=P(u^2+\Theta).
$$

Let

$$
\bar a(\lambda,t)
=
\mathbb E[\ddot\lambda_i\mid\lambda_i=\lambda]
$$

be the one-point conditional spacing acceleration.

The exact smooth moment equations are

$$
\boxed{
\partial_tP+\partial_\lambda(Pu)=0
}
$$

and

$$
\boxed{
\partial_t(Pu)
+
\partial_\lambda\left[P(u^2+\Theta)\right]
=
P\bar a.
}
$$

These are **EXACT / IDENTITY** statements for the smooth moment representation wherever the required moments exist and boundary terms in velocity vanish.

## 2. Exact equation for the logarithmic slope of $P$

Using the continuity equation to simplify the momentum equation gives

$$
\boxed{
D_tu
=
\bar a
-
\frac1P\partial_\lambda(P\Theta),
}
$$

where

$$
D_tu
=
\partial_tu+u\partial_\lambda u.
$$

Expanding the last term yields

$$
\bar a-D_tu
=
\partial_\lambda\Theta
+
\Theta\partial_\lambda\ln P.
$$

Therefore, wherever $P>0$ and $\Theta>0$,

$$
\boxed{
\partial_\lambda\ln P
=
\frac{\bar a-D_tu}{\Theta}
-
\partial_\lambda\ln\Theta.
}
$$

This is the central result of the milestone. It directly constrains the shape of $P$ from the exact transport fields.

## 3. Exact instantaneous functional representation

Integrating the log-slope at fixed time gives

$$
\boxed{
P(\lambda,t)
=
\frac{C(t)}{\Theta(\lambda,t)}
\exp\left[
\int_{\lambda_*}^{\lambda}
\frac{\bar a(s,t)-D_tu(s,t)}{\Theta(s,t)}\,ds
\right].
}
$$

$C(t)$ is fixed by normalization,

$$
\int_0^\infty P(\lambda,t)\,d\lambda=1.
$$

This is an **EXACT INSTANTANEOUS SHAPE REPRESENTATION under the stated smoothness and positivity conditions**. It is not yet a closed solution because the three fields

$$
\Theta(\lambda,t),
\qquad
D_tu(\lambda,t),
\qquad
\bar a(\lambda,t)
$$

must themselves be obtained from the mechanics.

The important change in viewpoint is that the shape of $P$ is not selected from a named distribution family. Its logarithmic slope is the balance between conditional acceleration, mean-flow material acceleration, and the gradient of conditional velocity variance.

## 4. Full nonlinear LJ enters through the conditional acceleration

For an interior spacing the exact microscopic equation is

$$
\ddot\lambda_i
=
\phi'(\lambda_{i+1})
-2\phi'(\lambda_i)
+
\phi'(\lambda_{i-1}).
$$

Therefore

$$
\boxed{
\bar a(\lambda,t)
=
m_+(\lambda,t)
+m_-(\lambda,t)
-2\phi'(\lambda),
}
$$

where

$$
m_+(\lambda,t)
=
\mathbb E[\phi'(\lambda_{i+1})\mid\lambda_i=\lambda],
$$

and

$$
m_-(\lambda,t)
=
\mathbb E[\phi'(\lambda_{i-1})\mid\lambda_i=\lambda].
$$

Equivalently, with ordered neighboring-spacing joint densities,

$$
m_+(\lambda,t)
=
\frac{
\int\phi'(\lambda')P_2^+(\lambda,\lambda',t)\,d\lambda'
}{P(\lambda,t)},
$$

and analogously for $m_-$. No neighbor-independence assumption is made.

Substitution into the exact shape law gives

$$
\boxed{
\partial_\lambda\ln P
=
\frac{
m_++m_- -2\phi'(\lambda)-D_tu
}{\Theta}
-
\partial_\lambda\ln\Theta.
}
$$

This is currently the most direct governing-equation clue for the function form of $P$.

## 5. Why the earlier exponential-LJ closure appeared

The previous large-$M$ closure had the form

$$
P\propto\exp[-\alpha\lambda-\beta\phi(\lambda)].
$$

The exact shape equation now shows one mechanical route by which such a form can appear as a restrictive special case.

Assume only for this diagnostic special case that:

1. $\Theta$ is independent of $\lambda$ at the considered instant;
2. $D_tu=0$;
3. the conditional left/right neighbor-force means are independent of the central spacing, so $m_++m_-=2f_N(t)$.

Then

$$
\partial_\lambda\ln P
=
\frac{2f_N-2\phi'(\lambda)}{\Theta},
$$

and therefore

$$
\boxed{
P(\lambda,t)
\propto
\exp\left[
\frac{2f_N(t)}{\Theta(t)}\lambda
-
\frac{2}{\Theta(t)}\phi(\lambda)
\right].
}
$$

Thus the old exponential-LJ family is recovered without maximum entropy, but only after imposing strong special conditions. The previously measured strong spatial correlations directly violate the spirit of the neighbor-independence condition, explaining why that closure need not reproduce the driven distribution.

This special case is **NOT** promoted to the active global model.

## 6. What the exact shape law says about tails and peaks

At any smooth point with $\Theta>0$,

$$
\operatorname{sign}(\partial_\lambda P)
=
\operatorname{sign}
\left[
\bar a-D_tu-\partial_\lambda\Theta
\right]
$$

because $P/\Theta>0$.

A local stationary point of the density therefore satisfies

$$
\boxed{
\bar a-D_tu
=
\partial_\lambda\Theta.
}
$$

This gives a mechanics-based peak condition without assuming symmetry or a named density.

For the tensile tail, its logarithmic decay or growth rate is exactly controlled by

$$
\boxed{
\partial_\lambda\ln P
=
\frac{m_++m_- -2\phi'(\lambda)-D_tu}{\Theta}
-
\partial_\lambda\ln\Theta.
}
$$

Hence the crack-relevant tail cannot be determined from the LJ potential alone. It also requires the conditional velocity spread and neighboring-spacing force statistics. This is a structural conclusion, not a fitting choice.

## 7. Degenerate case $\Theta=0$

The shape formula divides by $\Theta$ and therefore does not apply where the conditional velocity distribution is monokinetic,

$$
\Theta(\lambda,t)=0.
$$

This is not a numerical nuisance. It marks a genuinely different transport regime. In a monokinetic region the phase-space density collapses onto a velocity graph and $P$ must be evolved through the first-order continuity equation instead. Broadening from an initially singular or multi-stream state therefore cannot be inferred from the smooth shape formula alone.

## 8. Current conclusion

The governing equations now constrain the full-support one-point density to the instantaneous form

$$
\boxed{
P
=
\Theta^{-1}
\times
\exp\left[
\int
\frac{\text{conditional LJ acceleration}-\text{mean-flow material acceleration}}
{\text{conditional velocity variance}}
\,d\lambda
\right]
\times
\text{normalization}.
}
$$

The remaining problem is no longer "choose the form of $P$." It is to determine the mechanically generated fields $\Theta$, $D_tu$, and $m_\pm$ with the minimum additional state information.

## 9. Next falsification target

The next numerical step should record spacing and spacing velocity together in the deterministic chain, estimate

$$
u(\lambda,t),
\qquad
\Theta(\lambda,t),
\qquad
\bar a(\lambda,t),
$$

and test whether the measured histogram log-slope satisfies the exact shape identity after grid/bin refinement. This tests the representation itself before any closure for $P_2$ is introduced.

---

# 한국어 번역 — $P(\lambda,t)$의 정확한 순간 함수형 제약

## 범위

이번 milestone의 최우선 목표는 one-point normalized layer-spacing density $P(\lambda,t)$의 함수형이다. 유도에는 full nonlinear 1D layer-LJ dynamics에서 얻는 정확한 one-point moment equation만 사용한다. LJ force의 Taylor 전개, finite harmonic 가정, Gaussian/Weibull family, 경험적 damping, fatigue-damage variable은 사용하지 않는다.

finite-$M$ empirical density는 여전히 정확히

$$
P_M(\lambda,t)=\frac1M\sum_i\delta[\lambda-\lambda_i(t)]
$$

이다.

아래의 smooth $P(\lambda,t)$는 이 empirical measure의 continuum representation이며, finite represented chain 자체가 문자 그대로 smooth density를 가진다는 주장이 아니다.

## 1. 정확한 phase-space 1차·2차 moment

$F(\lambda,v,t)$를 spacing-velocity phase-space density라고 하자. 다음을 정의한다.

$$
P(\lambda,t)=\int F(\lambda,v,t)\,dv,
$$

$$
u(\lambda,t)=\mathbb E[v\mid\lambda],
$$

그리고

$$
\Theta(\lambda,t)=\operatorname{Var}(v\mid\lambda).
$$

그러면

$$
J=Pu
$$

이고 second velocity-moment density는

$$
K=P(u^2+\Theta)
$$

이다.

또

$$
\bar a(\lambda,t)
=
\mathbb E[\ddot\lambda_i\mid\lambda_i=\lambda]
$$

를 one-point conditional spacing acceleration으로 정의한다.

정확한 smooth moment equation은

$$
\boxed{
\partial_tP+\partial_\lambda(Pu)=0
}
$$

및

$$
\boxed{
\partial_t(Pu)
+
\partial_\lambda\left[P(u^2+\Theta)\right]
=
P\bar a
}
$$

이다.

필요한 moment가 존재하고 velocity 방향 boundary term이 사라지는 smooth moment representation에서는 이 식들은 **EXACT / IDENTITY**다.

## 2. $P$의 logarithmic slope에 대한 정확한 식

continuity equation을 이용해 momentum equation을 정리하면

$$
\boxed{
D_tu
=
\bar a
-
\frac1P\partial_\lambda(P\Theta)
}
$$

를 얻는다. 여기서

$$
D_tu
=
\partial_tu+u\partial_\lambda u
$$

이다.

마지막 항을 전개하면

$$
\bar a-D_tu
=
\partial_\lambda\Theta
+
\Theta\partial_\lambda\ln P
$$

이므로 $P>0$, $\Theta>0$인 곳에서

$$
\boxed{
\partial_\lambda\ln P
=
\frac{\bar a-D_tu}{\Theta}
-
\partial_\lambda\ln\Theta
}
$$

가 성립한다.

이 식이 이번 milestone의 핵심 결과다. exact transport field로부터 $P$의 shape를 직접 제약한다.

## 3. 정확한 순간 함수형 표현

고정된 시간에서 log-slope를 적분하면

$$
\boxed{
P(\lambda,t)
=
\frac{C(t)}{\Theta(\lambda,t)}
\exp\left[
\int_{\lambda_*}^{\lambda}
\frac{\bar a(s,t)-D_tu(s,t)}{\Theta(s,t)}\,ds
\right]
}
$$

를 얻는다.

$C(t)$는

$$
\int_0^\infty P(\lambda,t)\,d\lambda=1
$$

이라는 normalization으로 정한다.

이 식은 명시한 smoothness와 positivity 조건 아래 **EXACT INSTANTANEOUS SHAPE REPRESENTATION**이다. 다만

$$
\Theta(\lambda,t),
\qquad
D_tu(\lambda,t),
\qquad
\bar a(\lambda,t)
$$

세 field를 mechanics에서 별도로 얻어야 하므로 아직 closed solution은 아니다.

중요한 관점 변화는 $P$의 shape를 이름 붙은 distribution family에서 고르는 것이 아니라, conditional acceleration, mean-flow material acceleration, conditional velocity variance gradient의 balance가 $P$의 logarithmic slope를 결정한다는 점이다.

## 4. Full nonlinear LJ는 conditional acceleration을 통해 직접 들어간다

interior spacing의 정확한 microscopic equation은

$$
\ddot\lambda_i
=
\phi'(\lambda_{i+1})
-2\phi'(\lambda_i)
+
\phi'(\lambda_{i-1})
$$

이다.

따라서

$$
\boxed{
\bar a(\lambda,t)
=
m_+(\lambda,t)
+m_-(\lambda,t)
-2\phi'(\lambda)
}
$$

이고

$$
m_+(\lambda,t)
=
\mathbb E[\phi'(\lambda_{i+1})\mid\lambda_i=\lambda],
$$

$$
m_-(\lambda,t)
=
\mathbb E[\phi'(\lambda_{i-1})\mid\lambda_i=\lambda]
$$

이다.

ordered neighboring-spacing joint density를 쓰면

$$
m_+(\lambda,t)
=
\frac{
\int\phi'(\lambda')P_2^+(\lambda,\lambda',t)\,d\lambda'
}{P(\lambda,t)}
$$

이고 $m_-$도 같은 방식으로 쓸 수 있다. neighbor independence나 left/right symmetry는 가정하지 않는다.

이를 exact shape law에 대입하면

$$
\boxed{
\partial_\lambda\ln P
=
\frac{
m_++m_- -2\phi'(\lambda)-D_tu
}{\Theta}
-
\partial_\lambda\ln\Theta
}
$$

가 된다.

현재로서는 이것이 지배방정식에서 직접 얻은 $P$ 함수형에 대한 가장 직접적인 단서다.

## 5. 기존 exponential-LJ closure가 왜 나타났는가

이전에 사용했던 large-$M$ closure는

$$
P\propto\exp[-\alpha\lambda-\beta\phi(\lambda)]
$$

꼴이었다.

이번 exact shape equation을 이용하면 이 형식이 제한적인 특수조건 아래 mechanics에서 어떻게 나타날 수 있는지 알 수 있다.

오직 이 diagnostic special case에서만 다음을 가정하자.

1. 고려하는 순간에 $\Theta$가 $\lambda$에 무관하다.
2. $D_tu=0$이다.
3. left/right conditional neighbor-force mean이 central spacing에 무관해서 $m_++m_-=2f_N(t)$이다.

그러면

$$
\partial_\lambda\ln P
=
\frac{2f_N-2\phi'(\lambda)}{\Theta}
$$

이고 따라서

$$
\boxed{
P(\lambda,t)
\propto
\exp\left[
\frac{2f_N(t)}{\Theta(t)}\lambda
-
\frac{2}{\Theta(t)}\phi(\lambda)
\right]
}
$$

가 된다.

즉 maximum entropy를 쓰지 않고도 old exponential-LJ family를 복원할 수 있지만, 강한 특수조건이 필요하다. 앞에서 측정된 강한 spatial correlation은 neighbor-independence 조건의 취지와 직접 충돌하므로 기존 closure가 driven distribution을 정확히 못 맞춘 이유도 자연스럽게 설명된다.

이 special case는 **active global model로 승격하지 않는다.**

## 6. Exact shape law가 tail과 peak에 대해 말해주는 것

$\Theta>0$인 smooth point에서는 $P/\Theta>0$이므로

$$
\operatorname{sign}(\partial_\lambda P)
=
\operatorname{sign}
\left[
\bar a-D_tu-\partial_\lambda\Theta
\right]
$$

이다.

따라서 density의 local stationary point는

$$
\boxed{
\bar a-D_tu
=
\partial_\lambda\Theta
}
$$

를 만족한다.

즉 symmetry나 특정 density를 가정하지 않고 mechanics 기반 peak condition을 얻는다.

또 tensile tail의 logarithmic decay/growth rate는 정확히

$$
\boxed{
\partial_\lambda\ln P
=
\frac{m_++m_- -2\phi'(\lambda)-D_tu}{\Theta}
-
\partial_\lambda\ln\Theta
}
$$

에 의해 정해진다.

따라서 crack-relevant tail은 LJ potential만으로 정할 수 없다. conditional velocity spread와 neighboring-spacing force statistics도 필요하다. 이는 fitting 선택이 아니라 구조적 결론이다.

## 7. $\Theta=0$인 퇴화 경우

shape formula는 $\Theta$로 나누므로 conditional velocity distribution이 monokinetic인

$$
\Theta(\lambda,t)=0
$$

구간에는 적용되지 않는다.

이는 단순 수치 문제가 아니라 실제로 다른 transport regime이다. monokinetic region에서는 phase-space density가 velocity graph 위로 collapse하며 $P$는 first-order continuity equation으로 진화시켜야 한다. 따라서 initially singular 또는 multi-stream state의 broadening을 smooth shape formula 하나만으로 설명할 수는 없다.

## 8. 현재 결론

지배방정식으로부터 full-support one-point density는 순간적으로

$$
\boxed{
P
=
\Theta^{-1}
\times
\exp\left[
\int
\frac{\text{conditional LJ acceleration}-\text{mean-flow material acceleration}}
{\text{conditional velocity variance}}
\,d\lambda
\right]
\times
\text{normalization}
}
$$

형태로 제한된다.

이제 남은 문제는 "$P$의 형식을 고르는 것"이 아니다. $\Theta$, $D_tu$, $m_\pm$를 최소한의 추가 state information으로 mechanics에서 결정하는 것이다.

## 9. 다음 반증시험

다음 numerical step에서는 deterministic chain에서 spacing과 spacing velocity를 동시에 기록하고

$$
u(\lambda,t),
\qquad
\Theta(\lambda,t),
\qquad
\bar a(\lambda,t)
$$

를 추정한 뒤 measured histogram의 log-slope가 exact shape identity를 grid/bin refinement 후에도 만족하는지 확인해야 한다. 이 검증을 먼저 한 다음에야 $P_2$ closure를 도입하는 것이 맞다.
