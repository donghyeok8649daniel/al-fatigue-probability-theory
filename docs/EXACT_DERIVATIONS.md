# Exact Derivations

This file records only identities or equations that are exact under an explicitly stated microscopic model.

## Thermodynamic-limit spacing density

For local spacing variables $a_i(t)$, define

$$
P_N(a,t)=\frac1N\sum_{i=1}^N\delta(a-a_i(t)),
\qquad
P(a,t)=\lim_{N\to\infty}P_N(a,t).
$$

## Exact continuity equation in spacing space

Differentiating the empirical density under smooth deterministic dynamics gives

$$
\partial_tP+\partial_a(Pv)=0,
$$

where

$$
v(a,t)=\langle \dot a_i\mid a_i=a\rangle.
$$

This equation is exact as a kinematic identity. It is not, by itself, a closed evolution law because $v(a,t)$ depends on unresolved microscopic information.

## Moment identities

For any sufficiently regular moment,

$$
\frac{d}{dt}\langle a^n\rangle
=n\langle a^{n-1}v\rangle.
$$

In particular,

$$
\dot{\bar a}=\langle v\rangle,
$$

and

$$
\frac{d}{dt}\operatorname{Var}(a)=2\operatorname{Cov}(a,v).
$$

## Phase-space lift

Because atomistic mechanics is second order, introduce a joint spacing-velocity density

$$
F(a,c,t),\qquad c=\dot a,
$$

with

$$
P(a,t)=\int F(a,c,t)\,dc.
$$

The projected Liouville form is

$$
\partial_tF+\partial_a(cF)+\partial_c(AF)=0,
$$

where

$$
A(a,c,t)=\langle \ddot a_i\mid a_i=a,\dot a_i=c\rangle.
$$

## Nearest-neighbor chain relation

For $V=\sum_i v(a_i)$ and $a_i=x_{i+1}-x_i$, Newton's equations imply

$$
m\ddot a_i=v'(a_{i+1})-2v'(a_i)+v'(a_{i-1}).
$$

Therefore a one-point density $P(a,t)$ is generally not exactly closed; neighboring-spacing correlations enter the acceleration.

## Exact pair-distance energy hierarchy

The distance to the $k$-th neighbor is

$$
R_i^{(k)}=a_i+a_{i+1}+\cdots+a_{i+k-1}.
$$

If $P_k(r,t)$ is the exact density of $R_i^{(k)}$, then the configurational pair-potential energy per atom is

$$
\mathcal U(t)=\sum_{k=1}^{\infty}\int_0^\infty v(r)P_k(r,t)\,dr.
$$

For a uniform lattice $P_k(r)=\delta(r-ka)$, this reduces to

$$
U(a)=\sum_{k=1}^{\infty}v(ka).
$$

---

# 한국어 번역 — 정확한 유도

이 파일에는 명시적으로 정의된 미시모델 아래에서 정확하게 성립하는 항등식 또는 방정식만 기록한다.

## 열역학적 극한의 원자간격 밀도

국부 원자간격 변수 $a_i(t)$에 대해 다음을 정의한다.

$$
P_N(a,t)=\frac1N\sum_{i=1}^N\delta(a-a_i(t)),
\qquad
P(a,t)=\lim_{N\to\infty}P_N(a,t).
$$

## spacing space에서의 정확한 연속방정식

매끄러운 결정론적 동역학 아래에서 empirical density를 시간미분하면

$$
\partial_tP+\partial_a(Pv)=0
$$

을 얻는다. 여기서

$$
v(a,t)=\langle \dot a_i\mid a_i=a\rangle
$$

이다.

이 식은 운동학적 항등식으로서 정확하다. 다만 $v(a,t)$가 해소되지 않은 미시정보에 의존하므로 이 식 하나만으로는 닫힌 evolution law가 아니다.

## moment 항등식

충분히 정칙한 임의의 moment에 대해

$$
\frac{d}{dt}\langle a^n\rangle
=n\langle a^{n-1}v\rangle
$$

이 정확히 성립한다.

특히

$$
\dot{\bar a}=\langle v\rangle
$$

이고,

$$
\frac{d}{dt}\operatorname{Var}(a)=2\operatorname{Cov}(a,v)
$$

이다.

## 위상공간으로의 확장

원자역학은 2계 시간미분 방정식이므로 spacing과 spacing velocity의 결합밀도를

$$
F(a,c,t),\qquad c=\dot a
$$

로 도입한다. 원래의 spacing density는

$$
P(a,t)=\int F(a,c,t)\,dc
$$

라는 marginal이다.

projected Liouville 형태는

$$
\partial_tF+\partial_a(cF)+\partial_c(AF)=0
$$

이며,

$$
A(a,c,t)=\langle \ddot a_i\mid a_i=a,\dot a_i=c\rangle
$$

이다.

## 최근접 이웃 사슬 관계

$V=\sum_i v(a_i)$ 및 $a_i=x_{i+1}-x_i$인 최근접 이웃 사슬에서는 Newton 방정식으로부터

$$
m\ddot a_i=v'(a_{i+1})-2v'(a_i)+v'(a_{i-1})
$$

을 얻는다.

따라서 one-point density $P(a,t)$는 일반적으로 정확히 닫히지 않는다. 가속도에 인접 spacing의 상관관계가 들어가기 때문이다.

## 정확한 pair-distance energy hierarchy

$i$번째 원자에서 $k$번째 이웃까지의 거리는

$$
R_i^{(k)}=a_i+a_{i+1}+\cdots+a_{i+k-1}
$$

이다.

$P_k(r,t)$를 $R_i^{(k)}$의 정확한 밀도라고 하면 원자당 configurational pair-potential energy는

$$
\mathcal U(t)=\sum_{k=1}^{\infty}\int_0^\infty v(r)P_k(r,t)\,dr
$$

이다.

균일격자에서 $P_k(r)=\delta(r-ka)$이므로

$$
U(a)=\sum_{k=1}^{\infty}v(ka)
$$

가 복원된다.
