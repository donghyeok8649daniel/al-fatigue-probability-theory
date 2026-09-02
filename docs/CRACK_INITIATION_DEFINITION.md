# Operational definition of crack initiation — active 1D mainline

## 1. Scope

The active paper stops at **initiation / local stability loss**. It does not
model propagation of an already formed macroscopic crack.

The active microscopic mainline is the 1D normal generalized-LJ chain. The
current operational instability threshold is therefore derived from the same
normalized interaction that drives the chain.

No Smoluchowski closure is required for the definition below.

---

## 2. Local mechanical initiation threshold

For the normalized generalized-LJ energy

$$
\phi(\lambda)
=
\frac{\lambda^{-m}}{m(m-n)}
-
\frac{\lambda^{-n}}{n(m-n)},
$$

the local tangent stiffness is $\phi''(\lambda)$. The operational initiation
stretch is

$$
\boxed{
\phi''(\lambda_c)=0.
}
$$

Hence

$$
\boxed{
\lambda_c
=\left(\frac{m+1}{n+1}\right)^{1/(m-n)}.
}
$$

For the current $m=12.19$, $n=6$ baseline,

$$
\boxed{\lambda_c\approx1.1077715386.}
$$

This criterion means that the ideal intact local bond/layer description has
lost positive tangent stiffness. It does **not** claim that a macroscopic free
surface has already fully formed.

---

## 3. Instantaneous tail is not cumulative initiation

For a nonabsorbing spacing density $P(\lambda,t)$, define the instantaneous
critical tail

$$
\boxed{
Q_c(t)=\int_{\lambda_c}^{\infty}P(\lambda,t)\,d\lambda.
}
$$

This mass can return below $\lambda_c$ in a reflecting/nonabsorbing calculation.
Therefore

$$
\boxed{
Q_c(t)\neq F_{\rm ci}(t)
}
$$

in general, where $F_{\rm ci}$ denotes cumulative first passage.

---

## 4. Exact finite empirical first passage

For each represented spacing trajectory $\lambda_i(t)$ define

$$
\boxed{
\tau_i^c
=\inf\{t\ge0:\lambda_i(t)\ge\lambda_c\}.
}
$$

Once this time has occurred, that trajectory is counted as initiated even if a
nonabsorbing continuation of the mechanics would later return below
$\lambda_c$.

Define

$$
\chi_i(t)=\mathbf1_{\{t<\tau_i^c\}}.
$$

The local survivor fraction is

$$
\boxed{
S_M(t)=\frac1M\sum_i\chi_i(t),
}
$$

and the cumulative local first-passage fraction is

$$
\boxed{
F_{{\rm ci},M}^{\rm local}(t)=1-S_M(t).
}
$$

In the distributional sense,

$$
\boxed{
-\dot S_M(t)
=\frac1M\sum_i\delta(t-\tau_i^c).
}
$$

These are exact definitions for the finite empirical spacing population.

---

## 5. Smooth phase-space absorbing formulation

Let $F_b(\lambda,c,t)$ denote the intact/survivor phase-space subdensity on

$$
0<\lambda<\lambda_c,
\qquad c=\dot\lambda.
$$

The interior projected transport equation remains

$$
\partial_tF_b+\partial_\lambda(cF_b)+\partial_c(AF_b)=0.
$$

Because the underlying dynamics is second order, an absorbing boundary is most
naturally stated as **no inflow from the failed side**. At the right boundary,

$$
\boxed{
F_b(\lambda_c,c,t)=0
\quad\text{for incoming }c<0.
}
$$

Outgoing states with $c>0$ are allowed to leave. The escape flux is

$$
\boxed{
j_{\rm esc}(t)
=\int_0^{\infty}cF_b(\lambda_c^-,c,t)\,dc\ge0.
}
$$

Assuming no loss through the lower-spacing boundary,

$$
\boxed{
S(t)
=\int_0^{\lambda_c}\int_{-\infty}^{\infty}
F_b(\lambda,c,t)\,dc\,d\lambda
}
$$

obeys

$$
\boxed{
\dot S(t)=-j_{\rm esc}(t).
}
$$

Therefore

$$
\boxed{
F_{\rm ci}^{\rm local}(t)=1-S(t)
}
$$

and, for $S>0$,

$$
\boxed{
h(t)=\frac{j_{\rm esc}(t)}{S(t)}=-\frac{d}{dt}\ln S(t).
}
$$

This is first passage of the mechanically generated phase-space state, not an
empirical fatigue-damage law.

---

## 6. Survivor-conditioned observables

The intact spacing marginal

$$
P_b(\lambda,t)=\int F_b(\lambda,c,t)\,dc
$$

is a subdensity with

$$
\int_0^{\lambda_c}P_b\,d\lambda=S(t),
$$

not one. If a normalized conditional survivor distribution is needed, define

$$
\boxed{
\widehat P_b(\lambda,t)=\frac{P_b(\lambda,t)}{S(t)}.
}
$$

Then, for example,

$$
\boxed{
\bar\lambda_{\rm surv}
=\frac1S\int_0^{\lambda_c}\lambda P_b\,d\lambda
}
$$

and

$$
\boxed{
\bar U_{\rm surv}
=\frac{U_{\rm ref}}{S}
\int_0^{\lambda_c}[\phi(\lambda)-\phi(1)]P_b\,d\lambda.
}
$$

Unnormalized moments over $P_b$ must not be silently called conditional means.

---

## 7. Local first-passage fraction versus specimen probability

The empirical measure $P_M$ used by the deterministic chain is, first of all, a
spatial counting distribution across represented spacings. Therefore

$$
1-S_M(t)
$$

is the fraction of represented local spacings that have experienced first
passage.

For one deterministic chain realization, the specimen first-initiation time is

$$
\boxed{
\tau_{\rm spec}^c=\min_i\tau_i^c.
}
$$

A specimen-level probability requires an ensemble of microscopic/specimen
realizations $\omega$:

$$
\boxed{
S_{\rm spec}(t)
=\Pr_\omega\left[\min_i\tau_i^c(\omega)>t\right].
}
$$

In general,

$$
\boxed{
1-S_M(t)
\neq
\Pr(\tau_{\rm spec}^c\le t).
}
$$

No independence product over atoms, statistical cells, or FEM elements is
assumed. The specimen-scale probability bridge remains an open calibration /
correlation problem.

---

## 8. Relation to G4

The active normalization/survival equation is therefore

$$
\boxed{
\int_0^{\infty}P(\lambda,t)\,d\lambda=1
}
$$

for a nonabsorbing intact calculation, or

$$
\boxed{
S(t)=\int_0^{\lambda_c}P_b(\lambda,t)\,d\lambda\le1
}
$$

for the first-passage survivor subdensity.

The cumulative local initiation fraction is $1-S$.

This document does not promote the earlier overdamped/Smoluchowski one-cycle
operator to the active mainline. Such an operator can only be reintroduced if
its additional dynamical assumptions are independently justified.

---

## 한국어 요약

현재 균열개시는 같은 1D generalized-LJ chain에서

$$
\phi''(\lambda_c)=0
$$

이 되는 국소 접선강성 상실점의 **첫 도달**로 정의한다.

비흡수 계산의 순간 tail

$$
\int_{\lambda_c}^{\infty}P\,d\lambda
$$

은 다시 돌아올 수 있으므로 누적 균열개시확률이 아니다.

2차 동역학의 흡수조건은 phase space에서 실패영역으로 나가는 $c>0$ flux는
허용하고, 실패영역에서 다시 들어오는 $c<0$ inflow를 차단하는 방식으로
쓴다. 그러면

$$
\dot S=-j_{\rm esc},
\qquad
F_{\rm ci}^{\rm local}=1-S,
\qquad
h=j_{\rm esc}/S
$$

가 된다.

다만 현재 $P_M$은 우선 한 시편 내부 spacing들의 공간 empirical
distribution이므로 $1-S$를 specimen-to-specimen 균열확률로 바로 동일시하지
않는다. specimen probability에는 별도의 realization ensemble과 spatial
correlation/statistical-length bridge가 필요하다.
