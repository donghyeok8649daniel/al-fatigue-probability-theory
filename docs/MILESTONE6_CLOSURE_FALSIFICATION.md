# Milestone 6 — Direct Falsification of the Two-Moment 1D Layer-LJ Distribution Closure

## Scope

The active theory remains strictly one-dimensional and normal-only. The represented microscopic coordinate is the normal separation between material layers, and the effective layer interaction is the calibrated generalized Lennard-Jones model.

The current large-system closure is

$$
\boxed{
p_\lambda(\lambda,t)
=
\frac{1}{Z(t)}
\exp\left[-\alpha(t)\lambda-\beta(t)\psi(\lambda)\right].
}
$$

This is a **CONTROLLED APPROXIMATION** derived from an equal-base-measure fixed-length/fixed-configurational-energy ensemble followed by a large-$M$ saddle-point reduction. It is not an exact law of driven deterministic fatigue dynamics.

## 1. Falsification protocol

For a deterministic spacing snapshot $\{\lambda_i(t)\}_{i=1}^{M}$, measure only

$$
\mu_{\rm sim}(t)=\frac1M\sum_i\lambda_i(t)
$$

and

$$
\mathcal E_{\rm sim}(t)=\frac1M\sum_i\psi(\lambda_i(t)).
$$

Then solve

$$
\mu_{\rm closure}=\mu_{\rm sim},
\qquad
\mathcal E_{\rm closure}=\mathcal E_{\rm sim}
$$

for $\alpha(t)$ and $\beta(t)$. No histogram parameter is fitted.

## 2. Tested deterministic states

Two phase-locked snapshots of the existing 32-atom 1D chain were used.

The slower case uses

$$
F_a^*=0.03,
\qquad
\omega^*=0.01
$$

and is sampled at approximately $t=10T$.

The second case uses

$$
F_a^*=0.03,
\qquad
\omega^*=0.02
$$

and is sampled at $t=2T$, before the previously observed first local $\lambda_c$ crossing near $2.25T$.

These are reduced-model dimensionless frequencies, not 20 Hz fatigue predictions.

## 3. Numerical result

| Quantity | slow $t=10T$ | dynamic $t=2T$ |
|---|---:|---:|
| $\mu_{\rm sim}$ | $0.98483617$ | $0.98856118$ |
| $\mathcal E_{\rm sim}$ | $1.50102\times10^{-4}$ | $8.70461\times10^{-5}$ |
| $\operatorname{Var}_{\rm sim}(\lambda)$ | $3.33394\times10^{-5}$ | $2.57645\times10^{-5}$ |
| $\operatorname{Var}_{\rm closure}(\lambda)$ | $3.20035\times10^{-5}$ | $2.56370\times10^{-5}$ |
| variance relative error | $4.01\%$ | $0.495\%$ |
| empirical skewness | $1.06221$ | $0.555922$ |
| closure skewness | $0.044083$ | $0.426772$ |
| empirical $Q_c$ | $0$ | $0$ |
| closure $Q_c$ | $6.75\times10^{-51}$ | $1.32\times10^{-63}$ |
| Kolmogorov distance | $0.15888$ | $0.14584$ |

The quoted comparison is a **NUMERICAL DIAGNOSTIC**, not an exact theorem.

## 4. What worked

The two-moment closure reproduces the spacing variance surprisingly well in these two near-equilibrium snapshots. However, this is not strong independent validation.

Near $\lambda=1$,

$$
\psi(\lambda)=\frac12(\lambda-1)^2+O((\lambda-1)^3),
$$

because

$$
\psi(1)=0,
\qquad
\psi'(1)=0,
\qquad
\psi''(1)=1.
$$

Therefore matching mean configurational energy near equilibrium already constrains the variance strongly.

## 5. What failed

The full empirical distribution is not reproduced exactly. For the slower case,

$$
\gamma_{1,\rm sim}\approx1.06221,
\qquad
\gamma_{1,\rm closure}\approx0.04408,
$$

and

$$
\boxed{D_{\rm KS}\approx0.159.}
$$

For the second case,

$$
\boxed{D_{\rm KS}\approx0.146.}
$$

Thus the statement

$$
\boxed{
\mu(t),\mathcal E(t)
\text{ uniquely determine the driven one-point spacing distribution}
}
$$

is not supported by this direct finite-chain test.

## 6. Why the missing state is not yet uniquely identified

The deterministic sample contains only $M=31$ spacings, while the closure is a large-$M$ saddle-point marginal of an equal-measure ensemble. The mismatch can therefore arise from at least three mechanisms:

1. finite-$M$ corrections;
2. failure of the equal-measure/ergodic assumption because the driven chain retains coherent spatial structure;
3. a genuinely missing one-point state variable such as a third central moment.

It would be incorrect to add a fitted third-moment multiplier immediately and call the problem solved.

## 7. Exact third-moment kinematic identity

Define

$$
\boxed{
m_3(t)=\int[\lambda-\mu(t)]^3p_\lambda(\lambda,t)\,d\lambda.
}
$$

For

$$
\partial_t p_\lambda+\partial_\lambda(p_\lambda v_\lambda)=0
$$

with vanishing boundary fluxes,

$$
\boxed{
\dot m_3(t)
=
3\int[\lambda-\mu(t)]^2[v_\lambda(\lambda,t)-\dot\mu(t)]p_\lambda(\lambda,t)\,d\lambda
}
$$

or equivalently

$$
\boxed{
\dot m_3
=
3\operatorname{Cov}\left((\lambda-\mu)^2,v_\lambda\right).
}
$$

This is an **EXACT / IDENTITY** under the stated continuity equation and boundary conditions.

## 8. Next decisive test

Before enlarging the closure, run

$$
\boxed{M\uparrow\quad\Longrightarrow\quad D_{\rm KS}(M)\ ?}
$$

under the same 1D layer-LJ protocol.

- If $D_{\rm KS}\to0$, the large-$M$ closure remains plausible and the current discrepancy is mainly finite-size/self-averaging error.
- If $D_{\rm KS}$ stays finite, the equal-measure closure is missing dynamical information.
- If the mismatch tracks $m_3(t)$, a mechanics-derived third-moment state becomes a strong next candidate.

---

# 한국어 번역 — 2-Moment 1D Layer-LJ 분포 Closure의 직접 반증시험

## 범위

활성 이론은 계속 엄격하게 1차원 수직변형만 다룬다. 표현하는 미시좌표는 material layer 사이의 수직간격이고, layer 간 유효상호작용은 calibration된 generalized Lennard-Jones model이다.

현재 large-system closure는

$$
\boxed{
p_\lambda(\lambda,t)
=
\frac{1}{Z(t)}
\exp\left[-\alpha(t)\lambda-\beta(t)\psi(\lambda)\right]
}
$$

이다.

이 식은 fixed-length/fixed-configurational-energy manifold에 equal base measure를 둔 뒤 large-$M$ saddle-point reduction으로 얻은 **CONTROLLED APPROXIMATION**이다. driven deterministic fatigue dynamics의 exact law가 아니다.

## 1. 반증 프로토콜

deterministic spacing snapshot $\{\lambda_i(t)\}_{i=1}^{M}$에서 오직

$$
\mu_{\rm sim}(t)=\frac1M\sum_i\lambda_i(t)
$$

와

$$
\mathcal E_{\rm sim}(t)=\frac1M\sum_i\psi(\lambda_i(t))
$$

만 측정한다.

그 다음

$$
\mu_{\rm closure}=\mu_{\rm sim},
\qquad
\mathcal E_{\rm closure}=\mathcal E_{\rm sim}
$$

을 만족하도록 $\alpha(t)$와 $\beta(t)$를 푼다. histogram parameter는 전혀 fitting하지 않는다.

## 2. 시험한 deterministic state

기존 32-atom 1D chain의 phase-locked snapshot 두 개를 사용했다.

느린 case는

$$
F_a^*=0.03,
\qquad
\omega^*=0.01
$$

이고 약 $t=10T$에서 sample했다.

두 번째 case는

$$
F_a^*=0.03,
\qquad
\omega^*=0.02
$$

이고 $t=2T$에서 sample했다. 기존 첫 local $\lambda_c$ crossing 약 $2.25T$보다 이전이다.

이 값들은 reduced-model dimensionless frequency이며 20 Hz fatigue prediction이 아니다.

## 3. 수치결과

| Quantity | slow $t=10T$ | dynamic $t=2T$ |
|---|---:|---:|
| $\mu_{\rm sim}$ | $0.98483617$ | $0.98856118$ |
| $\mathcal E_{\rm sim}$ | $1.50102\times10^{-4}$ | $8.70461\times10^{-5}$ |
| $\operatorname{Var}_{\rm sim}(\lambda)$ | $3.33394\times10^{-5}$ | $2.57645\times10^{-5}$ |
| $\operatorname{Var}_{\rm closure}(\lambda)$ | $3.20035\times10^{-5}$ | $2.56370\times10^{-5}$ |
| variance relative error | $4.01\%$ | $0.495\%$ |
| empirical skewness | $1.06221$ | $0.555922$ |
| closure skewness | $0.044083$ | $0.426772$ |
| empirical $Q_c$ | $0$ | $0$ |
| closure $Q_c$ | $6.75\times10^{-51}$ | $1.32\times10^{-63}$ |
| Kolmogorov distance | $0.15888$ | $0.14584$ |

위 비교는 exact theorem이 아니라 **NUMERICAL DIAGNOSTIC**이다.

## 4. 맞은 부분

두 snapshot에서 two-moment closure가 spacing variance를 상당히 잘 재현했다. 하지만 이는 강한 독립 validation이 아니다.

$\lambda=1$ 근처에서

$$
\psi(\lambda)=\frac12(\lambda-1)^2+O((\lambda-1)^3)
$$

이고

$$
\psi(1)=0,
\qquad
\psi'(1)=0,
\qquad
\psi''(1)=1
$$

이므로 equilibrium 근처에서 mean configurational energy를 맞추는 것 자체가 variance를 강하게 제약한다.

## 5. 실패한 부분

전체 empirical distribution은 정확히 재현되지 않았다. 느린 case에서

$$
\gamma_{1,\rm sim}\approx1.06221,
\qquad
\gamma_{1,\rm closure}\approx0.04408
$$

이고

$$
\boxed{D_{\rm KS}\approx0.159}
$$

이다.

두 번째 case에서도

$$
\boxed{D_{\rm KS}\approx0.146}
$$

이다.

따라서

$$
\boxed{
\mu(t),\mathcal E(t)
\text{만으로 driven one-point spacing distribution이 유일하게 정해진다}
}
$$

라는 주장은 이번 finite-chain 직접시험에서 지지되지 않는다.

## 6. missing state를 아직 하나로 확정할 수 없는 이유

현재 deterministic sample에는 $M=31$개의 spacing만 있고 closure는 equal-measure ensemble의 large-$M$ saddle-point marginal이다. 따라서 mismatch는 최소한 다음 세 원인에서 생길 수 있다.

1. finite-$M$ correction;
2. driven chain이 coherent spatial structure를 유지해 equal-measure/ergodic assumption이 깨지는 경우;
3. third central moment 같은 one-point state variable이 실제로 추가로 필요한 경우.

따라서 fitted third-moment multiplier를 바로 추가하고 문제가 풀렸다고 하면 안 된다.

## 7. 정확한 third-moment kinematic identity

$$
\boxed{
m_3(t)=\int[\lambda-\mu(t)]^3p_\lambda(\lambda,t)\,d\lambda
}
$$

를 정의한다.

$$
\partial_t p_\lambda+\partial_\lambda(p_\lambda v_\lambda)=0
$$

및 vanishing boundary flux 아래에서

$$
\boxed{
\dot m_3(t)
=
3\int[\lambda-\mu(t)]^2[v_\lambda(\lambda,t)-\dot\mu(t)]p_\lambda(\lambda,t)\,d\lambda
}
$$

이고 동등하게

$$
\boxed{
\dot m_3
=
3\operatorname{Cov}\left((\lambda-\mu)^2,v_\lambda\right)
}
$$

이다.

이는 stated continuity equation과 boundary condition 아래 **EXACT / IDENTITY**다.

## 8. 다음 결정적 시험

closure를 확장하기 전에 동일한 1D layer-LJ protocol에서

$$
\boxed{M\uparrow\quad\Longrightarrow\quad D_{\rm KS}(M)\ ?}
$$

를 본다.

- $D_{\rm KS}\to0$이면 현재 mismatch는 주로 finite-size/self-averaging error이고 large-$M$ closure는 여전히 가능성이 있다.
- $D_{\rm KS}$가 유한하게 남으면 equal-measure closure가 dynamical information을 놓치고 있는 것이다.
- mismatch가 $m_3(t)$와 함께 움직이면 mechanics-derived third-moment state가 강한 다음 후보가 된다.
