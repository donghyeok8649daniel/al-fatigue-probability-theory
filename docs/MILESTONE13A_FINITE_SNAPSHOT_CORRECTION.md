# Milestone 13A — Finite-Snapshot Correction for the 1D Statistical Length

## Status

This note corrects a subtle but important point in Milestone 13.

The exact finite-$M$ variance formula

$$
\tau_M
=
1+2\sum_{k=1}^{M-1}
\left(1-\frac{k}{M}\right)\rho_k
$$

belongs to the **true correlation sequence of a second-order stationary stochastic process**.

It must not be evaluated by naively inserting every lag of the open-chain correlation diagnostic computed from one deterministic snapshot centered by its own sample mean.

## 1. Existing open-chain snapshot correlation

For one finite snapshot define

$$
d_i=\lambda_i-\bar\lambda,
\qquad
\sum_{i=1}^M d_i=0.
$$

The existing spatial-correlation module uses

$$
C_0
=\frac1M\sum_{i=1}^M d_i^2
$$

and for $k>0$

$$
C_k
=\frac1{M-k}
\sum_{i=1}^{M-k}d_i d_{i+k},
$$

with

$$
\widehat\rho_k=\frac{C_k}{C_0}.
$$

This is a useful spatial-ordering diagnostic. Because the numerator and denominator use different finite-sample supports, $\widehat\rho_k$ is not required to behave exactly like a population correlation coefficient at every large lag.

## 2. Exact all-lag zero-sum identity for a sample-mean-centered snapshot

Multiply the open-chain covariance by its finite-length weight:

$$
\left(1-\frac{k}{M}\right)C_k
=\frac1M\sum_{i=1}^{M-k}d_i d_{i+k}.
$$

Then

$$
C_0
+2\sum_{k=1}^{M-1}
\left(1-\frac{k}{M}\right)C_k
$$

contains every diagonal product $d_i^2$ once and every off-diagonal product $d_i d_j$ twice. Therefore

$$
\boxed{
C_0
+2\sum_{k=1}^{M-1}
\left(1-\frac{k}{M}\right)C_k
=
\frac1M\left(\sum_{i=1}^M d_i\right)^2
=0.
}
$$

If $C_0>0$, division by $C_0$ gives

$$
\boxed{
1
+2\sum_{k=1}^{M-1}
\left(1-\frac{k}{M}\right)\widehat\rho_k
=0.
}
$$

This is an **EXACT FINITE-SNAPSHOT IDENTITY** caused by centering the same finite sample by its own mean.

Hence an all-lag plug-in estimate would return zero by construction and cannot estimate the population correlation factor.

## 3. Positive-window estimator

For a single deterministic snapshot, define $K_0$ as the last strictly positive lag before the first non-positive empirical correlation. Then use

$$
\boxed{
\widehat\tau_M^{(+)}
=
1+2\sum_{k=1}^{K_0}
\left(1-\frac{k}{M}\right)
\widehat\rho_k.
}
$$

The corresponding diagnostics are

$$
\boxed{
\widehat M_{\rm eff}^{(+)}
=\frac{M}{\widehat\tau_M^{(+)}},
}
$$

and

$$
\boxed{
\widehat\ell_{\rm stat}^{(2,+)}
=a_0\widehat\tau_M^{(+)}.
}
$$

These are **ESTIMATOR / DIAGNOSTIC** quantities, not exact material constants.

The first-positive-lobe rule is used because it removes the forced long-lag cancellation from the finite centered sample without introducing a fitted correlation length. Its physical credibility must still be checked by represented-system-size convergence.

## 4. Interpretation rule

The hierarchy is now:

- true stationary $\rho_k$ available: use the exact $\tau_M$ identity;
- one deterministic finite snapshot: use the separately labeled positive-window estimator;
- no convergence with increasing represented size: do not call the result a material correlation length.

If

$$
\widehat\ell_{\rm stat}^{(2,+)}(M)
\propto M,
$$

then the current protocol is dominated by system-scale coherent motion. An arbitrary cutoff must not be imposed to manufacture a local material cell.

## 5. Relation to event dependence

Even a converged second-moment statistical length does not prove independence of crack-initiation events. Complete identity, partial dependence, and full factorization remain distinct.

For independent blocks of size $b$ whose variables are exactly identical within each block,

$$
P(\text{any event})
=
1-(1-q)^{M/b}.
$$

A future crack-tail length must therefore be checked using exceedance/first-passage clustering rather than identified automatically with $\ell_{\rm stat}^{(2)}$.

---

# 한국어 번역 — 1D 통계 특성길이의 finite-snapshot 정정

## 상태

이 문서는 Milestone 13의 미묘하지만 중요한 점을 정정한다.

정확한 finite-$M$ 평균분산 식에 들어가는

$$
\tau_M
=
1+2\sum_{k=1}^{M-1}
\left(1-\frac{k}{M}\right)\rho_k
$$

는 **second-order stationary stochastic process의 true correlation sequence**에 대한 식이다.

하나의 deterministic finite snapshot을 자기 자신의 sample mean으로 center해서 계산한 open-chain correlation의 모든 lag를 이 식에 그대로 넣으면 안 된다.

## 1. 현재 open-chain snapshot correlation

finite snapshot 하나에 대해

$$
d_i=\lambda_i-\bar\lambda,
\qquad
\sum_{i=1}^M d_i=0
$$

라고 하자.

기존 spatial-correlation module은

$$
C_0
=\frac1M\sum_{i=1}^M d_i^2
$$

및 $k>0$에 대해

$$
C_k
=\frac1{M-k}
\sum_{i=1}^{M-k}d_i d_{i+k}
$$

를 사용하고

$$
\widehat\rho_k=\frac{C_k}{C_0}
$$

를 정의한다.

이 양은 spatial ordering을 보는 데 유용한 diagnostic이다. 하지만 numerator와 denominator의 finite-sample support가 다르므로 큰 lag에서 population correlation coefficient와 완전히 같은 성질을 강제할 수는 없다.

## 2. sample-mean-centered finite snapshot의 정확한 all-lag zero-sum identity

open-chain covariance에 finite-length weight를 곱하면

$$
\left(1-\frac{k}{M}\right)C_k
=\frac1M\sum_{i=1}^{M-k}d_i d_{i+k}
$$

이다.

따라서

$$
C_0
+2\sum_{k=1}^{M-1}
\left(1-\frac{k}{M}\right)C_k
$$

에는 모든 diagonal product $d_i^2$가 한 번, 모든 off-diagonal product $d_i d_j$가 두 번 들어간다.

그러므로

$$
\boxed{
C_0
+2\sum_{k=1}^{M-1}
\left(1-\frac{k}{M}\right)C_k
=
\frac1M\left(\sum_{i=1}^M d_i\right)^2
=0
}
$$

이다.

$C_0>0$이면

$$
\boxed{
1
+2\sum_{k=1}^{M-1}
\left(1-\frac{k}{M}\right)\widehat\rho_k
=0
}
$$

이다.

이것은 자기 sample mean으로 centering한 데서 생기는 **EXACT FINITE-SNAPSHOT IDENTITY**다.

따라서 모든 lag를 그대로 넣는 plug-in estimator는 구조적으로 0이 되어 population correlation factor를 추정할 수 없다.

## 3. positive-window estimator

하나의 deterministic snapshot에서는 첫 non-positive empirical correlation 직전의 마지막 positive lag를 $K_0$라고 정의한다.

그때

$$
\boxed{
\widehat\tau_M^{(+)}
=
1+2\sum_{k=1}^{K_0}
\left(1-\frac{k}{M}\right)
\widehat\rho_k
}
$$

를 사용한다.

이에 대응하는 diagnostic은

$$
\boxed{
\widehat M_{\rm eff}^{(+)}
=\frac{M}{\widehat\tau_M^{(+)}}
}
$$

및

$$
\boxed{
\widehat\ell_{\rm stat}^{(2,+)}
=a_0\widehat\tau_M^{(+)}
}
$$

이다.

이들은 **ESTIMATOR / DIAGNOSTIC**이지 exact material constant가 아니다.

first-positive-lobe rule은 finite centered sample에서 강제로 생기는 long-lag cancellation을 제거하면서 fitted correlation length를 추가하지 않기 때문에 사용한다. 그래도 represented system size를 늘렸을 때 수렴하는지 확인해야 물리적 신뢰성이 생긴다.

## 4. 해석 규칙

이제 다음처럼 구분한다.

- true stationary $\rho_k$가 있으면 exact $\tau_M$ 식 사용;
- deterministic finite snapshot 하나만 있으면 별도로 표시한 positive-window estimator 사용;
- represented size를 늘려도 수렴하지 않으면 material correlation length라고 부르지 않음.

만약

$$
\widehat\ell_{\rm stat}^{(2,+)}(M)
\propto M
$$

이면 현재 protocol은 system-scale coherent motion이 지배하는 것이다. local material cell을 만들기 위해 임의 cutoff를 넣지 않는다.

## 5. event dependence와의 관계

수렴한 second-moment statistical length가 생겨도 crack-initiation event의 independence가 자동으로 증명되지는 않는다. 완전 동일, 부분 종속, full factorization은 계속 서로 다른 개념이다.

block size가 $b$이고 block 내부 변수들이 완전히 동일하며 block끼리 독립이면

$$
P(\text{any event})
=
1-(1-q)^{M/b}
$$

이다.

따라서 미래의 crack-tail 특성길이는 $\ell_{\rm stat}^{(2)}$와 자동으로 같다고 두지 않고 exceedance/first-passage clustering으로 따로 확인해야 한다.
