# 1D Statistical-Cell Correlation Result

## Classification

**NUMERICAL RESULT / DIAGNOSTIC** for the dynamically matched deterministic 1D layer-LJ protocol. The reported statistical lengths use the finite-snapshot first-positive-lobe estimator from `docs/MILESTONE13A_FINITE_SNAPSHOT_CORRECTION.md`; they are not yet material constants.

## Protocol

The represented spacing counts are

$$
M=31,63,127,255,
$$

with dynamic similarity

$$
\omega M=0.62.
$$

The open-chain correlation profile is evaluated at the same phase-locked cycle used in the existing spatial-correlation sweep.

For a finite snapshot,

$$
\widehat\tau_M^{(+)}
=
1+2\sum_{k=1}^{K_0}
\left(1-\frac{k}{M}\right)\widehat\rho_k,
$$

where $K_0$ is the last positive lag before the first non-positive empirical correlation.

Then

$$
\widehat M_{\rm eff}^{(+)}
=\frac{M}{\widehat\tau_M^{(+)}},
$$

and in reduced length units

$$
\frac{\widehat\ell_{\rm stat}^{(2,+)}}{a_0}
=\widehat\tau_M^{(+)}.
$$

## Numerical result

| $M$ | $\widehat\tau_M^{(+)}$ | $\widehat M_{\rm eff}^{(+)}$ | $\widehat\ell_{\rm stat}^{(2,+)}/a_0$ | $\widehat M_{\rm eff}^{(+)}/M$ |
|---:|---:|---:|---:|---:|
| 31 | 10.5836 | 2.9291 | 10.5836 | 0.09449 |
| 63 | 21.1011 | 2.9856 | 21.1011 | 0.04739 |
| 127 | 41.9340 | 3.0286 | 41.9340 | 0.02385 |
| 255 | 83.4891 | 3.0543 | 83.4891 | 0.01198 |

The first zero crossing simultaneously remains near

$$
\frac{k_0}{M}
\approx0.35.
$$

## Interpretation

The effective independent count remains approximately

$$
\boxed{\widehat M_{\rm eff}^{(+)}\approx3}
$$

over an eightfold increase in represented system size.

Meanwhile

$$
\widehat\ell_{\rm stat}^{(2,+)}
$$

grows approximately in proportion to $M$. Numerically,

$$
\frac{\widehat\ell_{\rm stat}^{(2,+)}}{Ma_0}
\approx
0.341,\;0.335,\;0.330,\;0.327.
$$

Therefore the tested protocol does **not** show convergence to a local, system-size-independent material correlation length.

The stronger interpretation is that the current dynamically matched deterministic chain retains system-scale coherent structure. In a variance-equivalent sense it behaves roughly like only three independent axial probability blocks, even as the represented chain is enlarged.

This is not evidence that each block is internally exactly identical. Complete identical dependence would require

$$
\mathbb E[(X-Y)^2]=0.
$$

The present result establishes strong partial dependence and a small variance-equivalent independent count, not exact block identity.

## Consequence for later probability aggregation

A numerical mesh cell must not be equated with an independent probability cell.

For variance-based aggregation, the current protocol gives a system-scale effective cell length rather than a converged local one.

For crack-tail or first-passage aggregation, a separate event-clustering length is still required. The second-moment length above cannot by itself justify

$$
1-(1-q)^N.
$$

Full event independence requires joint factorization, while an exact identical-block model with block size $b$ gives

$$
1-(1-q)^{M/b}.
$$

## Current decision

Do not assign a fixed local 1D statistical mini-cell length from the present sweep. First determine whether the system-scale coherence is caused by the chosen boundary/loading protocol or survives under a physically appropriate 1D state ensemble. Continue the active research in one dimension.

---

# 한국어 번역 — 1D 통계셀 상관 결과

## 분류

동적으로 matched된 deterministic 1D layer-LJ protocol에 대한 **NUMERICAL RESULT / DIAGNOSTIC**이다. 아래 통계 특성길이는 `docs/MILESTONE13A_FINITE_SNAPSHOT_CORRECTION.md`의 finite-snapshot first-positive-lobe estimator를 사용하며 아직 material constant가 아니다.

## Protocol

represented spacing 수는

$$
M=31,63,127,255
$$

이고 dynamic similarity는

$$
\omega M=0.62
$$

이다.

기존 spatial-correlation sweep과 같은 phase-locked cycle에서 open-chain correlation profile을 계산한다.

finite snapshot에 대해

$$
\widehat\tau_M^{(+)}
=
1+2\sum_{k=1}^{K_0}
\left(1-\frac{k}{M}\right)\widehat\rho_k
$$

를 사용한다. 여기서 $K_0$는 첫 non-positive empirical correlation 직전의 마지막 positive lag다.

그리고

$$
\widehat M_{\rm eff}^{(+)}
=\frac{M}{\widehat\tau_M^{(+)}}
$$

이며 reduced length 단위에서는

$$
\frac{\widehat\ell_{\rm stat}^{(2,+)}}{a_0}
=\widehat\tau_M^{(+)}
$$

이다.

## 수치결과

| $M$ | $\widehat\tau_M^{(+)}$ | $\widehat M_{\rm eff}^{(+)}$ | $\widehat\ell_{\rm stat}^{(2,+)}/a_0$ | $\widehat M_{\rm eff}^{(+)}/M$ |
|---:|---:|---:|---:|---:|
| 31 | 10.5836 | 2.9291 | 10.5836 | 0.09449 |
| 63 | 21.1011 | 2.9856 | 21.1011 | 0.04739 |
| 127 | 41.9340 | 3.0286 | 41.9340 | 0.02385 |
| 255 | 83.4891 | 3.0543 | 83.4891 | 0.01198 |

first zero crossing도 동시에

$$
\frac{k_0}{M}\approx0.35
$$

근처에 유지된다.

## 해석

represented system size를 8배 키워도 유효 독립개수는

$$
\boxed{
\widehat M_{\rm eff}^{(+)}\approx3
}
$$

정도로 거의 유지된다.

반면

$$
\widehat\ell_{\rm stat}^{(2,+)}
$$

는 거의 $M$에 비례해서 증가한다. 실제로

$$
\frac{\widehat\ell_{\rm stat}^{(2,+)}}{Ma_0}
\approx
0.341,\;0.335,\;0.330,\;0.327
$$

이다.

따라서 tested protocol에서는 system size와 무관한 local material correlation length로 수렴하는 현상이 보이지 않는다.

더 강한 해석은 현재 dynamically matched deterministic chain에 system-scale coherent structure가 남아 있다는 것이다. variance-equivalent 관점에서는 represented chain을 키워도 대략 3개의 독립 axial probability block처럼 행동한다.

다만 이것이 각 block 내부가 정확히 같은 확률변수라는 뜻은 아니다. 완전 동일 종속은

$$
\mathbb E[(X-Y)^2]=0
$$

을 요구한다.

현재 결과는 강한 부분 종속과 작은 variance-equivalent independent count를 보인 것이지 exact block identity를 증명한 것은 아니다.

## 이후 probability aggregation에 대한 의미

numerical mesh cell을 independent probability cell과 동일시하면 안 된다.

variance-based aggregation에 대해서는 현재 protocol이 converged local length가 아니라 system-scale effective cell length를 준다.

crack-tail 또는 first-passage aggregation에는 별도의 event-clustering length가 필요하다. 위 second-moment length만으로

$$
1-(1-q)^N
$$

을 정당화할 수 없다.

full event independence는 joint factorization을 요구하고, block size가 $b$인 exact identical-block model에서는

$$
1-(1-q)^{M/b}
$$

가 된다.

## 현재 결정

현재 sweep만으로 fixed local 1D statistical mini-cell length를 지정하지 않는다. 먼저 system-scale coherence가 현재 boundary/loading protocol 때문에 생긴 것인지, 물리적으로 적절한 1D state ensemble에서도 유지되는지 확인한다. 연구범위는 계속 1D로 유지한다.
