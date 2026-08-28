# Al Fatigue Probability Theory

Mechanics-first research framework for fatigue crack initiation under **one-dimensional normal cyclic loading** in high-purity / single-crystal aluminum.

## Active scope

The active derivation is deliberately restricted to a one-dimensional stack of represented material layers. The microscopic reduced coordinate is the normal spacing

$$
a_i(t)>0,
$$

or the normalized spacing

$$
\lambda_i(t)=a_i(t)/a_0.
$$

The effective normal interaction between layers is represented by the calibrated generalized Lennard-Jones model. Three-dimensional FCC and shear work remain archived under `libraries/` and are not part of the active derivation.

Physical time $t$ is fundamental. Fatigue cycle count is not an independent state variable.

## Calibrated 1D layer-LJ mechanics

The normalized layer energy is

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

The calibration gives

$$
\phi'(1)=0,
\qquad
\phi''(1)=1.
$$

The interior spacing equation is

$$
\boxed{
\ddot\lambda_i
=
\phi'(\lambda_{i+1})
-2\phi'(\lambda_i)
+\phi'(\lambda_{i-1}).
}
$$

The current idealized tangent-instability stretch is

$$
\boxed{
\lambda_c
=
\left(\frac{m+1}{n+1}\right)^{1/(m-n)}
\approx1.1077715386.
}
$$

## One-point probability state

The normalized one-point spacing density is

$$
\int_0^\infty p_\lambda(\lambda,t)\,d\lambda=1.
$$

Its mean and shifted configurational energy are

$$
\mu(t)=\int_0^\infty\lambda p_\lambda(\lambda,t)\,d\lambda,
$$

$$
\psi(\lambda)=\phi(\lambda)-\phi(1),
$$

and

$$
\mathcal E(t)=\int_0^\infty\psi(\lambda)p_\lambda(\lambda,t)\,d\lambda.
$$

The exact energy-feasibility work shows that normalization, mean, and energy alone cannot force a tensile tail because arbitrarily large energy can mathematically be hidden in the LJ compression branch as $\lambda\to0^+$. Any exact safe-energy ceiling therefore requires an independently justified compression-side constraint.

## Earlier two-moment distribution closure

A fixed-length/fixed-configurational-energy equal-base-measure assumption plus a large-$M$ saddle-point reduction gave the controlled approximation

$$
\boxed{
p_\lambda(\lambda,t)
=Z^{-1}\exp[-\alpha(t)\lambda-\beta(t)\psi(\lambda)].
}
$$

The multipliers are determined by $\mu(t)$ and $\mathcal E(t)$ rather than fitted to histograms.

Direct deterministic tests showed that this closure can reproduce near-equilibrium variance reasonably well but does not reproduce the full driven distribution. In the tested 32-node states the Kolmogorov distance is about $0.15$, and the slower case shows a strong skewness mismatch. Therefore $\mu$ and $\mathcal E$ are not sufficient to determine the driven one-point distribution.

## Spatial-correlation result

The deterministic chain retains strong spatial ordering. Define

$$
C_k(t)
=
\frac{1}{M-k}
\sum_{i=1}^{M-k}
[\lambda_i-\mu][\lambda_{i+k}-\mu],
$$

and

$$
\rho_k=C_k/C_0.
$$

In the dynamically matched sweep with $\omega M=0.62$, the nearest-neighbor correlation rises from about $0.933$ at $M=31$ to about $0.991$ at $M=255$. The profile approximately collapses when plotted against $k/M$, with the first zero crossing near $0.35M$.

A one-point density is exactly invariant under permutation of the layer labels while $C_k$ is not. Therefore

$$
\boxed{p_\lambda(\lambda,t)\text{ cannot encode the complete spatial mechanical state}.}
$$

## New governing-equation clue for the form of $P$

The one-point density is also the exact spatial push-forward of the deterministic spacing field. With continuum layer label $\xi\in[0,1]$ and spacing field $\Lambda(\xi,t)$,

$$
\boxed{
p_\lambda(\lambda,t)
=
\int_0^1\delta[\lambda-\Lambda(\xi,t)]\,d\xi.
}
$$

For a piecewise monotone field,

$$
\boxed{
p_\lambda(\lambda,t)
=
\sum_{\xi_j:\Lambda(\xi_j,t)=\lambda}
\frac{1}{|\partial_\xi\Lambda(\xi_j,t)|}.
}
$$

This provides a direct clue for the functional form of $P$: its shape is the kinematic image of the mechanically generated spatial waveform rather than an independently selected probability family.

Linearization around $\lambda=1$ gives

$$
\ddot u_i=u_{i+1}-2u_i+u_{i-1},
$$

with dispersion

$$
\boxed{\omega_q^2=4\sin^2(q/2).}
$$

A single coherent linear mode therefore generates the exact spatial push-forward

$$
\boxed{
p_{\rm 1mode}(\lambda)
=
\frac{\mathbf 1_{|\lambda-\mu|<|A|}}
{\pi\sqrt{A^2-(\lambda-\mu)^2}}.
}
$$

This arcsine form is not a fitted distribution; it follows from $\Lambda=\mu+A\cos\vartheta$ with uniform spatial phase.

The calibrated LJ force is strongly nonlinear:

$$
\boxed{
\phi'(1+u)
=u-10.595u^2+62.97935u^3+O(u^4).
}
$$

Thus the base mechanics itself generates harmonic distortion, skewness, and potentially multimodal one-point densities.

A first-plus-second-harmonic waveform

$$
\Lambda=\mu+A\cos\vartheta+B\cos2\vartheta
$$

has

$$
\operatorname{Var}(\lambda)=\frac{A^2+B^2}{2},
$$

$$
\mu_3=\frac34A^2B,
$$

and the exact bound

$$
\boxed{|\gamma_1|\le\sqrt{2/3}\approx0.816497.}
$$

The slower deterministic snapshot has $\gamma_1\approx1.062$, so even two spatial harmonics are insufficient for that state. A single-mode arcsine reference also gives larger KS error than the earlier exponential closure. These are negative results that point toward richer mechanically generated mode content, not toward fitting another named probability family.

## Current research direction

The current hierarchy is

$$
\boxed{
\text{1D layer-LJ governing equation}
\rightarrow
\Lambda(\xi,t)\text{ / neighboring-spacing pair state}
\rightarrow
p_\lambda(\lambda,t)\text{ by push-forward}
\rightarrow
Q_c(t)
\rightarrow
\text{normal-opening first passage}.
}
$$

The next target is to derive the mechanically excited spatial-mode or pair-state evolution directly from the 1D governing equation and determine the smallest representation that reproduces both the one-point distribution $P_1$ and the observed $C_k$ without empirical relaxation laws or fitted probability families.

## Active files

- `theory/normal_lj_chain.py` — conservative 1D layer-LJ chain
- `theory/normal_lj_energy_feasibility.py` — exact energy-feasibility bounds under stated compression constraints
- `theory/normal_lj_distribution.py` — earlier two-moment large-$M$ closure
- `theory/normal_lj_spatial_correlation.py` — exact finite-chain correlation diagnostics
- `theory/normal_lj_pushforward.py` — spatial push-forward and mode-derived distribution identities
- `simulations/run_normal_lj_spatial_correlation.py` — dynamically matched correlation sweep
- `simulations/run_normal_lj_pushforward_clue.py` — single-mode falsification diagnostic
- `docs/MILESTONE8_SPATIAL_CORRELATION.md` — spatial-ordering result
- `docs/MILESTONE9_GOVERNING_EQUATION_PUSHFORWARD.md` — governing-equation clue for the form of $P$
- `results/data/result_manifest.json` — current machine-readable research state

## Reproduce active results

```bash
python -m pip install -r requirements.txt
python -m simulations.generate_results
python -m unittest discover -s tests
```

## Research labels

Important statements are classified as one of:

- **EXACT / IDENTITY**
- **DEFINITION**
- **ASSUMPTION**
- **CONTROLLED APPROXIMATION**
- **EMPIRICAL INPUT / PHYSICAL CONSTRAINT**
- **NUMERICAL RESULT / DIAGNOSTIC**

No fitted Gaussian/Weibull distribution, cycle-dependent LJ parameter, empirical damage variable, or fitted correlation length is accepted as a mechanics derivation.

---

# 한국어 번역

고순도 또는 단결정 알루미늄의 **1차원 수직 반복하중** 아래 피로 균열개시를 mechanics-first 방식으로 전개하는 연구 저장소다.

## 활성 범위

활성 derivation은 represented material layer의 1차원 stack으로 의도적으로 제한한다. microscopic reduced coordinate는 수직 layer spacing

$$
a_i(t)>0
$$

또는 normalized spacing

$$
\lambda_i(t)=a_i(t)/a_0
$$

이다.

layer 사이의 유효 normal interaction은 calibration된 generalized Lennard-Jones model로 표현한다. 3차원 FCC와 shear 연구는 `libraries/` 아래 archive로 유지하며 active derivation에 포함하지 않는다.

물리적 시간 $t$가 근본 evolution coordinate이며 fatigue cycle count는 독립 상태변수가 아니다.

## Calibration된 1D layer-LJ mechanics

normalized layer energy는

$$
\phi(\lambda)
=
\frac{\lambda^{-m}}{m(m-n)}
-
\frac{\lambda^{-n}}{n(m-n)},
$$

이고

$$
m=12.19,
\qquad
n=6
$$

이다.

calibration은

$$
\phi'(1)=0,
\qquad
\phi''(1)=1
$$

을 준다.

interior spacing 지배방정식은

$$
\boxed{
\ddot\lambda_i
=
\phi'(\lambda_{i+1})
-2\phi'(\lambda_i)
+\phi'(\lambda_{i-1})
}
$$

이다.

현재 이상화된 tangent-instability stretch는

$$
\boxed{
\lambda_c
=
\left(\frac{m+1}{n+1}\right)^{1/(m-n)}
\approx1.1077715386
}
$$

이다.

## One-point probability state

normalized one-point spacing density는

$$
\int_0^\infty p_\lambda(\lambda,t)\,d\lambda=1
$$

을 만족한다.

평균과 shifted configurational energy는

$$
\mu(t)=\int_0^\infty\lambda p_\lambda(\lambda,t)\,d\lambda,
$$

$$
\psi(\lambda)=\phi(\lambda)-\phi(1),
$$

그리고

$$
\mathcal E(t)=\int_0^\infty\psi(\lambda)p_\lambda(\lambda,t)\,d\lambda
$$

이다.

exact energy-feasibility 연구에서 normalization, mean, energy만으로는 tensile tail을 강제할 수 없다는 것이 나왔다. LJ compression branch가 $\lambda\to0^+$에서 수학적으로 임의로 큰 energy를 담을 수 있기 때문이다. 따라서 exact safe-energy ceiling을 쓰려면 독립적으로 정당화된 compression-side constraint가 필요하다.

## 이전 two-moment distribution closure

fixed-length/fixed-configurational-energy equal-base-measure assumption과 large-$M$ saddle-point reduction으로

$$
\boxed{
p_\lambda(\lambda,t)
=Z^{-1}\exp[-\alpha(t)\lambda-\beta(t)\psi(\lambda)]
}
$$

이라는 controlled approximation을 얻었다.

$\alpha$, $\beta$는 histogram fitting이 아니라 $\mu(t)$와 $\mathcal E(t)$로 결정된다.

하지만 deterministic 1D layer-LJ 직접시험에서 이 closure는 near-equilibrium variance는 상당히 잘 맞춰도 driven distribution 전체 shape를 재현하지 못했다. tested 32-node state에서 Kolmogorov distance는 약 $0.15$이고, 느린 case에서는 skewness mismatch가 크게 남았다. 따라서 $\mu$와 $\mathcal E$만으로 driven one-point distribution을 정할 수는 없다.

## Spatial-correlation 결과

deterministic chain은 강한 spatial ordering을 유지한다. 다음을 정의한다.

$$
C_k(t)
=
\frac{1}{M-k}
\sum_{i=1}^{M-k}
[\lambda_i-\mu][\lambda_{i+k}-\mu],
$$

$$
\rho_k=C_k/C_0.
$$

$\omega M=0.62$로 dynamic similarity를 맞춘 sweep에서 nearest-neighbor correlation은 $M=31$의 약 $0.933$에서 $M=255$의 약 $0.991$까지 증가했다. $k/M$에 대해 그리면 profile이 거의 collapse하며 첫 zero crossing은 약 $0.35M$에 있다.

one-point density는 layer label을 permutation해도 정확히 불변이지만 $C_k$는 그렇지 않다. 따라서

$$
\boxed{p_\lambda(\lambda,t)\text{만으로 complete spatial mechanical state를 담을 수 없다}}
$$

는 구조적 결과를 얻었다.

## $P$의 형식에 대한 새로운 지배방정식 단서

one-point density는 deterministic spacing field의 정확한 spatial push-forward이기도 하다. continuum layer label $\xi\in[0,1]$과 spacing field $\Lambda(\xi,t)$를 쓰면

$$
\boxed{
p_\lambda(\lambda,t)
=
\int_0^1\delta[\lambda-\Lambda(\xi,t)]\,d\xi
}
$$

이다.

field가 piecewise monotone이면

$$
\boxed{
p_\lambda(\lambda,t)
=
\sum_{\xi_j:\Lambda(\xi_j,t)=\lambda}
\frac{1}{|\partial_\xi\Lambda(\xi_j,t)|}
}
$$

이다.

즉 $P$의 함수형식은 독립적으로 고른 probability family가 아니라 mechanically generated spatial waveform의 kinematic image라는 직접적인 단서를 얻는다.

$\lambda=1$ 주변에서 선형화하면

$$
\ddot u_i=u_{i+1}-2u_i+u_{i-1}
$$

이고 dispersion은

$$
\boxed{\omega_q^2=4\sin^2(q/2)}
$$

이다.

하나의 coherent linear mode는 따라서 정확히

$$
\boxed{
p_{\rm 1mode}(\lambda)
=
\frac{\mathbf 1_{|\lambda-\mu|<|A|}}
{\pi\sqrt{A^2-(\lambda-\mu)^2}}
}
$$

이라는 spatial push-forward를 만든다.

이 arcsine form은 fitted distribution이 아니라 $\Lambda=\mu+A\cos\vartheta$와 uniform spatial phase에서 직접 나온다.

calibration된 LJ force는 강하게 nonlinear하다.

$$
\boxed{
\phi'(1+u)
=u-10.595u^2+62.97935u^3+O(u^4)
}
$$

이다.

따라서 base mechanics 자체가 harmonic distortion, skewness, 그리고 잠재적인 multimodal one-point density를 만든다.

first+second harmonic waveform

$$
\Lambda=\mu+A\cos\vartheta+B\cos2\vartheta
$$

에서는

$$
\operatorname{Var}(\lambda)=\frac{A^2+B^2}{2},
$$

$$
\mu_3=\frac34A^2B,
$$

그리고 정확히

$$
\boxed{|\gamma_1|\le\sqrt{2/3}\approx0.816497}
$$

이다.

느린 deterministic snapshot의 $\gamma_1\approx1.062$는 이 bound도 넘기 때문에 two spatial harmonic만으로도 충분하지 않다. single-mode arcsine reference 역시 이전 exponential closure보다 KS error가 더 컸다. 이는 실패결과이며, 새로운 named probability family를 fitting해야 한다는 뜻이 아니라 mechanically generated mode content가 더 풍부하다는 뜻이다.

## 현재 연구방향

현재 hierarchy는

$$
\boxed{
\text{1D layer-LJ governing equation}
\rightarrow
\Lambda(\xi,t)\text{ / neighboring-spacing pair state}
\rightarrow
p_\lambda(\lambda,t)\text{ by push-forward}
\rightarrow
Q_c(t)
\rightarrow
\text{normal-opening first passage}
}
$$

이다.

다음 목표는 1D 지배방정식에서 mechanically excited spatial-mode 또는 pair-state evolution을 직접 유도하고 empirical relaxation law나 fitted probability family 없이 one-point distribution $P_1$과 관찰된 $C_k$를 동시에 재현하는 최소 representation을 찾는 것이다.

## 활성 파일

- `theory/normal_lj_chain.py` — conservative 1D layer-LJ chain
- `theory/normal_lj_energy_feasibility.py` — stated compression constraint 아래 exact energy-feasibility bound
- `theory/normal_lj_distribution.py` — 이전 two-moment large-$M$ closure
- `theory/normal_lj_spatial_correlation.py` — exact finite-chain correlation diagnostic
- `theory/normal_lj_pushforward.py` — spatial push-forward 및 mode-derived distribution identity
- `simulations/run_normal_lj_spatial_correlation.py` — dynamically matched correlation sweep
- `simulations/run_normal_lj_pushforward_clue.py` — single-mode falsification diagnostic
- `docs/MILESTONE8_SPATIAL_CORRELATION.md` — spatial-ordering result
- `docs/MILESTONE9_GOVERNING_EQUATION_PUSHFORWARD.md` — $P$ 형식에 대한 지배방정식 단서
- `results/data/result_manifest.json` — 현재 machine-readable research state

## 활성 결과 재현

```bash
python -m pip install -r requirements.txt
python -m simulations.generate_results
python -m unittest discover -s tests
```

## 연구 분류 라벨

중요한 statement는 다음 중 하나로 분류한다.

- **EXACT / IDENTITY**
- **DEFINITION**
- **ASSUMPTION**
- **CONTROLLED APPROXIMATION**
- **EMPIRICAL INPUT / PHYSICAL CONSTRAINT**
- **NUMERICAL RESULT / DIAGNOSTIC**

fitted Gaussian/Weibull distribution, cycle-dependent LJ parameter, empirical damage variable, fitted correlation length은 mechanics derivation으로 인정하지 않는다.
