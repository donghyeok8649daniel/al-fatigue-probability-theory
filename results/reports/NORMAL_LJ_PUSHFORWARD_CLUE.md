# Normal-LJ Governing-Equation Push-Forward Result

## Result

The one-point normalized layer-spacing density is not an independently chosen statistical object. For a deterministic continuum spacing field $\Lambda(\xi,t)$ with uniformly weighted layer label $\xi\in[0,1]$,

$$
\boxed{
p_\lambda(\lambda,t)
=
\int_0^1\delta[\lambda-\Lambda(\xi,t)]\,d\xi.
}
$$

For a piecewise monotone field this becomes

$$
\boxed{
p_\lambda(\lambda,t)
=
\sum_{\xi_j:\Lambda(\xi_j,t)=\lambda}
\frac{1}{|\partial_\xi\Lambda(\xi_j,t)|}.
}
$$

This is an exact spatial push-forward representation and provides a direct governing-equation clue for the form of $P$.

Linearizing the calibrated 1D layer-LJ equation about $\lambda=1$ gives

$$
\ddot u_i=u_{i+1}-2u_i+u_{i-1},
$$

with

$$
\omega_q^2=4\sin^2(q/2).
$$

A single coherent mode therefore produces an arcsine spacing density rather than a Gaussian or Weibull form.

However, direct comparison with the existing deterministic snapshots shows that the single-mode law is insufficient. The Kolmogorov distances are approximately $0.260$ and $0.197$, compared with about $0.159$ and $0.146$ for the previous mean-energy exponential closure.

The LJ force expansion is strongly nonlinear:

$$
\phi'(1+u)
=u-10.595u^2+62.97935u^3+O(u^4).
$$

Hence nonlinear harmonic distortion is already present in the base governing equation. For a first-plus-second-harmonic spacing waveform,

$$
\Lambda=\mu+A\cos\vartheta+B\cos2\vartheta,
$$

the skewness obeys the exact bound

$$
|\gamma_1|\le\sqrt{2/3}\approx0.816497.
$$

The slow deterministic snapshot has $\gamma_1\approx1.062$, so even a two-harmonic waveform is insufficient for that state.

## Scientific meaning

The current evidence points away from searching for a new named probability distribution. The form of $P$ should instead be obtained by resolving the mechanically generated spatial waveform or mode content and then applying the push-forward formula.

The next task is to connect this push-forward structure with the neighboring-spacing pair hierarchy so that the same reduced description can reproduce both the one-point density and the observed spatial correlations.

This is a reduced-model structural result and numerical diagnostic, not an aluminum fatigue-life prediction.

---

# 한국어 번역 — Normal-LJ 지배방정식 Push-Forward 결과

## 결과

one-point normalized layer-spacing density는 독립적으로 선택하는 statistical object가 아니다. uniformly weighted layer label $\xi\in[0,1]$과 deterministic continuum spacing field $\Lambda(\xi,t)$가 주어지면

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

가 된다.

이는 정확한 spatial push-forward representation이며 $P$의 형식에 대해 지배방정식에서 직접 얻는 단서다.

calibrated 1D layer-LJ 방정식을 $\lambda=1$ 주변에서 선형화하면

$$
\ddot u_i=u_{i+1}-2u_i+u_{i-1}
$$

이고

$$
\omega_q^2=4\sin^2(q/2)
$$

이다.

따라서 하나의 coherent mode는 Gaussian이나 Weibull이 아니라 arcsine spacing density를 만든다.

하지만 기존 deterministic snapshot과 직접 비교하면 single-mode law는 충분하지 않다. Kolmogorov distance는 약 $0.260$, $0.197$이고 이전 mean-energy exponential closure의 약 $0.159$, $0.146$보다 크다.

LJ force expansion은 강한 nonlinear term을 가진다.

$$
\phi'(1+u)
=u-10.595u^2+62.97935u^3+O(u^4)
$$

이다.

따라서 base governing equation 자체에 nonlinear harmonic distortion이 이미 들어 있다. first+second harmonic spacing waveform

$$
\Lambda=\mu+A\cos\vartheta+B\cos2\vartheta
$$

에 대해서는 skewness가 정확히

$$
|\gamma_1|\le\sqrt{2/3}\approx0.816497
$$

을 만족한다.

slow deterministic snapshot은 $\gamma_1\approx1.062$이므로 이 상태는 two-harmonic waveform으로도 충분히 표현되지 않는다.

## 과학적 의미

현재 증거는 새로운 named probability distribution을 찾는 방향보다 mechanically generated spatial waveform 또는 mode content를 해석하고 그 결과에 push-forward formula를 적용해 $P$를 얻는 방향을 지지한다.

다음 작업은 이 push-forward structure와 neighboring-spacing pair hierarchy를 연결하여 같은 reduced description이 one-point density와 관찰된 spatial correlation을 동시에 재현하게 하는 것이다.

이는 reduced-model structural result와 numerical diagnostic이며 aluminum fatigue-life prediction은 아니다.
