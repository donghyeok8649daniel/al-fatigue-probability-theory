# Failed / Rejected Approaches

Keep failed approaches here instead of silently deleting them. A failed model can still be useful as a null test.

## Reversible one-coordinate LJ model

Model:

$$
\sigma(t)\leftrightarrow a(t),\qquad U(a)\text{ single-valued and conservative.}
$$

Result:

$$
\oint \sigma\,d\epsilon=0.
$$

Reason for rejection as a fatigue model: no internal irreversible state evolution. It remains useful as a reversible baseline/unit test.

## Prescribed Weibull-like spacing density

Earlier work prescribed a Weibull-type density, later multiplied by an oscillatory factor. This is not accepted as a foundational evolution law because the distribution was imposed rather than derived from mechanics.

## Arbitrary stochastic kernel / Kramers rates

Transition kernels, barrier-crossing rates, or damping constants must not be introduced solely to create hysteresis or fatigue accumulation. Such models may be used only after a derivation or controlled coarse-graining argument.

## Instantaneous tail probability as crack probability

$$
Q_c(t)=\int_{a_c}^{\infty}P(a,t)\,da
$$

is an instantaneous unstable fraction, not automatically a cumulative crack-initiation probability. A first-passage or absorbing-boundary formulation is required for initiation probability.

---

# 한국어 번역 — 실패했거나 배제한 접근

실패한 접근을 조용히 삭제하지 말고 이 파일에 남긴다. 실패한 모델도 null test 또는 반증 기준으로 유용할 수 있다.

## 가역적인 단일좌표 LJ 모델

모델은

$$
\sigma(t)\leftrightarrow a(t),\qquad U(a)\text{가 단일값이며 보존적}
$$

인 형태다.

결과는

$$
\oint \sigma\,d\epsilon=0
$$

이다.

피로모델로 배제하는 이유는 내부의 비가역적 구조진화가 없기 때문이다. 다만 완전히 가역적인 baseline 및 unit test로는 계속 유용하다.

## 미리 정해 놓은 Weibull 형태의 spacing density

초기 연구에서는 Weibull 형태의 density를 먼저 가정하고 이후 oscillatory factor를 곱했다. 이 분포는 역학으로부터 유도된 것이 아니라 외부에서 부여된 것이므로 이론의 기초 evolution law로 인정하지 않는다.

## 임의의 stochastic kernel / Kramers rate

히스테리시스나 피로누적을 만들어내기 위해 transition kernel, barrier-crossing rate, damping constant를 임의로 넣어서는 안 된다. 이런 모델은 미시역학으로부터의 유도 또는 controlled coarse-graining 논증이 먼저 존재할 때에만 사용할 수 있다.

## 순간 tail probability를 crack probability로 해석하는 접근

$$
Q_c(t)=\int_{a_c}^{\infty}P(a,t)\,da
$$

는 그 순간 unstable 영역에 존재하는 분율이다. 이것이 자동으로 누적 crack-initiation probability가 되는 것은 아니다. 균열개시 확률을 정의하려면 first-passage 또는 absorbing-boundary 정식화가 필요하다.
