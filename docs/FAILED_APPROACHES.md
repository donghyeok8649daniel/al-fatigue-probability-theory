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

## Retaining hysteresis work as pair-potential energy

Adding $H_k$ or an arbitrary fraction of it to $U(a)$ is rejected. $U$ is a
state function, whereas loop work is path dependent. In the isothermal
Smoluchowski reduction, current dissipation is transferred to eliminated bath
coordinates and the surviving conditional state becomes periodic. The
irreversible active observable is escaped probability. Persistent stored
energy would require a separately derived slow microscopic state; it cannot be
created by changing the same LJ potential after each cycle.

## Finite-harmonic form as a global spacing-distribution model

A single or finite set of spatial harmonics was explored as a push-forward diagnostic. This route is not accepted as the active global form of $P(\lambda,t)$ because a finite harmonic ansatz inserts spatial periodicity that is not guaranteed by the actual driven boundary-value problem.

The harmonic calculation remains useful only as a historical/local falsification diagnostic. It must not be promoted to the global nonequilibrium distribution law.

## Taylor expansion about equilibrium as a full-support distribution derivation

Expanding the generalized-LJ force about $\lambda=1$ is a local approximation. The fatigue theory must reason about the entire occupied support, including tensile and compression tails, so a local Taylor expansion is not used to derive the active full distribution.

The original nonlinear generalized-LJ force is retained in the active transport/hierarchy derivation. Local Taylor coefficients may be kept only as explicitly local diagnostics, not as the governing law for the full $P(\lambda,t)$.

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

## Hysteresis work를 pair-potential energy로 저장하는 접근

$H_k$ 또는 그 임의 분율을 $U(a)$에 더하는 방식은 배제한다. $U$는
상태함수이고 loop work는 경로의존량이다. isothermal Smoluchowski 축약에서
current dissipation은 생략한 bath 좌표로 전달되고, 생존조건부 상태는
주기상태가 된다. 현재 비가역량은 유출확률이다. 지속적인 저장에너지를
표현하려면 별도로 유도한 느린 미시상태가 필요하며 같은 LJ potential을
cycle마다 바꿔서 만들 수 없다.

## 유한 harmonic 형식을 전역 spacing distribution 모델로 사용하는 접근

single 또는 finite spatial harmonic을 push-forward diagnostic으로 시험했지만 이를 $P(\lambda,t)$의 활성 전역 형식으로 사용하지 않는다. 유한 harmonic ansatz는 실제 driven boundary-value problem에서 보장되지 않은 공간 주기성을 미리 넣기 때문이다.

harmonic 계산은 과거의 국소 진단/반증시험으로만 보존하며 global nonequilibrium distribution law로 승격하지 않는다.

## 평형점 Taylor 전개로 전체 support 분포를 유도하는 접근

generalized-LJ force를 $\lambda=1$ 주변에서 Taylor 전개하는 것은 국소 근사다. 현재 피로이론은 tensile/compression tail을 포함한 전체 occupied support를 다뤄야 하므로 local Taylor expansion으로 활성 전체 distribution을 유도하지 않는다.

활성 transport/hierarchy derivation에서는 원래의 nonlinear generalized-LJ force를 그대로 유지한다. Taylor coefficient는 필요하면 명시적인 국소 diagnostic으로만 보존하고 full $P(\lambda,t)$의 governing law로 사용하지 않는다.
