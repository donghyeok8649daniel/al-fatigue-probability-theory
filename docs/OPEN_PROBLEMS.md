# Open Problems — Active Normal-Deformation Mainline

## Milestone 1 — Normal microscopic state and exact transport

The active state variable is the normal-spacing density

$$
P(a,t)=\lim_{N\to\infty}\frac1N\sum_i\delta(a-a_i(t)).
$$

The exact kinematic transport equation is

$$
\partial_tP+\partial_a(Pv)=0.
$$

The unresolved problem is closure: derive $v(a,t)$, or the required enlarged state, from normal microscopic mechanics without hiding memory inside an empirical constitutive law.

## Milestone 2 — Normal cyclic hysteresis and secular evolution

For a prescribed normal stress

$$
\sigma_n(t)=\sigma_m+\sigma_a\sin\omega t,
$$

obtain a mechanics-derived loading/unloading difference and then the stronger fatigue condition

$$
\boxed{P_{N+1}(a)\neq P_N(a)}.
$$

A closed periodic loop with no cycle-state drift is internal friction, not fatigue accumulation.

The current generalized-LJ perfect-chain 100 MPa null test correctly shows essentially reversible behavior. This result must not be tuned away.

## Milestone 3 — Time-scale bridge to laboratory fatigue frequencies

The direct atomic normal-chain dynamics has an atomic time scale. The current high-frequency resonance experiments therefore do not constitute 20 Hz fatigue predictions.

The central problem is to derive

$$
\text{fast atomic LJ dynamics}
\rightarrow
\text{projected memory / coarse state}
\rightarrow
\text{slow evolution of }P(a,t)
$$

without inserting an arbitrary relaxation time.

Candidate exact starting points include Liouville projection, conditional propagators, and memory kernels derived from eliminated microscopic variables.

## Milestone 4 — Normal-opening crack initiation

The active crack-initiation picture is normal opening or normal bond instability. For the generalized LJ baseline, a local idealized instability occurs at

$$
\phi''(\lambda_c)=0,
\qquad
\lambda_c\approx1.10777154.
$$

Instantaneous tail occupancy

$$
Q_c(t)=\int_{a_c}^{\infty}P(a,t)\,da
$$

must be distinguished from first-passage initiation.

A preferred formulation is

$$
\tau_c=\inf\{t:\text{a mechanically defined normal-opening instability occurs}\}.
$$

## Central closure questions

- Is $P(a,t)$ sufficient?
- Is the phase-space lift $F(a,c,t)$ required?
- Which neighbor-spacing correlations are essential?
- Can the exact projected memory be reduced at laboratory frequencies without empirical damping?
- What physically generates cycle-to-cycle broadening or tail growth under 100 MPa-class loading?

## Falsification tests

Any active normal model must satisfy at least:

1. zero loading gives zero artificial fatigue accumulation;
2. reversible conservative limits recover reversible behavior;
3. probability normalization is preserved;
4. density remains non-negative;
5. dimensions are consistent;
6. energy balance is satisfied;
7. the uniform-lattice limit recovers the LJ baseline;
8. $P_{N+1}\neq P_N$ is not numerical diffusion;
9. a 20 Hz claim must use a physically derived time-scale bridge;
10. LJ parameters remain fixed unless the microscopic interaction model itself is explicitly changed.

## Auxiliary shear library

Earlier Rubin-chain, non-affine slip, and gamma-surface studies are preserved under `libraries/shear/`. They remain useful methodological references but are outside the active normal-deformation mainline.

---

# 한국어 번역 — 활성 수직변형 Mainline의 미해결 문제

## 마일스톤 1 — 수직 미시상태와 정확한 수송

활성 상태변수는 수직 원자간격 밀도

$$
P(a,t)=\lim_{N\to\infty}\frac1N\sum_i\delta(a-a_i(t))
$$

이다.

정확한 운동학적 수송식은

$$
\partial_tP+\partial_a(Pv)=0
$$

이다.

아직 닫히지 않은 문제는 $v(a,t)$ 또는 필요한 확장상태를 수직 미시역학에서 유도하는 것이다. memory를 경험적 구성식에 숨기면 안 된다.

## 마일스톤 2 — 수직 반복 히스테리시스와 secular evolution

수직 반복응력

$$
\sigma_n(t)=\sigma_m+\sigma_a\sin\omega t
$$

아래에서 loading/unloading 차이를 역학으로부터 만들고, 더 강한 피로조건

$$
\boxed{P_{N+1}(a)\neq P_N(a)}
$$

를 얻어야 한다.

cycle-state drift가 없는 닫힌 주기루프는 internal friction일 수 있지만 피로누적은 아니다.

현재 generalized-LJ 완전사슬의 100 MPa null test는 거의 가역적인 응답을 올바르게 보여준다. 이 결과를 tuning으로 없애면 안 된다.

## 마일스톤 3 — 실험 피로주파수로의 시간척도 연결

직접적인 atomic normal-chain dynamics는 원자 시간척도를 가진다. 따라서 현재의 고주파 resonance 실험은 20 Hz 피로예측이 아니다.

핵심 문제는

$$
\text{빠른 atomic LJ dynamics}
\rightarrow
\text{projected memory / coarse state}
\rightarrow
P(a,t)\text{의 느린 진화}
$$

를 임의의 relaxation time 없이 유도하는 것이다.

출발점 후보는 Liouville projection, conditional propagator, 제거된 미시변수에서 직접 유도한 memory kernel이다.

## 마일스톤 4 — normal-opening 균열개시

활성 crack-initiation 그림은 수직 opening 또는 normal bond instability다. generalized LJ baseline에서는

$$
\phi''(\lambda_c)=0,
\qquad
\lambda_c\approx1.10777154
$$

에서 이상화된 국부 instability가 발생한다.

순간적인 tail occupancy

$$
Q_c(t)=\int_{a_c}^{\infty}P(a,t)\,da
$$

와 first-passage initiation을 구분해야 한다.

선호하는 정식화는

$$
\tau_c=\inf\{t:\text{역학적으로 정의된 normal-opening instability 발생}\}
$$

이다.

## 핵심 closure 질문

- $P(a,t)$만으로 충분한가?
- phase-space lift $F(a,c,t)$가 필요한가?
- 어떤 neighbor-spacing correlation이 필수인가?
- 경험적 damping 없이 실험주파수에서 exact projected memory를 축약할 수 있는가?
- 100 MPa급 loading에서 cycle-to-cycle broadening 또는 tail growth를 만드는 실제 물리는 무엇인가?

## 반증 테스트

활성 normal model은 최소한 다음을 만족해야 한다.

1. zero loading에서 인공적인 피로누적이 없어야 한다.
2. 가역 보존계 한계에서는 가역응답이 복원되어야 한다.
3. 확률 정규화가 보존되어야 한다.
4. density는 음수가 되면 안 된다.
5. 차원적으로 일관되어야 한다.
6. 에너지 수지를 만족해야 한다.
7. uniform-lattice limit에서 LJ baseline이 복원되어야 한다.
8. $P_{N+1}\neq P_N$가 numerical diffusion 때문이면 안 된다.
9. 20 Hz를 주장하려면 물리적으로 유도된 시간척도 연결이 있어야 한다.
10. microscopic interaction model 자체를 명시적으로 바꾸는 경우가 아니라면 LJ parameter는 고정되어야 한다.

## 보조 전단 라이브러리

기존 Rubin-chain, non-affine slip, gamma-surface 연구는 `libraries/shear/` 아래 보존한다. 방법론적 참고로는 유용하지만 활성 normal-deformation mainline 바깥이다.
