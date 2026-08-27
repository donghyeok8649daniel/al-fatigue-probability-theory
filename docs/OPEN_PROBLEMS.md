# Open Problems

## Main research direction

The principal problem is cyclic **normal deformation** under normal stress. The central state remains the spacing distribution $P(a,t)$.

The intended mechanics chain is

$$
\boxed{
\sigma_n(t)
\rightarrow
\{a_i(t)\}
\rightarrow
P(a,t)
\rightarrow
\text{normal hysteresis}
\rightarrow
P_{N+1}(a)\neq P_N(a)
\rightarrow
\text{normal-opening instability / crack initiation}.
}
$$

Shear/slip coordinates are auxiliary unless they are later shown to be mathematically necessary for closure of the normal-deformation problem.

## Milestone 1 — Mechanics-derived normal hysteresis

Given

$$
\sigma_n(t)=\sigma_m+\sigma_a\sin(\omega t),
$$

derive a branch-dependent normal strain/spacing response from microscopic mechanics while satisfying probability conservation and energy balance.

Success condition:

$$
\boxed{
A_H=\oint\sigma_n\,d\epsilon_n>0
}
$$

without inserting an empirical hysteresis law.

### Current status

The Rubin-chain calculation proves the general mechanism that a conservative microscopic system can show reduced hysteresis when energy is transferred into unresolved propagating modes.

However, this proof currently uses a generic resolved coordinate $Q$. The mainline task is to derive the analogous effect directly for the normal-spacing coordinate $a$ or its exact phase-space/correlation extension.

## Milestone 2 — Secular normal-spacing evolution

A periodic hysteresis loop alone is not fatigue. The main success criterion is

$$
\boxed{
P_{N+1}(a)\neq P_N(a)
}
$$

at identical cycle phase, generated from microscopic normal-deformation mechanics.

The existing slip-bath result $s_{N+1}\neq s_N$ is an auxiliary existence proof only. It does not solve this Milestone.

### Highest-priority questions

1. Can anharmonic normal lattice dynamics generate a secular change in $P(a,t)$ under cyclic normal loading?
2. What spacing correlations are necessary for closure?
3. Does a phase-space density $F(a,c,t)$ suffice, or is a higher pair/joint hierarchy required?
4. Can a free surface or geometric normal stress concentration produce localized opening without introducing an empirical damage law?
5. How does finite temperature alter the initial phase-space ensemble and tail transport?
6. Can all observed cycle-to-cycle drift survive time-step, system-size, and boundary-condition refinement?

## Milestone 3 — Crack initiation by normal opening

The main crack-initiation route is treated as a normal-opening stability or first-passage problem.

The idealized one-coordinate criterion

$$
U''(a_c)=0
$$

is a useful baseline, not a complete crack criterion.

A cumulative initiation theory should distinguish

$$
Q_c(t)=\int_{a_c}^{\infty}P(a,t)\,da
$$

from the first-passage probability

$$
F_{\rm ci}(t)=\Pr(\tau_c\le t).
$$

## Central closure problem

Determine the minimum state required for a mechanically closed normal-deformation model:

- $P(a,t)$ only?
- phase-space density $F(a,c,t)$?
- neighboring-spacing joint densities?
- full pair-distance hierarchy $P_k(r,t)$?
- a surface-normal structural coordinate in addition to $a$?

The preferred state is the smallest one that follows from microscopic mechanics without hiding memory or irreversibility in fitted constitutive laws.

## Energy-model problem

Use a fixed microscopic potential and let the **state** evolve rather than changing the potential with damage.

The principal analytic baseline is the generalized Lennard-Jones pair law

$$
v(r)=\varepsilon_{\rm LJ}
\left[
\left(\frac{\sigma_{\rm LJ}}{r}\right)^m
-
\left(\frac{\sigma_{\rm LJ}}{r}\right)^n
\right].
$$

A major open question is how far the resulting pair-potential mechanics can be pushed before many-body corrections become quantitatively necessary for Al. Such corrections should be introduced as validation upgrades, not as hidden phenomenological fatigue laws.

## Numerical falsification tests

Any mainline normal-deformation model should pass at least:

1. zero loading $\rightarrow$ zero hysteresis and zero secular drift;
2. perfectly reversible affine limit $\rightarrow P_{N+1}=P_N$;
3. probability normalization preserved;
4. non-negative density;
5. dimensional consistency;
6. energy balance;
7. uniform-lattice limit recovers the baseline lattice energy;
8. secular $P$ evolution is not created by numerical diffusion;
9. convergence with time step and system size;
10. no arbitrary barrier reduction or damage-dependent LJ parameters.

---

# 한국어 번역 — 미해결 문제

## 주 연구방향

핵심 문제는 수직응력 아래의 반복 **수직변형**이다. 중심 상태는 계속 spacing distribution $P(a,t)$이다.

목표 역학 chain은

$$
\boxed{
\sigma_n(t)
\rightarrow
\{a_i(t)\}
\rightarrow
P(a,t)
\rightarrow
\text{수직 히스테리시스}
\rightarrow
P_{N+1}(a)\neq P_N(a)
\rightarrow
\text{수직 opening instability / crack initiation}
}
$$

이다.

전단/slip 좌표는 향후 수직변형 문제의 closure에 수학적으로 필요하다고 확인되지 않는 한 보조적인 위치에 둔다.

## Milestone 1 — 역학에서 유도되는 수직 히스테리시스

다음 수직응력이 주어졌다고 하자.

$$
\sigma_n(t)=\sigma_m+\sigma_a\sin(\omega t)
$$

확률보존과 에너지수지를 만족하면서 microscopic mechanics에서 loading/unloading branch가 다른 수직 strain/spacing response를 유도한다.

성공조건은

$$
\boxed{
A_H=\oint\sigma_n\,d\epsilon_n>0
}
$$

이며 empirical hysteresis law를 삽입해서는 안 된다.

### 현재 상태

Rubin-chain 계산은 보존적인 미시계에서도 에너지가 해소되지 않은 propagating mode로 전달되면 축약 히스테리시스가 생길 수 있다는 일반 메커니즘을 증명했다.

다만 현재 이 원리증명은 generic resolved coordinate $Q$를 사용한다. 메인 과제는 동일한 효과를 normal-spacing coordinate $a$ 또는 그 정확한 phase-space/correlation 확장에서 직접 유도하는 것이다.

## Milestone 2 — cycle별 normal-spacing 누적진화

주기 히스테리시스만으로는 피로가 아니다. 메인 성공조건은 동일한 cycle phase에서

$$
\boxed{
P_{N+1}(a)\neq P_N(a)
}
$$

가 microscopic normal-deformation mechanics에서 발생하는 것이다.

기존 slip-bath의 $s_{N+1}\neq s_N$ 결과는 가능성을 보여주는 보조 existence proof일 뿐 이 Milestone을 해결한 것은 아니다.

### 최우선 질문

1. anharmonic normal lattice dynamics가 반복 수직하중에서 $P(a,t)$의 secular change를 만들 수 있는가?
2. closure에 어떤 spacing correlation이 필요한가?
3. phase-space density $F(a,c,t)$만으로 충분한가, 아니면 더 높은 joint/pair hierarchy가 필요한가?
4. 자유표면 또는 수직방향 형상 응력집중이 empirical damage law 없이 국부 opening을 만들 수 있는가?
5. 유한온도가 초기 phase-space ensemble과 tail transport를 어떻게 바꾸는가?
6. 관찰된 cycle drift가 time-step, system-size, boundary-condition refinement 뒤에도 유지되는가?

## Milestone 3 — 수직 opening에 의한 균열개시

주 crack-initiation route는 수직 opening의 안정성 상실 또는 first-passage 문제로 취급한다.

이상화된 단일좌표 기준

$$
U''(a_c)=0
$$

은 유용한 baseline이지만 완전한 crack criterion은 아니다.

누적 균열개시 이론에서는

$$
Q_c(t)=\int_{a_c}^{\infty}P(a,t)\,da
$$

와 first-passage probability

$$
F_{\rm ci}(t)=\Pr(\tau_c\le t)
$$

를 구분해야 한다.

## 핵심 closure 문제

기계적으로 닫힌 수직변형 모델에 필요한 최소상태를 결정한다.

- $P(a,t)$만으로 충분한가?
- phase-space density $F(a,c,t)$가 필요한가?
- 인접 spacing joint density가 필요한가?
- 전체 pair-distance hierarchy $P_k(r,t)$가 필요한가?
- $a$ 외에 surface-normal structural coordinate가 필요한가?

목표는 memory 또는 irreversibility를 fitting된 구성식에 숨기지 않으면서 microscopic mechanics에서 따라오는 가장 작은 상태를 찾는 것이다.

## 에너지모델 문제

고정된 microscopic potential을 사용하고 damage에 따라 potential을 바꾸는 대신 **state가 진화하도록** 한다.

주 해석 baseline은 generalized Lennard-Jones pair law다.

$$
v(r)=\varepsilon_{\rm LJ}
\left[
\left(\frac{\sigma_{\rm LJ}}{r}\right)^m
-
\left(\frac{\sigma_{\rm LJ}}{r}\right)^n
\right].
$$

중요한 open question은 Al에 대해 many-body correction이 정량적으로 필요해지기 전까지 pair-potential mechanics를 어디까지 밀어붙일 수 있는가이다. many-body correction은 숨은 피로 phenomenology가 아니라 validation upgrade로 도입해야 한다.

## 수치적 반증 테스트

메인 수직변형 모델은 최소한 다음을 통과해야 한다.

1. 외력 0 $\rightarrow$ 히스테리시스 0, secular drift 0;
2. 완전 가역 affine limit $\rightarrow P_{N+1}=P_N$;
3. 확률 정규화 보존;
4. density 비음수;
5. 차원 일관성;
6. 에너지수지;
7. 균일격자 limit에서 baseline lattice energy 복원;
8. secular $P$ evolution이 numerical diffusion 때문에 생기지 않음;
9. time step과 system size에 대해 convergence;
10. arbitrary barrier reduction 및 damage-dependent LJ parameter 금지.
