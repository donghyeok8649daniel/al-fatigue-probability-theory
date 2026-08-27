# Open Problems

## Milestone 1 — Mechanics-derived hysteresis

Given a prescribed cyclic stress

$$
\sigma(t)=\sigma_m+\sigma_a\sin(\omega t),
$$

derive an internal evolution law from microscopic mechanics such that loading and unloading do not simply retrace the same structural state, while conserving probability and satisfying the correct energy balance.

Success condition:

$$
A_H=\oint \sigma\,d\epsilon>0
$$

without inserting an empirical hysteresis law.

## Milestone 2 — Secular fatigue accumulation

A periodic hysteresis loop alone is internal friction, not fatigue. Derive conditions for

$$
P_{N+1}\neq P_N
$$

or for an equivalent slow internal state to evolve cycle by cycle.

## Milestone 3 — Crack initiation

Formulate initiation as a mechanical stability loss or first-passage event. Candidate formulations include an absorbing boundary in an enlarged state space or a distribution-level stability criterion.

## Central closure problem

Determine the minimum state required for a mechanically closed model:

- $P(a,t)$ only?
- phase-space density $F(a,c,t)$?
- spacing-correlation hierarchy?
- joint structural state such as $P(a,s,t)$, where $s$ is a non-affine/slip coordinate?

The preferred answer is the smallest state that can be derived from mechanics without hiding essential memory or irreversibility in fitted constitutive terms.

## Numerical falsification tests

Any candidate model should pass at least:

1. zero loading → zero hysteresis;
2. perfectly reversible conservative limit → zero hysteresis;
3. probability normalization preserved;
4. non-negative density;
5. dimensional consistency;
6. energy balance;
7. uniform-lattice limit recovers the baseline lattice energy;
8. nonzero fatigue accumulation is not created by numerical diffusion alone.

---

# 한국어 번역 — 미해결 문제

## 마일스톤 1 — 역학으로부터 유도되는 히스테리시스

다음과 같은 주기 응력이 주어졌다고 하자.

$$
\sigma(t)=\sigma_m+\sigma_a\sin(\omega t)
$$

미시적 역학으로부터 내부 상태의 진화 법칙을 유도한다. 이 진화 법칙은 하중과 제하 과정이 단순히 동일한 구조 상태를 역으로 되짚지 않도록 해야 하며, 동시에 확률 보존과 올바른 에너지 수지를 만족해야 한다.

성공 조건은

$$
A_H=\oint \sigma\,d\epsilon>0
$$

이며, 이를 위해 경험적인 히스테리시스 법칙을 별도로 삽입해서는 안 된다.

## 마일스톤 2 — 장기적 피로 누적

주기적으로 완전히 반복되는 히스테리시스 루프만 존재한다면 그것은 내부 마찰일 뿐, 그 자체로 피로 누적을 의미하지 않는다. 따라서 다음 조건이 성립하는 역학적 조건을 유도해야 한다.

$$
P_{N+1}\neq P_N
$$

또는 이와 동등하게, 어떤 느린 내부 상태가 cycle마다 누적적으로 진화해야 한다.

## 마일스톤 3 — 균열 개시

균열 개시를 기계적 안정성 상실 또는 first-passage 사건으로 정식화한다. 후보 접근으로는 확장된 상태공간에서의 흡수 경계 조건이나 분포 수준의 안정성 기준이 있다.

## 핵심 closure 문제

기계적으로 닫힌 모델을 구성하기 위해 필요한 최소 상태를 결정한다.

- $P(a,t)$만으로 충분한가?
- 위상공간 밀도 $F(a,c,t)$가 필요한가?
- 원자간격 상관관계의 hierarchy가 필요한가?
- 비아핀/non-affine 또는 slip 좌표 $s$를 포함하는 $P(a,s,t)$와 같은 결합 구조 상태가 필요한가?

목표는 필수적인 memory 또는 irreversibility를 fitting된 구성식 안에 숨기지 않으면서, 미시적 역학으로부터 유도 가능한 가장 작은 상태 표현을 찾는 것이다.

## 수치적 반증 테스트

어떤 후보 모델이든 최소한 다음 테스트를 통과해야 한다.

1. 외력이 0이면 히스테리시스도 0이어야 한다.
2. 완전히 가역적인 보존계 한계에서는 히스테리시스가 0이어야 한다.
3. 확률 정규화가 보존되어야 한다.
4. 확률밀도는 음수가 되어서는 안 된다.
5. 모든 주요 식은 차원적으로 일관되어야 한다.
6. 에너지 수지를 만족해야 한다.
7. 균일 격자 한계에서 기존의 기준 lattice energy가 복원되어야 한다.
8. 0이 아닌 피로 누적이 단순한 numerical diffusion 때문에 인공적으로 발생해서는 안 된다.
