# Milestone 2 Research Note — Experimental constraints from pure Al cyclic deformation

**Purpose:** use experimental observations only as constraints on the state variables that a mechanics-derived theory must be able to reproduce. These observations are **not** inserted as empirical evolution laws.

## 1. Why this matters

Milestone 1 can produce a nonzero reduced hysteresis loop from conservative microscopic mechanics through energy transfer into hidden modes. That alone does not produce fatigue accumulation.

Milestone 2 requires a slow structural change such that

$$
P_{N+1}\neq P_N
$$

or, more generally, a joint structural state $\mathcal P_N$ changes irreversibly from cycle to cycle.

The key question is therefore: what structural information must the reduced state retain for pure Al?

## 2. Experimental constraint A — early dislocation-cell development in 99.99% pure Al single crystals

Watanabe, Yamazaki, and Koga studied 99.99%-purity Al single crystals under fully reversed tension-compression fatigue with single-slip-oriented specimens.

At room temperature they reported:

- no simple stress saturation;
- initial hardening;
- softening for orientations with active cross-slip;
- secondary hardening;
- dislocation cell structures forming and developing from the early stage of fatigue.

At 77 K the cyclic response changed substantially and showed a saturation/plateau behavior instead.

**Constraint on theory:** a realistic reduced model for room-temperature Al must allow a slow structural state to evolve even after a hysteresis loop already exists. A purely periodic internal-friction state is insufficient.

Reference:

- C. Watanabe, S. Yamazaki, N. Koga, “Effects of cross-slip activity on low-cycle fatigue behavior and dislocation structure in pure aluminum single crystals with single-slip orientation,” *Materials Science and Engineering A* 815 (2021) 141221. DOI: https://doi.org/10.1016/j.msea.2021.141221

## 3. Experimental constraint B — cell and deformation-band structures depend on crystal orientation

Experiments on cyclically deformed Al single crystals and bi-crystals found:

- elongated dislocation cells along primary $\{111\}$ slip traces;
- deformation bands associated with crystallographic orientation;
- local misorientation evolving during cycling;
- crystal orientation strongly affecting local stress and the resulting dislocation structures.

**Constraint on theory:** a scalar bond-length density $P(a,t)$ may be an important marginal, but it may not uniquely distinguish structures that have similar spacing statistics while having different shear/slip topology or orientation-dependent internal organization.

Reference:

- “Low cycle fatigue in aluminum single and bi-crystals: On the influence of crystal orientation,” *Materials Science and Engineering A* 668 (2016) 166–179. DOI: https://doi.org/10.1016/j.msea.2016.05.054

## 4. Experimental constraint C — crack initiation is not controlled by only one slip system

A study of 99.99% pure Al single crystals found that fatigue crack-initiation life depends not only on the primary slip system but also on secondary and higher-order slip systems.

**Constraint on theory:** a final crack-initiation model should not collapse all crystallographic shear information into a single scalar resolved shear coordinate unless the reduction can be rigorously justified.

Reference:

- “Effect of crystal orientation on fatigue crack initiation life in pure aluminum single crystals,” *International Journal of Fatigue* 156 (2022) 106661. DOI: https://doi.org/10.1016/j.ijfatigue.2021.106661

## 5. Implication for the state hierarchy

The present evidence suggests the following hierarchy should be tested rather than assumed:

### Level 1 — spacing marginal

$$
P(a,t).
$$

Keep this as the central marginal because it directly represents the population of local spacing states and connects naturally to bond-energy distributions.

### Level 2 — mechanically derived non-affine/shear coordinate

Introduce a coordinate $s^\alpha$ for slip system $\alpha$ only if it is defined directly from atomistic displacement/disregistry rather than by an empirical plastic-strain rule.

Candidate joint density:

$$
P(a,s^\alpha,t).
$$

For multiple active slip systems,

$$
P\!\left(a,s^1,s^2,\ldots,t\right).
$$

The original spacing state remains the exact marginal

$$
P(a,t)
=
\int P\!\left(a,\mathbf s,t\right)d\mathbf s.
$$

### Level 3 — correlation / topology information if necessary

If two microscopic configurations have the same one-point $P(a,t)$ and the same local shear coordinates but different dislocation-cell topology, a further correlation or defect-topology descriptor may be unavoidable.

The theory should add such a variable only after constructing an explicit counterexample showing that the lower-dimensional state is not closed.

## 6. Proposed mechanics-first route for Milestone 2

Do **not** write an empirical equation such as

$$
\frac{d\rho}{dN}=f(\rho,\sigma_a)
$$

and fit it.

Instead:

1. define local non-affine/slip coordinates directly from atomic positions;
2. derive their equations of motion from the same interatomic potential used for the spacing dynamics;
3. form the joint thermodynamic-limit density;
4. derive the exact continuity / Liouville hierarchy;
5. identify the first point where closure fails;
6. only then introduce a controlled approximation;
7. test whether one-cycle evolution gives

$$
\mathcal P_{N+1}\neq\mathcal P_N.
$$

The first successful secular evolution should be checked against the experimentally observed qualitative sequence of hardening / softening / secondary hardening **without fitting that sequence into the equations**.

## 7. Immediate theoretical target

The next concrete derivation should define a local atomistic disregistry coordinate across an FCC $\{111\}$ plane,

$$
s^\alpha
=
\left(\overline{\mathbf u}^{+}-\overline{\mathbf u}^{-}\right)\cdot\hat{\mathbf b}^{\alpha},
$$

where $\overline{\mathbf u}^{+}$ and $\overline{\mathbf u}^{-}$ are coarse displacements on the two sides of the selected slip plane and $\hat{\mathbf b}^{\alpha}$ is the corresponding slip direction.

Then derive, rather than postulate, the coupled state evolution of

$$
P(a,s^\alpha,t).
$$

This is currently the most direct bridge from the original spacing-distribution theory to the experimentally observed slow cyclic structural evolution of pure Al.

---

# 한국어 번역 — 순 Al 반복변형 실험이 Milestone 2에 주는 제약

**목적:** 실험 결과를 경험식으로 넣는 것이 아니라, 우리가 만드는 미시역학 이론이 반드시 재현해야 하는 물리적 제약으로만 사용한다.

## 1. 왜 중요한가

Milestone 1에서는 숨은 격자모드로 에너지가 전달되면서 경험적 damping 없이도 관심 좌표에 히스테리시스가 생길 수 있음을 확인했다.

하지만 피로 누적이 되려면

$$
P_{N+1}\neq P_N
$$

또는 더 일반적으로 joint structural state가 cycle마다 변해야 한다.

따라서 핵심은 **순 Al의 느린 구조진화를 표현하기 위해 최소한 어떤 상태변수가 필요한가**이다.

## 2. 제약 A — 99.99% 순 Al 단결정에서 초기부터 전위 cell이 발달

99.99% 순 Al 단결정의 완전 반전 tension-compression 실험에서는 실온에서 단순한 stress saturation이 나타나지 않았고, orientation에 따라

- 초기 hardening,
- softening,
- secondary hardening

이 나타났다. 동시에 전위 cell 구조가 피로 초기부터 형성되고 발달했다.

77 K에서는 거동이 달라져 saturation/plateau가 나타났다.

**이론에 대한 제약:** 실온 Al에서는 히스테리시스가 이미 존재한 뒤에도 느린 구조상태가 계속 변할 수 있어야 한다. 완전히 주기적인 internal-friction 상태만으로는 부족하다.

## 3. 제약 B — 구조진화는 결정방향과 slip에 민감

Al 단결정/쌍결정 반복변형 실험에서는

- primary $\{111\}$ slip trace를 따라 길게 형성되는 전위 cell,
- deformation band,
- cycle에 따른 local misorientation,
- orientation에 따른 국부응력 및 구조 차이

가 관찰됐다.

**이론에 대한 제약:** $P(a,t)$는 중요한 중심 marginal이지만, 동일하거나 유사한 bond-spacing 통계를 가지면서 shear/slip topology가 다른 두 상태를 완전히 구별하지 못할 가능성이 있다.

## 4. 제약 C — 균열개시는 primary slip 하나만으로 정해지지 않음

99.99% 순 Al 단결정의 crack-initiation 연구에서는 primary slip뿐 아니라 2차, 3차 이상의 slip system도 균열개시 수명에 영향을 주는 것으로 보고됐다.

**이론에 대한 제약:** 최종 crack-initiation 이론에서는 여러 결정학적 shear 자유도를 하나의 scalar에 무조건 압축하면 안 된다. 압축한다면 그 reduction이 역학적으로 정당화되어야 한다.

## 5. 현재 추천 상태 hierarchy

중심 변수는 계속

$$
P(a,t)
$$

로 유지한다.

다만 closure가 불가능하다는 것이 확인되면 원자좌표에서 직접 정의되는 non-affine/slip coordinate $s^\alpha$를 최소한으로 추가한다.

$$
P(a,s^\alpha,t)
$$

또는 여러 slip system이 필요하면

$$
P\!\left(a,s^1,s^2,\ldots,t\right).
$$

원래의 spacing distribution은 항상 marginal로 남는다.

$$
P(a,t)=\int P(a,\mathbf s,t)d\mathbf s.
$$

그마저 closure가 안 될 때에만 correlation 또는 defect-topology 정보를 추가한다.

## 6. Milestone 2의 역학 우선 경로

전위밀도 같은 변수를 먼저 놓고

$$
\frac{d\rho}{dN}=f(\rho,\sigma_a)
$$

형태를 fitting하는 방식은 사용하지 않는다.

대신

1. 원자좌표에서 local slip/disregistry를 정의하고,
2. 같은 interatomic potential로 운동방정식을 유도하고,
3. joint density를 정의하고,
4. exact continuity/Liouville hierarchy를 만들고,
5. closure가 깨지는 지점을 찾고,
6. 그 다음에만 controlled approximation을 사용한다.

최종적으로 한 cycle map이

$$
\mathcal P_{N+1}\neq\mathcal P_N
$$

을 실제 역학 결과로 만들어야 한다.

## 7. 바로 다음 유도 목표

FCC $\{111\}$ slip plane 양쪽 원자들의 평균 상대변위에서

$$
s^\alpha
=
\left(\overline{\mathbf u}^{+}-\overline{\mathbf u}^{-}\right)\cdot\hat{\mathbf b}^{\alpha}
$$

를 정의하고, 이 $s^\alpha$와 spacing $a$의 결합 분포

$$
P(a,s^\alpha,t)
$$

의 진화를 미시역학에서 직접 유도하는 것이 현재 가장 좋은 다음 단계다.
