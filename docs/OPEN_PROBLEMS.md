# Open Problems — Active 1D Normal-LJ Mainline

## 1. Continuous-time probability state

The active state is

$$
P(a,t),
$$

not a cycle-indexed family.

The exact kinematic identity is

$$
\boxed{
\partial_tP+\partial_a(Pv)=0.
}
$$

The theory may use finite empirical densities for numerical work, but any finite index counts represented spacings or atoms, not fatigue cycles.

## 2. Mean and energy as reduced observables

Define

$$
\mu(t)=\int\lambda P(\lambda,t)\,d\lambda
$$

and

$$
\mathcal E(t)=\int\psi(\lambda)P(\lambda,t)\,d\lambda.
$$

The current hypothesis is that the mean may remain nearly conserved while the distribution broadens and stores configurational energy.

The exact identity

$$
\mathcal E(t)-\psi(\mu(t))
=
\int D_\psi(\lambda\mid\mu(t))P(\lambda,t)\,d\lambda
$$

makes this statement precise inside the convex LJ region.

## 3. Exact obstruction: reverse compression

Normalization, mean, and energy are not enough to force a tensile tail.

Because the generalized-LJ repulsion diverges at small spacing, a crack-free measure can store arbitrarily large energy through sufficiently strong compression while remaining below $\lambda_c$ on the tensile side.

This is now an exact falsification result, not an open question.

## 4. Highest-priority open problem: derive the compression bound

The energy-feasibility theorem becomes finite only after imposing a physically justified lower support bound

$$
\boxed{
\lambda\ge\lambda_L(t)>0.
}
$$

The highest-priority problem is therefore

$$
\boxed{
\text{derive or independently constrain }\lambda_L(t)
\text{ from 1D normal-LJ mechanics.}
}
$$

The bound must not be selected to obtain a desired fatigue life.

Candidate one-dimensional routes include:

1. finite-chain total-energy accessibility;
2. imposed-force and boundary-condition bounds;
3. a rigorously defined maximum reverse-compression work;
4. experimentally measured minimum normal spacing or strain bounds if available;
5. a stronger integral bound on the compression-side energy instead of a hard support bound.

## 5. Crack-free energy ceiling

Once $\lambda_L(t)$ is established, the exact maximum crack-free energy at mean $\mu(t)$ is

$$
\boxed{
\mathcal E_{\rm safe}^{\max}(t)
=
\frac{\lambda_c-\mu(t)}{\lambda_c-\lambda_L(t)}\psi(\lambda_L(t))
+
\frac{\mu(t)-\lambda_L(t)}{\lambda_c-\lambda_L(t)}\psi(\lambda_c).
}
$$

The associated energy margin is

$$
M_E(t)=\mathcal E_{\rm safe}^{\max}(t)-\mathcal E(t).
$$

The next validation problem is to compute $\mu(t)$, $\mathcal E(t)$, and the mechanically derived compression constraint in the same 1D simulation and verify whether the energy margin behaves consistently under null and non-null loading.

## 6. Continuous-time initiation

Define

$$
\boxed{
\tau_E
=
\inf\{t\ge0:M_E(t)<0\}.
}
$$

If $\lambda_L(t)$ is a true hard lower bound, then $M_E<0$ implies that some probability mass must leave the safe interval through the tensile side.

The instantaneous tail is

$$
Q_c(t)=\int_{\lambda_c}^{\infty}P(\lambda,t)\,d\lambda.
$$

A later task is to connect the energy-feasibility crossing to a physically precise crack-initiation event for a finite specimen or representative region.

## 7. Energy input from cyclic loading

A separate unresolved problem is the relation

$$
\boxed{
\sigma(t)
\rightarrow
W(t)
\rightarrow
\mathcal E(t).
}
$$

Not all external work is retained as configurational energy. Kinetic energy and recoverable elastic work must be separated by an exact energy balance.

No arbitrary retained-energy fraction is allowed.

## 8. 1D-only closure questions

- Is $P(a,t)$ plus $\mathcal E(t)$ sufficient for the crack-feasibility bound even if it is not sufficient for full dynamics?
- Can $\lambda_L(t)$ be obtained from total 1D LJ energy and finite-system constraints?
- Is a hard support bound too strong, and would a bound on compression-side Bregman energy be a more natural third condition?
- What finite representative length is required for a crack-initiation statement?
- How should kinetic energy be excluded from $\mathcal E(t)$ without losing exact energy accounting?
- Can the measured external normal work define a rigorous upper/lower bound on the configurational-energy trajectory?

## 9. Required falsification tests

Any active model must satisfy:

1. zero loading produces no artificial energy accumulation;
2. reversible conservative limits remain reversible;
3. probability normalization is preserved;
4. the mean and configurational energy are computed directly from the same $P(a,t)$;
5. dimensions are consistent;
6. full work-energy balance is checked;
7. the LJ potential parameters remain fixed;
8. no named probability family is inserted for convenience;
9. the energy-feasibility theorem is never used without explicitly stating the compression constraint;
10. illustrative values of $\lambda_L$ are never reported as Al material constants;
11. any claimed $\tau_E$ is a physical-time first passage, not a fitted cycle-life law;
12. 3D FCC and shear archives are not imported into the active default calculation.

---

# 한국어 번역 — 활성 1D Normal-LJ Mainline의 미해결 문제

## 1. 연속시간 확률상태

활성 상태는 cycle-indexed family가 아니라

$$
P(a,t)
$$

이다.

정확한 kinematic identity는

$$
\boxed{
\partial_tP+\partial_a(Pv)=0
}
$$

이다.

numerical work에서 finite empirical density를 사용할 수 있지만 finite index는 represented spacing 또는 atom 수를 뜻하며 fatigue cycle count가 아니다.

## 2. 축약 observable로서 평균과 에너지

$$
\mu(t)=\int\lambda P(\lambda,t)\,d\lambda
$$

및

$$
\mathcal E(t)=\int\psi(\lambda)P(\lambda,t)\,d\lambda
$$

를 정의한다.

현재 가설은 평균이 거의 보존되는 동안 distribution이 넓어지고 configurational energy를 저장할 수 있다는 것이다.

convex LJ region에서

$$
\mathcal E(t)-\psi(\mu(t))
=
\int D_\psi(\lambda\mid\mu(t))P(\lambda,t)\,d\lambda
$$

라는 exact identity가 이 주장을 수학적으로 명확하게 만든다.

## 3. 정확한 장애물: reverse compression

정규화, 평균, 에너지만으로 tensile tail을 강제할 수 없다.

generalized-LJ repulsion은 small spacing에서 발산하므로 tensile side에서 $\lambda_c$ 아래에 머무르면서도 충분히 강한 compression으로 임의로 큰 energy를 저장할 수 있다.

이것은 이제 open question이 아니라 exact falsification result다.

## 4. 가장 우선적인 open problem: compression bound 유도

energy-feasibility theorem이 유한한 ceiling을 가지려면 물리적으로 정당화된 lower support bound

$$
\boxed{
\lambda\ge\lambda_L(t)>0
}
$$

이 필요하다.

따라서 가장 중요한 문제는

$$
\boxed{
\text{1D normal-LJ mechanics로부터 }\lambda_L(t)
\text{를 유도하거나 독립적으로 제약하는 것}
}
$$

이다.

원하는 fatigue life를 만들기 위해 이 bound를 선택하면 안 된다.

가능한 1차원 경로는 다음과 같다.

1. finite-chain total-energy accessibility;
2. imposed-force 및 boundary-condition bound;
3. 엄밀하게 정의된 maximum reverse-compression work;
4. 가능하다면 직접 측정된 minimum normal spacing 또는 strain bound;
5. hard support bound 대신 compression-side energy에 대한 더 강한 integral bound.

## 5. crack-free energy ceiling

$\lambda_L(t)$가 확보되면 평균 $\mu(t)$에서 정확한 maximum crack-free energy는

$$
\boxed{
\mathcal E_{\rm safe}^{\max}(t)
=
\frac{\lambda_c-\mu(t)}{\lambda_c-\lambda_L(t)}\psi(\lambda_L(t))
+
\frac{\mu(t)-\lambda_L(t)}{\lambda_c-\lambda_L(t)}\psi(\lambda_c)
}
$$

이다.

energy margin은

$$
M_E(t)=\mathcal E_{\rm safe}^{\max}(t)-\mathcal E(t)
$$

이다.

다음 validation 문제는 같은 1D simulation에서 $\mu(t)$, $\mathcal E(t)$, mechanics-derived compression constraint를 함께 계산하고 null/non-null loading에서 energy margin이 일관되게 거동하는지 확인하는 것이다.

## 6. 연속시간 initiation

$$
\boxed{
\tau_E
=
\inf\{t\ge0:M_E(t)<0\}
}
$$

로 정의한다.

$\lambda_L(t)$가 true hard lower bound라면 $M_E<0$는 일부 probability mass가 safe interval을 tensile side로 빠져나가야 한다는 뜻이다.

순간 tail은

$$
Q_c(t)=\int_{\lambda_c}^{\infty}P(\lambda,t)\,d\lambda
$$

이다.

향후에는 energy-feasibility crossing을 finite specimen 또는 representative region의 물리적으로 정확한 crack-initiation event와 연결해야 한다.

## 7. cyclic loading으로부터의 energy input

별도의 미해결 문제는

$$
\boxed{
\sigma(t)
\rightarrow
W(t)
\rightarrow
\mathcal E(t)
}
$$

관계다.

외부 work 전체가 configurational energy로 남는 것은 아니다. kinetic energy와 recoverable elastic work를 exact energy balance로 분리해야 한다.

임의의 retained-energy fraction은 허용하지 않는다.

## 8. 1D-only closure 질문

- $P(a,t)$가 full dynamics에는 부족해도 $\mathcal E(t)$와 함께 crack-feasibility bound에는 충분한가?
- total 1D LJ energy와 finite-system constraint로 $\lambda_L(t)$를 얻을 수 있는가?
- hard support bound가 지나치게 강하다면 compression-side Bregman energy bound가 더 자연스러운 세 번째 조건인가?
- crack-initiation statement에 필요한 finite representative length는 얼마인가?
- exact energy accounting을 유지하면서 kinetic energy를 $\mathcal E(t)$에서 어떻게 분리할 것인가?
- measured external normal work로 configurational-energy trajectory의 rigorous upper/lower bound를 만들 수 있는가?

## 9. 필수 falsification test

모든 활성 model은 다음을 만족해야 한다.

1. zero loading에서 artificial energy accumulation이 없어야 한다.
2. reversible conservative limit는 reversible해야 한다.
3. probability normalization이 보존되어야 한다.
4. mean과 configurational energy는 동일한 $P(a,t)$에서 직접 계산해야 한다.
5. dimension이 일치해야 한다.
6. full work-energy balance를 확인해야 한다.
7. LJ potential parameter는 고정되어야 한다.
8. 편의를 위한 named probability family를 넣지 않는다.
9. compression constraint를 명시하지 않고 energy-feasibility theorem을 사용하지 않는다.
10. illustrative $\lambda_L$를 Al material constant로 보고하지 않는다.
11. 어떤 $\tau_E$ 주장도 fitted cycle-life law가 아니라 physical-time first passage여야 한다.
12. 3D FCC 및 shear archive를 active default calculation에서 import하지 않는다.
