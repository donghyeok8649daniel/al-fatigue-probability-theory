# Open Problems — Active Normal-Deformation Mainline

## Milestone 1 — Normal microscopic state and exact transport

The active state variable is the normal-spacing density

$$
P(a,t)=\lim_{N\to\infty}\frac1N\sum_i\delta(a-a_i(t)).
$$

The exact kinematic transport equation is

$$
\boxed{\partial_tP+\partial_a(Pv)=0.}
$$

The unresolved problem is closure: derive $v(a,t)$, or the minimum enlarged state required to determine it, from normal microscopic mechanics without hiding memory inside an empirical constitutive law.

## Milestone 2 — Normal cyclic hysteresis and secular evolution

For prescribed normal stress

$$
\sigma_n(t)=\sigma_m+\sigma_a\sin\omega t,
$$

derive a mechanics-generated loading/unloading difference and then the stronger fatigue condition

$$
\boxed{P_{N+1}(a)\neq P_N(a)}.
$$

A closed periodic loop with no cycle-state drift is internal friction, not fatigue accumulation.

The current generalized-LJ perfect-chain 100 MPa null test correctly shows essentially reversible behavior. This result must not be tuned away.

Key questions are:

1. Can anharmonic normal lattice dynamics generate a secular change in $P(a,t)$ under experimentally relevant cyclic normal loading?
2. Which neighbor-spacing correlations are necessary for closure?
3. Is the phase-space density $F(a,c,t)$ sufficient, or is a higher correlation hierarchy required?
4. Can a free surface or mechanically derived normal stress concentration produce localized opening without an empirical damage variable?
5. Does any observed cycle-to-cycle drift survive time-step, system-size, and boundary-condition refinement?

## Milestone 3 — Three-dimensional FCC normal pair-lattice validation

### Current status: partially achieved

The active 3D homogeneous FCC pair model is

$$
\boxed{
U(\mathbf F)=\frac12\sum_{\mathbf R\ne0}v(|\mathbf F\mathbf R|).
}
$$

For [001] normal loading,

$$
\mathbf F=\operatorname{diag}(\lambda_t,\lambda_t,\lambda_n),
$$

with $\lambda_t$ relaxed at fixed $\lambda_n$.

Using the inherited generalized-LJ exponents

$$
m=12.19,\qquad n=6,
$$

and calibrating only the LJ energy scale to the directional experimental modulus

$$
E_{[001]}\approx62.7024\ \mathrm{GPa},
$$

the model predicts

$$
C_{11}^{\rm LJ}\approx107.169\ \mathrm{GPa},
$$

$$
C_{12}^{\rm LJ}\approx61.180\ \mathrm{GPa}.
$$

These closely match the external normal-elastic reference values $107$ GPa and $61$ GPa.

The unfitted relaxed [001] ideal engineering strength is approximately

$$
\boxed{9.045\ \mathrm{GPa}},
$$

which is in the same scale as a first-principles [001] reference near $10.63$ GPa.

This is a strong positive result for the active normal-mechanics backbone.

### Exact limitation exposed by the same calculation

A cubic central pair potential at zero pressure obeys the Cauchy relation

$$
\boxed{C_{12}=C_{44}.}
$$

The numerical FCC sum reproduces this relation to numerical accuracy, predicting $C_{44}\approx61.180$ GPa, whereas real Al is approximately $29$ GPa.

Therefore the fixed central pair model is not a complete quantitative three-dimensional Al potential. This limitation is structural, not numerical.

## Milestone 4 — Absolute cohesive-energy consistency

This is now the highest-priority energy-model problem before quantitative thermal probability is introduced.

With the same fixed generalized-LJ shape:

- fitting the experimental cohesive energy $E_{\rm coh}\approx3.43$ eV/atom gives
  $$
  E_{[001]}^{\rm LJ}\approx220.466\ \mathrm{GPa},
  $$
  which is far too stiff;
- fitting the correct normal modulus gives
  $$
  E_{\rm coh}^{\rm LJ}\approx0.976\ \mathrm{eV/atom},
  $$
  which is far too small.

Thus the current pair law cannot simultaneously reproduce the tangent normal mechanics and the absolute separation-energy scale.

The central problem is

$$
\boxed{
\text{retain the successful normal LJ mechanics}
+\text{derive the minimum physically justified cohesive/many-body correction}.
}
$$

Requirements for any correction:

1. its microscopic physical origin must be explicit;
2. it must enter through an energy or Hamiltonian, not through a fitted fatigue-damage law;
3. it must not vary with cycle count merely to imitate degradation;
4. it must preserve the already successful normal small-strain limit unless the new physics requires a calculable correction;
5. the number of independent material inputs must be minimized and each must be classified;
6. the correction must be tested against cohesion and normal traction before any fatigue-life calculation.

Candidate directions may include an explicitly derived electronic-density / many-body cohesive term, but no specific functional form is accepted yet.

## Milestone 5 — Finite-temperature distribution and time-scale bridge

Thermal activation must not be added quantitatively until Milestone 4 fixes or bounds the relevant absolute energy landscape.

A probability containing

$$
\exp\!\left(-\frac{\Delta U}{k_B T}\right)
$$

is exponentially sensitive to $\Delta U$. Therefore an effective LJ potential that has the correct normal stiffness but the wrong cohesive-energy scale cannot be used directly for a quantitative escape rate.

After the energy problem is controlled, the next task is to derive

$$
\boxed{
\text{fast atomic dynamics}
\rightarrow
\text{exact projected memory / coarse state}
\rightarrow
\text{slow evolution of }P(a,t)
}
$$

without inserting an arbitrary relaxation time.

Exact starting points include microscopic Liouville dynamics, conditional propagators, and projection-derived memory kernels. Any Markovian reduction must be justified from this level rather than postulated.

Key questions include:

1. What is the correct finite-temperature microscopic state: $P(a,t)$ plus temperature, phase-space density, or a larger correlation hierarchy?
2. Which part of the thermal energy changes the normal-spacing distribution and which part remains kinetic?
3. Can a controlled separation of atomic and laboratory time scales be obtained from the microscopic spectrum or projected memory kernel?
4. Can a 20 Hz cycle map be derived without inserting an empirical relaxation time?

## Milestone 6 — Normal-opening crack initiation

The active crack-initiation picture is normal opening or normal mechanical stability loss.

The old reduced 1D LJ criterion

$$
\phi''(\lambda_c)=0,
\qquad
\lambda_c\approx1.10777154
$$

remains a useful reduced-model diagnostic, but the 3D FCC model shows that the ultimate normal stability condition should ultimately be defined from the relevant multidimensional energy Hessian or loss of a stable equilibrium branch.

Instantaneous tail occupancy

$$
Q_c(t)=\int_{a_c}^{\infty}P(a,t)\,da
$$

must be distinguished from first-passage initiation.

The preferred cumulative formulation is

$$
\tau_c=\inf\{t:\text{a mechanically defined normal-opening instability occurs}\},
$$

with

$$
F_{\rm ci}(t)=\Pr(\tau_c\le t).
$$

## Central closure questions

- Is $P(a,t)$ sufficient?
- Is the phase-space lift $F(a,c,t)$ required?
- Which neighbor-spacing joint densities are essential?
- Is the full pair-distance hierarchy $P_k(r,t)$ needed?
- What is the minimum extra energy variable or density variable required by the failure of the central pair potential?
- Can projected memory be reduced at laboratory frequencies without empirical damping?
- What physical mechanism produces cycle-to-cycle broadening or tail growth under 100 MPa-class normal loading?

## Falsification tests

Any active normal model must satisfy at least:

1. zero loading gives zero artificial fatigue accumulation;
2. reversible conservative limits recover reversible behavior;
3. probability normalization is preserved;
4. density remains non-negative;
5. dimensions are consistent;
6. energy balance is satisfied;
7. the uniform-lattice limit recovers the stated microscopic energy model;
8. $P_{N+1}\neq P_N$ is not numerical diffusion;
9. a 20 Hz claim uses a physically derived time-scale bridge;
10. microscopic potential parameters remain fixed unless the microscopic interaction model itself is explicitly changed;
11. lattice-sum and finite-difference errors are much smaller than any claimed material discrepancy;
12. thermal first-passage predictions are not made from an unvalidated absolute barrier-energy scale.

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
\boxed{\partial_tP+\partial_a(Pv)=0}
$$

이다.

아직 닫히지 않은 문제는 $v(a,t)$ 또는 이를 결정하기 위해 필요한 최소 확장상태를 수직 미시역학에서 유도하는 것이다. memory를 경험적 구성식 안에 숨기면 안 된다.

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

핵심 질문은 다음과 같다.

1. anharmonic normal lattice dynamics가 실험적으로 의미 있는 수직 반복하중에서 $P(a,t)$의 secular change를 만들 수 있는가?
2. closure를 위해 어떤 neighbor-spacing correlation이 필요한가?
3. phase-space density $F(a,c,t)$로 충분한가, 아니면 더 높은 correlation hierarchy가 필요한가?
4. 경험적 damage variable 없이 자유표면 또는 역학적으로 유도된 normal stress concentration이 국부 opening을 만들 수 있는가?
5. 관찰된 cycle-to-cycle drift가 time-step, system-size, boundary-condition refinement 후에도 남는가?

## 마일스톤 3 — 3차원 FCC normal pair-lattice 검증

### 현재 상태: 부분 달성

활성 3D homogeneous FCC pair model은

$$
\boxed{
U(\mathbf F)=\frac12\sum_{\mathbf R\ne0}v(|\mathbf F\mathbf R|)
}
$$

이다.

[001] normal loading에서는

$$
\mathbf F=\operatorname{diag}(\lambda_t,\lambda_t,\lambda_n)
$$

을 사용하고 고정 $\lambda_n$에서 $\lambda_t$를 relaxation한다.

기존 generalized-LJ exponent

$$
m=12.19,\qquad n=6
$$

을 유지하고 LJ energy scale 하나만 방향성 experimental modulus

$$
E_{[001]}\approx62.7024\ \mathrm{GPa}
$$

에 맞추면

$$
C_{11}^{\rm LJ}\approx107.169\ \mathrm{GPa},
$$

$$
C_{12}^{\rm LJ}\approx61.180\ \mathrm{GPa}
$$

가 나온다.

외부 normal-elastic reference $107$ GPa, $61$ GPa와 매우 가깝다.

peak를 fitting하지 않은 relaxed [001] ideal engineering strength는 약

$$
\boxed{9.045\ \mathrm{GPa}}
$$

이고, first-principles [001] reference 약 $10.63$ GPa와 같은 scale에 있다.

이는 active normal-mechanics backbone에 상당히 긍정적인 결과다.

### 같은 계산에서 드러난 정확한 한계

zero pressure의 cubic central pair potential은 Cauchy relation

$$
\boxed{C_{12}=C_{44}}
$$

을 만족한다.

수치 FCC sum도 이를 numerical accuracy까지 복원하여 $C_{44}\approx61.180$ GPa를 예측하지만 실제 Al은 약 $29$ GPa다.

따라서 fixed central pair model은 complete quantitative 3D Al potential이 아니다. 이 한계는 numerical artifact가 아니라 model class의 구조적 한계다.

## 마일스톤 4 — Absolute cohesive-energy consistency

이제 quantitative thermal probability를 넣기 전에 가장 우선해야 할 energy-model 문제다.

같은 fixed generalized-LJ shape에서

- experimental cohesive energy $E_{\rm coh}\approx3.43$ eV/atom을 맞추면
  $$
  E_{[001]}^{\rm LJ}\approx220.466\ \mathrm{GPa}
  $$
  로 지나치게 단단해지고,
- 올바른 normal modulus를 맞추면
  $$
  E_{\rm coh}^{\rm LJ}\approx0.976\ \mathrm{eV/atom}
  $$
  으로 separation-energy scale이 지나치게 작아진다.

따라서 현재 pair law 하나로 tangent normal mechanics와 absolute separation-energy scale을 동시에 재현할 수 없다.

중심 문제는

$$
\boxed{
\text{성공적인 normal LJ mechanics 유지}
+\text{최소한의 물리적으로 정당화된 cohesive/many-body correction 유도}
}
$$

이다.

어떤 correction이든 다음 조건을 만족해야 한다.

1. microscopic physical origin이 명확해야 한다.
2. fitted fatigue-damage law가 아니라 energy 또는 Hamiltonian을 통해 들어와야 한다.
3. degradation을 흉내내기 위해 cycle에 따라 변하면 안 된다.
4. 새로운 physics가 계산 가능한 correction을 요구하지 않는 한 이미 잘 맞는 normal small-strain limit를 유지해야 한다.
5. independent material input 수를 최소화하고 각각을 분류해야 한다.
6. fatigue-life 계산 전에 cohesion과 normal traction에 대해 먼저 검증해야 한다.

explicitly derived electronic-density / many-body cohesive term은 후보가 될 수 있지만 아직 특정 functional form은 채택하지 않는다.

## 마일스톤 5 — Finite-temperature distribution과 time-scale bridge

Milestone 4에서 relevant absolute energy landscape를 수정하거나 오차범위를 제어하기 전에는 quantitative thermal activation을 넣지 않는다.

$$
\exp\!\left(-\frac{\Delta U}{k_BT}\right)
$$

같은 probability는 $\Delta U$에 지수적으로 민감하다. normal stiffness는 맞지만 cohesive-energy scale이 틀린 effective LJ potential을 quantitative escape rate에 바로 쓰면 안 된다.

energy 문제가 제어된 뒤 다음 과제는

$$
\boxed{
\text{빠른 atomic dynamics}
\rightarrow
\text{exact projected memory / coarse state}
\rightarrow
P(a,t)\text{의 느린 진화}
}
$$

를 임의의 relaxation time 없이 유도하는 것이다.

출발점은 microscopic Liouville dynamics, conditional propagator, projection-derived memory kernel이다. Markovian reduction은 이 수준에서 정당화되어야 한다.

핵심 질문은 다음과 같다.

1. 올바른 finite-temperature microscopic state는 $P(a,t)$와 temperature인가, phase-space density인가, 아니면 더 큰 correlation hierarchy인가?
2. thermal energy 중 어떤 부분이 normal-spacing distribution을 바꾸고 어떤 부분이 kinetic energy로 남는가?
3. microscopic spectrum 또는 projected memory kernel에서 atomic/laboratory time-scale separation을 controlled하게 유도할 수 있는가?
4. empirical relaxation time 없이 20 Hz cycle map을 유도할 수 있는가?

## 마일스톤 6 — Normal-opening crack initiation

활성 crack-initiation 그림은 수직 opening 또는 normal mechanical stability loss다.

기존 reduced 1D LJ criterion

$$
\phi''(\lambda_c)=0,
\qquad
\lambda_c\approx1.10777154
$$

는 유용한 reduced-model diagnostic으로 남지만, 3D FCC 결과를 고려하면 최종 normal stability condition은 관련 multidimensional energy Hessian 또는 stable-equilibrium branch의 소멸로 정의해야 한다.

순간적인 tail occupancy

$$
Q_c(t)=\int_{a_c}^{\infty}P(a,t)\,da
$$

와 first-passage initiation을 구분해야 한다.

선호하는 누적 정식화는

$$
\tau_c=\inf\{t:\text{역학적으로 정의된 normal-opening instability 발생}\}
$$

이며

$$
F_{\rm ci}(t)=\Pr(\tau_c\le t)
$$

이다.

## 핵심 closure 질문

- $P(a,t)$만으로 충분한가?
- phase-space lift $F(a,c,t)$가 필요한가?
- 어떤 neighbor-spacing joint density가 필수인가?
- full pair-distance hierarchy $P_k(r,t)$가 필요한가?
- central pair potential의 실패 때문에 필요한 최소 extra energy variable 또는 density variable은 무엇인가?
- 경험적 damping 없이 실험주파수에서 projected memory를 축약할 수 있는가?
- 100 MPa급 수직하중에서 cycle-to-cycle broadening 또는 tail growth를 만드는 실제 물리는 무엇인가?

## 반증 테스트

활성 normal model은 최소한 다음을 만족해야 한다.

1. zero loading에서 인공적인 피로누적이 없어야 한다.
2. 가역 보존계 한계에서 가역응답이 복원되어야 한다.
3. 확률 정규화가 보존되어야 한다.
4. density는 음수가 되면 안 된다.
5. 차원적으로 일관되어야 한다.
6. 에너지 수지를 만족해야 한다.
7. uniform-lattice limit에서 명시된 microscopic energy model이 복원되어야 한다.
8. $P_{N+1}\neq P_N$가 numerical diffusion 때문이면 안 된다.
9. 20 Hz를 주장하려면 물리적으로 유도된 time-scale bridge가 있어야 한다.
10. microscopic interaction model 자체를 명시적으로 바꾸는 경우가 아니라면 potential parameter는 고정되어야 한다.
11. lattice-sum 및 finite-difference error는 주장하는 material discrepancy보다 훨씬 작아야 한다.
12. 검증되지 않은 absolute barrier-energy scale로 thermal first-passage prediction을 만들면 안 된다.
