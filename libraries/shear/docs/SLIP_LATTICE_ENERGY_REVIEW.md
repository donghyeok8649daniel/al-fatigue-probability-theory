# Audit of the imported slip-lattice energy note (historical source audit)

> **2026-09-02 supersession:** the exact two-row expression is now $W(d,s)$
> inside the active unweighted sum $U_0=\sum_{k\ge1}W(ka,s)$. The earlier
> separated-branch conclusion below is preserved as audit history. EAM/DFT
> remains future validation only.

## Decision

The shifted Epstein--Hurwitz lattice sum and its Poisson--Bessel representation
are retained as an **EXACT / IDENTITY** for the stated two-row geometry. The
corrected scalar registry energy and unwrapped probability dynamics are now an
optional active ideal-slip branch. They remain separate from the primary
normal-tensile solver. Unsupported mixed patch energy, automatic
irreversibility, and quantitative aluminum-plasticity claims remain rejected.

Source under review:
`research/source/slip_lattice_energy_mn_K_derivation_KR_v3_23pages.pdf`.

## What is accepted

1. **EXACT / IDENTITY -- finite-chain counting.**
   The formula $E_N=\sum_{k=1}^{N-1}(N-k)v(ka)$ and the factor-one-half
   explanation for a bidirectional bulk chain are correct.
2. **EXACT / IDENTITY -- two-row cross interaction per repeat.**
   For one upper reference atom and an infinite commensurate lower row,
   $r_p^2=a^2+(pb+s)^2$ and the sum over $p\in\mathbb Z$ are correct. By
   translational invariance this is cross-row energy per upper atom or per
   row repeat, not the total energy of two infinite rows.
3. **EXACT / IDENTITY -- special-function transformation.**
   Mellin transformation, Poisson summation, and the resulting modified
   Bessel-$K$ series have the correct powers, factor four, and convergence
   condition $q>1$. The integral representation agrees with
   [NIST DLMF 10.32](https://dlmf.nist.gov/10.32).
4. **EXACT / IDENTITY -- periodicity and registry force.**
   $W(a,s+b)=W(a,s)$ and differentiation with respect to $s$ gives the
   conservative per-repeat registry force. Division by a physically defined
   interfacial area per repeat converts it to traction.
5. **CONTROLLED APPROXIMATION -- fixed-load overdamped first passage.**
   The double-integral MFPT formula is correct for a one-dimensional,
   constant-mobility Smoluchowski process with the stated reflecting and
   absorbing boundaries. Its inverse becomes a time-local rate only after a
   quasi-stationary/renewal separation has been demonstrated.
6. **EXACT / IDENTITY -- EAM structure.**
   The standard EAM ordering, first summing neighbor density and then applying
   the nonlinear embedding function, is correct. The original method is
   described by [Daw and Baskes (1984)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.29.6443).

## Corrections governing the active reduced branch

| ID | Source claim | Finding | Conservative correction |
|---|---|---|---|
| C1 | $\varepsilon_{\rm LJ}$ is changed to the pair-well depth by multiplying the pair law by $C_{m,n}$. | Algebraically valid, but it conflicts with the active repository convention $v=\varepsilon_c[(\sigma/r)^m-(\sigma/r)^n]$. | Keep $\varepsilon_c$ as the coefficient used by the active model. If a well-depth parameter $\varepsilon_w$ is desired, introduce a different symbol and state $\varepsilon_c=C_{m,n}\varepsilon_w$. Never call both quantities $\varepsilon_{\rm LJ}$. |
| C2 | $W$ is called the energy of two infinite rows. | The total energy of two infinite rows is infinite. The displayed sum fixes one upper atom. | Call $W$ the cross-row energy per upper atom/per commensurate repeat. |
| C3 | An arbitrary chosen $s_0$ is called a perfect registry and $\Delta W\ge0$ is implied. | A central-force two-row model does not guarantee that a chosen $s_0$, especially $s_0=0$, minimizes $W$ at the specified $a$. | Define $s_0(a)$ by a verified local or global minimum, or retain the signed difference without calling it a fault energy. |
| C4 | $N_AU_\infty(a)+A_0\gamma(a,s)$ automatically avoids double counting. | $U_\infty$ is a collinear homogeneous-chain energy, whereas $W$ is a cross-row geometry. No common atomistic Hamiltonian or disjoint pair partition was supplied. | Do not add them. Keep the active normal energy and the archived two-row interface energy as alternative geometries. A future coupled model must begin from one FCC half-space Hamiltonian and define its bulk/interface subtraction explicitly. |
| C5 | $N_A=A_0/A_{\rm at}$ is used without an activated crystal plane. | $A_0$ is the mechanical representative area. $A_{\rm at}$ depends on orientation and interface construction; neither equals a correlation area or FEM element area. | Keep $N_A$, $A_0$, and $A_{\rm at}$ independent symbols until a crystallographic representative patch is declared and counted. |
| C6 | $\int_{a_b}^\infty P\,da$ is proposed as an intact tail while $a_b$ is also absorbing. | An absorbing intact density has no mass outside its boundary, and the dead-load full-domain density is not normalizable. | In a reflecting diagnostic use $S^{-1}\int_{a_c}^{a_b}\rho\,da$. In an absorbing model use outgoing flux, survival, hazard, and cumulative initiation. The active project currently absorbs at $a_c$ by its stated operational definition. |
| C7 | $\operatorname{Var}[a]=\iint(a-\bar a)^2P$ is used when $\int P=S<1$. | This is not the variance of the surviving population because $\bar a$ is then an unnormalized raw moment. | Use $M_k=\iint a^k\rho$; conditional mean $M_1/S$; conditional variance $M_2/S-(M_1/S)^2$, including $\dot S$ terms. |
| C8 | The projected GLE is written with independent scalar kernels and the bare lattice potential. | A projection generally produces a matrix, possibly coordinate-dependent memory kernel, cross-correlated noise, and a potential of mean force. Diagonal constant kernels are extra assumptions. | Label diagonal constant mass/friction and the identification of the PMF with the proposed energy as explicit assumptions to be tested by MD. Matrix memory kernels are standard in Mori--Zwanzig reductions; see the [GLE/PMF derivation](https://refubium.fu-berlin.de/bitstream/handle/fub188/35751/PhysRevE.105.054138.pdf?sequence=1). |
| C9 | $\omega_Lm/\Gamma\ll1$ and $\omega_L\tau_K\ll1$ are presented as sufficient overdamped tests. | They compare fast scales only with the loading period. Momentum and memory must also relax faster than the resolved intrawell/configurational evolution; local curvature matters. | Test $m/\Gamma\ll\tau_q$, $\tau_K\ll\tau_q$, and loading-scale separation. Near an instability, recheck the separation rather than assuming it. |
| C10 | Every change in the unwrapped well index $z$ is called irreversible plastic slip. | The periodic conservative landscape permits forward and backward jumps. Under symmetric potential and symmetric zero-mean driving, net residual drift is not guaranteed. | Treat $z$ as a signed net-translation counter. Plastic residual strain exists only if the final distribution of $z$ remains shifted after unloading/relaxation. Transition rates must satisfy local detailed balance when derived from a thermal bath. |
| C11 | A periodic well alone is said to create permanent dissipative plasticity. | Without a bath, emitted phonons, defects, or another sink, the landscape is conservative. With Langevin friction, loop work is transferred to the eliminated bath. | State the energy destination. A two-row Langevin model can show ideal thermally activated slip and bath-mediated hysteresis; it does not by itself derive dislocation plastic dissipation. |
| C12 | An exponential EAM density kernel is proposed generically. | An EAM parameterization supplies a mutually calibrated pair term, density function, embedding function, cutoff, and units. An isolated invented density kernel is not an aluminum EAM. | Remove the proposed kernel. Import a complete validated Al potential and compute the GSF surface from its stated geometry. DFT/EAM GSF results for Al are compared by [Lu et al. (2000)](https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.62.3099/fulltext). |
| C13 | One two-row scalar registry is described as aluminum plasticity. | It has no FCC half-space, dislocation line, Burgers-vector topology, multiplication, backstress, or hardening. | Label it an ideal single-registry mechanism model. Quantitative single-crystal Al plasticity remains unresolved. |
| C14 | The pasted TeX is treated as a buildable bilingual source. | Its Korean text is mojibake and two `\input` files are absent. | Preserve the PDF as source evidence and rewrite corrected UTF-8 TeX. Do not copy the corrupted translation. |

## Reproducible numerical checks

The archive tests do not calibrate aluminum. They use nondimensional values to
check identities and expose assumptions.

- For $m=12.19$, $n=6$, the conversion is
  $C_{m,n}=3.914858619766638$ and $r_e/\sigma=1.121331023181714$.
  Substitution gives $v(r_e)/\varepsilon_w=-1$ to floating-point precision.
- For $q=2$, the independent closed form

  $$
  \sum_{p\in\mathbb Z}{1\over(p+\delta)^2+\eta^2}
  ={\pi\over\eta}{\sinh(2\pi\eta)\over
  \cosh(2\pi\eta)-\cos(2\pi\delta)}
  $$

  verifies the shifted direct sum. Raising the symmetric numerical half-width
  from 100 to 1000 reduces the observed truncation error by approximately a
  factor of ten, as expected for the $q=2$ tail.
- In the illustrative case $\sigma/b=1$, the evaluated two-row energy has its
  lower sampled registry at $\delta=1/2$ for $a/b=1.0$, but at $\delta=0$ for
  $a/b=1.5$. This is a theorem-check example, not an Al parameter result. It
  directly falsifies treating one unverified fixed $s_0$ as the minimum for
  every separation.

The checks are implemented in `libraries/shear/tests/test_two_row_registry.py`.

## Energy and irreversibility interpretation

A fixed multiwell energy makes metastable states possible, but does not by
itself dissipate energy. For a closed conservative model, cycle work equals the
change in kinetic plus potential energy. A positive steady loop requires a
specified sink or irreversible state change. In the archived Langevin model,
friction transfers energy to unresolved atomic/phonon modes and fluctuation--
dissipation fixes the associated noise. If both bath transfer and irreversible
microstructure are excluded, the steady cycle loop must vanish unless escape
or fracture occurs.

The unwrapped registry can be useful without overclaiming it:

$$
s=zb+\widetilde s,
\qquad z\in\mathbb Z,
\qquad 0\le\widetilde s<b.
$$

This is a kinematic bookkeeping identity. It becomes a plastic-strain model
only after the interface geometry, homogenization thickness, transition law,
and persistence of the shifted $z$ population have been validated.

## Repository integration decision

- **Primary active normal branch:** unchanged; uses only $a$ and $P(a,t)$.
- **Optional active exact mathematics:** the two-row direct sum and
  Poisson--Bessel identity drive `theory/registry_lattice.py`.
- **Optional active controlled kinetics:** the unwrapped scalar Smoluchowski
  registry dynamics in `theory/registry_plasticity.py` use an explicit
  isothermal bath and constant-mobility overdamped assumption.
- **Rejected as active physics:** the mixed patch energy, automatic
  irreversibility of $z$, invented EAM kernel, and quantitative-Al claim.
- **Required before quantitative aluminum use:** full FCC plane/direction,
  validated complete EAM or DFT GSF surface, MD-derived memory/mobility,
  homogenization thickness, hardening model, and comparison to single-crystal
  cyclic data.

---

# 반입 slip 격자에너지 노트 검토

## 판정

이동된 Epstein--Hurwitz 격자합과 Poisson--Bessel 표현은 명시된 두 원자열
기하에서 **EXACT / IDENTITY**로 보존한다. 소성, 대표영역 에너지, 확률
tail 및 비가역 전이에 관한 주장은 수정이 필요하다. 결과는 archive로
보존하며 normal tensile solver에서는 활성화하지 않는다.

검토 원본은
`research/source/slip_lattice_energy_mn_K_derivation_KR_v3_23pages.pdf`다.

## 채택하는 내용

1. **EXACT / IDENTITY -- 유한 원자열 계수법.**
   $E_N=\sum_{k=1}^{N-1}(N-k)v(ka)$와 bulk 양방향 합의 $1/2$ 설명은
   정확하다.
2. **EXACT / IDENTITY -- 반복단위당 두 원자열 교차상호작용.**
   위쪽 기준원자 하나와 commensurate한 무한 아래 원자열에 대해
   $r_p^2=a^2+(pb+s)^2$와 $p\in\mathbb Z$ 합은 정확하다. 평행이동
   대칭성 때문에 이는 위쪽 원자 하나 또는 원자열 반복단위 하나당
   교차에너지이며, 두 무한 원자열의 무한한 총에너지가 아니다.
3. **EXACT / IDENTITY -- 특수함수 변환.**
   Mellin 변환, Poisson 합공식 및 수정 Bessel-$K$ 급수의 거듭제곱,
   계수 4, 수렴조건 $q>1$은 정확하다. 적분표현은
   [NIST DLMF 10.32](https://dlmf.nist.gov/10.32)와 일치한다.
4. **EXACT / IDENTITY -- 주기성과 registry force.**
   $W(a,s+b)=W(a,s)$이며 $s$ 미분은 반복단위당 보존적 registry force를
   준다. 물리적으로 정의된 반복단위 계면면적으로 나눌 때만 traction이
   된다.
5. **CONTROLLED APPROXIMATION -- 고정하중 overdamped 최초통과.**
   이중적분 MFPT 식은 명시한 반사·흡수경계를 가진 1차원 constant-
   mobility Smoluchowski 과정에서 정확하다. 그 역수를 시간국소 rate로
   쓰려면 준정상상태 또는 renewal 시간척도 분리가 추가로 입증돼야 한다.
6. **EXACT / IDENTITY -- EAM 구조.**
   이웃 전자밀도를 먼저 합하고 그다음 비선형 embedding function을
   적용하는 표준 EAM 순서는 정확하다. 원래 방법은
   [Daw와 Baskes (1984)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.29.6443)에
   제시돼 있다.

## 활성화 전에 반드시 고칠 사항

| ID | 원본 주장 | 판정 | 보수적 교정 |
|---|---|---|---|
| C1 | $C_{m,n}$을 곱해 $\varepsilon_{\rm LJ}$를 pair well depth로 바꾼다. | 대수적으로 가능하지만 활성 저장소의 $v=\varepsilon_c[(\sigma/r)^m-(\sigma/r)^n]$ convention과 충돌한다. | 활성 모델의 계수는 $\varepsilon_c$로 유지한다. well depth $\varepsilon_w$가 필요하면 다른 기호를 쓰고 $\varepsilon_c=C_{m,n}\varepsilon_w$를 명시한다. 두 값을 같은 $\varepsilon_{\rm LJ}$로 부르지 않는다. |
| C2 | $W$를 두 무한 원자열의 에너지라고 부른다. | 두 무한 원자열의 총에너지는 무한하며 식은 위쪽 원자 하나를 고정한다. | $W$를 위쪽 원자 하나 또는 commensurate 반복단위당 교차에너지라고 부른다. |
| C3 | 임의로 정한 $s_0$를 완전 registry라 부르고 $\Delta W\ge0$를 암시한다. | 중심력 두 원자열은 특정 $a$에서 임의의 $s_0$, 특히 $s_0=0$이 최소임을 보장하지 않는다. | 검증된 최소점으로 $s_0(a)$를 정하거나, 부호 있는 차이로만 두고 fault energy라고 부르지 않는다. |
| C4 | $N_AU_\infty(a)+A_0\gamma(a,s)$가 자동으로 중복계산을 피한다. | $U_\infty$는 동일선상 homogeneous chain이고 $W$는 교차 원자열 기하다. 공통 Hamiltonian과 서로 겹치지 않는 pair 분할이 없다. | 둘을 더하지 않는다. 활성 normal energy와 archive 두 원자열 계면에너지를 대안 기하로 분리한다. 미래 결합모델은 하나의 FCC half-space Hamiltonian에서 시작해 bulk/interface subtraction을 정의해야 한다. |
| C5 | 활성 결정면 없이 $N_A=A_0/A_{\rm at}$를 쓴다. | $A_0$는 mechanical representative area이고 $A_{\rm at}$는 방향과 계면기하에 의존한다. 둘은 correlation area나 FEM element area도 아니다. | 결정학적 patch의 실제 원자수를 세기 전까지 $N_A$, $A_0$, $A_{\rm at}$를 독립 기호로 유지한다. |
| C6 | $a_b$를 흡수경계로 두면서 $\int_{a_b}^\infty P\,da$를 intact tail로 쓴다. | 흡수된 intact density는 경계 밖 질량이 없고 dead-load 전구간 density는 정규화되지 않는다. | 반사 진단에서는 $S^{-1}\int_{a_c}^{a_b}\rho\,da$를 쓴다. 흡수모델에서는 outgoing flux, survival, hazard, cumulative initiation을 쓴다. 현재 활성 연구는 명시한 operational definition에 따라 $a_c$에서 흡수한다. |
| C7 | $\int P=S<1$인데 $\operatorname{Var}[a]=\iint(a-\bar a)^2P$를 쓴다. | $\bar a$가 비정규 raw moment이므로 생존집단의 분산이 아니다. | $M_k=\iint a^k\rho$, 조건부 평균 $M_1/S$, 조건부 분산 $M_2/S-(M_1/S)^2$를 쓰고 $\dot S$ 항을 포함한다. |
| C8 | projected GLE를 서로 독립인 scalar kernel과 bare lattice potential로 쓴다. | 일반적인 projection은 matrix 및 위치의존 memory, 교차상관 noise, potential of mean force를 만든다. 대각 constant kernel은 추가 가정이다. | 대각 constant mass/friction과 PMF를 제안 에너지와 동일시하는 선택을 MD로 검사할 ASSUMPTION으로 표시한다. Mori--Zwanzig의 matrix memory 예는 [GLE/PMF 유도](https://refubium.fu-berlin.de/bitstream/handle/fub188/35751/PhysRevE.105.054138.pdf?sequence=1)를 따른다. |
| C9 | $\omega_Lm/\Gamma\ll1$, $\omega_L\tau_K\ll1$을 충분한 overdamped 검사로 제시한다. | 빠른 시간척도를 loading period와만 비교한다. momentum과 memory가 intrawell/configurational 변화보다도 빨라야 하며 국소 곡률이 중요하다. | $m/\Gamma\ll\tau_q$, $\tau_K\ll\tau_q$와 loading 분리를 함께 검사하고 instability 근처에서 다시 검증한다. |
| C10 | unwrapped well index $z$의 모든 변화를 비가역 plastic slip이라 부른다. | 주기적 보존 에너지면에는 정방향과 역방향 jump가 모두 있다. 대칭 potential과 대칭 zero-mean 하중만으로 net residual drift는 보장되지 않는다. | $z$를 signed net-translation counter로 쓴다. unloading/relaxation 뒤 $z$ 분포가 이동된 채 남을 때만 잔류 소성변형이 있다. thermal bath에서 유도한 전이율은 local detailed balance를 만족해야 한다. |
| C11 | periodic well 하나가 영구 소산소성을 만든다고 한다. | bath, phonon 방출, 결함 또는 다른 sink가 없으면 에너지면은 보존적이다. Langevin friction이 있으면 loop work는 제거한 bath로 간다. | 에너지의 도착지를 명시한다. 두 원자열 Langevin 모델은 이상적 열활성 slip과 bath-mediated hysteresis를 보일 수 있지만 전위 소성소산을 유도하지는 않는다. |
| C12 | 일반적인 지수형 EAM density kernel을 제안한다. | EAM parameterization은 pair term, density function, embedding function, cutoff와 단위가 함께 보정돼 있다. 단독으로 만든 density kernel은 Al EAM이 아니다. | 제안 kernel을 제거한다. 검증된 완전한 Al potential을 가져와 명시한 기하에서 GSF를 계산한다. Al의 DFT/EAM GSF 비교는 [Lu 등 (2000)](https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.62.3099/fulltext)을 따른다. |
| C13 | 두 원자열 scalar registry 하나를 알루미늄 소성이라 부른다. | FCC half-space, 전위선, Burgers-vector topology, 증식, backstress, hardening이 없다. | ideal single-registry mechanism model로 표시한다. 단결정 Al의 정량 소성은 미해결로 둔다. |
| C14 | 붙여 넣은 TeX를 빌드 가능한 이중언어 원본으로 취급한다. | 한국어가 mojibake이고 두 `\input` 파일도 없다. | PDF를 원본 증거로 보존하고 정상 UTF-8 교정 TeX를 다시 쓴다. 깨진 번역은 복사하지 않는다. |

## 재현 가능한 수치 검산

archive 테스트는 알루미늄 보정이 아니다. 무차원 값을 사용해 항등식을
검사하고 가정을 드러낸다.

- $m=12.19$, $n=6$이면 변환계수는
  $C_{m,n}=3.914858619766638$, $r_e/\sigma=1.121331023181714$다.
  이를 대입하면 floating-point 정밀도에서
  $v(r_e)/\varepsilon_w=-1$이다.
- $q=2$에서는 독립적인 닫힌형

  $$
  \sum_{p\in\mathbb Z}{1\over(p+\delta)^2+\eta^2}
  ={\pi\over\eta}{\sinh(2\pi\eta)\over
  \cosh(2\pi\eta)-\cos(2\pi\delta)}
  $$

  으로 이동된 직접합을 검증했다. 대칭 수치 half-width를 100에서 1000으로
  늘리면 관찰된 truncation error가 약 10분의 1로 감소해 $q=2$ tail의
  예상과 일치한다.
- 설명용 $\sigma/b=1$에서 $a/b=1.0$이면 표본 registry 중
  $\delta=1/2$의 에너지가 낮지만, $a/b=1.5$이면 $\delta=0$이 더 낮다.
  이는 정리 검산용 예이며 Al 물성결과가 아니다. 검증되지 않은 하나의
  고정 $s_0$를 모든 separation의 최소점으로 취급할 수 없음을 직접
  보여준다.

검산은 `libraries/shear/tests/test_two_row_registry.py`에 구현했다.

## 에너지와 비가역성 해석

고정된 multiwell 에너지는 준안정상태를 만들지만 그 자체로 에너지를
소산하지 않는다. 닫힌 보존계에서는 한 cycle의 일이 운동에너지와
위치에너지 변화와 같다. 정상상태에서 양의 loop를 유지하려면 명시된
sink 또는 비가역 상태변화가 필요하다. archive Langevin 모델에서는
마찰이 제거된 원자·phonon mode로 에너지를 보내고 fluctuation--dissipation
관계가 그에 대응하는 noise를 정한다. bath 전달과 비가역 미세구조를 모두
제외하면 escape나 fracture 전 정상 cycle loop는 0이어야 한다.

unwrapped registry는 과장 없이 다음 kinematic bookkeeping으로 쓸 수 있다.

$$
s=zb+\widetilde s,
\qquad z\in\mathbb Z,
\qquad 0\le\widetilde s<b.
$$

이 식만으로는 소성법칙이 아니다. 계면기하, 균질화 두께, 전이법칙 및
이동된 $z$ population의 잔류성이 검증된 뒤에만 plastic strain 모델이 된다.

## 저장소 반영 판정

- **활성 normal branch:** 변경하지 않으며 $a$와 $P(a,t)$만 사용한다.
- **보관할 정확한 수학:** 두 원자열 직접합과 Poisson--Bessel 항등식은
  보존할 수 있다.
- **보관할 controlled kinetics:** Kramers/Smoluchowski registry dynamics는
  bath와 시간척도 가정을 명시할 때만 연구할 수 있다.
- **활성 물리에서 기각:** 혼합 patch energy, $z$의 자동 비가역성, 임의
  EAM kernel, 정량 Al 주장.
- **활성화 전 필수:** 완전한 FCC 면·방향, 검증된 전체 EAM 또는 DFT GSF,
  MD로 얻은 memory/mobility, 잔류 slip 시험, hardening 모델, 단결정 cyclic
  data와의 비교.
