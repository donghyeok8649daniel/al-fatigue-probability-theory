# Si 국소 좌표의 조건부 측도와 원자 동역학 — v3

## 범위

v2의 같은 순수 Si/original-SW/고정 grip/(111) shuffle 균열을 사용한다.
이번에는 조건부 Gaussian 자유에너지 후보와 주변 원자를 제거할 때 남는
기억 효과를 유도하고, **새 보존적 원자 운동 13개/합계13 ps**로 대조했다.
실제 결과는 [완료 보고서](../results/silicon_conditional_v3/COMPLETED_SUMMARY.md)다.

공통 확률 보존 구조를 유지한다. 기존 Al 에너지, a/s SG/PDE, Si UI 실행 차단,
이동도와 물리 Hz gate는 그대로다. 원자 질량은 아래 Newton 운동에만 쓰며
생산 overdamped 시간에 질량 기반 clock을 대입하지 않는다.

## 1. 선형 좌표의 조건부 측도

고정 경계를 제외한 Cartesian 원자 변위를 x, 국소 쌍의 실제 normal gap과
선택적 slip을 q=Cx로 쓴다. v2의 orthonormal bath B와 최소 lift D는

```text
x = B z + D q,       B^T B=I, B^T D=0, C B=0, C D=I.
D = C^T (C C^T)^(-1).
dx delta(Cx-q) = [det(C C^T)]^(-1/2) dz.
```

두 원자의 상대 gap 한 개는 D의 두 항이 -1/2,+1/2이므로 D^T D=1/2다.
두 직교 gap/slip도 각 대각이1/2다. 이 상수 coarea factor는 **같은 좌표군의
q별 차이에서만** 소거한다. bath B는 실제 Cartesian 길이의 정규직교 기저다.
온도는 총 에너지에 대한 kBT이며 원자수/전면 길이로 나누지 않는다.

추적한 하나의 조건부 최소 가지 z*(q) 근방에서

```text
U(q,z) ~= U_min(q) + (z-z*)^T A(q) (z-z*)/2,    A = B^T H B > 0
Delta F_G(q,T) = Delta U_min(q) + (kBT/2) log[det A(q)/det A(q_ref)].
```

`conditional_harmonic`은 Cholesky 분해로 양정성을 확인하고 log determinant를
계산한다. 음/영 고유값을 절댓값·clipping·regularization으로 보정하지 않는다.
이 결과는 단일 가지의 **고전 Gaussian 근사**다. 전 조건부 분포, 가지 간 합,
비조화 결합항, 재구성, quantum 진동 보정을 계산한 PMF와 구별한다.

0 K 정지점에 대입한 F_G 차이는 온도별 안장점이나 활성화 자유에너지가 아니다.
q에 대한 determinant의 변화가 정지점 위치도 옮길 수 있다. 이번17개 경로 표본은
곡선 진단이며 유한 온도의 안장점 탐색·q-grid 수렴을 대신하지 않는다.

### opening만 남기는 경우와 opening/slip을 남기는 경우

q=(a,s)의 bath를 먼저 제거한 곡률을 K_as라 하면, s도 적분하는 Gaussian에는
추가 요인이 필요하다. 상대 s와 정규직교 s/sqrt(2)의 차이까지 포함해

```text
log det A_a = log det A_as + log(2 K_ss)
K_a = K_aa - K_as^2/K_ss.
```

세 정지점에서 이 두 축약의 일치를 별도로 검사했다. 자유에너지와 기계적
minimum을 같은 것으로 놓거나, eliminated slip의 열적 적분을 빠뜨리지 않는다.

## 2. 정적 Schur 곡률과 동적 기억은 별개다

한 정지점에서 potential의 quadratic block을

```text
U = (q^T H_qq q + 2 z^T h q + z^T A z)/2
H_qq=D^T H D, h=B^T H D, A=B^T H B
K=H_qq-h^T A^(-1)h
```

로 둔다. 동일 원자 질량 m의 운동에너지는
`(m/2)[zdot^T zdot + qdot^T(D^T D)qdot]`다. 주변 원자 운동을 정확히 풀어
제거하면, conditional equilibrium 준비에 대해 다음 식을 얻는다.

```text
m_q qddot + K q + integral_0^t Gamma(t-u) qdot(u) du = eta(t)
m_q = m D^T D,       Omega^2=A/m
Gamma(t) = h^T A^(-1) cos(Omega t) h
<eta(t) eta(0)^T> = kBT Gamma(t).
```

A v_j=lambda_j v_j, c_j=v_j^T h이면
`Gamma(t)=sum_j c_j c_j^T/lambda_j cos(sqrt(lambda_j/m)t)`다.
이는 이 quadratic Hamiltonian에서 직접 얻은 식이며 Markov/white-noise 가정이 아니다.
각 주파수의 행렬 가중치는 PSD지만 **시간 영역의 Gamma(t)는 음수일 수 있다**.
음수를 버리거나 절댓값 적분으로 마찰계수를 만들지 않는다.

일반 초기조건에는

```text
eta(t) = -h^T[cos(Omega t)(z0+A^(-1)h q0)
                       + Omega^(-1) sin(Omega t) zdot0].
```

따라서 initial slip 항을 임의로 없애면 같은 초기화가 아니다. 이번 응답 비교는
`x0=(D-B A^(-1)h) delta_q`, 초기 속도0으로 준비해서 quadratic noise가0인 경우다.
실제 SW에서는 이것이 **선형 조건부 최소의 접선 준비**이며 유한 변위의 정확한
비선형 조건부 평형이라고 부르지 않는다. +/-진폭 대조로 비선형 오차를 분리했다.

정적 최소 경로의 질량 `m_path=m(D-B A^(-1)h)^T(D-B A^(-1)h)`도 계산한다.
이는 bath가 따라오는 저주파 전개의 관성항에 대응하지만, 일반 주파수의
Gamma(t)를 대체하지 않는다. 두 질량으로 만든 즉시 이완 대조군을 각각 표시했다.

유한 undamped cosine 합은 저주파 상수 마찰의 인증값을 주지 않는다. 실제 유한
온도의 비조화 산란, 큰 bath 극한과 coarse observation time은 별도로 검증해야 한다.
짧은 시간의 관성 진동은 더 긴 시간의 확산적 축약 전체를 부정하는 증거도 아니다.

## 3. 새 원자 운동과 단위 검증

모든 자유 원자는 기존 SpatialSW의 같은 실제 force로 velocity-Verlet 적분한다.
고정 원자는 위치/속도가 변하지 않으며, 그 site 에너지가 free atom force에
기여하는 부분은 계속 포함된다. 힘 clipping, bond 삭제, thermostat, 경험 감쇠가 없다.

기준 질량28.0855 amu와 명시된 eV/Angstrom/ps 환산을 사용한다.
`m[eV ps^2/Angstrom^2] = m[amu] * atomic_mass[kg] * 1e4 / electron_volt[J]`다.
실행 시 scipy.constants에서 읽은 수치를 running_summary에 저장했다.
이 질량은 자연 Si의 기준 선택이며 도핑/동위원소의 실제 분포를 적합한 값이 아니다.

- 같은1152원자/자유672원자와 고정 경계에서 새13실행, 각1 ps.
- 무변위 대조1개; +/-0.004, +/-0.002 Angstrom에 dt=1,0.5,0.25 fs 대조12개.
- potential은 per-site 기준 차분을 합하고 kinetic은 모든 자유 원자의 총량이다.
  에너지 잔차는 실제 매 step에서 검사한다. 관측 저장 간격은1 fs다.
- 2% 에너지 오차 한도는 실행의 숫자 오류를 거부하는 사전 한도일 뿐 물리 검증
  임계값이 아니다. 실제 최대값과 시간 간격에 따른 감소율을 별도로 보고한다.
- full Cartesian harmonic 진화는 독립 eigen expansion, 그 memory 식은 작은
  coupled ODE 및 convolution 적분으로 대조한다. 적분기에는 에너지2차수렴,
  시간역전, 고정경계 시험을 적용한다.

Richardson 조합 `(4 q_dt/2-q_dt)/3`은 O(dt^2) 오차를 추정·제거하는 후처리다.
새 trajectory나 독립 sample로 세지 않는다. 실제 coarse/fine 원시 이력을 보존한다.

## 4. 어떤 좌표를 더 남길 것인가

정적 Gaussian bath에서 일반화 힘 요동 `f=-h^T z`와 후보 관측값 `y=Lz`의
covariance를 계산한다. 설명되는 분산 비율은

```text
R^2 = [h^T A^(-1) L^T (L A^(-1) L^T)^(-1) L A^(-1) h]
      / [h^T A^(-1) h].
```

여기서는 kBT가 소거된다. 후보를 적합 매개변수로 보정하지 않으며 중복·singular
관측값은 거부한다. local slip, 같은 전면의 다른 gap, 다음 전면 gap을 각각 검사한다.
많은 gap을 넣은 경우는 표현력 진단이지 채택할 다차원 PDE 모델이 아니다.
이 비율은 정적 조건부 force variance이며 그만큼 실제 memory duration이나
동역학 오차가 줄었다는 보증은 없다.

## 5. 현재 판정과 다음 단계

Gaussian 측도·정적/동적 축약의 수학적 일치와 실제 작은 변위 동역학은 검증했다.
finite-T joint PMF, 온도별 마찰/Markov 축약, local-cell 확률법칙과 실제 파손 사건,
실측 재료·표면·도핑 에너지 보정은 미완료다.

다음은 조건부 finite-T joint 분포와 선택한 주변 변형 좌표의 응답을 구하고,
canonical 준비 후 보존적 동역학으로 memory의 온도·크기·시간 범위를 검사하는 단계다.
정적 분포를 준비하는 thermostat의 임의 시간상수를 물리 이동도로 사용하지 않는다.
단시간 관성이 필요한 경우 같은 확률 보존 구조에서 운동량/주변 좌표를 남기며,
생산 Smoluchowski/SG 경로로 이동하기 전에 coarse-time 축약 조건을 따로 검증한다.

## 참고한 1차 자료

- R. Zwanzig (1973), *Nonlinear generalized Langevin equations*, Journal of
  Statistical Physics9,215–220, [DOI:10.1007/BF01008729](https://link.springer.com/article/10.1007/BF01008729).
  Publisher abstract의 exact bath reduction/Markov approximation 구분을 확인했다.
  위 block harmonic 식은 이 문서에서 직접 유도한 것이며 논문 본문을 인용한 것이 아니다.
- [LAMMPS units](https://docs.lammps.org/stable/units.html),
  [velocity-Verlet / fix nve](https://docs.lammps.org/fix_nve.html): 원자 기준 단위·적분법.
  이번 새 MD는 자체 적분기다. 과거 v2의 LAMMPS run0 force 대조를 새 MD로 세지 않는다.
- [CIAAW Si standard weights](https://ciaaw.org/atomic-weights.htm),
  [historical values](https://www.ciaaw.org/historical-atomic-weights.htm): 현재 interval과
  명시적 기준 질량28.0855의 역사적 출처를 구별한다.

모두2026-09-21 확인. SW 매개변수는 v2의 hash-bound 원본을 유지한다.
