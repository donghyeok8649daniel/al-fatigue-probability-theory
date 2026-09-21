# Si 조건부 비선형 열분포와 동역학 v4

## 연구 범위

Original SW, 순수 Si, 같은 fixed-grip (111) shuffle 구조의 별도 연구다.
기존 Al/SG/PDE, 물리 Hz와 Si 실행 gate는 보존한다. 실제 결과의 완료 상태는
`results/silicon_thermal_v4`의 최종 보고서를 따른다. 이 문서의 수식 구현과
실제 Si 재료 검증은 다른 단계다.

## 1. 공통 측도와 유한 영역

선형 좌표에서 `R = R_ref + B x + D (q-q_ref)`, `B.T B=I`, `B.T D=0`.
q는 한 결합의 실제 normal gap, x는 모든 다른 원자의 orthonormal Cartesian
좌표다. 고정된 `Omega_R = {x : |x_i| < R}`에서

```text
Z_R(q,T) = constant_coarea * integral_Omega_R exp[-U(q,x)/(kBT)] dx
F_R = -kBT log Z_R
F_R'(q) = <partial_q U>_(q,T,R).
```

Box 중심과 B, R은 q에 따라 바꾸지 않는다. 그러므로 위 미분에서 움직이는
경계항이 없다. 모든 원자의 전체 에너지를 쓰며 kBT를 원자수로 나누지 않는다.
고정 grip만 있고 자유 원자의 공간을 무한히 열어 놓은 finite-range 원자계는
원자가 멀리 떨어져도 에너지가 유한해 무한부피 partition을 정규화할 수 없다.
유한 box를 명시하는 것은 임의로 실제 시편의 파손 법칙을 정하는 것이 아니다.
그 크기와 다른 가지의 탐색에 결과가 의존하는지 실제로 검사해야 한다.

한 국소 minimum의 H_bath=A=L L.T에 대해
`x=mu(q)+sqrt(kBT) L^-T u`, `u ~ N(0,I)`를 참고 분포로 둔다.
`Phi = [U-U_min]/kBT-u.u/2`라 쓰면 정확한 유한 영역의 비선형 target은
`exp(-Phi)N(0,I) 1_Omega`다. Gaussian importance weights의 overlap이 나쁘면
단일 지수 평균으로 자유에너지를 인증하지 않는다. Box 밖의 표본도 분모 N에는
남고 target weight만 0이어야 한다.

## 2. 두 독립 표본화 방식

pCN: `u' = sqrt(1-b²) u + b xi`, xi 표준정규. Gaussian reference를 정확히
보존하므로 수락률은 `min(1,exp(Phi(u)-Phi(u')))`. 제안/거절 횟수는 물리 시간이 아니다.

Harmonic-reference HMC는 보조 Hamiltonian
`H=(u.u+p.p)/2+Phi(u)`의 harmonic 부분을 정확한 회전으로 적분하고 Phi의
gradient half-kick과 대칭 합성한다. 마지막에 Metropolis correction을 적용한다.
보조 운동량은 제안마다 새 표준정규수로 준비한다. 허용 영역을 벗어난 경로는
전체 제안을 거절한다. 허용된 경로의 역경로도 허용되며, 거절 표본은 이력에 남긴다.
보조 질량/적분간격/제안 수를 물리 mobility나 fracture lifetime에 대입하지 않는다.

독립 seed 반복, 초기/후기 평균, block 크기, 원시 에너지/힘의 split-Rhat,
표본화 방식과 제안 설정을 비교한다. 같은 basin에서 여러 체인이 일치해도
다른 basin을 놓쳤다면 global equilibrium 인증이 아니다.

## 3. 유한 box에서 정확한 평균 힘 control

`v=A^-1 h`, `h=B.T H D`를 선택하되 다음 identity는 임의 고정 v에 성립한다.
각 face에서 그 normal 성분이 0인 `g_i(x)=v_i[1-(x_i/R)^2]`를 쓰면

```text
integral_Omega div(g exp(-U/kBT)) dx = 0
<g.grad_x U - kBT div(g)> = 0
force_control = partial_q U - [g.grad_x U - kBT div(g)]
div(g) = -2 sum_i v_i x_i/R².
```

무한 영역 identity를 유한 box에 그대로 적용하는 편향을 피했다. 실제 raw force와
controlled force를 둘 다 저장하고, 독립 적분에서 평균 불변을 검증한다. 체인이
전체 target을 탐색했다는 가정은 별도로 검사한다. Control로 줄어든 관측 분산을
기억 효과 감소나 물리 소산 감소로 해석하지 않는다.

## 4. 열분포 준비 뒤의 실제 원자 운동

저장된 비선형 conditional sample에 독립 Maxwell bath 속도를 부여한다.
q를 고정한 orthonormal x에서는 운동에너지가 `m xdot.xdot/2`다.
실제 Si 기준 질량으로 velocity-Verlet, 같은 box에서 elastic reflection을
사용하고 반사 횟수를 기록한다. 외부 thermostat/감쇠는 없다.
원자 MD의 ps는 이 Newton reference에만 해당한다. 생산 overdamped t0와 다르다.

`C_f(t)=<delta(partial_q U)(t) delta(partial_q U)(0)>`를 측정한다.
Harmonic bath에서는 C_f/kBT가 v3의 정확한 kernel이다. 비선형 finite-T에서는
이는 **고정-q force correlation 진단**이며 원래 q dynamics의 정확한 memory와
같다고 자동으로 선언하지 않는다. 관측 window·온도·시편크기·초기 준비를 비교한다.
진동 kernel의 음수를 abs/clip하여 양의 상수 마찰로 만들지 않는다.

## 5. 적합 없는 harmonic bath 압축

v3 kernel의 양의 spectral measure `sum w_j delta(lambda-lambda_j)`에
Lanczos/Gauss quadrature를 적용한다. 양의 node/weight에서

```text
Gamma_k(t) = sum_(j=1..k) w_j cos(sqrt(lambda_j/m) t)
U_k(q,y) = K q²/2 + sum_j lambda_j [y_j+sqrt(w_j/lambda_j)q]²/2.
```

이 counterterm은 static marginal K를 정확히 유지한다. 계수는 같은 Hessian의
projection에서 나오며 damping/경험 피로 매개변수를 fitting하지 않는다.
확장된 (q,y,p_q,p_y)의 Hamiltonian 확률 보존 식은 같은 에너지 기반 구조다.
하지만 rank를 작게 잘랐을 때의 실제 시간응답 오차는 별도의 검증 대상이다.
짧은 창에 맞춘 것을 장시간 closure로 승인하지 않는다.

## 6. 같은 간격에서 여러 구조가 존재할 때

고정 q에서 전체 bath를 이완한 여러 minimum이 확인되면 하나의 `U_min(q)`나
그 Hessian만으로 전체 `F_R(q,T)`를 결정할 수 없다. 여러 basin의 partition은
합산한 뒤 log를 취해야 한다. 서로 다른 reference Hessian은 Metropolis 제안의
preconditioner일 뿐, 같은 q/box의 실제 target `exp(-U/kBT)`를 바꾸지 않는다.
동일 minimum 주변 체인의 split-Rhat가 1에 가까워도 이 합산을 검증한 것은 아니다.

특히 force control의 zero-mean 증명은 **선언된 전체 box 분포**에 대한 것이다.
체인이 한 basin에만 갇혔을 때는 그 내부 경계에서 g의 flux가 0이라는 보장이 없다.
따라서 raw/control 차이, 여러 초기 basin, box 크기 검사 없이 그 적분을
인증된 basin 자유에너지라고 바꾸어 부르지 않는다. 이번 mean-force profile은
이 검사를 동반하는 sampling 진단이며 global PMF 인증은 닫혀 있다.

재배열 변위를 나타내는 정규화된 방향 w를 남기는 후보는 `u=w.T x`다.
Householder 직교변환으로 `x=B_u y+w u`를 구성한다. q는 상대 원자 간격이므로
질량 m/2, u와 y의 질량은 m이며 coarea 상수도 좌표 정의에 맞게 유지한다.
분율 `u/path_length`는 그림의 진행 표시일 뿐 새로운 에너지나 온도 단위가 아니다.

원래 Cartesian bath box를 그대로 옮기면 y의 허용 단면은
`|B_u y+w u|_infinity<R`가 된다. 이 영역은 u에 따라 움직인다. 따라서 u 방향
free-energy 미분에는 경계항이 필요하다. 별도의 고정 y-box를 도입하고 원래
ensemble과 같다고 선언하지 않는다. 현재 collective 구현은 정적 제약 이완과
기하학 검사이며, 이 두 좌표의 전체 finite-T PMF 구현/검증은 후속 작업이다.

## 7. 경로·안장점·전이율의 구분

제약 최소를 한 방향으로 추적하다가 다른 basin으로 뛰면 곡선의 최고값을 연속
전이 경로의 장벽으로 사용할 수 없다. 정지점의 전체 **고정-q bath Hessian**을
검사하며 음의 고유값이 하나인 index-one saddle와 양쪽 downhill minimum의
연결성을 따로 확인한다. 고정 q의 index와 q까지 풀었을 때의 index는 다르다.

독립 원자 좌표의 NEB는 실제 SW 에너지/힘으로 경로를 찾는다. 경로 스프링,
FIRE의 가상 시간, 이동 제한은 수치 최적화 설정이며 물질 에너지/강도/동역학에
더하지 않는다. 최고 image를 climbing image로 바꾸기 전에 일반 band의 힘 잔차를
검사한다. 반복 횟수 소진, 경로 전체의 수렴, polish된 정지점의 수렴을 구별한다.
국소 saddle 하나를 찾았다고 전역 최저 경로나 실제 전이율을 확정하지 않는다.

## 8. 조화 양자 보정의 민감도

안정한 조건부 oscillator들의 동일한 개수에 대해서만

```text
F_bath,Q = sum_j [hbar omega_j/2 + kBT log(1-exp(-hbar omega_j/kBT))]
Delta F_Q = Delta U + F_bath,Q(q)-F_bath,Q(q_ref)
```

를 계산한다. 고온 극한은 같은 차원의 classical determinant 차이로 수렴한다.
0 K에서는 zero-point energy 차이만 남는다. 이 계산은 q를 classical 외부
좌표로 고정한 국소 조화 민감도다. 터널링, 비조화 양자 운동, 전자/도핑 효과,
실제 Si 재료 보정을 포함하지 않으며 finite-T 정지점도 새로 찾은 것이 아니다.

## 9. 지배방정식과의 연결 범위

에너지와 확률 보존이라는 공통 출발점은 유지한다. 검증된 extended harmonic
bath에는 Hamiltonian Liouville 방정식을 쓸 수 있고, q만 남기면 memory와
그에 대응하는 열적 힘이 나타난다. 온도별 PMF와 memory의 Markov 한계가
검증되기 전에는 이 연구의 원자 질량 시간을 생산 SG 방정식의 M/t0로 대입하지
않는다. 작은 memory rank, 하나의 q, 양의 상수 drag는 각각 검증할 축약 가정이다.
현재 결과는 이 가정을 시험한 것이며 기존 SG 확률 법칙을 폐기하거나 경험적
피로식·임의 파손 threshold를 추가한 것이 아니다.

## 10. 1차 참고 자료

- Cotter, Roberts, Stuart, White (2013), *MCMC Methods for Functions: Modifying
  Old Algorithms to Make Them Faster*, Statistical Science 28, 424–446,
  [arXiv:1202.0709](https://arxiv.org/abs/1202.0709),
  [DOI 10.1214/13-STS421](https://doi.org/10.1214/13-STS421).
  Gaussian reference를 보존하는 proposal의 원리. 이번 finite-dimensional target의
  normalization/경계 control은 위 식과 별도 수치 적분으로 검사한다.
- Duane, Kennedy, Pendleton, Roweth (1987), *Hybrid Monte Carlo*, Physics
  Letters B195, 216–222, [DOI](https://doi.org/10.1016/0370-2693(87)91197-X).
  Publisher abstract에서 Metropolis correction으로 적분 편향을 제거하는 원리를
  확인했다. 본문을 읽었다고 주장하지 않는다.
- Fernandez-Torre, Albaret, De Vita (2010), *Role of Surface Reconstructions
  in (111) Silicon Fracture*, PRL105,185502,
  [DOI](https://doi.org/10.1103/PhysRevLett.105.185502).
  QM/classical 연구가 다룬 실제 표면 재구성을 이번 SW 재배열과 동일시하지 않는다.
- Bernstein, Hess (2003), *Lattice Trapping Barriers to Brittle Fracture*,
  PRL91,025501, [DOI](https://doi.org/10.1103/PhysRevLett.91.025501).
  Empirical reference의 계산 일치만으로 실제 Si fracture를 인증하지 않는 근거.
- Henkelman, Jonsson (2000), *Improved tangent estimate in the nudged elastic
  band method for finding minimum energy paths and saddle points*, JCP113,9978,
  [primary PDF](https://henkelmanlab.org/pubs/henkelman00_9978.pdf),
  DOI 10.1063/1.1323224.
- Henkelman, Uberuaga, Jonsson (2000), *A climbing image nudged elastic band
  method for finding saddle points and minimum energy paths*, JCP113,9901,
  [primary PDF](https://hj.hi.is/papers/paperCI-NEB.pdf), DOI 10.1063/1.1329672.

확인일 2026-09-21. 결과의 primary source는 이 저장소의 실제 원시 계산이다.
