# Si 국소 균열 좌표: 이상적인 균일 분리 가정의 제거

## 목적과 실제 변경

`silicon_crack_research.py`는 기존 순수 Si original-SW 기준 에너지를 그대로
공간 원자 배열에 적용한다. 기존 강체 계면의 원자면 전체 동시 이동 제약을 제거하고,
초기 균열 끝의 원자들이 독립적으로 이완·분리·재결합할 수 있게 했다.
에너지 계수, 임계응력, 피로 법칙을 맞춰 강도를 낮추는 수정은 없다.

기존 `RigidInterface`는 대조 계산으로 보존한다. 새 계산은 연구 모듈이며
Si 생산 해석, 실측 웨이퍼 강도, 피로 수명, physical seconds/Hz 승인이 아니다.
실행 결과와 검증은 [계산 보고서](../results/silicon_local_crack_v2/COMPLETED_SUMMARY.md)에 있다.

## 1. 같은 원자 에너지, 더 넓은 상태공간

전체 원자 좌표를 R=(R_1,...,R_N), 고정된 외곽 원자를 R_B(lambda)라 한다.
lambda는 지정한 경계 변위장의 진폭이다. 자유 원자 집합 I의 에너지는

```text
U(R_I;lambda) = sum_i [1/2 sum_j phi(r_ij)
                       + sum_(j<k) h(r_ij,r_ik,cos(theta_jik))].
```

중심 i의 합에는 고정 원자도 포함한다. 자유 원자를 움직이면 그 원자와 접한
고정 원자의 3체 에너지도 변한다. 이를 누락하면 힘·Hessian·경계 반력이 틀린다.
균열은 실제 벌어진 좌표로 초기화한다. 끊어진 bond 목록, 제외된 pair/triplet,
opening에 따른 힘 clipping을 쓰지 않는다. 가까워진 원자 간 상호작용은 복구된다.

이웃 목록은 original SW의 실제 cutoff `a*sigma`와 별개인 skin을 둔다.
모든 원자의 이동이 skin/2 미만일 때만 재사용하고, 이 조건을 벗어나면 다시 만든다.
주기 전면의 이웃 이미지는 명시적으로 포함한다. 얇은 주기 cell에서 한 원자의
여러 이미지가 들어갈 수 있다는 점을 유지한다.

### 동일한 환경 모멘트 표현

각 중심에서 단위벡터 n_j, 원래 SW 반경 가중치 w_j를 사용해

```text
rho = sum_j w_j,       v = sum_j w_j n_j,
T = sum_j w_j n_j n_j^T,       d = sum_j w_j^2,
c = cos(theta0),       Q = T - rho I/3

U_3/(lambda_SW epsilon)
 = 1/2 T:T - c v.v + c^2 rho^2/2 - (1-c)^2 d/2
 = 1/2 Q:Q - c v.v + (1/6+c^2/2)rho^2 - (1-c)^2 d/2.
```

반경, 방향, 자기항 제거를 모두 미분한 힘을 쓴다. 직접 3체 합과 별도의
per-site forward 2-jet Hessian 및 외부 LAMMPS와 대조한다. 원래 SW의 유한
cutoff가 정의인 기준 계산이며, Al의 무한 LJ/Bessel 합을 cutoff로 바꾸지 않는다.
이 표현 일치가 Al의 현재 LJ/지수 radial family로 Si 물성을 적합했다는 뜻은 아니다.

## 2. 기하학·하중·정규화

직교 프레임은 x=[11-2]/sqrt(6), y=[111]/sqrt(3), z=[1-10]/sqrt(2)다.
y가 Si(111) shuffle 면의 법선, x가 균열 전진, z가 전면 방향이다.
12원자 직교 diamond cell의 반복으로 시편을 만들고 z만 주기적이다.
x/y 외곽 grip을 고정하고 내부 3성분 변위를 모두 허용한다.

초기 변위는 원자면 전체의 균일 opening이 아니다.

```text
w(x) = [1-tanh((x-x_tip)/ell)]/2
u_y(x,y;lambda) = lambda [(1-w(x))*y/H + w(x)*sign(y)/2].
```

뒤쪽에서는 두 반쪽이 벌어지고 앞쪽은 인장된다. 이 함수는 초기 배치와
명시적 외곽 Dirichlet 조건이며, 정확한 탄성 K_I 해나 실측 wafer 결함이 아니다.
자유 원자는 이 함수에 계속 묶어 두지 않고 에너지로 이완한다.

고정 변위에서 비교하는 에너지는 U 자체다. 외력 일을 다시 빼면 이중 계산이다.
경계 진폭에 대한 가상일은

```text
d U_min/d lambda = sum_(i in B) grad_i U . dR_B,i/dlambda.
```

측면 grip의 변위도 포함하므로 이 반력을 임의 면적으로 나눠 균일 원격응력이나
웨이퍼 파괴강도로 부르지 않는다. 실제 에너지 차분과 위 반력을 대조한다.

단위는 총 eV, Angstrom, eV/Angstrom이다. 총 에너지를 원자수로 나누거나,
단위 셀 장벽에 임의 전면 길이·메시 면적·A_c를 곱하지 않는다. 배경 에너지는
각 site에서 같은 기준 상태를 뺀 뒤 합해 extensive total의 상쇄 오차를 줄인다.

## 3. 국소 좌표와 조건부 이완

전면의 실제 원자 쌍 (i,j)을 골라

```text
a = (R_j-R_i).n,
s = (R_j-R_i).m,       n=y 방향, m=x 방향
q = C R
```

로 정의한다. `RelaxedCoordinates(..., components=(1,0))`는 두 좌표를 정확히
고정하고 나머지를 이완한다. 이번 경로 검사는 a만 고정하고 s와 다른 성분도
모두 이완했다. 원래 강체 registry 좌표나 full-period well index와 같다고
가정하지 않는다. 아직 생산 `n=floor(s/b+1/2)`로 연결하지 않는다.

R=R_ref+Bz+Dq의 선형 제약에서

```text
G_0(q) = U(R_ref+B z_*(q)+Dq),     B^T grad U = 0,
dG_0/dq = D^T grad U,
H_eff = D^T H D - D^T H B (B^T H B)^(-1) B^T H D.
```

여기서 G_0는 추적한 조건부 국소 최소 가지의 **0 K 에너지**다. 전역 최소나
finite-T PMF라고 부르지 않는다. B^T H B>0을 확인하고 독립 힘 차분 및 이완된
에너지·반력 차분으로 Hessian과 Schur 곡률을 확인한다.

균일 전면 대조군에서는 z 반복 사이 같은 원자의 변위만 묶는다. 실제 총
에너지를 그대로 비교한다. 이 대조군의 안장점이 제한된 공간에서 index 1이어도,
전체 원자 공간에 여러 음의 방향이 있으면 국소 전이 안장점으로 채택하지 않는다.
국소 후보에서는 전체 Hessian index 1, 잔여 양의 곡률, 힘 잔차, 불안정 고유벡터의
양방향 하강이 서로 다른 안정 상태에 연결되는지를 모두 확인한다.

## 4. 공통 확률방정식과 연결되는 조건

Al/Si에 별도 피로 법칙을 도입하지 않는다. 공간 좌표와 원자 에너지를 정한 뒤
같은 확률 흐름 구조를 사용한다.

```text
partial_t P = -div J,
J = -M(q) [P grad G(q;lambda,T) + k_B T grad P].
```

다만 위 식에 필요한 자유에너지는 실제 조건부 적분이다.

```text
G(q;lambda,T) = -k_B T log [ integral exp(-U(R;lambda)/k_B T)
                                      delta(CR-q) dR / reference_measure ].
```

선형 C의 coarea/Jacobian 인자는 q에 무관하므로 기준 척도에 포함할 수 있다.
조건부 단일 안정 가지의 저온 Gaussian 근사에서도
`G = G_0 + (k_B T/2) log(det H_zz(q)/det H_ref) + q-independent terms`가 남는다.
여러 가지가 차지하는 확률은 합산해야 하며, 가장 낮은 가지 하나만 고르는
minimization과 finite-T marginalization은 다르다.

또한 static Hessian의 양의 부호만으로 빠른 조건부 이완, Markov 성질,
overdamped 이동도 또는 물리 시간을 얻지 않는다. 원자 좌표의 이동도를 알고
있어도 투영한 q의 기억 효과와 변동을 검증해야 한다. 실제 균열 전면의 인접 결합
상태가 느리면 여러 국소 q를 함께 남겨야 한다. 이는 같은 이론의 상태공간 확장이다.
기존 N=1 SG가 그 고차원 generator까지 검증했다는 뜻은 아니다.

따라서 이번 계산에서는 `physical_activation_free_energy_eV`, 물리 이동도,
seconds/Hz 및 wafer strength를 제공하지 않는다. 0 K 정적 에너지 장벽은 별도
명칭으로 저장한다. 임의 attempt frequency나 열항으로 정적 한계를 숨기지 않는다.

## 5. 파손 사건과 남은 검증

특정 초기 bond의 normal gap이 SW cutoff를 넘는 것은 기하학적 관측량이다.
흡수 질량·균열 확률·불가역 손상으로 즉시 바꾸지 않는다. 다시 닫힐 수 있고,
국소 분리 후에도 전면 나머지는 남을 수 있다. 실제 wafer 파괴 사건에는 국소
전이들의 연결, 균열 성장/치유, 외곽 시편과의 일관된 역학 연결이 필요하다.

실측 결함·표면/산화층·도핑 및 Si 에너지 적합, 온도별 PMF/동역학, 실제 하중과
원자 영역의 경계 매칭은 별도 미완료다. 초기 결함은 초기 상태/기하학 자료이며
강도를 맞추는 임의 감소계수로 대체하지 않는다. 현재 산출값은 지정된 유한 strip의
순수 Si original-SW 기준값이며, 실제 웨이퍼 강도나 피로 수명으로 외삽하지 않는다.

## 출처

- Stillinger & Weber (1985), *Computer simulation of local order in condensed
  phases of silicon*, [DOI:10.1103/PhysRevB.31.5262](https://doi.org/10.1103/PhysRevB.31.5262).
  구현식은 [LAMMPS SW 문서](https://docs.lammps.org/pair_sw.html), 매개변수는 기존
  hash-bound `results/silicon_wafer_feasibility/source_Si.sw`를 사용한다.
- Kermode et al. (2015), *Low Speed Crack Propagation via Kink Formation and Advance
  on the Silicon (110) Cleavage Plane*,
  [DOI:10.1103/PhysRevLett.115.135501](https://doi.org/10.1103/PhysRevLett.115.135501).
  국소 전면 전진과 동시 전진을 구분하는 근거다. 이 논문의 modified-SW/DFT,
  (110) 기하학과 본 계산의 original-SW/(111)을 혼동하지 않으며 논문 수치를 재현한 것이 아니다.
