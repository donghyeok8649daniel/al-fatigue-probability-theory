# Si v5 — 원자 에너지부터 확인하는 재료 기반

## 목적과 위치

사용자의 “원자모델링부터” 지시에 따라, v4의 열분포 계산을 늘리기 전에
그 계산에 쓰인 **원자 에너지 자체**를 감사한다. 공통 흐름은 유지한다.

```text
원자 배열·원소·경계 조건
  → 원자 에너지 U와 힘, Hessian
  → 검증할 집단 좌표 q와 주변 원자 y
  → 조건부 자유에너지 G(q,T), 제거된 자유도의 기억
  → 검증된 확률 연산자
  → 개시·생존·구성 응답
```

이 단계는 첫 번째 화살표를 검증한다. v4의 여러 조건부 구조, 미완료
sampling과 경계 충돌은 그대로 남는다. 원자적 시간 보정이나 실제 웨이퍼
강도·피로 수명 인증을 이 정적 감사로 대신하지 않는다.

## 1. 같은 수학적 계약, 명시적으로 다른 에너지

원자 위치를 R, 셀 벡터를 행으로 갖는 행렬을 C라 한다. 현재 두 대조군은
모두 **순수 Si**이며 도펀트·산화막·표면 수소·전하 상태를 포함하지 않는다.

\[
U(R,C)=\sum_i u_i(\mathcal E_i),\qquad
\mathbf F_i=-\nabla_{R_i}U,\qquad H=\nabla_R^2U.
\]

- Original SW: 기존 삼체 angular 에너지를 포함한 v1–v4 기준.
- Ordinary Tersoff: 공식 SiC.tersoff의 **Si–Si–Si** 항, 1989 매개변수.
  1988 Si.tersoff, screened Tersoff, 재적합 후보와 서로 구분한다.
- GAP/DFT: 독립 재료 기준의 후보. 공개 DFT 파일을 대조에 쓴다고 해서
  GAP를 실행했거나 새 DFT 계산을 수행한 것이 되지 않는다.

단위는 Å, eV, eV/Å, eV/Å²다. 총 에너지와 원자당 에너지를 구분한다.
LAMMPS의 per-atom energy는 삼체 항의 분배 규약이 기존 center-site 분해와
다를 수 있다. 합계와 힘을 비교하며, 서로 다른 분해의 site 값을 같은
국소 물리량처럼 맞추지 않는다.

파일 SHA-256으로 매개변수 출처를 고정한다. 표준 Tersoff의 다른
parameterization이나 screened potential을 같은 이름으로 대체하지 않는다.
이 구현은 optional research adapter이며 생산 registry에서 참조하지 않는다.

## 2. 기하·힘·응력 검증

LAMMPS의 restricted triclinic cell B를 직교 회전 Q로 만든다.

\[
B=CQ,\quad R'=RQ,\quad F=F'Q^T,\quad
\sigma=Q\sigma'Q^T,\quad \det Q=1.
\]

왼손 셀·퇴화 셀을 조용히 뒤집지 않는다. 원자 ID 순서를 보존한다.
주기 셀이 cutoff의 두 배보다 작아도 모든 이미지가 계산되어야 한다.
최소 이미지 하나만 사용한 근사로 바꾸지 않는다.

수치 대조:

1. 에너지 방향 미분과 Cartesian 힘의 일치, 세 차분 간격.
2. 전체 회전·병진에서 에너지 보존 및 힘·응력 변환.
3. 주기 셀 복제에서 에너지 extensivity와 같은 힘.
4. \(dU/d\eta=V\sigma:D\), \(F_{deform}=I+\eta D\)의 응력 부호·단위.
5. SW의 독립 direct-triple 구현, Tersoff의 독립 ASE 구현 대조.
6. 에너지 최적화 종료 메시지와 실제 잔류 힘을 별도로 검사.

고정 원자에 걸린 힘도 계산한다. 고정된 위치는 하중 경계이고, 그 반력을
삭제해서 자유 시편의 평형인 것처럼 해석하지 않는다.

## 3. 결정 내부 이완과 탄성

Diamond에는 내부 sublattice 변위 u가 있다. Affine 변형만 가한 곡률과
내부 원자가 이완한 곡률은 다른 관측량이다.

\[
U=U_0+\tfrac12 V\epsilon^T C_{aff}\epsilon
 +\epsilon^T B u+\tfrac12u^TH_{uu}u,
\]
\[
u_*=-H_{uu}^{-1}B^T\epsilon,\qquad
C_{rel}=C_{aff}-V^{-1}BH_{uu}^{-1}B^T.
\]

이완 Hessian이 양수라는 조건이 필요하다. 두 원자 primitive cell에서
독립 내부 변위 세 개를 풀고 그 3×3 Hessian을 검사한다. 전단은 engineering
shear γ로 정의한다. 0 K 자체 격자 평형에서 세 strain 간격을 비교한다.
실온 문헌값을 0 K 계산에 무표기 혼합하여 calibration 성공을 선언하지 않는다.

이 차이는 후속 q 축약에서도 반복된다. 원자 전체를 affine하게 묶은
에너지, 주변 원자를 최소화한 0 K 에너지, finite-T 조건부 자유에너지는
같은 함수가 아니다.

## 4. 기존 재배열 구조의 대조

v4의 state00/04/01/02를 동일 원자 배열에서 각각 평가한다. 이어서
**같은 고정 grip과 같은 선택 bond gap**에서 각 모델로 다시 이완한다.

- 같은 배열의 U 차이: 원자 에너지 함수의 차이를 직접 검사한다.
- 각 모델의 이완 후 차이: 같은 구속에서 서로 다른 정상 상태를 비교한다.
- 서로 다른 모델의 총 절대 에너지 영점은 비교하지 않는다.
- 같은 grip이 각 모델의 같은 응력·에너지 방출률을 뜻하지 않는다.
- 양의 fixed-q bath Hessian은 full-q 안정성·전이 장벽 인증이 아니다.

처음 L-BFGS 종료 시 잔류 힘이 기준을 넘은 기록을 보존한다. 힘 방정식의
Newton 보정과 두 차분 간격의 Hessian으로 별도 검사하며 원래 결과를
성공으로 덮어쓰지 않는다.

## 5. 독립 재료 대조의 판정

공개 DFT의 **동일 원자 배열**에서 에너지와 힘을 평가한다. 힘에는 에너지
영점 문제가 없다. 에너지는 명시된 단일 reference를 기준으로 상대값을
비교하며, 구조 종류마다 별도 offset을 맞춰 오차를 숨기지 않는다.

힘 component RMSE와 vector norm RMS를 구분한다. 큰 구조가 작은 구조보다
더 많은 원자/성분을 갖는 weighting을 공개하고, 구조 종류별 결과도 남긴다.
결정 내부, 표면, 균열 끝, 결함을 전체 평균 하나에 묻지 않는다.

Cambridge 공개 자료는 GAP 학습 자료다. 이를 SW/Tersoff와 대조하는 것은
두 기존 potential의 새 평가이며 **GAP의 독립 held-out 검증은 아니다**.
프레임 수와 구분은 실제 archive에서 확인한다. 논문의 표와 같다고 가정하지
않는다. 준비 온도·교환상관 함수·완화 여부가 파일에 없으면 미표기로 남긴다.

채택에는 필요한 구성·에너지·힘·계면·변형 경로를 함께 만족해야 한다.
두 경험 potential의 불일치는 어느 하나의 정답, 모든 analytic family의
불가능, 또는 새로운 물리 보정의 성공을 뜻하지 않는다.

## 출처와 재현

- [Stillinger–Weber, 1985](https://doi.org/10.1103/PhysRevB.31.5262)
- [Tersoff, 1989](https://doi.org/10.1103/PhysRevB.39.5566),
  [공식 매개변수](https://raw.githubusercontent.com/lammps/lammps/stable_22Jul2025_update4/potentials/SiC.tersoff)
- [LAMMPS Tersoff 식](https://docs.lammps.org/pair_tersoff.html),
  [ASE 독립 구현](https://ase-lib.org/ase/calculators/tersoff.html)
- [Bartók et al., 2018](https://doi.org/10.1103/PhysRevX.8.041048),
  [원본 모델·DFT 자료](https://doi.org/10.17863/CAM.65004)

실제 실행·수치·제한·재현 명령은
`results/silicon_atomistic_v5/COMPLETED_SUMMARY.md`에서 확인한다.
