# Si 도핑 v8: 독립 검증과 온도 범위 수정

## 범위와 실제 수정

v7(4d1e3b9)의 수학·출처를 재감사했다. Jaakkola et al.의 Section III는
As1.7만 최고 측정 온도가 80 C라고 명시한다. v7은 모든 시편에 85 C를
허용했으므로 이 예외를 놓쳤다. 새 source CSV는 원래 표의 모든 값을
그대로 유지하고 시편별 온도 범위 열만 추가한다. 보간에는 양의 가중치를
갖는 시편들의 범위 교집합을 적용한다. 정확히 As2.5를 선택한 경우에는
가중치 0인 As1.7 때문에 85 C가 거부되지 않는다.

v7 결과는 역사적 기록으로 보존했다. As1.7의 85 C 행은 실측 범위 안의
결과로 해석하면 안 된다. v8은 이 행을 80 C로 다시 계산했다. 25 C 결과와
같은 Wsep를 가정한 K비(.989609~1)는 그대로다. 원래 v7 재현·해시 검사는
그 커밋의 코드에서 실행한다. 현재 코드와 과거 v7 manifest가 다른 것은
수정을 숨기지 않기 위한 정상적인 차이다.

출처: A. Jaakkola et al., *Determination of Doping and Temperature Dependent
Elastic Constants of Degenerately Doped Silicon from MEMS Resonators* (2014),
[arXiv:1401.1363, Section III, Tables I/III, Eq.8](https://arxiv.org/html/1401.1363).

## 독립 탄성 풀이

v6의 Airy compliance quartic과 독립적으로 변위 평형에서 출발했다.
고정된 (111)[11-2], front [1-10] plane-strain 프레임에서 cubic tensor를
성분별로 구성·회전한다. Fourier mode u=A exp[i k(x+p y)]에 대해:

    Q_ik=C_i1k1, R_ik=C_i1k2, T_ik=C_i2k2
    [Q + p(R+R^T) + p^2 T] A = 0
    B = (R^T+pT) A
    surface traction = -sigma_i2 = -i k B amplitude
    Y = i A B^(-1),   H = Re(Y_yy)/2

Im(p)>0인 두 모드의 Schur invariant subspace를 사용한다. 정확히 등방성인
중근도 고유벡터 두 개를 무리하게 분리하지 않고 처리한다. Y의 Hermitian
대칭·양의 에너지 및 불변 부분공간 잔차를 확인한다. 등방성 한계는
H=(1-nu^2)/E와 대조한다.

실측 fit 21조건과 c0 민감도 corner 56조건, 총77조건에서 Airy 풀이와
H의 최대 상대차는 2.89e-15다. 이는 같은 선형 탄성 경계 문제의 두 수치
표현이 일치함을 뜻한다. 문헌 탄성 오차나 실제 균열 재료 오차가 그만큼
작다는 뜻은 아니다. 기존 J contour42개의 검사도 별도로 유지한다.

## 전하 자유에너지의 독립 기준

두 개의 독립적인 spin-resolved orbital을 갖는 합성 전자계를 사용한다.
각 orbital은 비점유 또는 점유 가능하므로 N=0,1,2 canonical sector가 있다.
N=1 자유에너지에는 두 내부 상태의 entropy가 이미 포함된다.

    F0=0
    F1=-kT log[exp(-e1/kT)+exp(-e2/kT)]
    F2=e1+e2
    Z_grand = product_i [1+exp((mu-e_i)/kT)]

v7의 sector 합과 독립적인 Fermi product를27조건에서 비교했다.
Omega, gradient, Hessian, 평균 N 및 d<N>/dmu=Var(N)/kT가 일치한다.
Hessian 최대 절대차2.56e-15는 합성 기준계의 대수적 검사다.
실제 Si의 전자 상태, screening, charged surface 에너지를 보정한 것이 아니다.

### 고정 N, 고정 평균 N, 고정 mu

정확한 canonical N sector와 고정 평균 전하의 Legendre 조건은 유한 계에서
같지 않다. 평균 nbar를 유지하도록 mu(q)를 조정한
A(q,nbar)=Omega(q,mu(q))+mu(q)nbar는:

    Hess A = Hess Omega + (grad_q <N>)(grad_q <N>)^T / (d<N>/dmu)

를 따른다. 우변의 추가항은 양의 준정부호다. 고정 mu보다 항상 안정하다는
전체 물리적 선언은 아니다. 같은 q에서의 전하 제약에 따른 국소 곡률 차이다.
mu를 매 q에서 독립적으로 다시 풀어 central difference로 대조했고,
step 1e-3에서1e-5로 줄일 때 최대 오차8.95e-6에서8.96e-10으로 감소했다.
실제 exact-N와의 등치는 열역학적 한계 등 추가 조건을 필요로 한다.

### 정규화와 charge-sector 절단

두 orbital 에너지가0인 경우 sector degeneracy는1:2:1이다. N=2를 누락해도
남은 가중치 합은1이지만 평균 전하는1에서2/3으로 틀린다. kT=.025 eV에서
grand potential 오차는 .025 log(4/3)=.00719205 eV다. 정규화 검사만으로
sector 범위 수렴을 승인할 수 없다. 실제 자료에는 누락 sector 에너지나
수렴 대조가 필요하며, 현재 함수가 이를 자동으로 인증한다고 주장하지 않는다.

## 문헌 영률의 제약과 비식별성

Noda et al. (2023), DOI [10.1038/s41598-023-42676-z](https://www.nature.com/articles/s41598-023-42676-z),
Tables S3/S4의 12행은 0.1 GPa 반올림 구간 안에서 cubic 영률 항등식
4/E110=1/E100+3/E111과 양립한다. 이 검사는 표 내부의 대수적 호환성으로,
DFT의 실제 재료 정확도를 인증하거나 그림 판독을 대체하지 않는다.

S11과 Q=S11-S12-S44/2를 고정하고 S12를 바꾼 세 개의 안정 cubic tensor는
세 방향 영률이 완전히 같지만 H가 .00551361, .00568863, .00580822 GPa^-1로
다르다. 이는 세 영률만으로 full tensor나 crack coefficient를 유일하게
결정할 수 없다는 구성적 반례다. 임의 세 tensor를 Si 보정값으로 채택하지 않는다.

## 재현과 상태

    python -m results.silicon_wafer_feasibility.audit_doping_v8 --output NEW_DIR
    python -m results.silicon_wafer_feasibility.run_doping_audit --output NEW_DIR_CONTINUUM

새 출력 폴더를 사용한다. 기존 결과를 덮어쓰지 않는다. 숫자와 source SHA256,
두 독립 재실행의 일치, 테스트 결과는 results/silicon_doping_v8에 보존한다.
현재 stage는 출처 범위·수학 구현 검증이다. 새 DFT/MD/원자 에너지 계산0회,
도핑된 Wsep/균열 장벽/이동도 및 physical Hz 미보정이다. 기존 Al/Si 생산
모델이나 UI 소재 선택 의미를 변경하지 않았다.
