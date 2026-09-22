# Si v7 — 도핑 상태와 공통 확률모델의 연결

## 1. 이번에 실제로 구현한 범위

`silicon_doping_research.py`는 화학적 도펀트 수, 이온화된 도펀트 수,
전자·정공 수를 구분한다. 문헌의 도핑별 cubic 탄성을 v6 이방성 균열장에
연결한다. 같은 도펀트 배치에서 빠른 전하 평형을 소거하는 자유에너지의
gradient/Hessian도 구현했다. 이는 연구용 수치 연결부다.

순수 SW/Tersoff에 B/P/As/Sb 상호작용이나 전자 자유도를 추가한 것이 아니다.
새 DFT/MD/원자 에너지 계산, 도핑된 균열 장벽, 물리 이동도/Hz는 이번에 없다.
생산 Al, Si 실행 gate, UI의 입력 의미를 변경하지 않는다.

## 2. 세 종류의 농도

화학 농도 C_D는 재료 부피당 도펀트 원자 수다. 이온화된 비율 f_D가 알려져야
ND+ 또는 NA-를 계산한다. 전자 n과 정공 p는 별도 상태량이다.

    rho/e = p - n + sum_d z_d f_d C_d
    z_B=-1, z_P=z_As=z_Sb=+1

이는 위 네 원소의 단일 이온화 상태에 대한 bookkeeping이다. 중성 bulk에서
우변은 0이어야 하지만, 공간 전하·외부 보상 상태는 0일 필요가 없다.
`charge_balance_cm3`는 잔차를 반환하고 이를 임의로 없애지 않는다.
`ionized_fraction=None`은 미지이며 1로 바꾸지 않는다. 활성화/이온화는
열이력·온도·농도·보상 도핑을 포함한 별도 전자 상태 모델이나 측정이 필요하다.

예: P 2e19, B 1e19 cm^-3에서 각각 f=.75,.5이면 순 양이온은 1e19 cm^-3다.
전체 화학 도펀트 수 3e19를 순 전자 수 1e19로 바꾸어 원자 배치를 만들면 안 된다.
UI의 기존 농도는 화학적 원자 농도다. 문헌 탄성표의 carrier 농도에 자동 대입하지 않는다.

## 3. 유한 원자 상자가 표현할 수 있는 농도

Diamond cubic 관례 셀에는 8개 site가 있으므로:

    N_Si = 8 / a_lat^3
    x_D = C_D / N_Si
    V_material = N_sites a_lat^3 / 8
    E[N_D] = N_sites x_D

길이 a_lat를 Angstrom으로 주면 N_Si=8e24/a_lat^3 cm^-3다.
독립·균일 치환 배치라는 명시적 가설에서만:

    Pr(N_D=0) = (1-x_D)^N_sites
    Pr(N_D>=1) = -expm1(N_sites log1p(-x_D))
    Var(N_D) = N_sites x_D (1-x_D)

이 확률은 도펀트 존재 확률이며 균열 확률이 아니다. Monte Carlo 없이 정확히
계산한다. 응집·편석·침전·격자간 도펀트가 있으면 이 배치 가설이 달라진다.
여러 종을 함께 배치할 때는 동일 site의 배타성을 갖는 multinomial 구성이
필요하다. 현재 단일 종 진단을 종별 독립 배치 생성기로 쓰지 않는다.

5.431 Angstrom, 7600 sites에서는 한 원자가 6.57108e18 cm^-3다.
1e15 cm^-3는 평균 0.000152182개, 1e19는 1.521819개다. 작은 셀마다
한 원자를 강제로 넣거나 기대값을 반올림해 희박 도핑을 표현하면 농도가 바뀐다.
외부 전자 reservoir의 영향은 국소 상자 안 도펀트가 0개라는 이유로 사라지지 않는다.
희박 도핑은 국소 배치별 반응과 외부 전자 상태를 구분한 조건부 연구가 필요하다.

과잉 전자 수는 Delta N = n_excess V_material이다. 전자 추가를 양수, 제거를
음수로 사용한다. vacuum 포함 simulation box 부피를 사용하지 않는다.
8원자 셀과 1e19 cm^-3이면 Delta N=0.00160191이다. 전자 수의 분수 점유와
원자 종을 분수로 치환하는 것은 별개다.

## 4. 도핑별 탄성과 균열 하중

Jaakkola et al., arXiv:1401.1363의 Tables I/III를 입력 데이터로 보존했다.
7개 B/P/As 시편의 C11,C12,C44 및 온도 계수를 쓴다. 농도는 비저항에서
추정한 carrier 값이다. 표의 nominal label과 실제 평균·범위를 모두 남겼다.

    Cij(T) = Cij(25 C) [1 + aij*1e-6*(T-25) + bij*1e-9*(T-25)^2]

aij는 ppm/K, bij는 ppb/K^2. 지원 범위는 -40~85 C다.
농도 보간은 같은 원소의 Table I 평균 농도 사이에서만 명시적으로 요청한다.
선형 carrier-density 보간은 연구용 가정이며 측정된 constitutive law가 아니다.
외삽, Sb 대체, 화학 농도→carrier 환산을 자동 수행하지 않는다.
독립 cubic 안정성 조건 C44>0, C11-C12>0, C11+2C12>0를 확인한다.

v6와 같은 (111)[11-2], 전면 [1-10], plane-strain 균열장에서:

    K = Y sigma_infinity sqrt(pi c)
    G_release = 1000 H(Cij) K^2    [H: GPa^-1, K: MPa sqrt(m), G: J/m^2]

25 C, B0.6를 비교 기준으로 **Wsep가 같다고 가정한 민감도**는:

    K_G,i / K_G,B0.6 = sqrt(H_B0.6 / H_i)

B0.6는 이미 도핑된 시편이다. 이 비를 intrinsic 대비 변화나 실제 K_Ic로
보고하지 않는다. 실제 도핑 변화에는 Wsep와 전이 장벽도 함께 필요하다.
측정 탄성을 순수 SW의 원자 경계에 붙이면 서로 다른 재료의 hybrid가 되므로,
이번에는 continuum reference 계산까지만 연결하고 원자 statics를 재실행하지 않았다.
실제 결합에서는 동일 doped energy의 탄성 및 계면 에너지와 일관성을 확인해야 한다.

탄성 오차막대 .3/.1/.2 GPa의 8개 corner는 민감도 진단이다.
공분산이 없으므로 이를 joint CI 또는 파괴확률 불확도로 읽지 않는다.

또한 cubic 결정의 방향별 영률은

    1/E(n) = S11 - 2(S11-S12-S44/2)*(n1^2 n2^2+n2^2 n3^2+n3^2 n1^2)
    4/E110 = 1/E100 + 3/E111

를 만족한다. E100/E110/E111 세 값만으로 세 cubic Cij를 독립 식별할 수 없다.
Noda의 영률표를 역산해 유일한 탄성 tensor를 만들어 넣지 않았다.

## 5. 공통 자유에너지와 전하 조건

alpha를 고정된 화학 도펀트 배치, q를 균열/registry 집단 좌표라 하자.
원자 환경과 전자 상태를 포함한 같은 원자 에너지에서 조건부 자유에너지
F_alpha(q,N,T)를 정의한다. N은 중성 화학계 기준의 과잉 전자 수다.
다른 원소 조성을 비교할 때는 중성 기준 전자 수와 화학적 reference를 함께 명시한다.

### 고정 전자 수와 고정 화학 퍼텐셜

닫힌 전자계에서는 F_alpha(q,N,T)를 사용한다. 전자 reservoir와 충분히 빠르게
교환하는 계에서만 고정 mu의 grand potential로 줄일 수 있다:

    Omega_alpha(q,mu,T) = -kBT log sum_N exp[-(F_alpha(q,N,T)-mu N)/kBT]
    w_N = exp[-(F_alpha-mu N)/kBT] / Z
    grad Omega = <grad F>_w
    Hess Omega = <Hess F>_w - Cov_w(grad F,grad F)/kBT
    dOmega/dmu = -<N>_w

구현은 이 discrete charge-branch 식이다. 내부 전자 entropy는 각 F_N에 이미
포함돼야 하며 degeneracy를 다시 곱하지 않는다. 이를 실제 Si에 쓰려면 동일
q/alpha/N 기준의 데이터, 충분한 charge-state 범위, electrostatic correction,
reservoir 조건, 전하 평형 시간 검증이 필요하다. 현재 검증은 합성된 smooth
branch에서의 수학적 derivative 검사이며 실제 Si 전하 상태를 계산한 것이 아니다.
이 합의 N은 서로 다른 정수 전자수 sector다. Fractional-charge periodic DFT의
샘플 점들을 독립적인 정수 전하 상태처럼 합하면 sampling 간격에 따라 가짜 entropy가
생기므로 코드가 거부한다. 연속 평균 전하의 DFT F(N)는 해당 ensemble의 Legendre
변환과 적절한 전기적 경계 검증이 따로 필요하다.

고정 carrier density n을 유지한 채 V(q)가 변하면 N(q)=n V(q)이므로
고정-N 미분과 다르다. dF/dq에는 (dF/dN)n dV/dq가 붙는다. reservoir work와
전기장 에너지를 일관되게 포함하지 않은 채 응력 보정항으로 추가하면 안 된다.
균일 background를 둔 periodic DFT와 자유표면/균열의 전기적 경계도 같지 않다.

### 움직이지 않는 도펀트와 빠른 전자를 구별

도펀트 원자가 관측 시간에 이동하지 않으면 alpha는 quenched 배치다.
각 배치에서 공통 확률밀도 법칙을 사용한다:

    partial_t P_alpha = -div J_alpha
    J_alpha = -M_alpha [P_alpha grad Omega_alpha + kBT grad P_alpha]
    P(q,t) = sum_alpha omega_alpha P_alpha(q,t)

omega_alpha는 제조/배치 통계다. 고정 도펀트 배치를 Boltzmann 합으로
anneal해 -kBT log sum_alpha exp(-F_alpha/kBT) 하나로 치환하면 원자가
확산·재배치할 수 있는 다른 물리를 가정한다. 혼합 density도 일반적으로
하나의 평균 potential과 평균 M의 Markov PDE를 따르지 않는다.
이는 소재별 새 피로 법칙을 만드는 것이 아니라 같은 법칙의 조건을 명확히 하는 것이다.
M_alpha, 시간 단위, fast-electronic 가정은 아직 보정·검증되지 않았다.

## 6. 문헌의 비교 범위와 남은 일

- Noda et al. (2023), DOI 10.1038/s41598-023-42676-z, 부록 S3~S5:
  HSE06 격자·LDA 영률/균일 [111] 인장. 도펀트 원자 없는 전하 주입이다.
  0~5e21 cm^-3의 기존 표를 재사용하며 새 DFT나 실제 균열 장벽으로 부르지 않는다.
- Chen et al. (2025), DOI 10.1063/5.0270955, Tables I/III:
  6종의 연마 Cz 웨이퍼, 각50개, ball-on-ring. 특성강도와 SE를 보존하며
  원시300개 관측값을 확보했다고 하지 않는다. 유의차 미검출은 동등성 증명이 아니다.
  Weibull fit은 문헌 근거일 뿐 공통 확률 generator에 넣지 않는다.
- Xu et al. (2026), DOI 10.1016/j.engfracmech.2025.111721:
  전하 주입 벽개 연구. 이번에는 초록/section snippet만 확인했으며
  정량적인 cleavage 곡선이나 장벽을 추출·적합하지 않았다.

다음 실제 원자 계산은 B/P 치환의 중성 chemical reference와 전하 주입 대조군을
구분하고, bulk/표면/균열 끝에서 같은 전자 조건의 에너지·힘을 비교해야 한다.
농도·위치·셀 크기·전기적 경계·표면 재구성 수렴을 통과한 뒤에만
Omega_alpha(q)와 장벽, 이동도/기억 kernel을 공통 PDE에 연결할 수 있다.
실행한 수치와 Git 상태는 `results/silicon_doping_v7/COMPLETED_SUMMARY.md` 참조.
