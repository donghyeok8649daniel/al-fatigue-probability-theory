# Si v9 — 도핑 배치, 전하 평형과 공통 확률 법칙

## 1. 이번 연구의 질문

Si의 도핑 영향을 농도에 따른 단일 강도 배율로 바꾸지 않는다. 같은 원자 에너지
→ 조건부 자유에너지 → 보존 확률 흐름의 구조에서 **화학 배치**, **전자 수의 제약**,
**전하 평형 속도**를 각각 확인한다. v9는 그 축약의 수학 검증과 실제 공개 B 자료의
감사, 독립 중성 ML 원자모델 검사를 수행한다. 실제 Si 도핑 균열 장벽·확률·Hz의
보정 완료를 뜻하지 않는다. 최종 수치/실패/테스트는
[완료 요약](../results/silicon_doping_v9/COMPLETED_SUMMARY.md)을 따른다.

## 2. 화학 배치와 전하 상태는 다른 변수다

- α: 도펀트 원소, 치환/격자간 위치, 군집, 주변 결함의 **화학 배치**.
- N: 같은 α에서의 정수 초과 전자 수. 양의 N은 전자가 더 많다는 규약이다.
  여기서 N은 기존 확률모델의 차원/인터페이스 수를 뜻하는 N=1 표기와 별개다.
- q: 공통 기계 좌표. 현재 생산 N=1의 q=(a,s)와 실제 Si collective coordinate의
  대응은 별도로 검증해야 한다.
- F_{αN}(q,T): 빠른 원자/전자 내부 자유도를 정리한 canonical 자유에너지.
  같은 내부 엔트로피를 다시 degeneracy로 곱하지 않는다.

α가 실험 시간 동안 움직이지 않으면 α별 확률 법칙을 풀고 준비된 배치 분포로
관측량을 평균한다. 형성에너지가 낮다는 이유만으로 α에 Boltzmann 가중을 주면
도펀트가 계속 재배열되는 별도 가정을 넣게 된다. 농도는 α를 완전히 결정하지 않는다.

## 3. 같은 확률 법칙에 전하 상태를 보존하기

고정 α와 전자 저장고 μ에서 Φ_N=F_N−μN으로 놓으면

```text
∂t P_N = −div J_N + Σ_M K_NM P_M − κ_N P_N
J_N = −M_N [P_N grad Φ_N + kT grad P_N]
```

K는 열별 합이 0인 전하 전이 생성자이고 off-diagonal rate는 음수가 아니다.
각 M_N은 대칭 양의 이동도이며 실제 값은 외부 kinetic 검증이 필요하다.
κ_N은 검증된 first-passage 축약이 있을 때만 쓸 수 있는 sink다. 이번 계산의 κ는
**수학 검사만을 위한 합성 함수**이며 Si 균열률이 아니다.

`silicon_charge_dynamics.py`는 기존 SG Bernoulli 함수로 1차원 reflecting 공간을
이산화한다. Cell probability mass를 시간에 따라 결정론적으로 진화시킨다.
이 검증 모듈은 production registry/PDE를 바꾸거나 Monte Carlo로 대체하지 않는다.

## 4. 빠른 전하 평형에서만 단일 자유에너지로 축약된다

```text
Ω(q) = −kT log Σ_N exp[−Φ_N(q)/kT]
w_N(q) = exp[−(Φ_N−Ω)/kT]
P_N ≈ w_N P,       Mbar = Σ_N w_N M_N
kT grad w_N = −w_N(grad Φ_N−grad Ω)
Σ_N J_N = −Mbar [P grad Ω + kT grad P]
κbar = Σ_N w_N κ_N
```

이 결과는 같은 공통 확률 구조다. 단, 전하 평형이 기계 좌표의 진화보다 충분히
빠르며 사용한 sector가 완전하고 μ·열적 조건이 적절해야 한다. Ω만 계산했다고
그 시간척도 분리가 증명되지는 않는다. 위치에 따라 달라지는 Mbar는 위 **flux**
형식으로 사용한다. Ito drift를 별도로 쓸 때 열적 이동도 기울기를 빠뜨리지 않는다.

### 유한 격자와 연속 극한

S는 전하별 확률을 더하는 연산자, T는 w_N으로 되돌리는 연산자다. ST=I이고,
이번 합성 전하 생성자는 λ(TS−I)다. 같은 유한 격자의 정확한 fast limit는

```text
L_grid = S L_spatial T.
```

이를 Ω와 Mbar를 넣은 SG 행렬과 그대로 같다고 두면 안 된다. 두 행렬은 같은
Gibbs 평형을 보존하지만 유한 격자의 전이율은 다르다. v9는 λ와 공간 간격을
별도 축으로 바꾸고 두 해가 접근하는지 검사한다. 초기 전하가 평형인 경우와
한 sector에만 있는 경우도 구별한다. 흡수 counter를 붙여 생존+흡수 질량을 검사한다.

## 5. 전하가 유한한 속도로 움직이면 기억항이 남는다

j>0에 대해 y_j=P_j−w_j P를 정의한다. P_joint=T P+U y이고, 고정된 계수에서는

```text
P' = A P + B y,       y' = C P + D y
P'(t) = A P(t) + B exp(Dt)y(0)
        + ∫_0^t B exp[D(t−s)] C P(s) ds.
```

Laplace 변환의 정확한 진화 연산자는 A+B(zI−D)⁻¹C다. 이 항을 상수 이동도로
부를 수 없다. `charge_memory_blocks`는 숨은 전하 좌표의 가역 선형 변환을 만들며,
2/3sector의 전체 연립식과 복소 resolvent를 직접 대조한다. 시간에 따라 w가
변하는 경우에는 좌표변환의 시간미분 항이 추가되므로 이 고정계수 공식을 그대로
사용하지 않는다. 관측 고유시간이나 전도도만으로 생산 a/s의 M을 정하지 않는다.

화학 배치가 고정된 α 혼합도 일반적으로

```text
Σ_α p_α exp(L_α t) ≠ exp[(Σ_α p_α L_α)t]
```

이다. 중간 시점에 숨은 배치를 새로 추첨한 것처럼 평균 연산자를 반복 적용하면
원래 준비된 시편의 기억을 지운다. v9의 deterministic 반례가 이를 정량화한다.

## 6. 고립 웨이퍼와 전자 저장고의 전하 제약

μ가 외부에서 고정된 grand canonical 계, **전체 평균** 전자 수가 고정된
Legendre 계, **정확한 전체** 전자 수가 고정된 canonical 계는 유한 계에서 다르다.
전체 전하 고정은 각 mesh/cell의 평균 전하를 따로 고정하라는 뜻도 아니다.

서로 독립인 전자 site를 가정해도 정확한 전체 N을 고정하면

(이 예제의 site는0/1점유를 갖는 스핀 분해 전자 궤도 상태이며 Si 원자 site가 아니다.)

```text
Z_N(q) = Σ_{n_1+...+n_m=N} exp[−Σ_i n_i ε_i(q_i)/kT]
∂i F_N = <n_i> ε_i'
∂ij F_N = δij <n_i> ε_i'' − Cov(n_i,n_j) ε_i' ε_j'/kT.
```

전하 제약이 site 사이의 공분산과 공간 교차 곡률을 만든다. 고정 평균 N의
Legendre Hessian은 grand canonical Hessian에 양의 rank-one 보정을 더하지만,
그것이 위 exact-N 식과 유한 크기에서 일치하는 것은 아니다. v9는 2/4/8/16site,
3온도에서 정확한 조합 합과 독립 유한차분을 비교한다. 전자 상태 공간의 정확한
합이며 Monte Carlo 궤적 집계가 아니다.

실제 웨이퍼에서는 도펀트 이온, 전자·정공, 접촉/전극, Coulomb 에너지와 Poisson
방정식의 경계조건도 맞춰야 한다. 이 합성 무상호작용 검사에는 Coulomb 항이 없다.
고립 웨이퍼 전체의 전하 중성과 모든 위치의 엄격한 국소 전하 중성을 혼동하지 않는다.
이 구현을 self-consistent semiconductor device solver라고 부르지 않는다.

## 7. 공개 B 자료에서 확인한 것

출처: Dhamotharan et al., *Electronic and vibrational properties of interstitial
clusters in degenerately boron-doped silicon*, J. Phys.: Condens. Matter 37 (2025)
465901, [DOI](https://doi.org/10.1088/1361-648X/ae191e),
[공개 데이터](https://doi.org/10.15128/r1mw22v5500), CC BY 4.0.

- CASTEP 입력 구조 7개. 64개 host site에 Si/B를 배치했고 실제 원자 수는64–67개다.
  10.862 Å cube에서 B1개는 약7.8e20/cm³로 저농도 웨이퍼의 직접 셀이 아니다.
- `.cell`은 입력 좌표이며 이완된 DFT 출력이 아니다. 전하·pseudopotential·최종
  force/energy/Hessian 출력이 없으므로 논문의 DFT 재현에 필요한 파일이 모두
  있는 것은 아니다. 모든 공개 입력은 FIX_ALL_CELL TRUE다. 논문에 설명된 일부
  가변 셀 계산과 입력 archive의 관계는 별도 확인이 필요하다.
- 같은 neutral 3B 조성의 논문 형성에너지 차이는 configuration1 대비1.42/4.51eV다.
  이는 같은 기준이 상쇄되는 정적 상대에너지이며 이동/균열 활성화 장벽이 아니다.
- B acceptor⁻와 neutral interstitial을 같은 화학 배치의 전하 sector로 합치지 않는다.
  원자 조성과 배치도 달라진다. 대전 셀의 Fermi reference/보정도 따로 검증해야 한다.
- Hall 원자료는 298K의 반복 측정과10–320K 온도 자료를 포함한다. 반복 측정4회를
  독립 웨이퍼4개로 세지 않는다. Hall 이동도는 전하 수송량이며 기계 M_a/M_s가 아니다.
- SIMS에는 B11/Si29 검출 counts와 시간만 있다. sensitivity factor가 없어 화학 농도,
  이온화율, cluster 비율을 산출할 수 없다. 원본의 #REF! 세 셀은 그대로 보존하고
  raw row로 계산한 요약을 별도 파일에 저장했다.

## 8. 독립 ML 원자모델의 역할과 한계

[MACE-MP-0b3 medium](https://github.com/ACEsuit/mace-foundations/releases/tag/mace_mp_0b3)
(MIT, MPTrj/PBE+U 기반, 정확한 model hash 고정)을 optional isolated 환경에서 실행한다.
새 pretrained potential 평가와 에너지 최소화이며 새 DFT, potential 적합, MD가 아니다.

1. 에너지 미분↔힘·응력, 병진/원자순서와 bulk 평형을 검사한다.
2. B 입력7개를 source의 고정 셀에서 각각 이완하고 단계·잔류 힘·실패를 보존한다.
   대전 acceptor의 문헌 에너지는 독립 charge 입력이 없는 이 모델과 수치 비교하지 않는다.
3. 같은 neutral 3B 조성끼리 상대에너지·구조·국소 안정성을 검사한다. LDA와 PBE 기반의
   방법 차이, 다른 국소 최소로 재배열될 가능성을 숨기지 않는다.
4. 기존 v5와 같은 Cambridge Si 원자배열2,475개에서 새 에너지/힘을 계산해 SW/Tersoff와
   비교한다. 원래 GAP 학습 자료이며 MACE의 학습 중복은 독립 감사하지 않았다.
   PW91/PBE/미표기 source group과 에너지 기준을 섞지 않는다.

MACE 전체 자동미분 Hessian 시도는 메모리 증가로 중단했다. 실패 기록을 남기고
좌표별 중앙 힘 차분으로 바꿨다. 잠재에너지와 구조를 바꾼 것이 아니다.
Translation을 명시적으로 제거한 Cartesian Hessian과 질량 가중 Gamma 모드를 쓴다.
Gamma supercell 모드는 full phonon spectrum이나 finite-T PMF/물리 clock이 아니다.
논문은 Gamma–X 분산을 계산하고 제시된 결함 modes가 nondispersive라고 기술한다.
그 선택된 높은3NB modes를 비교하는 것과 전체 저주파 mode의 안정성 검사는 다르다.
실제로 이번 MACE의3B configuration2는 힘 허용오차를 통과했지만 음의 Cartesian
곡률을 보였다. 원래 상태를 보존하고 음의 mode를 독립 차분으로 확인한 뒤 양쪽
방향으로 재이완한다. 최종 최소점 검증과 원래 stationary state 비교를 구별한다.
이 결과만으로 공개되지 않은 LDA 전체 Hessian의 부호를 추측하지 않는다.

## 9. 실제 Si 모델로 넘어가는 조건

추가 [1차 문헌 감사](../results/silicon_doping_v9/sources/dopant_fracture_source_screening.json)는
같은 B도 균열계에 따라 trapping을 늘리거나 진행을 쉽게 할 수 있음을 보여 준다
([Liu2022](https://doi.org/10.1016/j.microrel.2022.114653)). B 한 원자의 재구성 효과가
균열 통과 시간과 경쟁하는 QM/MM·실험 결과도 있다
([Kermode2013](https://doi.org/10.1038/ncomms3441)). 이 Newtonian crack speed를
overdamped 확률 clock으로 대입하지 않는다. P 농도와 미세보의 피로수명 증가를
보고한 실험([Tao2013](https://doi.org/10.1016/j.microrel.2013.07.074))은 추가 target
후보지만, 확보하지 못한 농도/수명 원자료를 만들거나 그 추세를 피로한도 증명으로
사용하지 않는다. Si-P 2NN-MEAM2016은 논문 후보만 확인했고 source-bound parameter
파일을 확보하지 못했으므로 실행하지 않았다.

현재 부족한 입력은 α별 전하 에너지·힘·계면/균열 경로, 그에 맞는 전하 제약,
finite-T 조건부 자유에너지, 기계 및 전하 kinetic scale이다. 도핑 농도를 알고
있는 것만으로 이 입력들이 정해지지 않는다. source/수치/재료/kinetic gate를
구분하고, 어느 항목이 실패하더라도 임의 강도 배율이나 M 조정으로 맞추지 않는다.
기존 Al/SW/Tersoff/SG와 생산 Si 비활성·physical-Hz 차단은 유지한다.
