# 고체의 kBT 적용 가능성: 현재 Si 시편의 열항 감사

2026-09-29 · silicon-wafer-research · 출발 커밋 4bfe1e1.

## 판정

고체의 강한 결합은 kBT를 없애는 근거가 아니다. 구속된 상태공간과 전체 결합
Hessian을 사용한 고전적 평형 계산은 수학적으로 일관된다. 그러나 현재의 국소
조화 근사와 상수 이동도만으로 실제 Si의 열분포·균열 개시 확률·수명을 계산할
물리 근거는 확보하지 못했다. 아래의 통과는 구현 검증이며 재료 채택이 아니다.

검사 대상은 기존 MP-0b3의 360원자 순수 Si 시편이다. 216원자는 자유롭고
144원자는 고정되며, 도펀트·산화층·미리 넣은 내부 균열은 없다. 처음 8%, 복귀 8%,
10% 고정변위의 648차원 Hessian을 재사용했다. 과거 SW 동역학은 별도 대조 자료다.
생산 Al/UI/에너지/초·Hz 설정과 v13 원자료·보고서는 바꾸지 않았다.

## 1. 구속과 에너지 정규화

허용된 직교 Cartesian 좌표 z에서 U=U0+g·z+z^T H z/2로 쓸 때
평균 이동은 -H^-1 g, 고전적 공분산은 kBT H^-1이다. 전체 원자 좌표에 대한
공분산은 B(kBT H^-1)B^T이며 고정 원자의 열적 변위 분산은 정확히 0이었다.
자유원자들을 독립 스프링으로 취급한 것이 아니다. 약한 모드와 상관된 움직임을
포함하므로 H의 대각 성분만 역수로 취한 결과와 다를 수 있다.

300 K의 평균 조화 위치에너지는 648 kBT/2 = 8.37604793 eV다. 이 총 진동 에너지가
하나의 결합 또는 균열 좌표에 집중된다는 뜻은 아니다. 단일 모드 열에너지와
시편 총 진동에너지를 섞지 않는다. 기존 v13의 원자 RMS를 8e-16 A 이내로 재현했다.

동일한 q에서 총에너지 F를 원자당 f=F/N으로 바꿀 때는 열항도 kBT/N,
이동도는 N M으로 함께 표현해야 같은 generator다. 실제 SG 구현으로 N=2,216,360
합성 대조를 수행했고 최대 상대 차이는 7.33e-16이었다. 열항을 그대로 둔 대조군은
generator와 평형 분포가 달라졌다. 이 환산은 실제 온도를 바꾼다는 뜻이 아니다.
이 대조군의 f는 합성 에너지이며 실제 Si의 자유에너지가 아니다.

## 2. 조화 양자 민감도

같은 양의 고유값 lambda와 명시한 Si 기준 질량 28.0855 amu를 사용했다.

    omega = sqrt(lambda / m)
    C_classical = kBT / lambda
    C_quantum = hbar / (2 m omega) coth(hbar omega / (2 kBT))

| 300 K 상태 | 고전 원자 RMS (A) | 양자 조화 RMS (A) | 모드별 최대 분산 비율 | hbar omega > kBT 모드 |
|---|---:|---:|---:|---:|
| 처음 8% | 0.235088 | 0.238052 | 1.36487 | 358/648 |
| 복귀 8% | 0.237623 | 0.240553 | 1.36960 | 359/648 |
| 10% | 0.249861 | 0.252652 | 1.36454 | 348/648 |

전체 원자 RMS 보정은 약 1.1~1.3%지만 일부 모드의 분산은 약 37% 달라진다.
평균 변위가 비슷하다고 모든 모드의 고전 근사가 정확하다고 선언할 수 없다.
10% 오차와 hbar omega/kBT=1은 보고용 구간이며 임의의 물리 승인 기준이 아니다.

같은 grip의 복귀8%-처음8% 진동 자유에너지 차이는 300 K에서 고전 0.262814 eV,
양자 조화 0.296690 eV다. 차이 0.0338759 eV는 1.31 kBT 정도다. 이는 두 국소
우물의 진동 보정 차이이며 활성화 장벽·전체 PMF·균열 확률이 아니다.
이전 고전 국소 Delta F=-9.456409 eV는 이 근사에서 -9.422533 eV가 되지만,
이 값으로 실제 우물 점유율이나 전이율을 채택하지 않는다.

동일한 고전 퍼텐셜에 양자 조화 통계를 적용한 민감도 검사다. 실제 Si의 양자
진동·전자구조·비조화 PMF를 검증한 것은 아니다. 649차원 힘제어 grip 좌표에는
하중 장치의 관성이 별도로 정해지지 않아 이 원자 질량 진동 비교를 적용하지 않았다.

## 3. 실제 전차원 에너지 평가

64개 예정 변위는 두 상태·두 가정 온도(50/300 K)에서 같은 648차원 정규 난수를
사용했다. 각 상태·온도에 8개 난수벡터의 양/음 대조 16점이다. 상태·온도 간
난수는 공통이고 대칭 쌍은 독립 표본이 아니다. 이는 고전 조화 Gaussian의 제안점이며
실제 MACE 열평형 표본이나 물리 MD 궤적이 아니다. 변위를 clipping하지 않았다.

처음 8% /300 K의 16점에서 실제 에너지와 국소 조화 예측의 차이는
-8.54156~+5.19373 kBT였다. 상대 Boltzmann 가중치 exp(-DeltaU/kBT)의 집중도 지표
(sum w)^2/sum(w^2)는 16점 중 2.04279였다. 이 수치는 통계적으로 검증된 독립
표본 수가 아니다. 작은 대조만으로도 단순 조화 제안분포의 가중치가 심하게
치우칠 수 있음을 보여주며, 이를 그대로 실제 열분포로 사용하는 것은 부적절하다.
표본의 basin membership은 확인하지 않았으므로 근방의 비조화성과 basin 이탈의
기여를 여기서 분리하지 않는다. 전이 장벽이나 자유에너지는 추정하지 않았다.

첫 실행은 300.358초 한도에서 37/64점 + 출발 상태 대조 2회, 총 39회 평가 후
시간 한도로 종료했다. partial_replay.json의 427개 독립 산술·해시 대조를 통과했다.
full_dimensional_probes 폴더에 당시 원자료·partial 상태·코드 snapshot을 보존했다.
당시 summary의 independent_draws는 실제로 완성 대칭쌍 수(len(rows)//2)다.
후속 코드에서는 normal_vectors_used와 complete_antithetic_pairs를 구분했다.
최종 재개 결과와 실제 추가 호출 수는 full_dimensional_completed/summary.json 및
COMPLETION.md를 따른다. 이미 계산한 37개 좌표·힘·에너지는 재평가하지 않는다.

## 4. 관성·기억과 이동도

현재 Hessian과 기준 질량에서 얻는 조화 진동 주기는 약 0.073~2.30 ps다.
평균 내부 신장 좌표의 무감쇠 조화 상관은 진동하며 음수가 된다. 이 좌표는
자유원자 위/아래 절반의 평균 z 차이이며 검증된 균열 좌표가 아니다.
이 진동 시간은 생산 확률모델의 t0나 물리 Hz 변환값이 아니다. 감쇠를 넣지 않은
조화 모형의 결과만으로 실제 긴 시간 Markov 근사의 불가능성을 증명하지 않는다.

별도로 과거 SW /1152원자 /고정-gap /300 K /40 ps 4개 기록을 원자료에서 다시
계산했다. 선택한 조건의 모든 seed를 포함했고, 새 MD는 0이다. FFT 상관을 직접
내적의 5개 lag와 대조했으며 최대 차이는 1.56e-15다. 9개 cutoff 적분의 기존
결과와 최대 차이는 2.15e-16 eV ps/A2다.

| 적분 cutoff | 평균 힘 상관 적분 (eV ps/A2) |
|---|---:|
| 0.02 ps | 0.0298303 |
| 0.5 ps | -0.000487147 |
| 1 ps | 0.00600706 |
| 10 ps | 0.00128321 |

10 ps의 seed별 값은 -0.000675032~0.00276252로 0 양쪽에 있다. 이 범위는
신뢰구간이 아니다. 부호를 abs/clip하거나 양의 cutoff만 골라 상수 마찰을 만들지
않는다. 이 nonlinear frozen-q correlation은 정확한 Mori memory kernel과도
동일시할 수 없다. geometry·퍼텐셜이 다른 SW 결과를 현재 MP-0b3 시편의
동역학 보정으로 전이하지 않는다. 현재 시편의 유한온도 상관·외력 응답은 미확보다.

## 5. 실제 검증 및 재현

10개 새 구현 시험이 통과했다(6.65초). 직접 Gaussian 적분, 양자 oscillator의
유한 partition 합, 고온/저온 극한, 좌표 회전 불변성, 고정 grip 및 잘못된
곡률/온도/질량 거부를 검사했다. 관련 코드는 solver_v1/test_silicon_thermal_validity.py다.
추가 검산과 실제 호출 수는 COMPLETION.md 및 replay receipt를 확인한다.

저장소 root에서 Python, NumPy, SciPy, Matplotlib, pytest 환경을 사용한다.
실제 새 에너지 평가에는 검증된 SHA의 기존 MACE MP-0b3와 ASE/PyTorch가 추가로 필요하다.
모든 출력은 기존 결과와 다른 새 경로를 지정한다.

    python -m results.silicon_wafer_feasibility.audit_thermal_validity_v14 --root . --output NEW_AUDIT
    python -m pytest -q -p no:cacheprovider solver_v1/test_silicon_thermal_validity.py
    python -m results.silicon_wafer_feasibility.probe_thermal_cloud_v14 --root . --model MODEL_FILE --output NEW_PROBES --max-seconds 600
    python -m results.silicon_wafer_feasibility.replay_thermal_validity_v14 --root . --result results/silicon_thermal_validity_v14 --cloud results/silicon_thermal_validity_v14/full_dimensional_completed --output NEW_REPLAY.json

부분 재개는 --resume PARTIAL_DIR --replay VERIFIED_PREFIX_RECEIPT를 사용한다.
재개 폴더는 새 경로여야 한다. 입력·모델·원자료 해시 및 기존 제안 좌표가 달라지면
중단한다. 원자료 검산은 포텐셜 재호출을 하지 않는다.

## 다음 물리 검증

같은 MP-0b3 시편과 명시된 균열 후보 좌표에서, 구속과 유한 영역을 보존한
비선형 분포 샘플링·basin 혼합·열평형 공분산부터 검사해야 한다. 그 다음 충분한
시간과 반복을 가진 보존적/물리적으로 정당한 동역학 및 약한 conjugate 외력 대조로
기억·관성·이동도를 검사한다. thermostat 마찰을 재료의 실제 마찰로 채택하지 않는다.
MC는 열평형 준비에만 쓸 수 있고 MC step을 물리 시간이나 피로 확률로 세지 않는다.

## 근거 문헌

- D. S. Kim et al., PNAS (2018), DOI 10.1073/pnas.1707745115:
  https://arxiv.org/abs/1610.08737 — Si의 양자 핵운동·비조화성이 중요할 수 있다는
  실험/계산 근거다. 본 MP-0b3 결과가 검증됐다는 뜻은 아니다.
- N. Di Pasquale et al., arXiv:2011.00996v3:
  https://arxiv.org/abs/2011.00996 — 축약한 기억항과 Markov 근사 검증.
- Phys. Rev. B 83, 174116, Appendix A:
  https://link.aps.org/accepted/10.1103/PhysRevB.83.174116 — 조화 진동자의 유한온도 변위 분포.
