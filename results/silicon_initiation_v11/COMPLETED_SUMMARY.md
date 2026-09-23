# Si 균열 개시 v11 - 이번 계산 묶음의 완료·실패·한계

2026-09-24. 사용자 요청은 초기 균열을 넣지 않은 시편의 첫 균열 개시다.
예약된 원자 계산과 후처리는07:42KST까지 모두 종료했다. 계산 종료와 재료/개시
보정 성공은 다르다. 추가 agent0, 새 DFT0, MD0, 재학습0. Al/UI/생산 Hz는 보존했다.

## 핵심 결과

1. 제공 열산화층 논문의 기하·응력 분담27대조를 계산했다. 논문 FEM은 측정
   파단강도를 하중으로 사용하며 독립 강도 예측이 아니다. 200nm 막의 core 치수와
   산화막 제거 행은0.20um 차이로 미해결이다. 역산 결함을 독립 입력으로 재사용하지 않았다.
2. 공개 QE/PBE1159상태와 CP2K/PBE1013상태에 새 MACE 평가를 실행했다.
   CP2K Si/O/H군의 산소이웃1~3개 Si 힘 RMSE는 약0.572eV/A이며 거리 규약3개에서
   큰 오차가 유지된다. 동일조성 에너지 순서의 역전도 있다. 모델은 재료 미채택이다.
3. 균열 seed 없는 bare360Si 시편의0~8%5상태와 별도10% 후속 이완이 수렴했다.
   8->10%에서 큰 재배열/응력 하락이 있지만 첫 균열 판정은 못 했다.
4. 10%에서8%로 되돌린 정적 이완도 수렴했다. 같은8% grip에서 에너지 차이
   -9.719223eV, 응력8.624053->2.233825GPa, 원자 matching RMS1.054057A다.
   정적 하중 이력 의존성이며 실제 소성/피로소산/균열 개시 인증은 아니다.
5. 축력0/100MPa를 별도로 풀고 힘 허용오차를50배 엄격하게 재검사했다.
   전체 지지잔차/면적은 각각0.070630/0.063927MPa로 감소했다. 축력0 길이는
   초기 고정길이보다2.78964% 짧다. 변위0에서의3.184GPa 반력을 무하중으로 부르지 않았다.
6. 원인별 첫 통과, committor/MFPT, 고정 화학배치 혼합을 공통 확률 구조로 구현했다.
   시간단위 변경 때 확률 누출을 놓치던 보존 검사를 수정했다. 관련41PASS(19.53s),
   합성 결과3파일 byte 재현. Si의 실제 전이율/수명/physical time 보정은 아니다.

## 실제 실행 상태

| 결과 폴더 | 실제 상태 | 시간 또는 범위 |
|---|---|---|
| thermal_oxide_audit | 완료 | 원문5그룹/기하·잔류응력27대조 |
| oxide_mace_mtpu_1159_v2 | 완료 및 독립 replay |1159상태/17413원자,1023.68s |
| oxidized_surfaces_max320 | 완료 및 독립 replay |1013상태/179857원자,5363.40s |
| surface_environments / surface_environment_summary_v2 | 완료 | 3종 거리 규약, 모든 원자 오차 보존 |
| surface_energy_differences | 완료 |14개 동일조성·관찰그룹, offset 적합0 |
| intact_prism_360_linesearch_v2 | 부분 종료 |1237force calls/10806.22s;0~8%5개 수렴,10% 원래 값은 미수렴 |
| prism_10pct_continuation | 완료 및 독립 replay |364force calls/2250.20s,197steps |
| prism_unload_10_to_8 | 완료 및 독립 replay |242force calls/1926.65s,131steps |
| force_controlled_prism_360 | 완료 및 독립 replay |102force calls/937.54s |
| force_controlled_prism_360_tight | 완료 및 독립 replay |129force calls/990.07s |
| prism_curvature | 미수렴 종료 |519force calls/4500.02s, 저장/독립검사mode0 |
| prism_curvature_tight_k64 | 미수렴 종료 |495force calls/3601.82s, 저장/독립검사mode0 |
| first_passage_math_scale_guard_replay | 완료 | 합성3파일 기존 결과와 byte-exact |

원래14개 인장 계획 전체를 완료한 것이 아니다. 원래 미수렴10%는 그대로 보존하고
같은 grip에서 이어간 별도 계산의 수렴 결과를 추가했다. 곡률 큐가 정상 종료됐다고
곡률 수렴으로 읽지 않는다. 두 곡률 실패 모두 보존했고 안정성은 미확정이다.
정밀100MPa의 H=U-F*Delta gradient 최대2.93206e-6eV/A, raw 재검증 차이1.11e-16이다.

## 검증의 범위

- QE/CP2K의 독립 parser로 원자료 배열을 대조했다. raw 차이0, CP2K 프레임 지표
  차이 최대3.35e-14다. MACE 사전훈련과의 중복은 미감사이므로 held-out 인증은 아니다.
- Force-only/표준 계산, work에 정확히 대응하는 힘, 고정 grip, 내부 힘/토크,
  합성 확률 보존/격자수렴을 검사했다. 남는 지지 반력 잔차는0으로 숨기지 않았다.
- 초기 LBFGS 정체, extxyz 자릿수 읽기 검사, left-handed cell parser 실패와 원시
  prefix를 보존했다. parser는 수정 후 전체 자료를 새로 평가했다.
- 동일 Si 원자의 순서 교환을 허용한 비교도 구조 차이를 없애지 못했으나,
  이것만으로 실제 물리 소성/균열 또는 특정 상전이를 분류하지 않는다.
- 도펀트1개/360자리 시편은1.35703e20cm^-3다. 1e18cm^-3의 균일 독립 치환 예시에서
  적어도1개 존재 확률0.7342%이며 개시확률/활성전하/편석의 측정값이 아니다.

## 원자료와 재현

- RESULTS_AND_LIMITS.md: 상세 수치와 해석 한계.
- REPRODUCE.md: 실제 실행/독립 replay/최종 Git 바이트 검증 명령.
- validation.json: 실제41PASS 결과와 시험한 구현/테스트7파일 SHA256.
- source_manifest.json: 공개 원자료 commit/DOI/SHA, 모델 가중치 SHA, 코드 의존성.
- artifact_manifest.json: 종료 후의 결과·실패·checkpoint 파일 SHA256.
- verify_initiation_bundle_v11.py: working 파일과 지정 Git commit/index의 실제 blob을
  각각 대조한다. 검증 결과의 실제 커밋은 최종 PDF와 저장소 인계를 따른다.
- 사용자 제공 논문 PDF, 모델 가중치, 개인 경로가 든 로그/PID는 Git에 복제하지 않는다.

## 아직 달성하지 못한 것

실제 최초 균열의 물리적 상태와 핵생성 경로, 해당 Si/O/도펀트 환경의 장벽,
유한온도 자유에너지, 축약 동역학의 정당성, 물리 이동도/초/Hz/수명.
GPa를 MPa로 임의 축소하지 않았고 성장 장벽을 개시 장벽으로 바꾸어 부르지 않았다.
