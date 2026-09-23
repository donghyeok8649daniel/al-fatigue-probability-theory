# v11 연구 재현

이 디렉터리는 작업 중이다. 완료/진행/실패를 `WORKING_STATUS.md`와 각 run의
summary/result/progress JSON으로 확인한다. 원자 개시 장벽·물리 확률·수명은 미보정이다.
`complete`는 해당 수치 실험 완료만 뜻한다.

## 입력과 실행 환경

- 중성 MACE-MP-0b3-medium SHA256:
  `2f2be696351ac9e94fbe01cdfb6f017679acdbd2db7645209ef55fec9826b012`.
- 공개 QE/PBE Si/O 자료: [MTPu](https://gitlab.com/Kazongogit/MTPu),
  commit `4ed842eea9c0b894310fee50e9bbca21404e356e`,
  `datasets/Unified_training_set2_1159.cfg`.
- 공개 CP2K Si/O/H 자료: [MLFF-SiOx](https://github.com/lukas-cvitkovich/MLFF-SiOx),
  commit `b61806d3a097e21d7880db39dcc077573e2c3c92`, `training_data.xyz.zip`.
- Python/NumPy/SciPy/ASE/MACE/Torch 버전과 입력 SHA/크기는
  `prepared_inputs_validation.json`에 기록했다.
- CPU float64, Torch2 threads, BLAS/OMP1 thread. 실제 DFT/MD 실행이나 재학습은 없다.
- 가중치, 공개 원자료 원본, 사용자 제공 PDF를 이 연구 기록에 복제하지 않는다.
  경로는 실행자가 지정한다. 원자료 SHA가 다르면 runner가 거부한다.

## 수학/원문 표 감사

저장소 root에서 다음을 실행한다. 이미 있는 결과를 덮어쓰지 않도록 새 output을 쓴다.

```text
python -m pytest solver_v1/test_silicon_initiation_research.py solver_v1/test_silicon_initiation_probability.py solver_v1/test_silicon_charge_dynamics.py -q
python results/silicon_wafer_feasibility/run_initiation_math_v11.py --output NEW_MATH_OUTPUT
python results/silicon_wafer_feasibility/audit_thermal_oxide_v11.py --help
```

41개 관련 테스트 및 합성 계산 반복 일치가 검증됐다. 이는 Si rate 또는 Hz 검증이 아니다.
generator의 열별 보존 검사는 각 열 전이율 크기에 상대적으로 적용한다. 단위 크기의
절대 허용오차는 느린 전이 또는 시간단위 변경에서 확률 손실을 숨길 수 있어 제거했다.
Thermal oxide source JSON의 값은 논문 표 전사이며 모든 기하 규약을 분리한다.

## 공개 DFT 원자료 대조

```text
python results/silicon_wafer_feasibility/run_oxide_mace_audit_v11.py --help
python results/silicon_wafer_feasibility/validate_oxide_replay_v11.py --help
python results/silicon_wafer_feasibility/run_oxidized_surface_audit_v11.py --help
python results/silicon_wafer_feasibility/validate_oxidized_surface_replay_v11.py --help
python results/silicon_wafer_feasibility/analyze_surface_environments_v11.py --help
python results/silicon_wafer_feasibility/summarize_surface_environments_v11.py --results results/silicon_initiation_v11 --output NEW_SUMMARY_OUTPUT
python results/silicon_wafer_feasibility/audit_surface_energy_differences_v11.py --run results/silicon_initiation_v11/oxidized_surfaces_max320 --output NEW_ENERGY_OUTPUT
```

QE 자료는1159개 전부, CP2K 자료는 사전 지정한 원자수 N<=320만 평가한다.
CP2K 제외453개는 source_inventory.csv에 남긴다. 힘 오차로 제외하지 않는다.
비교는 모델 독립 개발 자료와의 비교이며 MACE 사전훈련의 held-out 인증은 아니다.
QE/CP2K의 pseudopotential energy offset이 다르므로 원시 절대 에너지 RMSE로
서로 비교하지 않는다. 같은 조성 내 명시된 기준 상태의 에너지 차만 비교한다.
지정된 화학 coordination은 거리 cutoff에 의존하며 산화수/전하/유일한 계면 라벨이 아니다.

## 균열 seed 없는 유한 시편

```text
python results/silicon_wafer_feasibility/run_intact_tension_v11.py --help
python results/silicon_wafer_feasibility/run_force_controlled_prism_v11.py --help
python results/silicon_wafer_feasibility/check_prism_curvature_v11.py --help
python results/silicon_wafer_feasibility/audit_force_controlled_prism_v11.py --run results/silicon_initiation_v11/force_controlled_prism_360 --geometry results/silicon_initiation_v11/intact_prism_360_linesearch_v2/geometry.npz --output NEW_FORCE_REPLAY_OUTPUT
```

원래 LBFGS의 미수렴 checkpoint를 보존한 뒤 line-search로 이어간 실제 실행은
`intact_prism_360_linesearch_v2/protocol.json`과 restart SHA로 식별한다.
원래 격자를 처음부터 다른 optimizer로 이완한 경우와 같은 경로라고 가정하지 않는다.
원자 구조/힘은 full precision NPZ가 기준이며 extxyz는 시각화용이다.
변위0 상태의 반력은 표면 응력 때문에 0이 아닐 수 있다. 힘제어는 U-F*Delta를 푼다.
외부 grip 좌표는 정확히 보존하고 자유원자 잔여힘·반력 평형·내부 힘/토크를 함께 읽는다.
전체 graph 단절/응력 최대/단일 결합 길이는 최초 균열의 독립적인 정의가 아니다.

### 실제 인장 중단과 힘 수렴 대조

`intact_prism_360_linesearch_v2`는10800초 예산 도달로 중단됐다.0~8%5개 상태는
수렴했으나10%6번째 저장 상태는 미수렴이다. 저장 상태 수를 수렴 상태 수로 세지 않는다.
다음 정밀 계산은 동일한 축력0 시작 구조에서 힘 허용오차만50배 엄격하게 검사한다.
기존 결과를 덮어쓰지 않고 새 output을 사용한다.

```text
python results/silicon_wafer_feasibility/run_force_controlled_prism_v11.py --state results/silicon_initiation_v11/force_controlled_prism_360/state_000/raw.npz --geometry results/silicon_initiation_v11/intact_prism_360_linesearch_v2/geometry.npz --model MODEL_FILE --output NEW_TIGHT_FORCE_OUTPUT --stresses 0 .1 --fmax .00001 --max-steps 350 --max-seconds 2400
python results/silicon_wafer_feasibility/analyze_prism_rearrangement_v11.py --help
python results/silicon_wafer_feasibility/audit_force_ensemble_curvature_v11.py --curvature results/silicon_initiation_v11/prism_curvature --state results/silicon_initiation_v11/force_controlled_prism_360/state_001 --output NEW_CURVATURE_AUDIT_OUTPUT
python results/silicon_wafer_feasibility/compare_force_tolerances_v11.py --baseline results/silicon_initiation_v11/force_controlled_prism_360 --refined results/silicon_initiation_v11/force_controlled_prism_360_tight --baseline-replay results/silicon_initiation_v11/force_control_replay --refined-replay results/silicon_initiation_v11/force_control_replay_tight --output NEW_TOLERANCE_COMPARISON_OUTPUT
python results/silicon_wafer_feasibility/resume_prism_relaxation_v11.py --state results/silicon_initiation_v11/intact_prism_360_linesearch_v2/state_005/raw.npz --geometry results/silicon_initiation_v11/intact_prism_360_linesearch_v2/geometry.npz --model MODEL_FILE --output NEW_CONTINUATION_OUTPUT --fmax .001 --max-steps 600 --max-seconds 3000
python results/silicon_wafer_feasibility/audit_prism_continuation_v11.py --parent results/silicon_initiation_v11/intact_prism_360_linesearch_v2/state_005 --continuation results/silicon_initiation_v11/prism_10pct_continuation --previous results/silicon_initiation_v11/intact_prism_360_linesearch_v2/state_004 --geometry results/silicon_initiation_v11/intact_prism_360_linesearch_v2/geometry.npz --output NEW_CONTINUATION_AUDIT_OUTPUT
python results/silicon_wafer_feasibility/check_prism_curvature_v11.py --state results/silicon_initiation_v11/force_controlled_prism_360_tight/state_001/raw.npz --geometry results/silicon_initiation_v11/intact_prism_360_linesearch_v2/geometry.npz --model MODEL_FILE --output NEW_TIGHT_K64_CURVATURE_OUTPUT --ensemble force --ncv 64 --maxiter 24 --max-calls 850 --max-seconds 3600
```

원자 재배열은 수렴 상태만 사용하고 전체 자유원자의 affine+translation 최소제곱 잔차다.
국소 에너지/소성변형률/균열 detector가 아니다. 곡률 runner의 U gradient와
odd energy slope에는 dead load F를 빼야 힘제어 H=U-F*Delta의 정상점을 검사한다.
선형 외력의 일은 Hessian과 even energy curvature를 바꾸지 않는다.
정밀 상태/Krylov64 대조는 이완 상태와 Krylov 크기가 함께 달라졌다. 기존Krylov16
실패와의 차이를 어느 하나의 변화 때문이라고 단정하는 단일변수 대조는 아니다.
원래 초기 protocol의 `zero-load relaxation/grips` 문구는 실제로 변위0을 뜻했던
초기 명칭이다. 그 상태의 실제 반력은3.184GPa이며 축력0으로 읽지 않는다.

## 추가 재현 및 무결성

### 같은 grip으로 되돌리는 정적 대조

```text
python results/silicon_wafer_feasibility/run_prism_unload_v11.py --state results/silicon_initiation_v11/prism_10pct_continuation/raw.npz --result results/silicon_initiation_v11/prism_10pct_continuation/summary.json --geometry results/silicon_initiation_v11/intact_prism_360_linesearch_v2/geometry.npz --model MODEL_FILE --output NEW_UNLOAD_OUTPUT --target-strain .08 --fmax .001 --max-steps 600 --max-seconds 3000
python results/silicon_wafer_feasibility/audit_prism_unload_v11.py --baseline results/silicon_initiation_v11/intact_prism_360_linesearch_v2/state_004 --unload results/silicon_initiation_v11/prism_unload_10_to_8 --geometry results/silicon_initiation_v11/intact_prism_360_linesearch_v2/geometry.npz --output NEW_UNLOAD_REPLAY_OUTPUT
python results/silicon_wafer_feasibility/plot_prism_unload_v11.py --results results/silicon_initiation_v11 --output NEW_RETURN_FIGURE_OUTPUT
```

수렴한10% 상태의 에너지/힘을 먼저 재계산한 뒤 원래 reference의 축방향 비율로
자유원자 predictor를 이동하고8% grip을 고정한다. predictor와 최종 상태를 각각
저장한다. 이전8%와 labelled 원자 차이 및 동일 Si 원자 순서에 무관한 최소 제곱
거리 assignment를 비교한다. 정적 이력 의존성 대조이며 실제 피로주기·소산·
committor 또는 균열 개시의 판정이 아니다. 미수렴이면 상태와 잔여힘을 그대로 읽는다.

### 최종 원자료/소스와 Git 저장 바이트

```text
python results/silicon_wafer_feasibility/verify_initiation_bundle_v11.py --repo . --revision HEAD
```

artifact_manifest.json은 모든 관련 계산이 끝나거나 중단 상태로 보존된 뒤 작성한다.
위 검사는 working files와 HEAD의 실제 blob을 각각 SHA256으로 대조하고, 테스트한
7개 소스가 현재 source manifest와 같은지도 확인한다. --revision INDEX는 commit
전에 staged blob을 검사한다. 이 무결성 검사는 물리 검증이나 테스트 재실행이 아니다.

기존 공용 Python 소스4개는 Windows working tree의 CRLF와 Git의 LF가 다르다.
source manifest는 실행한 파일의 SHA를 보존하고, 줄바꿈만 LF로 바꾼 git_lf의
별도 bytes/SHA도 기록한다. 검증기는 실제 변환의 동일성을 검사한 뒤 Git blob과
대조한다. 기존 공용 소스의 물리 내용은 바꾸지 않았다. 결과 파일은 줄바꿈 변환을
금지하여 NPZ/CSV/JSON/문서를 포함한 실제 바이트가 그대로 일치해야 한다.
다른 운영체제의 fresh checkout이 이미 LF를 쓰면, 이4개 소스에 한해서 명시된
git_lf hash를 working tree에서도 허용하고 해당 파일 목록을 별도로 출력한다.

## 보존한 실패

- `initial_lbfgs_attempt`: 원래 optimizer 소스와 미수렴 실행.
- `intact_prism_360_linesearch/failed_restart_source.py`: extxyz 저장 자릿수 때문에
  검사가 실패한 소스. 허용 저장 오차 검사 뒤 prescribed grip을 정확 복원하도록 수정.
- `failed_parser_v1`: 유효한 left-handed cell을 배제한 초기 parser.
  v2는 부호가 아닌 cell nonsingularity를 검사하고 동등 기저 대조를 수행했다.
- raw log에는 로컬 경로가 있으므로 Git에 넣지 않는다. 최종 검증 기록은 경로를
  제거한 별도 텍스트/JSON으로 보존한다. 실패를 성공 run으로 합산하지 않는다.
