# Si 균열 개시 v11 — 진행 기록

시작 2026-09-24 02:25:46 KST, 사용자 요청 약7시간.
마감 2026-09-24 09:25:46 KST. 자동 이어가기 si-7도 같은 마감이다.

## 범위 정정

사용자의 주목표는 균열 개시다. v2의 crack seed 및 v3/v4 그 끝점의
열운동, v6 K 경계는 성장 관련 보조 연구다. v10 강체 계면 분리도
균열 없는 상태의 개시 장벽이 아니다. 그 자료와 실패 기록은 보존한다.
검증된 Si 개시 확률/장벽/수명은 현재 없다.

## 이번 실행 순서

1. Tsuchiya 외, Engineering Fracture Mechanics163(2016)523–532,
   DOI10.1016/j.engfracmech.2015.08.029의 표/기하/응력 환산 감사.
   제공 PDF SHA256 3f1b9270ca5266cce1f71c4a08a94de801b56ff8dc990f019b3e13c6c37cdf06.
   비공개 경로/PDF 자체를 Git에 넣지 않는다.
2. 균열을 넣지 않은 bare Si 유한 인장 시편의 실제 원자 기준 계산.
   제조 시편 전체·산화막·유한온도·실험 개시 강도를 검증한 것으로 읽지 않는다.
3. 산화막 모델의 독립 DFT 근거와 표면/계면 에너지 적합 가능성 확인.
4. 첫 통과/원인별 흡수 확률의 공통 구조를 정리하고 필요한 수학 검증.
5. 실제 수행, 미완료, 실패 및 검증 상태를 정리·정상 Git 게시.

## 중요 판정 원칙

- 단일 결합 거리 또는 힘 최대값만으로 균열 개시를 선언하지 않는다.
- 초기 상태에는 내부 면 절단, void, crack seed를 넣지 않는다.
  외부 자유표면은 정의하며 표면 재구성은 허용한다.
- 논문 FEM은 측정 평균강도를 하중으로 썼다. 독립 강도 예측이 아니다.
- 논문 역산 초기 결함 크기를 같은 강도의 독립 예측에 재사용하지 않는다.
- 산화막0.46GPa는 참고문헌에서 가져온 입력이지 이 시편의 측정값이 아니다.
- 추가agent0, Al/UI/생산Si/Hz gate 불변. 중단 prefix도 보존한다.

## 최신 계산 상태 - 2026-09-24 07:42KST까지 모든 예약 계산 종료

- 10% 후속 이완 및 독립 원자료 재검증 완료. fmax9.0183e-4eV/A,
  응력3.617024GPa. 큰 원자 재배열이며 첫 균열 판정은 미완료다.
- 정밀100MPa/Krylov64 곡률은07:39KST 495force calls/3601.82s로 예산 종료했다.
  Lanczos 미수렴/저장mode0/독립 검사mode0. 예약 후처리도 완료했고 gradient
  raw 대조차1.11e-16eV/A, H gradient 최대2.93206e-6eV/A다. 안정성 인증은 아니다.
  `postchecks_queue_status.json`은complete이며 재실행하지 않는다.
- 10->8% 정적 감하중도07:42KST 완료:242force calls/1926.65s,131steps,
  fmax9.35926e-4eV/A, 응력2.233825GPa. 이전8%와 에너지차-9.719223eV,
  동일 Si 원자 순서에 무관한 위치 RMS1.054057A. 원시 관측량/grip 재검증 차0.
  `unload_audit_status.json`도complete다. 실제 소성/개시/피로소산의 판정은 아니다.
- 기존 큐, 정밀힘, continuation, 곡률 후처리, 감하중 후처리는 모두complete다.
  07:46KST 이전에 해당 부모 프로세스가 종료한 것을 확인했다. 다시 실행하지 않는다.
- 감하중 비교 그림은 `prism_unload_figure_v2`를 시각검수했다. 첫 그림의 겹친
  라벨을 고친 것이며 수치/물리 입력 변경은 없다. 보고서 builder는16쪽이다.
  PDF의 실제 작성 시각/사용 커밋/전체 검수 상태는 최종 전달본과 원래 폴더 인계를 따른다.
- 재현 설명에 감하중 명령과 최종 Git blob 검증 명령을 추가했다.
  `verify_initiation_bundle_v11.py`는 실제 working/index/commit 바이트를 검사한다.
  최종 manifest와 함께 재실행할 수 있다. 파일 무결성을 물리 보정으로 읽지 않는다.
- v11 소스32개 AST 구문 검사 통과/CRLF0. 새 원자료/코드 개인 절대경로
  검색은 일치0이었다. 물리 검증이나 관련41PASS의 재실행은 아니다.
- 상세 최종 결과는 `COMPLETED_SUMMARY.md`, `RESULTS_AND_LIMITS.md`를 먼저 읽는다.
  이후 파일 무결성/Git/PDF 확인은 계산의 새 실행이 아니다. 새 장시간 계산을 시작하지 않는다.

## 아래는 각 시점에 기록한 실행 이력

### 2026-09-24 05:52 KST 기록

시작 기준 local/origin 736bcf5494871b483dfcc4f20fee54e509099ff6 일치/Si clean 확인.
v11은 아직 커밋/게시 전이다. 아래의 수학 검증과 실제 원자 계산을 구별한다.

### 완료

- `thermal_oxide_audit`: 원문 5그룹, 기하 3규약과 잔류응력 3해석의
  27조건을 재계산했다. 1D 힘 잔차 최대 1.42e-14 GPa um^2.
  원문 외형치와 막두께/면 규약에서 Si 응력은 원문 FEM값과 최대 0.609% 차이이나,
  측정 파단강도를 입력한 계산이다. 독립 강도 예측이 아니다.
  200nm 조건의 core 치수와 산화막 제거 행에 0.20um 차이가 남는다.
- `oxide_mace_mtpu_1159_v2`: 공개 QE/PBE 기준 1159상태/17413원자의
  중성 MACE 힘 비교 완료. SiO2 조성 RMSE 0.15823eV/A,
  기타 Si/O 조성 0.34272eV/A. Si-only 전체는 큰 오차 상태를 포함해
  1.27379eV/A이며 제거하지 않았다. 형태 레이블 없는 자료를 계면으로 부르지 않는다.
  실제 새 DFT/MD/fit은 0. MACE 훈련자료와의 중복은 미감사.
- `oxide_replay`: 독립 고정열 parser로 1159상태의 원본/저장 배열을
  대조하고 모든 프레임/그룹 오차를 재계산했다. 최대 차이 0.
- `oxidized_surfaces_max320`: CP2K Si/O/H 사전선택1013구조/179857원자
  평가 완료, 원래1466구조 중 제외453개 ID도 보존했다. 새 force-only1013회,
  표준평가9대조 E/F차0. Si/O/H vacuum군423구조/106148원자의 힘 성분
  RMSE0.328090eV/A, Si/H0.274488, Si/O(no-vacuum label)0.185805.
  원자1~2개 small계의 큰 오차도 별도 표로 보존하며 제거하지 않았다.
- `oxidized_surface_replay`: 별도 shlex/고정열 parser로 원자료 배열 대조차0,
  프레임지표 차이최대3.35e-14, 그룹지표 차이0. 원자료는 공개 DFT이며
  이번에 새 DFT/MD/학습을 실행한 것은 아니다. 사전학습 중복도 미감사다.
- `surface_environments`와 `surface_environment_summary_v2`:3종 이웃거리
  규약에서 모든 원자의 오차를 보존해 화학적 이웃 구성별로 분해했다.
  중앙 Si-O2.0A 규약에서 산소이웃1~3개 Si의 RMSE0.571980eV/A,
  원자비중7.173%/오차제곱합비중21.801%. Si-O1.8/2.2A에서도
  RMSE0.573512/0.569989로 큰 오차가 유지된다. 산화수/전하/개시위치의
  증거로 해석하지 않는다. 그룹 RMSE 재집계 차이최대1.78e-15.
  산화Si/O/H군의 균일 힘오차 성분은 전체 오차제곱합의0.0006553%뿐이다.
  원자료 net force가 전체 모델 오차를 설명하지 않으며 원래 힘을 수정하지 않았다.
  두 그림을 렌더 검수했고 첫 그림 축라벨 줄바꿈만 수정했다(v2가 최종 그림).
- `force_controlled_prism_360`: 축력0/100MPa 두 상태 모두 수렴, 새 MACE
  force-only102회/937.54s. 축력0 길이는 초기 grip간격 대비2.78965% 짧다.
  100MPa에서 이 기준 대비변형0.0890519%; 두 상태 secant112.438GPa는
  rigid-grip 유한시편 응답이지 bulk 영률/실제웨이퍼 검증값이 아니다.
- `force_control_replay`: raw좌표/힘/gradient/enthalpy 독립재계산 최대차0,
  prescribed grip 좌표차0. 축력목표 오차는0/100MPa에서 -0.113225/+0.015066MPa.
  전체 지지반력 합을 면적으로 나눈 잔차는 각각5.04941/0.822802MPa로 남는다.
  자유원자 잔여힘 합과 반력 사이 등식은1.58e-15eV/A 이내로 확인했다.
  축력 수렴과 완전한 지지힘/토크 평형을 혼동하지 않는다. 그림 시각검수 완료.
- `surface_energy_differences`: 동일조성/관찰그룹14개에서 상대에너지 raw
  재검증 차이0, 쌍 제곱합-분산 항등식 상대차1.62e-13. Si/O/H423상태의
  동일조성25955쌍 차이 RMSE0.020277eV/atom, DFT차이nonzero25954쌍 중
  2107쌍의 에너지 차이 부호가 반대다. 쌍은 독립 표본/반응경로/개시 장벽이 아니다.
  `MATERIAL_AUDIT_INTERPRETATION.md`에 원자료 범위·국소 환경·에너지 판정을 정리했다.
- 균열 개시의 cause별 first-passage collector/committor/MFPT 구현,
  고정 화학배치 혼합 검증을 추가해 기존 charge dynamics 포함 관련 테스트
  **39 passed (6.90s)**. 직전38PASS 결과도 raw 로그로 보존.
  합성 확률 예제 mass 잔차 2.00e-15, 확산 first-passage 격자오차 비
  4.17/4.04. 실제 Si basin/rate/physical clock을 보정한 결과가 아니다.
- `intact_prism_360_linesearch_v2/initial_force_validation.json`:
  force-only/표준 MACE 힘 최대차이 2.22e-16eV/A.
  자유원자 및 grip-work 유한차분은 간격 반감 때 약 4배 오차 감소.
- 같은 run의 state_000: 360원자/216자유원자의 초기 고정길이 이완 수렴,
  fmax 9.15e-4eV/A. **변위 0은 힘 0이 아니다.** 표면 이완 후 이 길이의
  명목 축응력은 3.18415GPa다. 실제 무하중 기준과 별도이며 실험 개시강도가 아니다.
- `prism_snapshot_audit_0406`: state_000의 raw힘/에너지로 경계·반력·토크를
  독립 재계산, 저장 관측량 최대차이0/처방 grip 좌표차이0. 반력과 자유원자
  잔여힘의 벡터 등식 오차1.07e-15eV/A, torque등식4.08e-14eV.
  지지 반력 자체의 합은 잔여 자유힘 때문에 축력의0.498%이며 0으로 숨기지 않았다.
  2.8A 거리기준 기존 이웃쌍 손실0/새 이웃쌍28. cutoff3.1/3.4/3.7A도
  전체 연결은 유지된다. 이 진단이 실제 최초 균열 detector는 아니다.
  실제 원자구조4패널 그림을 렌더·시각 검사했다.
- `prism_snapshot_audit_0545`:0/2/4/6/8%5상태의 raw관측량/처방 grip
  좌표차0. 네 거리기준에서 초기 격자 이웃쌍 손실0, 8%의2.8A기준 새쌍38개.
  모든 원자는 연결되어 있으나 이것만으로 개시 여부를 판정하지 않는다.
  6->8%의 deltaE4.46014eV와 coarse trapezoid work6.20311eV의 차이는
  -1.74297eV다. 하중 간격·구조 재배열·이완이 섞인 진단이며 소산 또는
  개시 장벽으로 부르지 않는다. 8% 구조 그림을 시각 검수했다.
- `doping_occupancy`: nominal360자리/부피7.36904e-21cm^3에서
  도펀트1개/시편은1.35703e20cm^-3다. 독립·균일 치환 예시에서
  1e18cm^-3일 때 적어도1개 있을 확률0.7342%. 측정된 배치분포/개시확률은 아니다.
  Binomial 독립 식 대조 차이5.55e-17. 고정 배치의 조건부 생존 혼합을
  평균 에너지/평균률 대입과 구분하도록 이론과 해석해 테스트에 추가했다.
- `validation.json`: 실제39PASS 로그와 해당 구현/테스트7파일SHA를 기록했다.
- `INITIATION_SOURCE_SCOPE.md`: 개시/표면/경계조건 1차 문헌6개를 검토했고
  원문/초록 확인 범위를 구분했다. 성장 선단 장벽을 개시에 재사용하지 않는다.
  native-oxide2023의 p=(W0-Wfail)/W0는 폭 감소 지표이며 우리 흡수 확률과
  다르다는 점을 원문17번째 PDF페이지 이미지에서 확인했다. 추가 재료 보정값 채택0.

### 06:30KST 최신 검증 및 실행 상태

### 07:20KST 후속 결과와 마지막 실행

- 10% continuation 완료:07:06KST,197steps/364force calls/2250.20s.
  fmax9.01826e-4eV/A, 명목응력3.617024GPa, E-1740.795500572eV.
  고정10%에서 에너지 변화-3.577201711eV는 정적 optimizer의 감소이며 물리 소산/장벽이 아니다.
- `prism_10pct_replay` 완료: 원시 관측량/grip 차0, 힘/토크 등식오차
  1.85e-15eV/A /3.61e-14eV.8% 대비 affine 제거RMS0.954638A/max2.343656A.
  2.8A 초기 이웃쌍 손실31/생성98, 네 거리규약에서 전체 연결 유지.
  그림 시각검수 완료. 큰 재배열의 진단이며 첫 균열 또는 물리 소성 인증이 아니다.
- `prism_unload_10_to_8` 시작07:09KST/부모PID19076. 수렴10%에서8%로 되돌린
  prescribed grip과 기준 위치의 affine predictor를 사용했다. 새로운 섭동/seed없음.
  fmax.001/maxsteps600/max_seconds3000. 종료 예상08:00전후, 남은 마감 시간 확보.
  source raw energy차0/force차2.22e-16/grip차0으로 시작을 검증했다.
  운영 정보는 비추적 unload_status.json. 지금 정밀Krylov64와 함께 총2개 실행 중.
- 감하중 원자료 대조도 이미 예약: 후처리 부모PID55496,
  비추적 unload_audit_status.json / postcheck_unload.ps1. 완료 뒤
  `prism_unload_replay`에 동일8% 경계의 에너지·힘·Si 원자 matching을 검사한다.
  원자번호 직접 비교와 자유Si permutation에 대한 최소제곱 거리 비교를 구분한다.
  이 후처리를 직접 중복 실행하지 않는다.
- `RESULTS_AND_LIMITS.md`에 완료된 결과를 결론 중심으로 정리했다. 현재 중간본이며
  정밀곡률과 감하중 결과, 최종무결성/Git은 실제 종료 뒤 갱신해야 한다.
- PDF builder를15쪽으로 확장했다.10%완료/재배열과 감하중 대조용 보충 페이지를
  추가했고 아직 이15쪽 버전은 생성·렌더 검수 전이다. 마커는 다시 실행하지 않는다.
  정밀곡률 후처리가 완료되면 --curvature-audit curvature_force_ensemble_audit_tight_k64를
  명시한다. 원래 곡률 실패를 계속 보존해 설명한다.
- 결과별 .gitattributes로 SHA-bound 원자료 바이트를 보존하도록 준비했다.
  최종 manifest 생성 전 적용하고 실제 index/blob과 파일 바이트의 일치를 검증한다.

### 이전 단계 기록(06:30~06:47)

- 06:38KST `prism_curvature_tight_k64` 별도 수치 대조 시작. 이번에는 fmax1e-5로
  이완한100MPa state를 사용하고 Krylov 공간을16->64로 확장했다.
  maxiter24/max_calls850/max_seconds3600, force ensemble649좌표, step2e-4A.
  시작 부모PID55808, 비추적 curvature_tight_k64_status.json에 운영 정보가 있다.
  원래 느슨한 상태의 곡률 실패는 보존한다. 새 검사가 같은 상태의 엄밀한 단일변수
  수렴 대조는 아니다(원자 이완 상태와 Krylov 크기가 함께 바뀜). 중복 실행 금지.
  완료 후 새 `curvature_force_ensemble_audit_tight_k64` 폴더로 force-ensemble
  gradient/energy 해석을 검증한다. 기존 audit 폴더를 덮어쓰지 않는다.
- `audit_prism_continuation_v11.py` 준비:10% 후속 계산이 끝난 후 parent와
  동일한 grip을 확인하고 원시힘/에너지/이웃쌍/8%대비 기하 재배열을 독립 검사한다.
  아직 미실행이며10% continuation은 현재 진행 중이다. 고정 변위에서 optimizer가
  낮춘 에너지를 물리 소산/시간/장벽으로 부르지 않는다.
- 06:47KST 후처리 큐PID43328도 시작했다. 비추적 postchecks_queue_status.json과
  postchecks_0645.ps1을 확인한다. continuation 부모PID2652 종료 뒤
  `prism_10pct_replay`를 생성하고, 정밀곡률 부모PID55808 종료 뒤
  `curvature_force_ensemble_audit_tight_k64`를 생성한다. 두 후처리는 이미 예약됐으므로
  직접 중복 실행하지 않는다. 곡률 source 실패/미수렴도 그대로 기록한다.
  새 후처리 코드의 실제 실행 전 AST 구문 검사는 통과했으며 검증 성공을 의미하지 않는다.
- 06:40대 새 fetch 성공 후 Si local/origin736bcf5 일치. 아직 v11커밋 전.

- 정밀 `force_controlled_prism_360_tight` 및 독립 `force_control_replay_tight` 완료:
  129force calls/990.07s, 두 상태 모두 수렴. 비교는 `force_tolerance_comparison`.
  0MPa: 자유힘max2.546e-6eV/A, 축력오차-0.00199779MPa,
  전체지지잔차/면적5.049414->0.0706296MPa.
  100MPa: 자유힘max3.001e-6eV/A, 축력오차-0.000269569MPa,
  전체지지잔차/면적0.822802->0.0639267MPa. 지지토크도0.0221772->6.3747e-5eV.
  두 점 secant112.438129->112.162573GPa. 수렴 대조이며 재료/개시 승인이 아니다.
  기존 느슨한 상태의 곡률 실패를 새 정밀 상태의 곡률 판정으로 전이하지 않는다.
  refine_queue_status는complete이다. 현재 무거운 계산은10% continuation과
  위 정밀100MPa/Krylov64 곡률 대조의2개다.

- 확률 generator 보존 검사에서 단위에 의존하는 허용오차 결함을 재현·수정했다.
  예전 max(1,max|L|) 기준은 열 합이0이 아닌 L도1e-15배하면 통과시켰다.
  각 열의 L1 크기 기준으로 수정하고1e-15/1/1e15 시간단위 및 빠른/느린 열의
  혼합 대조를 추가했다. 관련 **41PASS(19.53s)**, 현재 validation.json에 실제
  소스SHA를 기록했다. 이전39PASS 기록은 validation_before_scale_guard.json에 보존.
  `first_passage_math_scale_guard_replay`의3개 결과는 기존 합성과 byte-exact.
- `prism_curvature`는06:28KST 4500.02s/519force calls에서 예산 종료했다.
  Lanczos 미수렴, 검사된 mode0개다. 안정성 통과가 아니며 반간격·독립 에너지
  mode 검사는 수렴 mode가 없어 실행되지 않았다.
  `curvature_force_ensemble_audit` 완료: 원래 U gradient max0.118909eV/A,
  H=U-F*Delta gradient max0.000137357eV/A. 원시 force-state와 차이4.44e-16.
  외력 일의 선형성에 따른 Hessian 동일성과 정상점 기울기의 차이를 구분했다.
- `prism_10pct_continuation`은06:28KST 시작. 원래 미수렴 state_005/raw.npz에서
  같은 고정 grip/10%변위 그대로3000초 한도 후속 이완이다. 이전 자료 덮어쓰기 없음.
  `resume_prism_relaxation_v11.py`의 source_replay로 시작 energy/force를 대조한다.
  운영 큐PID2652, 상태는 비추적 continuation_queue_status.json이다. 중복 실행 금지.
  실제 source_replay: 에너지차0/힘차2.22e-16eV/A/grip좌표차0.

- 06:12KST 변위제어 실행이10800초 예산으로 종료됐다. 총1237force calls,
  elapsed10806.22s. 저장6상태 중0~8%5개만 수렴했고10%는 미수렴이다.
  state_005의 fmax0.165302eV/A, 명목응력4.00165GPa는 중단 시점 값이다.
  응력 하락을 파손강도/개시로 해석하지 않는다. raw/checkpoint/result/progress를 보존했다.
- 06:12KST 별도 `force_controlled_prism_360_tight`가 시작돼06:29 완료됐다.
  기존 축력0 raw를 시작점으로0/100MPa, fmax1e-5eV/A(기존5e-4의1/50),
  max350steps/2400s의 실제 수렴 대조다. 기존 결과는 덮어쓰지 않는다.
  비추적 `refine_queue_status.json`과 `refine_force_after_displacement.ps1`이
  상태를 기록하고 `force_control_replay_tight`로 독립 재검증했다.
  큐PID58872가 기존 변위제어 종료/메모리 여유를 확인하고 시작했다. 중복 실행 금지.
- `prism_rearrangement_0610`: 수렴5상태의 자유원자216개에서3D affine+translation을
  제거한 기하 잔차를 계산했다.6->8% RMS0.224178A는4->6%0.034747A의6.45배다.
  초기 grip이웃36개/나머지180개의 RMS는0.134284/0.238118A이며,
  해당 구간 grip이웃의 잔차제곱합 비중은5.98%다. 경계조건 영향을 배제한 것은 아니다.
  강체 좌표변환 불변성 오차5.56e-14A, LS 직교성 잔차2.04e-15,
  독립 all-pair와 기존 neighbor count 일치. 그림 시각검수 완료.
  새 potential/DFT/MD0. 소성변형률/에너지 국소화/개시 detector가 아니다.

- `intact_prism_360_linesearch_v2`: 변위제어 0..26%/2%간격 총14상태 중
  5상태 수렴(0%,2%,4%,6%,8%);10% 미수렴 중단.8% 명목응력8.62405GPa,
  fmax8.25e-4eV/A. 실제 상태는 각 state/result.json과 progress.json 참조.
  초기 crack seed/내부 면 절단 없음. 그래프/응력 최대만으로 개시를 판정하지 않는다.
- `prism_curvature`: 자유원자648좌표와 grip연장1좌표의 힘제어 검사였으며
  위와 같이 예산 종료/미수렴이다. 후처리도 완료했으므로 다시 실행하지 않는다.
  원래 runner의 g0/projected_gradient/odd_slope는 U 기준이며 힘제어의 정상점은
  H=U-F*Delta의 gradient로 읽는다. 선형 외력은 Hessian을 바꾸지 않는다.
  command: `python results/silicon_wafer_feasibility/audit_force_ensemble_curvature_v11.py
  --curvature results/silicon_initiation_v11/prism_curvature
  --state results/silicon_initiation_v11/force_controlled_prism_360/state_001
  --output results/silicon_initiation_v11/curvature_force_ensemble_audit`.
  지지반력의 엄격한 힘 허용오차 대조도 완료했다. 기존 축력0의 전체
  지지잔차5.05MPa를 정확한 완전평형이라고 숨기지 않는다.
- **실행 큐가 이미 있으므로 중복 실행하지 않는다.** 표면 감사 -> 독립재검증 ->
  화학환경 분해 -> 힘제어0/100MPa까지 완료했고 국소 곡률은 미수렴 종료했다.
  중간 단계가 실패하면 prefix/실패를 보존한다. 무거운 계산은 동시에2개 이내다.
  큐 자체의 운영 상태는 원래 Al worktree의 비추적 `.cache/si-initiation-v11/queue_status.json`,
  실행 제어 스크립트는 같은 cache의 `queue_after_surface.ps1`이다.
  각 출력의 파일 존재와 complete를 직접 확인한다. 큐가 시작된 것과 계산 완료는 다르다.
- `prepared_inputs_validation.json`: 별도의 shlex/fixed-column reader로
  CP2K1466개 inventory/사전선택1013개/제외453개 일치 확인.
  합성 확률 재실행3개 파일은 최초 실행과 byte-exact.
- 곡률 검사는 고정 변위와 힘제어의 허용 변위를 구분한다. 힘제어에서는
  자유원자648좌표에 grip extension1좌표를 추가한다. 음의 방향은 불안정성
  증거가 될 수 있으나 제한된 Lanczos의 양의 Ritz값만으로 전역 최소를 인증하지 않는다.

### 실패/중단 기록

- `intact_prism_360_seed1101`: 원래 LBFGS 이완이 정체되어 해당 작업만 중단.
  checkpoint를 보존하고 line-search 재시작의 시작점으로 사용했다.
  초기 상태가 수렴한 실행으로 세지 않는다.
- `intact_prism_360_linesearch`: extxyz 8자리 저장의 4.64e-9A round-trip 차이가
  기존1e-9A 읽기 검사에 걸려 평가 전 실패. 1e-7A 이내를 확인한 뒤 원래 grip
  좌표를 정확히 복구했다. 힘 허용오차를 늘린 수정이 아니다.
- `oxide_mace_mtpu_1159`: 유효한 left-handed cell을 det>0 검사로 거부한 parser
  결함. nonsingular abs(det)>0로 수정하고 v2를 새로 전 실행했다.
  frame575의 동등 격자기저 대조 E차0/F차3.22e-15eV/A.

### 남은 판정

실제 Si의 핵생성 장벽, committed basin, 유한온도 개시확률/물리시간은 미확정.
현재 계산의 국소 힘 수렴은 Hessian 안정성이나 재료 검증을 대신하지 않는다.
완료 시 source/실행 manifest, 재현 설명, 최종 검사와 Git 상태를 별도 기록한다.
05:11의 새 fetch 후 local/origin736bcf5 일치. 아직 v11 커밋 전.

### PDF 준비

- 06:20 초안에 원자 재배열 페이지와10% 미수렴 자료를 추가해14쪽이 됐다.
  변경9/10/14쪽 렌더 시각검수 완료. 곡률 실패와 새41PASS는 다음 재생성에 반영한다.
  최신 힘 결과는 --force-replay force_control_replay_tight로 선택한다.
- 06:33 중간PDF 재생성:14쪽/사용source26개. 정밀힘/곡률미수렴/41PASS를
  반영하고 변경11/12/14쪽을 다시 렌더·시각검수했다. 아직 최종 전달본은 아니다.
  이후 builder는10% continuation의 진행 또는 완료를9쪽에 자동으로 추가하므로
  다음 생성 뒤 해당 페이지도 다시 검수한다. marker는 다시 실행하지 않는다.

- 새13쪽 보고서 builder `results/silicon_wafer_feasibility/build_initiation_report_v11.py`.
  출력은 원래 작업 폴더 `output/pdf/silicon_wafer_initiation_v11_2026-09-24.pdf`.
  05:43 중간 기록을13쪽 모두100dpi로 렌더해 시각 검수했다. 글리프/겹침 문제없음.
  아직 최종 산출물이 아니며 전달하지 않았다. 이후8% snapshot 및 수렴여부 열,
  원문 표의 Bare Si4.09GPa와 문헌 폭감소지표 설명을 반영하고 변경 페이지도 재검수했다.
- PDF skill의 create marker는 이미1회 성공했다. 다시 실행하지 않는다.
- builder인자: --results, --output, --font, --bold-font, --as-of 필수.
  --snapshot prism_snapshot_audit_0545, --force-replay force_control_replay.
  최종 --label/--revision은 실제 마감/Git 상태로 바꾸고 재렌더 검수해야 한다.
  PDF 옆 `.manifest.json`은 사용한 실제 파일SHA와 페이지 경계를 기록한다.
  실제 PDF13쪽/북마크13개/읽은 원자료21파일 SHA일치도 확인했다.
- 번들 Python/Poppler 경로는 기존 report 작업과 같다. dependency 경로 조회 도구가
  오래 응답하지 않아 중단했으나 이미 확인된 번들 런타임에서 생성/렌더는 성공했다.

### 최종 무결성 기록 준비

- `record_initiation_manifest_v11.py`: v11 runner/테스트의 정적 local import를
  재귀 추적해 소스33개와 원문/DFT/모델 출처를 source_manifest.json에 기록했다.
  외부 package 버전은 기존 실제 실행환경 기록을 사용하며 완전한 environment lock은 아니다.
- 05:55대 중간 artifact hash2282파일/20,811,440bytes를 생성하고 즉시 재해시했다.
  `interim_artifact_manifest_0555.json`은 그 시점 기록이다. 이후 진행 문서/코드와
  source_manifest가 바뀌므로 마지막 결과의 무결성 인증으로 사용하지 않는다.
- 최종 실행이 멈춘 것을 확인한 뒤 같은 runner를 --filename artifact_manifest.json
  --include-mutable --label ACTUAL_FINAL_STATUS로 재실행한다. source_manifest도 갱신된다.
  로그/PID/pycache/다른 manifest JSON은 제외하며, 중간에는 running/progress/checkpoint도
  제외한다. 최종은 checkpoint/실패 prefix를 포함해 검사해야 한다. 파일 변경 중이면 거부한다.
- 새로운 v11 텍스트/코드의 개인 경로 검색은 일치0, git diff --check도 오류0.
  Git 게시 전 source/manifest/실제 저장blob의 일치와 최종 원격 상태는 다시 검증해야 한다.
