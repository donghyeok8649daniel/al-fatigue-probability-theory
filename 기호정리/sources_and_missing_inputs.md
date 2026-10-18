# 공개 출처와 미정 입력의 원자료

## 확인해서 채운 기준 입력

[입력 JSON](input_evidence.json)의 source_reference_only 기록 20개는 공개 저장소에서 확인했습니다. native N=1 UI와 일반 ModelParams 기본값은 별도 프로파일입니다. Grid/step 기본값은 실제 사용 run의 선택을 대신하지 않습니다.

공개 kinetic 파일이 선언한 L0=2.863782463805517e-10 m, E0=1.602176634e-19 J와 같은 파일의 calibrated=false, temperature/mobility/t0/source null을 함께 보존했습니다. plane/finite-band kinetic 연구의 존재와 local production coordinate 보정 미완료를 구별합니다. TwoRowLJ의 실제 물리 척도나 실제 기하를 이 값으로 지정하지 않습니다.

## 실제 결합 입력 15개

| 항목 | 요구 단위 | 값 | 확정에 필요한 원자료 |
|---|---|---|---|
| model/parameter/energy representative | `식별자` | `null` | 같은 모형 ID·parameter fingerprint·per-atom/cell/interface counting과 채택 근거 |
| 실제 결합 모형 길이 척도 | `m` | `null` | 같은 좌표의 측정 길이·기준 lattice/상태·단위와 대응식 |
| 실제 결합 모형 에너지 척도 | `J` | `null` | 같은 대표계 에너지·thermal normalization·보정 근거 |
| 3D cell vectors/원자 위치·유한 집단 영역 | `m` | `null` | 3D primitive vectors 또는 원자 좌표, 경계·두께·partition·집단 지지영역 |
| 실제 재료 셀 부피 | `m^3` | `null` | 위 3D 자료와 동일한 대표계/기준 또는 현재 상태의 부피 측정 |
| 실제 interface 면적 | `m^2` | `null` | 계면의 두 접벡터 또는 동일 영역의 면적/법선·방위 |
| 실제 chi/strain-work mapping | `1` | `null` | 좌표 jump→strain의 orientation/기하, 공액 stress·strain·volume convention |
| 같은 방향의 탄성계수 | `Pa` | `null` | 같은 소재·온도·방위·상태의 modulus와 측정/보정 원자료 |
| 실제 운전/보정 온도 | `K` | `null` | 같은 specimen/MD record의 온도 로그 및 thermostat/ensemble |
| local a/s 이동도와 공통 t0 | `m^2/(J s); s` | `null` | 같은 집단좌표·PMF·열항의 signed response 또는 covariance/확산 원자료, 오차와 coarse-time 검증 |
| 실제 하중·시간·주기 단계 | `Pa; s; cycles` | `null` | controller raw trace, 센서 단위/교정, cycle-to-stage map·sampling clock |
| 실제 초기법칙·raw cohort·원점 | `density/mass` | `null` | 초기 준비·선별·support, 실제 grid와 quadrature, preload·초기 absorber 기록 |
| opening mask와 실측 균열 사건 대응 | `event` | `null` | 동일 stage의 원래 mask/geometry, 검사 last-negative/first-positive, 검출한계와 DCPD 교정 |
| 실제 상관 영역·집단 법칙 | `m^2; joint law` | `null` | 동일 하중/조직의 공간 위치·correlation/공동 사건 자료; 독립성 검증 |
| 실제 quadrature/runtime/source 오차 | `해당량 단위` | `null` | 같은 run의 grid/step/order/mask/source identity와 오차 근거 |

## 공개 원격 출처

아래 공개 GitHub 코드/config/문서만 이 공개본의 근거로 사용했습니다. 로컬 실험·메일 기록과 미공개 보조 이론/분석 기호는 포함하지 않습니다. 기호와 기준값을 정리하는 문서이며 실제 소재·시계·피로수명을 새로 보정한 결과가 아닙니다.

| ID | 자료 | 사용 범위 |
|---|---|---|
| S1 | [solver_v1/model.py](https://github.com/donghyeok8649daniel/al-fatigue-probability-theory/blob/c8ccff7553cbb2e66b0c37a76ee8a24808478824/solver_v1/model.py) | ModelParams 및 TwoRowLJ의 좌표·에너지·strain |
| S2 | [app/solver_adapter.py](https://github.com/donghyeok8649daniel/al-fatigue-probability-theory/blob/c8ccff7553cbb2e66b0c37a76ee8a24808478824/app/solver_adapter.py) | canonical_model_params:237-246의 N=1 UI 기준 |
| S3 | [solver_v1/probability_pde_2d.py](https://github.com/donghyeok8649daniel/al-fatigue-probability-theory/blob/c8ccff7553cbb2e66b0c37a76ee8a24808478824/solver_v1/probability_pde_2d.py) | density, cell_volume, SG face/generator, 초기화와 흡수 |
| S4 | [solver_v1/PHYSICAL_TIME_AND_MOBILITY.md](https://github.com/donghyeok8649daniel/al-fatigue-probability-theory/blob/c8ccff7553cbb2e66b0c37a76ee8a24808478824/solver_v1/PHYSICAL_TIME_AND_MOBILITY.md) | 단위 변환 및 물리 시계 미보정 상태 |
| S5 | [solver_v1/data/aluminum_kinetic_calibration.json](https://github.com/donghyeok8649daniel/al-fatigue-probability-theory/blob/c8ccff7553cbb2e66b0c37a76ee8a24808478824/solver_v1/data/aluminum_kinetic_calibration.json) | 선언된 단위 척도와 calibrated=false/null fields |
| S6 | [solver_v1/energy_model_registry.py](https://github.com/donghyeok8649daniel/al-fatigue-probability-theory/blob/c8ccff7553cbb2e66b0c37a76ee8a24808478824/solver_v1/energy_model_registry.py) | 모형 ID와 각 energy/geometry scale의 적용 범위 |
| S7 | [paper/fatigue_probability_main_v4.tex](https://github.com/donghyeok8649daniel/al-fatigue-probability-theory/blob/c8ccff7553cbb2e66b0c37a76ee8a24808478824/paper/fatigue_probability_main_v4.tex) | 일반 correlated 2N Smoluchowski 이론; 과거 N=3 screening 별도 |
| S8 | [AGENTS.md](https://github.com/donghyeok8649daniel/al-fatigue-probability-theory/blob/c8ccff7553cbb2e66b0c37a76ee8a24808478824/AGENTS.md) | N=1 역사적 a0/kappa, 기하·통계·시계 구별 |
| S9 | [README.md](https://github.com/donghyeok8649daniel/al-fatigue-probability-theory/blob/c8ccff7553cbb2e66b0c37a76ee8a24808478824/README.md) | 현행 PDE 경로와 역사적 MC 검사의 구별 |
