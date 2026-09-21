# 웨이퍼 연구 브랜치 안내

브랜치: **`silicon-wafer-research`**

공통 기반 `probability-pde-solver-v1`의 `5f878370`에서 분기했다.
Al과 Si를 같은 원자 에너지 → 자유에너지 → 확률 흐름의 구조로 연구한다.
웨이퍼 소재 UI, 공통 공간 해석 기반, Si 연구 코드와 원시 결과를 이 브랜치에 모았다.

## 현재 상태

| 항목 | 상태 | 코드·문서 |
|---|---|---|
| Al / Si 소재 선택, B/P/As/Sb·농도 설정과 저장 | 구현·검증 완료 | [소재 UI](app/MATERIAL_SELECTION.md), [materials.py](app/materials.py) |
| 공통 3D 메시·면하중·힘/토크 평형·프로젝트 저장 | 기존 검증 기반 포함 | [공간 역학](app/SOLID_MECHANICS.md), [메시·평형 검사](app/MESH_BALANCE_VERIFICATION.md) |
| Si 원자 에너지의 환경 표현·내부 원자 이완 | 정적 연구 검증 완료 | [환경 모델](solver_v1/silicon_environment_research.py), [검증 결과](results/silicon_unified_v1/COMPLETED_SUMMARY.md) |
| 균열 끝의 국소 분리·재결합·순차 전이 | 정적 연구 구현·검증 완료 | [공간 모델](solver_v1/silicon_crack_research.py), [국소 균열 결과](results/silicon_local_crack_v2/COMPLETED_SUMMARY.md) |
| 조건부 Gaussian 후보·주변 원자 기억·새 단시간 원자 운동 | 연구 구현·수치 대조 완료 | [v3 이론](solver_v1/SILICON_CONDITIONAL_DYNAMICS_V3.md), [새13 ps 결과](results/silicon_conditional_v3/COMPLETED_SUMMARY.md) |
| 비선형 열분포·재배열·새1.74 ns 원자 운동 | 부분 완료; 전체 PMF·크기 수렴 실패/미완료 보존 | [v4 이론](solver_v1/SILICON_FINITE_T_ENSEMBLES_V4.md), [완료/미완료 결과](results/silicon_thermal_v4/COMPLETED_SUMMARY.md) |
| 원자 에너지·내부 이완 탄성·공개 DFT 대조 | 2,475구조×2모델 평가 완료; SW/Tersoff 재료 미채택 | [v5 이론](solver_v1/SILICON_ATOMISTIC_FOUNDATION_V5.md), [실제 계산](results/silicon_atomistic_v5/COMPLETED_SUMMARY.md) |
| 실측 Si 물성·표면·도핑 에너지 적합 | 미완료 | [탐색 근거](results/silicon_wafer_feasibility/RESEARCH_MEMO.md) |
| finite-T PMF·이동도·실제 수명·Hz | 미보정 | [공통 이론](results/silicon_wafer_feasibility/UNIFIED_FATIGUE_THEORY.md), [국소 축약 조건](solver_v1/SILICON_LOCAL_CRACK_V2.md) |
| 생산 UI의 실제 Si 해석 | 비활성 | [실행 차단 계약](app/MATERIAL_SELECTION.md#실행-상태) |

도핑 농도는 도펀트 원자 수/cm³다. 캐리어 농도·비저항·실제 도핑 물성의 자동 환산이 아니다.
현재 공간 역학의 Al 해석 경로를 Si 예측으로 사용하지 않는다.

## 연구를 읽는 순서

1. [공통 확률이론과 Si 피로·피로한도](results/silicon_wafer_feasibility/UNIFIED_FATIGUE_THEORY.md)
2. [환경 에너지·좌표 축약 감사](results/silicon_unified_v1/COMPLETED_SUMMARY.md)
3. [국소 균열 좌표의 유도](solver_v1/SILICON_LOCAL_CRACK_V2.md)
4. [국소 균열 실제 계산과 검증](results/silicon_local_crack_v2/COMPLETED_SUMMARY.md)
5. [조건부 자유에너지 후보와 원자 동역학](results/silicon_conditional_v3/COMPLETED_SUMMARY.md)
6. [비선형 분포·구조 재배열·기억 검증과 중단 상태](results/silicon_thermal_v4/COMPLETED_SUMMARY.md)
7. [원자 에너지부터 다시 확인한 재료 기반](results/silicon_atomistic_v5/COMPLETED_SUMMARY.md)

기존 원자면 전체 동시 분리 계산은 대조군으로 남겨 둔다. 같은 최종 상태까지
전면의 네 결합을 순차 전이시키는 검증 경로의 최대 에너지는 1.075 eV,
동시 전이는 2.206 eV다. 이는 지정된 유한 시편·0 K·original-SW의 **에너지 장벽**이며
실측 웨이퍼 강도 또는 피로 수명이 아니다.

## 설치와 실행

```text
python -m pip install -r solver_v1/requirements.txt
python -m app.desktop_ui
```

`Pre / 사전 조건`에서 Si와 도핑 설정을 저장할 수 있다. Si 연구 계산은 아래 명령으로
별도 실행한다. 기존 연구 결과 폴더는 보존하고 새 출력 폴더를 지정한다.

```text
python -m results.silicon_wafer_feasibility.run_local_crack_audit --output .cache/wafer_reproduction
python -m results.silicon_wafer_feasibility.run_local_crack_sequence --calculation .cache/wafer_reproduction
python -m results.silicon_wafer_feasibility.validate_local_crack_reduction --calculation .cache/wafer_reproduction
```

선택적 외부 LAMMPS가 설치되어 있으면 같은 매개변수의 독립 에너지·힘 대조도 실행할 수 있다.

```text
python -m results.silicon_wafer_feasibility.check_spatial_sw_lammps --calculation .cache/wafer_reproduction --output .cache/wafer_reproduction/lammps_validation.json
```

위 계산은 MD, 열적 전이율, 물리 시간 보정이 아니다. UI에는 LAMMPS 설치가 필요 없다.

## 검증과 재현 자료

Si 전용 테스트의 매개변수는 저장소에 포함돼 있다. 전체 `solver_v1` 회귀에는
별도의 Al 참조 파일도 필요하다. 기존 다운로드 명령이 NIST 원본을 받아 SHA-256을
검사하고 `.cache/al-reference/Al99.eam.alloy`에 저장한다. 다른 내용의 기존 파일은
덮어쓰지 않는다. 이 파일은 시험용 참조이며 Si 에너지나 생산 PDE의 입력이 아니다.

```text
python -m pip install pytest
python -m solver_v1.reference_eam_targets --download
python -m pytest solver_v1/test_silicon_environment_research.py solver_v1/test_silicon_crack_research.py app/test_materials.py app/test_project_file.py -q
python -m pytest solver_v1 -q
python -m pytest app -q
python -m app.desktop_ui --smoke
```

- [국소 연구의 원시 검증 기록](results/silicon_local_crack_v2/validation_summary.json)
- [이 브랜치 분리 후의 검증](results/wafer_branch_setup/validation.json)
- [검증 요약·전체 테스트 목록·파일 해시](results/wafer_branch_setup/README.md):
  최종 solver 827 PASS, UI 189 PASS + 3 subtests, desktop smoke 통과.
- JSON/CSV/NPZ는 결과·원자 좌표·재현 입력이며, PNG/SVG는 실제 계산의 그림이다.
- `.gitattributes`가 hash-bound Si 소스와 매개변수의 LF를 유지한다.
  Windows checkout의 줄바꿈 변경 때문에 원본 매개변수 검증이 깨지지 않게 했다.
- 연구 결과에 남은 이전 브랜치·HEAD·작업 폴더 이름은 계산 당시의 이력이다.
  새 브랜치에서 재실행한 것처럼 과거 기록을 바꾸지 않았다.

## 다음 연구

1. 최신 “원자모델링부터” 범위: GAP/screened 에너지를 실제 실행하고,
   v5의 표면·균열 DFT 기준과 현재 재배열 최소·saddle을 독립 검증한다.
2. Glide/shuffle·표면·도핑의 필요한 물성을 검증한 뒤 조건부 분포·기억으로 이어간다.
3. v4의 profile56개/multibasin2개 체인, 전면8의 경계·dt 이슈는 미완료로 유지한다.
4. 검증된 좌표·PMF·이동도·파손 사건이 갖춰진 뒤 공통 확률 솔버에 연결한다.

이 브랜치에서는 정적 구현, 재료 적합, kinetic calibration, 생산 활성화 상태를
구분해 기록한다. 자세한 재개 위치는 [현재 인계](CURRENT_WORK_HANDOFF.md)를 따른다.
