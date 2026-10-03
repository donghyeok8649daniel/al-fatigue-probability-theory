# 실리콘 SW 베셀 에너지의 물성 보정

## 결과

베셀 표현을 유지한 순수 Si 정적 보정 연구다. 재료 채택, 최초 균열 개시,
유한 온도 자유에너지, 이동도, 초/Hz 및 피로수명은 검증하지 못했다.
Al/UI/생산 PDE/에너지/물리 시간 설정을 변경하지 않았다.

현재 후보는 `raw_runs/local_boundary_fit_results/selected_fit.json`의
profile 1, beta=0.3640142880725576, c0=-0.37357245686009144,
sigma=1.75 Å, gamma=1.5다. 학습 자료와 선언된 guard로 선택했다.
격자5.4610215037Å와 C11/C12/C44=153.28991/56.25009/72.17693GPa를
맞추며 원래 SW 벌크 에너지 guard도 통과했다. 전역 최적성 인증은 아니다.
PW91 Si(111) 47개를 손실에서 제외한 검사에서 힘 RMSE
1.0407327→1.0210356 eV/Å, 상대 E RMSE 0.2193716→0.1700248 eV/atom.
벌크 상대 E RMSE 0.0288054→0.0253361지만 벌크 힘은 악화했다.
외부 구조 일부와 진동수도 미보정이다. 이 파일의 `selected`는 계산의
선택 결과이며 재료 승인이나 생산 기본값이라는 뜻이 아니다.

추가 이웃16개/선호각9개/초기공동19개/각도경계9개는 모두 원래 SW 벌크
에너지 제한을 넘었다. 이 단계의 선택 파일은 미통과 후보 중 **진단용**이다.
후속형상으로전체양의각도경계30개중12개가통과했다. final_selection.json에
전체선택목록과학습손실·guard·정적anchor·수치파수gate를고정했다.
앞선angular_shape profile28은
C44가9.6%낮지만 학습·표면오차가 더작은 별도비교후보로 보존했다.
조사한 형상 범위의 결과를 전체 SW/베셀 계열 불가능 증명으로 읽지 않는다.

## 원자료의 역할

- 공개 GAP 구축용 원자료 2475구조/171815원자, 새 GAP·DFT·MD 실행 없음.
- PW91 dia489/surface_00110/surface_11012 = 511개 학습. 나머지1964개 제외.
- PW91/PBE/미상 XC를 분리. 에너지 기준은 고정 PW91 frame548의 per-atom E.
  열적 변형 구조이며 이완 기저상태·응집 에너지와 다르다.
- 알려진 0 K 공개 격자/탄성/EOS는 보정 기준. 탄성은 과거 유한 변형 적합이며,
  정확한 DFT zero-strain Hessian으로 인증된 값이 아니다.
- 별도 외부28구조/1020원자는 최종 형상 선택 후 검사했다. 손실·선택에 쓰지 않았다.
  앞선 후보에서도 이미 본 세트이므로 완전히 새 독립 검증은 아니다.
  정렬한 wrapped fraction/cell 지문에서 중복0. 모든 동등 셀·회전 검증은 아니다.
- 손실 척도는 선언된 가중치다. 측정 불확실성/신뢰구간이 아니다.

## 보존 파일

`selected_predictions.npz`는 12후보의 전체 에너지·원자별 힘과 원래 DFT 참조다.
`offsets`로 각 frame의 원자 구간을 찾는다. 구조 종류/XC/학습 여부는
`frame_metadata.json`, 예측의 원래 파일/SHA/진폭은 `prediction_cases.json`이다.

`raw_runs`에는 단계별 JSON/CSV, 당시 소스 스냅샷, 실패와 교정이 있다.
`design_gram_replay.json`은 반복 설계 배열의 quadratic Gram 통계로
학습 손실과 벌크 제한을 재생한다. 전체 raw design/features 배열은 기존
로컬 계산 캐시에 보존하고 `raw_cache_catalog.json`에 상대 위치·크기·SHA를
기록했다. 이 패키지는 그 큰 배열 전부를 포함하지 않는다. 새 구조 평가를
재현하려면 공개 원자 입력과 아래 계산 순서를 사용해야 한다.

## 실패·교정 결과를 읽는 순서

| 이전 자료 | 현재 판정 |
|---|---|
| angular_continuum_results | EOS sum/mean 정규화 오류. corrected_results가 권위 있는 결과. |
| density_results/coordination_results | 압축에서 EOS 불변 가정이 깨짐. eos_reaudit_results의 새864평가를 함께 읽음. 선택nu0 유지. |
| angular_shape_results/summary.json | 마지막 NumPy bool JSON 저장만 실패. 모든 raw 완료. recovery_results로 독립 복구. 이전 원 실행 시간 미상. |
| elastic_anchor_results | 실제 H 호출118. 별도 counter audit 확인. 현재 소스는 counter 이름 충돌 수정. |
| interface_audit_results | 최초 GPa 표시10배 과소. raw E/grad/H/root와 분리 일은 유효. recovery_results의 SI 교정 하중 사용. |
| independent_elastic_results | 첫 최소화 force gate 실패. force_root_results의 실제 force-root 결과와 구별. |

계면의 강체 법선 하중은 eV/Å³를 GPa로 바꾼다. 160.2176634와
표면 에너지 eV/Å²→J/m²의16.02176634를 혼동하지 않는다.
계면 전체 평면의 국소 안정성 상실·분리 일을 실제 균열 개시로 부르지 않는다.

## 계산을 하지 않는 독립 검증

Python과 NumPy만으로 패키지의 해시, 원자 예측 집계, 손실, SI 단위를 검사한다.
현재 디렉터리를 저장소 루트로 둔다.

```powershell
python results/silicon_sw_material_fit/verify_material_package.py
```

기존 독립 LAMMPS 예측과도 비교하려면:

```powershell
python results/silicon_sw_material_fit/verify_material_package.py --lammps-control results/silicon_atomistic_v5/dft_predictions.npz
```

이 검증은 새로운 구조·DFT·MD 평가가 아니다. 당시 source hash와 현재 소스를
동일하다고 가정하지 않는다. 과거 소스는 스냅샷의 정확한 바이트로 검증한다.

## 새 계산을 다시 하는 순서

최종 패키지의 모든 파일은 `package_manifest.json`에 실제 Git blob 바이트로
기록했다. Git에서 복원한 패키지는 아래 명령으로 파일 누락과 SHA-256을 검사한다.
이 검사는 새 물리 계산을 실행하지 않는다.

```text
python results/silicon_sw_material_fit/verify_package_manifest.py
```

`git_byte_test_results.xml`과 `tests_receipt.json`은 실제 Git 소스로 실행한 시험의
관측 기록이다. 원래 JUnit의 개인 호스트명만 공개 파일에서 제거했다.
원자료, 수치, 소스 해시는 바꾸지 않았다. PDF 검수 기록은 `report_qa.json`이다.

`build_material_package.py`는 기존 계산 캐시에서 원자료·예측·Gram 통계를
새 폴더로 다시 모으는 도구다. 시험을 다시 실행하거나 최종 보고서 검수를
자동으로 완료하는 도구는 아니다. 시험 기록·독립 검증·그림·보고서는 각각의
실행과 검수가 필요하다. 최종 배포 파일의 정확한 바이트 재현은 Git과 manifest를 따른다.

표면 힘만 맞춘 계수가 벌크 조건을 벗어나는지는 `anchor_force_tradeoff.json`에
기록했다. 고정 형상의 저장된 3×3 행렬을 재생한 검사이며 새 원자 계산이 아니다.
아래 명령의 `FRESH_JSON`은 아직 없는 새 출력 파일이어야 한다.

```text
python results/silicon_sw_material_fit/diagnose_anchor_tradeoff.py --output FRESH_JSON
```

`bulk_optical_diagnostic.json`은 저장된 상대 변위 Hessian과 Γ 광학 모드를
원자 두 개의 환산질량 m/2로 대조한다. 내부 이완 전후의 전단 곡률도 분리한다.
같은 질량을 가정해 문헌 진동수에서 환산한 곡률은 새 DFT 계산이 아니다.
재생은 `diagnose_bulk_optical.py --output FRESH_JSON`으로 수행한다.

`optical_census.json`에는 기존 양의 각도 경계 30개와 조건을 통과한 12개의
저장 곡률에서 환산한 Γ 광학 모드를 담았다. 후보 선택은 그대로 유지했다.
재생은 `diagnose_optical_census.py --output FRESH_JSON`으로 수행한다.
각 후보의 전체 밴드를 새로 계산한 작업이나 모든 형상의 불가능성 증명이 아니다.

Python, NumPy, SciPy, ASE, pytest가 필요하다. 실제 버전과 관측된 테스트 명령은
`tests_receipt.json`에 기록한다. 출력 폴더는 매번 새 경로여야 한다.
동일 결과 재현과 전체 최적화 재실행을 구별한다.

공개 archive DOI10.17863/CAM.65004의 `Si_PRX_GAP.zip`을 먼저 확보한다.
SHA256 `1d3efe976c53bcd2889c0556f2d600c26bf4603522b5c92fd8633dff7b7ab5c2`.
member `gp_iter6_sparse9k.xml.xyz` SHA256
`74ac5f387aa0149ec9ea8f3e275d1087d7287db68c15ba8bc4c28a450de67e2a`.
실행 시 `--archive`로 그 경로를 지정한다. 자동 다운로드에 의존하지 않는다.

아래 명령의 BASE/GUARD/STRESS/ANCHOR/ANGLE은 사용자가 정한 새 출력 폴더다.
`REF`는 이 패키지 `primary_sources/test-results`의 CASTEP bulk JSON이다.

```text
python -m solver_v1.run_silicon_sw_material_fit --archive ARCHIVE --output BASE
python -m solver_v1.run_silicon_sw_guarded_fit --archive ARCHIVE --baseline BASE --output GUARD
python -m solver_v1.run_silicon_sw_stress_fit --archive ARCHIVE --baseline BASE --output STRESS
python -m solver_v1.run_silicon_sw_elastic_anchor --archive ARCHIVE --baseline BASE --profiles STRESS --reference REF --output ANCHOR
python -m solver_v1.run_silicon_sw_angular_fit --archive ARCHIVE --baseline BASE --anchors ANCHOR --profiles STRESS --reference REF --output ANGLE
python -m solver_v1.run_silicon_sw_angular_shape --archive ARCHIVE --baseline BASE --prior ANGLE --reference REF --output SHAPE
python -m solver_v1.run_silicon_sw_angle_preference --reference REF --positive-boundary --sigmas 1.8 --gammas 1.2 1.5 1.8 --output JOINT_PREFLIGHT
python -m solver_v1.run_silicon_sw_extended_fit --archive ARCHIVE --baseline BASE --reference REF --preflight JOINT_PREFLIGHT --output JOINT_FIT
python -m solver_v1.run_silicon_sw_external_audit --candidate JOINT_FIT/selected_fit.json --dataset results/silicon_sw_material_fit/primary_sources/tests/testing_database_force_error/testing_database.xyz --archive ARCHIVE --baseline BASE --output EXTERNAL
```

정확한 인자는 각 모듈의 `--help` 및 당시 소스 스냅샷으로 확인한다.
과거 실행은 당시 source를 사용했다. 현재 source에는 알려진 구현 오류 교정이
포함되어 있어 실패 상태 자체를 재현하려면 해당 snapshot을 사용해야 한다.
Gram 재생은 이 전체 재적합을 대체했다고 보고하지 않는다.

### 최종 후보의 전체 형상 선택 재실행

새 디렉터리 `CACHE`와 위의 `BASE`, `ARCHIVE`, `REF`를 사용한다.
아래 네 행을 각각 preflight로 계산하고 같은 이름의 fit 결과를 만든다.

| preflight 이름 (fit은 preflight를 fit으로 바꿈) | --sigmas | --gammas |
|---|---|---|
| angular_boundary_preflight_results | 1.6 2.1 2.8 | .9 1.2 1.5 |
| refined_boundary_preflight_results | 1.8 | 1.2 1.5 1.8 |
| local_boundary_preflight_results | 1.75 1.85 1.9 | 1.35 1.5 1.65 |
| guard_boundary_preflight_results | 1.625 1.65 1.675 | 1.35 1.5 1.65 |

```text
python -m solver_v1.run_silicon_sw_angle_preference --reference REF --positive-boundary --sigmas SIGMAS --gammas GAMMAS --output CACHE/PREFLIGHT_NAME
python -m solver_v1.run_silicon_sw_extended_fit --archive ARCHIVE --baseline BASE --reference REF --preflight CACHE/PREFLIGHT_NAME --output CACHE/FIT_NAME
python results/silicon_sw_material_fit/select_final_shape.py --cache CACHE
```

첫 두 명령을 표의 네 행마다 실행한 뒤 selector를 한 번 실행한다.
전체30개/eligible12개의 같은 학습/EOS 손실로 선택하며, 제외 오차 파일을
읽지 않는다. 네 anchor와 벌크 guard, 수치 파수 gate를 통과해도
전역 shape 최적성이나 재료 채택을 뜻하지 않는다.

선택 이후에는 `run_silicon_sw_external_audit`, `run_silicon_sw_final_phonon`,
`run_silicon_sw_finite_strain --candidate-only`, `run_silicon_sw_interface_audit
--candidate-only`, `run_silicon_sw_final_convergence`를 해당 선택 파일에 실행한다.
실제 인자와 입력 SHA는 raw_runs의 source/input snapshot과 binding에 있다.

## 보고서와 그림 재생

그림에는 Matplotlib, 보고서에는 ReportLab 및 한국어 TTF가 필요하다.
각 도구는 `--package`로 이 디렉터리를 지정할 수 있다.

```text
python results/silicon_sw_material_fit/create_figures.py
python results/silicon_sw_material_fit/create_report.py --output report.pdf --font KOREAN_TTF --bold-font KOREAN_BOLD_TTF
```

보고서는 렌더링 후 모든 페이지의 수식·표·줄바꿈·그래프를 확인했다.
생성 성공만으로 렌더 검수가 완료되었다고 판단하지 않는다.

## 원 출처

- Bartók 외, PRX8,041048(2018), DOI10.1103/PhysRevX.8.041048.
- libAtoms silicon-testing-framework, DOI10.5281/zenodo.1250555,
  고정 tree fc252cb7d41df7e2bc672d614f3d76a40c9f2ecb.
- 과거 matscipy path revision1202c4e9. 실제 benchmark 환경의 commit pin은 미확보.
- primary_sources의 manifest에는 실제 다운로드 URL, Git blob 및 SHA256이 있다.
  외부 코드 전체를 이 패키지에 복제하지 않았다.
