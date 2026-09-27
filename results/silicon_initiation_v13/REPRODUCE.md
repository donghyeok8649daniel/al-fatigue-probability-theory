# 현재 Si 개시 연구의 재현과 원자료 검증

이 문서는 실행 방법이다. 실행 완료 여부는 WORKING_STATUS와 각 summary를
따른다. 기본 큐·전체QE 보완·새 무축력 상태의649차원 Hessian과 독립 재검증을
완료했다. 큰 산소 구조는 원래 시도와1회재시도 모두 모델 할당 전 실패했다.
최종 수치 집계는 execution_accounting.json, 실제 실행 이력은
execution_timeline.json, 판정과 한계는 COMPLETED_SUMMARY.md에 있다.
기준 브랜치는 silicon-wafer-research, 시작 revision은
c3a93f5d2e4df6e21a68d6080ed99cf6aee840bc다.

모든 명령은 Si 저장소 root에서 실행한다. MODEL, MPA_MODEL, CP2K_SOURCE는
로컬 파일 경로 자리표시자다. 새 실행은 NEW_* 출력 폴더를 사용한다.
실행 중인 큐와 완료 결과를 중복 시작하거나 덮어쓰지 않는다.

## 실제 8% Hessian의 안전한 재개

```text
python results/silicon_wafer_feasibility/check_prism_dense_hessian_v13.py --state results/silicon_initiation_v11/intact_prism_360_linesearch_v2/state_004/raw.npz --geometry results/silicon_initiation_v11/intact_prism_360_linesearch_v2/geometry.npz --model MODEL --output NEW_LOADING8 --ensemble displacement --check-modes 4 --max-seconds 6000 --resume-from results/silicon_initiation_v12/dense_loading8 --resume-sha256 3495625f164cb1557b914cdb0abc4a964af7e79a48f6dda77a4e9b81ad69f71e
```

기존 539/648행의 전체 열을 보존한다. 입력·모델·기저·기울기를 대조하고
0/269/538행을 새로 계산해 일치한 뒤 109행을 추가했다. 최종 행렬 뒤에는
낮은 4모드와 무작위 2방향의 두 간격 힘·에너지 차분을 따로 검사한다.
복귀 8%는 state를 prism_unload_10_to_8/raw.npz로, 10%는
prism_10pct_continuation/raw.npz로 바꾸며 두 resume 옵션은 사용하지 않는다.
모두 v11 결과 아래의 상태 파일이다. 기존 100 MPa 행렬은 새 계산 없이
v12에서 byte 동일하게 복사했으며 reuse_provenance.json에 근거가 있다.

displacement protocol의 external_nominal_stress_GPa=0은 사용하지 않은
force 옵션의 기본값이며 이 상태의 반력이0이라는 뜻이 아니다. 실제 basis는
1080×648이고 고정부 행은0이다. 공통 basis 설명 문자열에 나오는 grip
연장은 force ensemble의649번째 좌표에만 존재한다. 실제 배열·ensemble·
기저 생성 코드를 기준으로 해석하고 당시 실행 원자료를 바꾸지 않는다.

## 저장 행렬, 국소 근사와 경로 초기화 기하

```text
python results/silicon_wafer_feasibility/audit_dense_hessians_v12.py --results results/silicon_initiation_v13 --output NEW_MATRIX_REPLAY
python results/silicon_wafer_feasibility/analyze_local_harmonic_v12.py --results results/silicon_initiation_v13 --geometry results/silicon_initiation_v11/intact_prism_360_linesearch_v2/geometry.npz --output NEW_HARMONIC
python results/silicon_wafer_feasibility/audit_prism_chord_geometry_v13.py --before results/silicon_initiation_v11/intact_prism_360_linesearch_v2/state_004/raw.npz --after results/silicon_initiation_v11/prism_unload_10_to_8/raw.npz --geometry results/silicon_initiation_v11/intact_prism_360_linesearch_v2/geometry.npz --output NEW_CHORD
```

이 세 후처리는 새 원자 에너지/힘 평가를 하지 않는다. 모두 완료했고,
수치 결과와 한계는 각 summary에 있다. 조화 근사는 같은 측도에서의 한 국소 quadratic
well 근사다. 실제 유한온도 F나 활성화 장벽이 아니다. 기하 경로도 원자
번호를 잇는 직선 보간일 뿐이며 경로 좌표는 물리 시간이 아니다.

완료한 두8% 상태의 같은 Cartesian 방향 및 유한변위 외삽 비교:

```text
python results/silicon_wafer_feasibility/compare_same_grip_modes_v13.py --results results/silicon_initiation_v13 --geometry results/silicon_initiation_v11/intact_prism_360_linesearch_v2/geometry.npz --output NEW_MODE_COMPARISON
python results/silicon_wafer_feasibility/replay_same_grip_modes_v13.py --results results/silicon_initiation_v13 --comparison NEW_MODE_COMPARISON --output NEW_MODE_REPLAY.json
```

행렬 최저값의 비교와 같은 변형 방향의 곡률을 구별한다. 제곱 내적,
부분공간 주각도와 chord 투영은 기하학적 양이며 확률·경로·장벽이 아니다.
독립 대조73개는 배열 산술 검사이고 추가 pytest 또는 새 원자 계산이 아니다.

## Si/O 모델 비교와 추가 QE 감사

run_large_oxide_audit_v12.py의 --min-atoms 320 --max-atoms 1200
--require-oxygen 조건은 모델 할당 전12GB 메모리 기준에 걸려0구조에서 중단됐다.
계산 종료 뒤 idle12.53GB에서 새 폴더로1회재시도했으나, 라이브러리와 입력을
준비한 뒤에도 같은 사전 기준에 걸렸다. 두 시도 모두0평가이고 실패 원자료를
보존했다. 세 번째 시도를 하거나 기준을 낮추지 않았다.
compare_mpa_model_v12.py의 CP2K135/QE64 선택과 원자료 재검증은 완료했다.
실제 선택은 예측 전에
mpa_sampling_scope/source_raw_manifest.json 및 선택 JSON으로 고정했다.
동일조성 에너지 차이와 힘 성분 오차를 구분한다. 전체 공개 자료도 실제
웨이퍼 상태의 무작위 표본이라고 가정하지 않는다.

complete_mpa_qe_v13.py는 위 비교의199개 결과를 정확히 재사용하고
선택하지 않았던QE1095개와 표준 대조9개를 실제로 추가했다. 전체1294기록은
CP2K135+QE1159이며 독립 재검증까지 완료했다. 실행 명령은 다음과 같다.

```text
python results/silicon_wafer_feasibility/complete_mpa_qe_v13.py --cp2k results/silicon_initiation_v11/oxidized_surfaces_max320 --qe results/silicon_initiation_v11/oxide_mace_mtpu_1159_v2 --model MPA_MODEL --reuse results/silicon_initiation_v13/mpa_comparison --output NEW_FULL_QE --max-seconds 1200 --deadline-utc USER_APPROVED_UTC_DEADLINE
python results/silicon_wafer_feasibility/replay_mpa_qe_v13.py --cp2k results/silicon_initiation_v11/oxidized_surfaces_max320 --qe results/silicon_initiation_v11/oxide_mace_mtpu_1159_v2 --reuse results/silicon_initiation_v13/mpa_comparison --result NEW_FULL_QE --output NEW_FULL_QE_REPLAY
```

이번 실행의 종료 guard는2026-09-28T01:19:39Z였다. 재실행 때 지난 시각을
그대로 쓰지 않고 실제 승인된 종료 시각을 넣는다. 최종 데이터와 source hash를
검사하는 것만 필요하면 원자 계산을 반복하지 않고 아래 manifest 검증을 사용한다.

비선형 국소72점도 완료했다. 다음 독립 후처리와 힘 방향 분해를 실행했다.

```text
python results/silicon_wafer_feasibility/replay_anharmonic_v13.py --results results/silicon_initiation_v13 --probes results/silicon_initiation_v13/anharmonic_probes --output NEW_ANHARMONIC_REPLAY.json
python results/silicon_wafer_feasibility/summarize_anharmonic_scope_v13.py --results results/silicon_initiation_v13 --output NEW_ANHARMONIC_SCOPE
```

## 일정 축력 확장과 무축력 복귀

```text
python results/silicon_wafer_feasibility/replay_grip_augmented_v13.py --results results/silicon_initiation_v13 --augmented results/silicon_initiation_v13/grip_augmented --geometry results/silicon_initiation_v11/intact_prism_360_linesearch_v2/geometry.npz --output NEW_GRIP_REPLAY.json
python results/silicon_wafer_feasibility/run_force_controlled_prism_v11.py --state results/silicon_initiation_v11/prism_unload_10_to_8/raw.npz --geometry results/silicon_initiation_v11/intact_prism_360_linesearch_v2/geometry.npz --model MODEL --output NEW_ZERO_RETURN --stresses 0 --fmax 1e-5 --max-steps 600 --max-seconds 1800
python results/silicon_wafer_feasibility/audit_zero_force_return_v12.py --geometry results/silicon_initiation_v11/intact_prism_360_linesearch_v2/geometry.npz --baseline results/silicon_initiation_v11/force_controlled_prism_360_tight/state_000 --returned NEW_ZERO_RETURN/state_000 --output NEW_ZERO_REPLAY.json
```

grip 재검증은 새 모델0회다. 두 번째 명령은 새 에너지·힘 평가를 하며 실제147회로
수렴했다. 무축력 복귀가 수렴해도 전체 곡률 검사를 대신하지 않는다. 새649차원
계산은 아래와 같이 실행해 전체649행·6독립방향과 별도 원자료 재검증까지
완료했다. 최저 곡률0.0150298619131eV/Å², 음의 고유값0이다. 정확한
정지점·전역 안정성·실제 개시의 인증과 구별한다.

```text
python results/silicon_wafer_feasibility/check_prism_dense_hessian_v13.py --state results/silicon_initiation_v13/return_zero_force/state_000/raw.npz --geometry results/silicon_initiation_v11/intact_prism_360_linesearch_v2/geometry.npz --model MODEL --output NEW_ZERO_HESSIAN --ensemble force --stress 0 --check-modes 4 --max-seconds 5400 --deadline-utc USER_APPROVED_UTC_DEADLINE
python results/silicon_wafer_feasibility/audit_dense_hessians_v13.py --results results/silicon_initiation_v13 --states return_zero --output NEW_ZERO_MATRIX_REPLAY
```

두 번째 명령은 최종 dense_return_zero 폴더를 읽는 검증 명령이다. 새 재실행을
검증하려면 새 결과 root에 dense_return_zero라는 명칭으로 저장하고 그 root를
--results에 지정한다. 이름이 force100인 이전 자료와 바꿔치기하지 않는다.

## 시험 기록

validation.json은 이미 실행한 17+10+11개 시험을 서로 다른 시험 파일과
소스 SHA로 다시 결합한 기록이다. 이 결합은 새 pytest 실행이 아니다.
기존 v12의 63개 시험도 이번 시험 수에 합산하지 않는다. 세 원래 기록은
보존하며, 첫 resume 시험의 TEMP setup 실패와 성공한 재실행도 구분한다.
재실행이 필요한 경우 새 임시 폴더로 다음을 사용한다.

```text
python -m pytest -q -p no:cacheprovider --basetemp=NEW_TEMP solver_v1/test_hessian_resume.py solver_v1/test_silicon_anharmonic_replay_v13.py solver_v1/test_mpa_qe_replay_v13.py
```

구현 시험은 실제 원자 계산, 재료 적합성, 균열 개시 인증과 다르다.

## 새 계산과 재사용 집계

audit_execution_accounting_v13.py는 저장된 완료 수치와 재사용 명세를 읽어
새 forward 평가와 자동미분 행을 별도로 합산한다. 진행 중 단계는 null로
표시해 완료 소계에서 제외한다. 모델을 다시 실행하거나 새 pytest를 하지 않는다.
2026-09-28 08:10KST 초안 집계는 새 forward1707/AD1411이며, 아직 진행 중인
무축력 복귀 Hessian은 포함하지 않았다. 이는 최종 집계가 아니다.

```text
python results/silicon_wafer_feasibility/audit_execution_accounting_v13.py --results results/silicon_initiation_v13 --output NEW_ACCOUNTING.json
```

모든 계산의 완료·실패·부분 결과를 확인한 뒤 --final을 붙여
execution_accounting.json에 최종 기록한다. 산소 재시도가 있으면 새 결과 폴더도
함께 확인한다. 실행 도중 중단된 부분 원자료를 임의 0회로 처리하지 않도록,
지원하지 않는 부분 상태에서는 집계기가 중단하고 명시적인 감사가 필요하다.
재사용한 기존100MPa의 과거 호출수·시간은 새 합계에 넣지 않는다.
이 합계는 독립 시편 수나 물리 시간, 재료 검증 성공 횟수가 아니다.

최종 집계는 새 forward1756/AD2060이다. 두 메모리 실패의 새 모델 호출은
0이며 기존 재사용 행·예측은 새 합계에 넣지 않았다. 최종 집계기의 재실행은
모델 평가나 pytest를 다시 수행하지 않는다.

## 최종 보존과 독립 확인

모든 계산·문서 작성 프로세스가 종료하고 완료·부분·실패를 기록한 뒤에만:

```text
python results/silicon_wafer_feasibility/record_initiation_manifest_v13.py --repo . --results results/silicon_initiation_v13 --baseline-model MODEL --comparison-model MPA_MODEL --cp2k-source CP2K_SOURCE --label FINAL_LABEL --final
python results/silicon_wafer_feasibility/verify_initiation_bundle_v12.py --repo . --results results/silicon_initiation_v13
```

새 recorder는 v12 artifact와 그 upstream v11/v9 원자료의 기존 해시까지
다시 대조한다. 실제 실행한 코드와 준비만 한 코드가 함께 포함되므로
source_manifest를 실행 내역으로 읽지 않는다. 외부 모델/자료는 파일명·
SHA·공개 URL만 저장하고 개인 경로를 커밋하지 않는다.

기존 verifier는 공통 manifest schema를 읽는다. 선택적으로 --revision INDEX
또는 최종 commit을 붙이면 실제 Git blob도 검사한다. CRLF 실행 바이트와
기록된 LF 전송 표현을 구별한다. 원자료 해시 검증은 수치·물리 검증을
대신하지 않는다. 최종 manifest 뒤 결과 파일을 바꾸면 manifest를 갱신하고
다시 검증해야 한다.

## 현재 설명 보고서

PDF는 공통 확률 법칙과 도핑/산화·원자 시편의 의미부터 설명하는 새 초안이다.
버전 관리할 생성기는 build_current_model_report_v13.py다. 기존 설명의17개
원자료 해시를 확인하고 현재 결과를 읽는다. 필요한 초기 설명 JSON과 두 그림은
report_inputs/baseline에 포함했다. 개인 캐시를 재현 입력으로 요구하지 않는다.
ReportLab과 한국어 TTF 글꼴 두 개(본문·굵은체)가 필요하다. 글꼴 파일은 배포하지
않으며 실제 사용한 파일 해시는 생성되는 PDF manifest에 기록한다.

```text
python results/silicon_wafer_feasibility/build_current_model_report_v13.py --output NEW_TMP_PDF --font KOREAN_TTF --bold-font KOREAN_BOLD_TTF --as-of "YYYY-MM-DD HH:MM KST"
```

최종 생성에는 추가로 --final-addendum FINAL_ADDENDUM.json이 필요하다.
이 JSON의 sections는 버전 관리할 report_inputs/final_sections.json과 같아야
하며, revision은 원자료·소스 검증을 마친 실제40자리 Git commit이다.
reviewed_actual_outcomes/git_status_verified/manifest_verified는 실제 완료한
검증을 확인한 후에만 true로 기록한다. 준비된 계획이나 예상 커밋을 넣지 않는다.
--as-of는 당시 표시 시각으로 고정할 수 있지만, PDF 메타데이터 생성 시각까지
byte 동일하다고 가정하지 않는다.28쪽 초안의 캐시 방식/버전관리 방식 사이에서
추출 텍스트와 페이지 content stream이 모두 같음을 확인했다.
이후 완료된 무축력 안정성·실행 집계·실패를 담은3쪽을 추가해30쪽 초안을
생성하고 새3쪽을 시각검수했다. --preview-outcomes는 이 결과 쪽을 포함해
커밋 전 초안으로 확인하는 옵션이다. 최종본·Git 검증 완료를 주장하지 않는다.
실제 최종 PDF 생성 뒤에는 모든 쪽을 다시 렌더링하고 시각검수한다.

보고서의 최종 계산·실패·부분 범위를 확정하고 전쪽 렌더링·시각검수를 마친다.
과거 crack seed, 정적 분리 일, 국소 양의 곡률을 미리 균열이 없는 Si 시편의
첫 균열 개시로 바꾸지 않는다.
