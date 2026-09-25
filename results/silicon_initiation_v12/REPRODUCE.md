# v12 재현 및 검증

저장소 root에서 실행한다. 버전과모델/source SHA는`source_manifest.json`에
있다. 아래BASELINE_MODEL,MPA_MODEL,CP2K_SOURCE는받은파일의로컬경로다.
개인경로나가중치를Git에넣지않았다. 출력은새폴더를지정하고부분결과를보존한다.

## 저장된 결과 검증

```text
python results/silicon_wafer_feasibility/verify_initiation_bundle_v12.py --repo . --revision HEAD
```

테스트/원자계산 재실행이 아니다.현재소스·시험소스·결과·변하지않은v11입력의
SHA와Git blob을검증한다. Windows실행bytes와Git LF표현은별도결합한다.

## 실제 완료한 전체 Hessian

```text
python results/silicon_wafer_feasibility/check_prism_dense_hessian_v12.py --state results/silicon_initiation_v11/force_controlled_prism_360_tight/state_001/raw.npz --geometry results/silicon_initiation_v11/intact_prism_360_linesearch_v2/geometry.npz --model BASELINE_MODEL --output NEW_FORCE100_OUTPUT --ensemble force --stress 0.1 --max-seconds 6000
```

실제4행시험은같은조건의`--max-rows 4`였다.
`validate_hessian_probe_v12.py REPOSITORY BASELINE_MODEL NEW_JSON_OUTPUT`은
당시독립표준calculator검사의byte동일소스사본이다. 결과사본은
`hessian_probe/standard_calculator_replay.json`이며새실행으로집계하지않는다.

`audit_dense_hessians_v12.py --results RESULT_ROOT --output NEW_REPLAY`는
저장행렬·단위·독립차분기록을대조한다. 이번`dense_first_completed_replay`의
complete=false는다른3상태부재이며100MPa검증실패가아니다.

8%부분조건은v11`intact_prism_360_linesearch_v2/state_004/raw.npz`와displacement
ensemble이다.539행을보존했으며현재러너는resume를지원하지않는다.
`closing_records_source_snapshot.py`는실제종료기록수집코드의사본이다.
원래위치는작업cache였다. 과거전원이벤트는새실행에서똑같이만드는재현대상이
아니다.부분블록감사는저장행렬좌상539x539를대칭화해`numpy.linalg.eigvalsh`로
독립재계산할수있다.양의부분블록은전체안정성증거가아니다.

## 구조와 실제12fit

```text
python results/silicon_wafer_feasibility/analyze_prism_structure_v12.py --help
python results/silicon_wafer_feasibility/fit_oxide_radial_delta_v12.py --help
python results/silicon_wafer_feasibility/audit_total_energy_errors_v12.py --help
```

입력은v11`oxidized_surfaces_max320`과`oxide_mace_mtpu_1159_v2`다. 각family의
5조성제외+전체fit을실행했다.초기radial실행소스는
`oxide_radial_delta/runner_snapshot.py`에별도로있다.이후angular선택을추가한
최종소스와당시실행소스를구분한다.새MACE/DFT/MD평가는없었고두후보미채택이다.

## 확률·합성 정확해

```text
python results/silicon_wafer_feasibility/record_validation_v12.py --output NEW_VALIDATION_JSON
python results/silicon_wafer_feasibility/audit_charge_basin_limit_v12.py --help
python results/silicon_wafer_feasibility/audit_charge_boundary_continuum_v12.py --output NEW_CONTINUUM_OUTPUT
python results/silicon_wafer_feasibility/replay_charge_continuum_v12.py --result NEW_CONTINUUM_OUTPUT
```

실제통합pytest63PASS/5.16s/ASE경고2개.전하반례는baselineGit86360a0소스로
재현했다.v9감사재실행5파일이기존과byte동일했다.보존·정확해검증은
실제Si전이율/clock보정과다르다.

## 준비만 하고 미실행한 코드

grip추가,mode위치,국소Gaussian,비선형변위점,큰구조,MACE-MPA비교,
0축력복귀와각replay러너는이번에실행하지않았다.큰구조를재개할때
`run_large_oxide_audit_v12.py --min-atoms 320 --max-atoms 1200 --require-oxygen`
조건은69개다.옵션없는168개와다르며pureSi99개는보류했다.MPA199표본은
CP2K그룹+조성별최대12,QE그룹별최대16을원자료순서에서고르게선택한다.

## PDF

```text
python results/silicon_wafer_feasibility/build_initiation_report_v12.py --results results --output OUTPUT.pdf --font KOREAN_FONT.ttf --bold-font KOREAN_BOLD_FONT.ttf --as-of "YYYY-MM-DD HH:MM:SS" --revision COMMIT_SHA
```

실제결과·중단기록으로15쪽을생성한다.옆manifest는사용자료·폰트·builder·
revision을고정한다.모든페이지렌더링과시각검수는별도로수행한다.
