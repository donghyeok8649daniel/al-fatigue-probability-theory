# v9 실행과 검증

명령은 `silicon-wafer-research` 저장소 루트에서 실행한다. 아래 `python`은 해당
단계의 의존성을 가진 실행기를 뜻한다. 새 출력 폴더를 사용한다. 저장된 결과를
읽는 감사와 실제 원자모델 재실행은 다르다.

## 환경과 원본

- 수학/회귀: Python3.13.5, NumPy2.5.0, SciPy1.18.1. 전체 환경은 실제 기록 참고.
- MACE: 별도 venv, MACE0.3.16, torch2.12.1+cpu, e3nn0.4.4, ASE3.29.0,
  matscipy1.2.0. CPU/float64, torch2threads. 전역 패키지 변경 없이 기존 torch를
  공유한 격리 환경에 MACE 의존성을 설치했다.
- 정확한 패키지/배포 wheel hash와 모델은 `sources/mace_runtime_provenance.json`.
  이 기록은 실제 실행 환경이며 모든 플랫폼에서 설치된다는 보장은 아니다.
- MACE weight는 [공식 배포](https://github.com/ACEsuit/mace-foundations/releases/download/mace_mp_0b3/mace-mp-0b3-medium.model)에서
  받는다. SHA256 `2f2be696351ac9e94fbe01cdfb6f017679acdbd2db7645209ef55fec9826b012`.
  큰 weight 파일은 Git에 넣지 않았다. 아래 예시는 `.cache/mace-mp-0b3-medium.model`.
- Cambridge Si archive는 v5와 동일한 `Si_PRX_GAP.zip`이다. 원본 획득과 hash는
  [v5 자료](../silicon_atomistic_v5/COMPLETED_SUMMARY.md)를 따른다.
  아래 예시는 `.cache/si-atomistic-v5/Si_PRX_GAP.zip`.
- Durham B archive는 `sources/durham_boron_2025.zip`에 원본 그대로 포함했다.
  워크북을 읽는 단계에는 openpyxl이 필요하다.

## 수학 계산

각 스크립트의 `--help`로 정확한 옵션을 확인한다.

```text
python -m results.silicon_wafer_feasibility.audit_charge_dynamics_v9 --output .cache/v9-replay/dynamics
python -m results.silicon_wafer_feasibility.audit_charge_memory_v9 --output .cache/v9-replay/memory
python -m results.silicon_wafer_feasibility.audit_global_charge_v9 --output .cache/v9-replay/global_charge
python -m pytest solver_v1/test_silicon_charge_dynamics.py -q
```

기록된 별도 재실행은 수학3폴더와 B source audit의 JSON/CSV19개를 byte 대조했다.
그림의 binary identity나 무거운 MACE 전체 재실행을 뜻하지 않는다.

## 실제 pretrained potential 평가

```text
python -m results.silicon_wafer_feasibility.run_mace_dft_audit_v9 --model .cache/mace-mp-0b3-medium.model --archive .cache/si-atomistic-v5/Si_PRX_GAP.zip --output .cache/v9-new/mace_si_dft
python -m results.silicon_wafer_feasibility.validate_mace_dft_v9 --model .cache/mace-mp-0b3-medium.model --archive .cache/si-atomistic-v5/Si_PRX_GAP.zip --results .cache/v9-new/mace_si_dft --baseline results/silicon_atomistic_v5 --output .cache/v9-new/mace_si_dft_validation
```

첫 명령은 새2,476회 계산(2,475frames+energy anchor)이다. 두 번째 명령은 원시
배열 재집계, v5 source 배열 일치, 7개 독립 backend 대조와 회전/복제 대조다.
새 DFT나 재적합이 아니다. XC별로 다른 source를 유지한다.

B 초기 이완은 `run_mace_boron_v9 --help`를 따른다. 생성한 각 `*_final.extxyz`는
힘 허용오차를 통과한 stationary candidate이며 Hessian 전에는 최소점이라고
보장하지 않는다. 원본 extxyz 위치는8자리 소수 정밀도이며 Hessian은 이 저장
좌표를 다시 읽어 계산했다. 원본 reloaded force도 각각 기록한다.

```text
python -m results.silicon_wafer_feasibility.validate_mace_boron_v9 --model .cache/mace-mp-0b3-medium.model --states results/silicon_doping_v9/mace_boron --output .cache/v9-new/hessian --force-backend direct --cases b_silicon_64_interstitial 2boron_int_cluster 3boron_int_cluster_confi_1 3boron_int_cluster_confi_2 3boron_int_cluster_confi_3
python -m results.silicon_wafer_feasibility.validate_mace_curvature_modes_v9 --model .cache/mace-mp-0b3-medium.model --states results/silicon_doping_v9/mace_boron --hessians .cache/v9-new/hessian --output .cache/v9-new/hessian_refinement
```

`direct`는 같은 MACE의 힘만 미분하는 경로다. 표준 ASE 계산과의 일치를 먼저
검사하고 좌표당2회씩 순차 호출한다. 전체 Hessian을 벡터화한 첫 시도는 메모리
증가로 중단했으며 `failed_attempts/`에 보존했다. potential을 바꾸지 않았다.

음의 mode 추적은 `follow_boron_instability_v9 --help`를 따른다. 최초 line search
없는 LBFGS 시도는 중단했다. 완료된 재실행은 BFGSLineSearch로 ±방향을 모두
추적했다. 두 새 끝점의 Hessian은 별도로 계산했다. 최초 상태와 새 끝점을 덮어쓰지
않았다. 유한 셀 대조는 `run_mace_dilute_cells_v9 --repeats 3 4`로 실제217/218개,
513/514개 원자를 생성한다. 원본64개 host Si의 격자 일치를 실행 전에 검사한다.

가장 높은 진동수의 셀 크기 대조는 `validate_largest_mode_v9 --help`를 따른다.
64/216/512host의1B interstitial 상태를 순서대로 넣고 첫 상태의 전체 Hessian을
`--reference-hessian`으로 제공한다. 질량 가중 Hessian-vector product를 중앙 힘
차분으로 계산해 ARPACK의 largest eigenpair를 구한다. 동일 원소 질량, 고정 셀,
독립 차분 간격의 eigen-residual을 저장한다. 가장 낮은 mode나 전체 안정성을
검사하는 계산은 아니다. 큰 행렬의 전체 메모리 적재는 필요하지 않다.

## 벌크 탄성 보충

추가 벌크 탄성은 `audit_mace_elastic_v9 --model MODEL --bulk BULK_JSON --output NEW_FOLDER`로
실제 다시 계산한다. BULK_JSON은 `mace_boron/bulk_reference.json`이다.
`validate_mace_elastic_v9 --model MODEL --reference ELASTIC_SUMMARY --output NEW_FOLDER`는
서로 다른 변형 모드를 이용한 독립8상태 대조다. 최초 `mace_bulk_elastic`의37표준평가
기록은 캐시를 포함한 상태 조회 수였으며 `mace_bulk_elastic_counted`의 실제 호출 계측으로
구분했다. 두 실행의 CSV3개는 byte 일치한다. 두 폴더의 합을 하나의 실행 횟수로 부르지 않는다.

## 저장 결과만 감사하기

```text
python -m results.silicon_wafer_feasibility.verify_artifacts_v9 --root results/silicon_doping_v9
```

이 명령은 hash, source/case 수, 원시 Hessian 고유값, 판정, 완료 기록의 일관성을
읽는다. 원자 계산이나 pytest를 다시 실행하지 않는다. 전체 solver 회귀의 실제
명령·종료코드·reported time은 `verification/full_solver.json/.log`에 있다.
`--write`는 새 audit/manifest를 만드는 옵션이므로 배포된 결과를 검증할 때는 빼고 쓴다.

## 해석 제한

MACE의 원소 구분은 사용했지만 전자 수를 독립적으로 조절하지 않았다. B acceptor⁻
문헌 결과를 neutral MACE와 같은 전하라고 비교하지 않는다. 고정 셀 Gamma 안정성과
세 셀 크기의 에너지 차이는 전체 q/strain/finite-T 수렴이나 실제 수명 검증이 아니다.
검증된 Si free-energy branch, 균열 장벽, kinetic rate와 physical clock은 여전히 없다.
