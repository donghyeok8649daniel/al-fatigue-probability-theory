# v15 재현

Si 브랜치 루트에서 Python, NumPy, SciPy, pytest를 사용한다. 실제 원자 파일럿에는
기존 v14와 같은 ASE/PyTorch/MACE 환경과 SHA가 고정된 MP-0b3 모델이 추가로 필요하다.
개인 경로와 인증 정보는 저장하지 않는다. OMP/MKL/OpenBLAS 스레드는 1로 제한했다.

## 검증 (원자 모델 실행 없음)

```powershell
python -m pytest solver_v1/test_silicon_device_ensemble.py solver_v1/test_silicon_thermal_research.py -q -p no:cacheprovider
python -m results.silicon_wafer_feasibility.replay_device_ensemble_v15
python -m results.silicon_wafer_feasibility.seal_device_ensemble_v15 --revision HEAD
```

재검산기에 `--cad <original-file.cad.json>`을 주면 외부 CAD 파일 SHA도 확인한다.
그 경우 총 33검사, 생략 시 외부 파일 검사를 제외한 32검사다.
재검산기는 이번 v15의 preallocation 중단 기록을 대상으로 한다.
추후 별도 폴더에서 실제 chain을 실행한 결과는 별도 chain 재검산이 필요하다.

## 기존 Hessian 재분석

```powershell
python -m results.silicon_wafer_feasibility.audit_device_ensemble_v15 --cad <original-file.cad.json> --output <new-empty-audit-directory>
```

입력 파일 SHA는 `device_audit/input_manifest.json`에 있다. CAD 원본은 사용자 파일을
그대로 두고 필요한 사실과 전체 파일 SHA만 저장했다. 최신 CAD 소스는
`CAD_REPOSITORY_REFERENCE.json`의 commit 및 파일 SHA로 확인할 수 있다.

## 실제 비선형 원자 파일럿 (이번에는 모델 로드 전 중단)

```powershell
python -m results.silicon_wafer_feasibility.sample_nonlinear_ensemble_v15 --model <mace-mp-0b3-medium.model> --output <new-empty-pilot-directory> --max-seconds 1200
```

- 모델 SHA: `2f2be696351ac9e94fbe01cdfb6f017679acdbd2db7645209ef55fec9826b012`.
- 300 K, β=0.35, 2 chain, 각각 warmup 16 / 저장 48, 총 제안 128.
- 두 원자료 재평가 포함 최대 130회 모델 호출. 실행 전 메모리 6 GB 기준 유지.
- 다른 프로그램이나 기존 연구 프로세스를 종료하지 않는다.
- 기존 output이 있으면 거부하며, 중간 제안 원자료는 단계마다 저장한다.
- MC 단계는 물리 시간/실제 균열 횟수가 아니다. chain 완료도 혼합 인증이 아니다.
- 실행 완료 시 독립 난수/수용 이력 검산 및 표본 혼합·영역 검사를 이어야 한다.

`package_manifest.json`은 v15 결과·소스와 사용한 기존 원자료를 exact byte로 묶는다.
외부 CAD 전체 파일, 모델 가중치, 변경 가능한 원본 handoff는 패키지에 복사하지 않는다.
Windows에서 생성된 JSON의 CRLF도 해시 대상 원자료 바이트로 보존한다. 첫 Git
공백 검사는 CR을 줄 끝 공백으로 판정했다. JSON에만 cr-at-eol 규약을 명시해
해결했으며, 원자료 바이트나 소스 코드의 공백 검사를 완화/변경하지 않았다.
