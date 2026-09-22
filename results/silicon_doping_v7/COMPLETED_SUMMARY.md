# Si 도핑 연구 v7 — 탄성 연결과 전하 조건의 구현·검증

## 결론

도핑 원자 농도, 이온화된 도펀트, 전자·정공 농도를 분리하고,
**실측 도핑별 탄성을 기존 시편 균열 계산에 연결했다.** 같은 화학 배치에서
빠른 전하 평형을 소거하는 자유에너지의 gradient/Hessian도 구현·수학 검증했다.

실제 도핑된 균열 에너지·전이 장벽·물리 이동도 보정은 아직 미완료다.
새 DFT/MD/원자 에너지 평가는 0회다. 순수 SW/Tersoff를 도핑 재료로
바꾸거나, 생산 Si/Hz를 활성화하지 않았다. Al/UI 코드는 변경하지 않았다.

## 1. 새로 실행한 계산

| 항목 | 실행·결과 |
|---|---|
| B/P/As 문헌 탄성 | 7시편 × -40/25/85 C = 21개 균열장 |
| 독립 J contour 대조 | 각2개 반경/분할, 합계42개; 최대 상대차 1.77636e-15 |
| 탄성계수 오차 민감도 | 실온 7개×8 corner=56개; 통계적 joint CI 아님 |
| 화학 농도의 원자 수 환산 | 1904/7600 sites×29농도 = 58개, 반올림·MC 없음 |
| 전자 평형 축약 | 합성된 2차 energy branch에서 force/Hessian/전자수 미분 검사 |
| 문헌 데이터 보존 | 탄성7개, 전하 DFT12행(무도핑 중복 포함), 웨이퍼6개 조건 |

21개는 continuum 계산이다. 도핑된 원자 시편21개를 새로 이완했다는 뜻이 아니다.
J 대조의 작은 오차는 구현 일치이며 재료 오차가 1e-15라는 뜻이 아니다.

## 2. 정량적으로 드러난 내용

### 탄성 변화만 넣었을 때

(111)[11-2], 전면[1-10], plane strain, 25 C에서:

| 문헌 시편 | carrier 평균 /cm^3 | C11/C12/C44 (GPa) | H (GPa^-1) | 같은 Wsep 가정의 K 기준 비 |
|---|---:|---|---:|---:|
| B0.6 | 6.00e18 | 165.5 / 63.7 / 79.6 | 0.005688632 | 1 |
| B3 | 2.96e19 | 164.2 / 63.5 / 78.5 | 0.005751239 | 0.994542 |
| P7.5 | 7.47e19 | 161.4 / 66.1 / 78.5 | 0.005808723 | 0.989609 |

전체7시편의 비는 0.989609~1이다. **분리 에너지가 같다는 가정 아래 약1.04%의
변화**이며 도핑의 전체 효과나 실제 강도 변화가 아니다. 기준 B0.6도 도핑되어 있다.
v6의 순수 SW에서 나온 약110 MPa에 이 비를 곱해 도핑 강도라고 표시하지 않는다.
응력/균열 형상→K→J에는 실측 탄성을 연결했지만, 실측 C와 SW의 원자 에너지를
혼합한 atomistic 계산은 하지 않았다.

### 작은 원자 모델의 농도 해상도

실험 기준 격자 5.431 Angstrom을 사용하면 site density=4.99402e22 cm^-3다.
7600 sites에서 도펀트 한 개는 6.57108e18 cm^-3다.

| 화학 농도 /cm^3 | 평균 도펀트 수 | 0개일 확률(독립 균일 치환 가설) |
|---:|---:|---:|
| 1e15 | 0.000152182 | 0.999847830 |
| 1e19 | 1.521819 | 0.218281138 |

저농도를 작은 셀에 원자1개로 표현하면 실제 농도가 크게 바뀐다.
도펀트가 없는 국소 영역도 외부의 전자 reservoir 영향을 받을 수 있으므로
이 표를 도핑 효과가 없을 확률 또는 균열 확률로 읽으면 안 된다.

### 큰 변화가 나온 문헌은 다른 농도·조건이다

Noda 2023 부록 S5의 5e21 cm^-3에서 균일[111] 이상 인장강도는
무전하20.98 GPa 대비 전자11.19 GPa(-46.66%), 정공27.23 GPa(+29.79%)다.
기존 문헌 LDA 수치다. 실제 도펀트 원자가 없는 전하 주입 모델이며
이를 1e19대 화학 도핑으로 외삽하지 않았다.
같은 부록 S4의 정공 5e21에서 E111은 179.4→113.4 GPa로 낮아진다.
영률 감소와 이상 인장강도 증가가 동시에 나타나는 이 예는 탄성 변화만으로
파괴 장벽 변화의 방향까지 정할 수 없음을 보여 준다.

Chen 2025 Tables I/III의 실제 웨이퍼6조건, 각50개에서는 특성강도
2.07~2.12 GPa, 비교 p=.365~.929다. 논문에서 유의차를 검출하지 못했다는
결과를 보존한다. 같은 강도라는 증명, 임의 허용오차 내 동등성 검정,
원시300개 데이터 재분석을 했다는 뜻은 아니다. 이 Weibull 계수는 PDE에 넣지 않았다.

## 3. 공통 이론에 연결한 내용

화학적 도펀트 배치 alpha가 고정되어 있을 때 전하 자유도 N이 충분히 빠르게
평형에 도달하면 다음 grand potential을 사용한다:

    Omega_alpha(q,mu,T) = -kBT log sum_N exp[-(F_alpha(q,N,T)-mu N)/kBT]
    grad Omega = <grad F>
    Hess Omega = <Hess F> - Cov(grad F)/kBT

기울기만 평균하고 마지막 covariance 항을 빼먹으면 안정성과 장벽 곡률이 틀린다.
고정 전자 수/고정 밀도/고정 화학 퍼텐셜의 경계조건도 구분했다.

관측 시간 동안 도펀트 원자가 이동하지 않으면 배치별로 같은 확률 방정식을
풀고 제조 통계로 결과를 평균해야 한다. 도펀트 배치를 전자처럼 평형화하는
단일 Boltzmann 평균은 다른 물리를 가정한다. 실제 Si의 fast-charge 가정과
mobility는 아직 검증 전이다. 상세 유도는 [v7 이론](../../solver_v1/SILICON_DOPING_V7.md).

## 4. 검증 상태

- 첫 집중 검사: 55 PASS + 6 subtests, 1.80 s.
- 정수 전하 sector 검증을 보완한 최종 집중 검사: 57 PASS + 6 subtests, 1.22 s.
  전체 solver의 test collection 뒤 추가한 2개 검사는 이 최종 검사에 포함된다.
  분수-charge DFT 샘플을 독립 전하 상태처럼 합산하는 입력을 거부한다.
- 전체 solver: 893 PASS + 6 subtests, 1335.06 s. App: 189 PASS + 3 subtests,
  246.36 s. Desktop smoke: exit0, 1.879 s.
  `solver_tests.json`, `app_tests.json`, `desktop_smoke.json`에 실제 결과/시간/nodeid 저장.
  전체 solver collection 이후의 마지막 sector guard는 위 최종57개 검사로 확인했다.
- 영률은 독립 Kelvin tensor 역행렬과 비교했다. 전자 축약은 자유에너지의
  유한차분 gradient, full Hessian, dOmega/dmu=-N 및 공통 에너지 shift로 검사했다.
- 농도 보간은 원소별 측정 평균 농도 내부만 허용한다. 온도 외삽, Sb 대체,
  미지 이온화율의 자동1처리, site 초과 농도와 음수 입력은 거부한다.
- 그림 `doping_audit.png`는 직접 시각 검토했다. 원자료·결과·소스 해시는
  `artifact_manifest.json`에 남긴다.
- 별도 새 출력 폴더의 재실행에서 audit JSON과 두 CSV가 byte 단위로 일치했다.

## 5. 재현

저장소 root에서 NumPy/SciPy/Matplotlib 및 pytest 환경으로:

```text
python -m results.silicon_wafer_feasibility.run_doping_audit --output .cache/si-doping-reproduce
python -m pytest -q -p no:cacheprovider solver_v1/test_silicon_doping_research.py
```

기존 결과 파일이 있으면 덮어쓰지 않는다. 외부 LAMMPS/ASE/DFT 엔진은 필요 없다.
`sources/provenance.json`에 논문 제목·저자·연도·DOI/ID·방법·조건·표 위치를 기록했다.
Jaakkola 표 I/III는 PDF를 렌더해 직접 읽었다. 원문 text 추출의 숫자 깨짐을
전사에 사용하지 않았다. Noda 부록은 web PDF text로 읽었고 직접 다운로드는
HTML challenge여서 PDF로 취급하지 않았다. 정량 cleavage의 Xu2026 원자료는 미확보다.

## 6. Git과 다음 연구

작업 위치 `aft-silicon-wafer`, 브랜치 `silicon-wafer-research`.
시작 fresh fetch의 local/origin은 `7da9864b7893807b9fc41c5154c51c111b495ec9`로 일치했다.
이번 변경은 새 연구 module/test/runner/source/results 및 handoff/브랜치 안내다.
AGENTS.md와 기존 Al/UI/순수Si 에너지·원시 결과는 수정하지 않는다.
최종 commit/push 상태는 Git log와 최종 응답을 확인한다.

남은 작업은 도핑된 bulk/계면/균열의 동일 전자 조건 에너지·힘 자료,
치환 위치·표면·셀 크기 수렴, 실제 조건부 장벽과 이동도다.
v5의 순수 Si 재료 오차, v4 부분 실행, Si/물리Hz gate는 해결된 것으로 재해석하지 않는다.
