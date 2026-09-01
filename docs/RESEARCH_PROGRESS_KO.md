# 연구 진행상황 — 2026-09-02

## 현재 활성 이론

단결정의 반복 단축 인장 하에서 normal spacing `a`와 하나의 결정학적
unwrapped slip coordinate `s`를 동시에 쓰는 통합 미시 에너지를 도입했다.

```text
U0(a,s) = sum_{k>=1} W(k*a,s)
```

여기서 `W(d,s)`는 row--row kernel이다. local fatigue counting에서는 각
normal layer를 기준 layer와 한 번씩 세므로 `k W(k*a,s)`가 아니며, `s`는
모든 layer에서 같은 collective coordinate이므로 `W(k*a,k*s)`도 아니다.

## 수학적 진전

- `H_q=sum_k sum_p[(p+delta)^2+k^2 eta^2]^(-q/2)`를 정의했다.
- Mellin--Poisson--Bessel 유도를 보존하여 zeta zero mode와
  Bessel--Lambert series로 이루어진 정확한 `H_q` 식을 얻었다.
- full absolute energy의 수렴조건을 `m>n>2`로 확정했다.
- normal energy와 slip energy를 서로 더하지 않고 동일한 `U0`의 정확한
  항등분해로 정의했다.
- 12--6의 `K_5/2`, `K_11/2`를 독립적으로 전개해 polylog closure를 얻었다.
- direct double sum과 analytic representation의 상대오차는 시험점에서
  `q=6`은 최대 약 `1.12e-11`, `q=12`는 약 `2e-16` 이하였다.
- polylog와 Bessel--Lambert 계산은 기계정밀도로 일치했다.
- 전체 회귀검사 `167 passed`, 새 multilayer/registry 집중검사 `25 passed`를
  확인했다.
- 영문 본문과 한국어 요약/기호표를 포함한 11쪽 논문 PDF를 로컬 Tectonic
  0.17.0으로 빌드했다. LaTeX 오류, 누락 글자, overflow, 미해결 참조가 없다.
  산출물은 `output/pdf/slip_lattice_energy_derivation.pdf`에 있다.
- `.github/workflows/**`는 변경하지 않았고 지속 자동화도 추가하지 않았다.

## 확률ㆍ소성ㆍ균열 연결

underlying law는 `(a,s)` Smoluchowski continuity equation이다. 공식적인
“4 governing equations”는 평균거리, 평균 intrinsic energy, 누적
hysteresis dissipation, normalization/survival이다.

소성변형은 `s=s0+z*b+s_tilde`의 well population 이동으로 계산하며, 단순
phase lag가 아니라 unloading/relaxation 후 `Delta<z> != 0`을 요구한다.

균열개시는 임의 `a_c`나 `E_hyst>E_c`로 정의하지 않는다. 순간 인장 force에서
`partial_a U0=Q_a`, `partial_a^2 U0<0`인 외측 장벽을 흡수경계로 두고 그
상대 probability outflux를 누적해 `P_crack=1-S`로 정의한다.

## 아직 미완료/미확정

- 현재 registry transport demo는 `a/b`를 고정한다. 완전한 2D `(a,s)`
  finite-volume solver와 moving absorbing boundary의 수치결합은 후속 작업이다.
- `A0`, mobility/memory, `h_slip`, active slip system, dislocation hardening은
  물리적으로 미확정이다.
- 현재 결과는 dimensionless mechanism verification이며 알루미늄 수명예측이
  아니다.
- EAM/DFT는 향후 정량 검증용일 뿐 현재 governing potential이 아니다.

진행상황은 모든 이론ㆍ코드ㆍ검증ㆍ논문 변경 commit 전에 README와 이 문서에
항상 함께 갱신한다.
