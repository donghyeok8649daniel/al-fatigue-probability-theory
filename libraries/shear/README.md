# Shear / Auxiliary Mechanics Research Library

This directory preserves the earlier Rubin-chain, non-affine slip, gamma-surface, and shear-oriented proof-of-principle work.

## Status

This library is **preserved research history and auxiliary methodology**, not the active mainline of the repository.

The active project now focuses on

$$
\sigma_n(t)\rightarrow a_i(t)\rightarrow P(a,t)\rightarrow\text{normal-opening instability}.
$$

Nothing in this library should be imported by the default root-level normal simulations or tests.

## Why it is preserved

The work here still contains useful results:

- conservative hidden-mode dynamics can generate reduced hysteresis without fitted viscous damping;
- nonlinear non-affine coordinates can produce cycle-to-cycle state changes in a proof-of-principle model;
- gamma-surface and pure-Al cyclic-deformation literature provide constraints on shear mechanisms;
- the numerical energy-balance and falsification tests remain useful methodological references.

These results are not deleted because they may later be useful for comparison, coupled-mode studies, or failure-mode analysis.

## Structure

- `docs/` — Rubin, slip, gamma-surface, shear constraints, and the broad historical variable dictionary
- `theory/` — Rubin-chain and Hamiltonian slip-bath model code
- `simulations/` — shear/auxiliary simulation runners
- `tests/` — shear/auxiliary unit tests
- `results/` — historical data, figures, and result reports

## Import policy

The files are preserved close to their original form. If they are executed from this subdirectory in the future, imports may need to be adjusted to the library package layout. Do not modify the active normal mainline merely to keep old shear scripts executable.

---

# 한국어 번역 — 전단 / 보조 역학 연구 라이브러리

이 디렉터리는 기존 Rubin-chain, non-affine slip, gamma-surface, 전단 지향 proof-of-principle 연구를 보존한다.

## 상태

이 library는 **보존된 연구이력 및 보조 방법론**이며 현재 repository의 활성 mainline이 아니다.

활성 프로젝트는 현재

$$
\sigma_n(t)\rightarrow a_i(t)\rightarrow P(a,t)\rightarrow\text{normal-opening instability}
$$

에 집중한다.

이 library의 어떤 파일도 기본 root-level normal simulation이나 test에서 import하면 안 된다.

## 보존하는 이유

여기에는 여전히 유용한 결과가 있다.

- fitted viscous damping 없이 conservative hidden-mode dynamics만으로 reduced hysteresis가 가능하다는 결과;
- nonlinear non-affine coordinate가 proof-of-principle에서 cycle-to-cycle state change를 만들 수 있다는 결과;
- gamma-surface 및 pure-Al cyclic-deformation 문헌이 주는 shear-mechanism 제약;
- energy-balance와 falsification test 방법론.

향후 비교연구, coupled-mode 연구, failure-mode 분석에 사용할 수 있으므로 삭제하지 않는다.

## 구조

- `docs/` — Rubin, slip, gamma-surface, shear constraint, 과거 broad 변수사전
- `theory/` — Rubin-chain 및 Hamiltonian slip-bath 코드
- `simulations/` — shear/auxiliary simulation runner
- `tests/` — shear/auxiliary unit test
- `results/` — 과거 data, figure, result report

## Import 정책

파일은 가능한 한 원형에 가깝게 보존한다. 향후 이 subdirectory에서 다시 실행할 경우 library package layout에 맞게 import를 조정해야 할 수 있다. 과거 shear script의 실행성을 유지하기 위해 활성 normal mainline을 변경하지 않는다.
