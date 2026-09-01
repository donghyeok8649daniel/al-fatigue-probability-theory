# Imported research sources

## Corrected ideal-registry derivation

The project owner supplied the current corrected source on 2026-09-01. It has
replaced the earlier repository copy at:

`slip_lattice_energy_mn_K_derivation_KR_v3_23pages.pdf`

- Current PDF SHA-256:
  `42C3D5086CA203C76F3DC8213A1718B5121AA1273067738C5B478BCBF12D999D`
- Current PDF size: 125,691 bytes
- Previous PDF SHA-256, recoverable from Git history:
  `0E293BC2D8C33788CF89B9667630F22D094E267360F35C52B497385FC2DB8208`
- Status: corrected research source and provenance for the optional active
  ideal-registry branch

The accompanying `symbol_index_en.tex` has SHA-256
`D72A3C3CA43339467489D47FB57E475F7DA5014BE236EF2E370D86B527D11A93`.
A build-local copy, with the repository coefficient and area conventions
corrected, is stored at `libraries/shear/docs/symbol_index_en.tex`.

The source is not accepted blindly. The exact shifted Epstein--Hurwitz /
Poisson--Bessel identity, its derivative, and the unwrapped registry notation
are active after independent numerical checks. The repository continues to
correct coefficient/well-depth notation, keeps `A0`, `Ac`, and registry area
distinct, rejects unsupported addition of the collinear and two-row energies,
and makes no quantitative aluminum-plasticity claim.

See `docs/ACTIVE_IDEAL_REGISTRY_PLASTICITY.md` and
`libraries/shear/docs/SLIP_LATTICE_ENERGY_REVIEW.md`.

## 한국어 상태

2026-09-01에 프로젝트 소유자가 제공한 오류수정 PDF가 기존 저장소 PDF를
대체했다. 이전판은 Git 이력에서 복구할 수 있다. 정확한 Bessel 격자합과
unwrapped registry 좌표는 독립 검산 뒤 선택적 활성 소성 branch에 반영했다.
다만 `A0`, correlation area, registry interface area는 동일시하지 않으며,
서로 다른 기하의 normal-chain energy와 two-row energy도 단순히 더하지 않는다.
현재 결과는 정량적 알루미늄 소성 보정값이 아니라 이상적 단일 registry
mechanism이다.
