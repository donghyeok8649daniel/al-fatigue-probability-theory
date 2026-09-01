# Al Fatigue Probability Theory

This `main` branch is the **stable project index/release baseline**. Active development is separated by responsibility so theory, numerical implementation, manuscript work, and physical tester development do not overwrite each other.

## Active branches

| Branch | Responsibility | Source of truth |
|---|---|---|
| `theory-core` | Physical/mathematical model: generalized LJ lattice energy, probability evolution, plastic slip, four governing equations, crack-initiation definition | Governing equations and derivations |
| `numerical-fem` | Numerical realization: lattice sums, Smoluchowski solvers, FEM coupling, UI, verification, result generation | Executable numerical implementation |
| `paper-manuscript` | LaTeX manuscript, paper figures/tables, notation/index, appendices, PDF-oriented writing | Published presentation of validated theory/results |
| `fatigue-tester` | Physical fatigue tester: firmware, hardware, BOM, DCPD, actuator/control, telemetry bridge | Experimental machine implementation |
| `integration` | Cross-branch merge and interface validation point | Combined system integration |
| `fem-probability-coupling` | Legacy/migration branch retained for history while the split is completed | Do not use as the long-term ownership boundary |

## Dependency and merge direction

Do not routinely merge feature branches directly into one another. Reviewed work converges through `integration`:

```text
theory-core ───────┐
numerical-fem ─────┼──→ integration ──→ main (stable checkpoints only)
fatigue-tester ────┤
paper-manuscript ──┘
```

Technical dependency still flows primarily as:

```text
theory-core
    ↓ validated equations
numerical-fem
    ↓ verified numerical results
paper-manuscript
```

Experimental data follows:

```text
fatigue-tester
    ↓ timestamped telemetry
numerical-fem / UI
```

`main` should receive only reviewed, cross-module stable checkpoints from `integration`. Branch-specific development should not be performed directly on `main`.

## Main-branch policy

`main` intentionally does **not** own active theory code, FEM/UI implementation, manuscript sources, or fatigue-tester firmware/hardware. Those artifacts are preserved in the dedicated branches above.

Repository-level infrastructure such as `.gitignore` and existing GitHub configuration remains here. CI/workflow changes are handled separately and are not implicitly changed by branch-content reorganization.

## Shared interfaces

Cross-branch integration should use explicit versioned interfaces rather than reaching into another branch's internal implementation. The main shared contract to freeze next is the tester/UI telemetry and command schema, including timestamp, cycle/phase, reference/measured force, displacement/strain, temperature, DCPD, actuator command, and fault flags.

---

# 한국어

`main`은 앞으로 **안정된 프로젝트 진입점/릴리스 기준선**으로만 사용한다.

- `theory-core`: 이론 및 지배방정식
- `numerical-fem`: 수치해석, FEM, UI, 검증코드
- `paper-manuscript`: 논문 LaTeX/PDF/그림/표
- `fatigue-tester`: 실제 시험기 펌웨어/하드웨어/BOM/통신
- `integration`: 브랜치 간 통합 및 인터페이스 검증

개발 브랜치끼리 직접 merge를 남발하지 않고 `integration`에서 합친 뒤, 검증된 안정 checkpoint만 `main`에 반영한다. `main`에서 직접 기능 개발하지 않는다.
