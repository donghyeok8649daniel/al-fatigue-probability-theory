# Repository Style Guide

## GitHub Markdown math

All mathematical expressions in repository Markdown files must use GitHub-compatible math delimiters.

### Inline math

Use single-dollar delimiters:

```markdown
$P(a,t)$
```

### Display math

Use double-dollar delimiters:

```markdown
$$
\partial_t P + \partial_a(Pv)=0
$$
```

Do not use `\[` and `\]` in Markdown files.

## Mandatory bilingual Markdown rule

**Every `.md` file in this repository, including files inside research libraries, must contain a complete Korean translation.**

Required order:

1. English technical version;
2. horizontal rule `---`;
3. heading `# 한국어 번역` or equivalent;
4. Korean translation of the complete technical content.

A short Korean summary is not a substitute for the full translation.

## Active normal-mainline rule

The repository root is reserved for the **normal-deformation / normal-opening** research path.

New active theory, simulation, test, result, and documentation files should support

$$
\sigma_n(t)
\rightarrow
\{a_i(t)\}
\rightarrow
P(a,t)
\rightarrow
\text{normal memory / secular evolution}
\rightarrow
\text{normal-opening instability}.
$$

A non-normal mechanism should not be promoted into the active mainline unless the normal problem itself shows that it is mathematically necessary.

## Auxiliary-library rule

Historical or side-path work may be preserved under `libraries/` without being part of the active mainline.

Current preserved auxiliary library:

- `libraries/shear/` — Rubin-chain, non-affine slip, gamma-surface, shear-oriented simulations, tests, data, and reports.

The default root-level normal simulations and tests must not import or depend on this library.

## Variable dictionary rule

Every new active symbol must be documented in the same change that introduces it.

Use:

- `docs/VARIABLE_DEFINITIONS_NORMAL_LJ.md` for spacing, probability, 1D normal dynamics, time-scale, and first-passage variables;
- `docs/VARIABLE_DEFINITIONS_FCC_NORMAL_LJ.md` for FCC geometry, deformation-gradient, lattice-sum, directional-elasticity, and FCC calibration variables;
- `firmware/VARIABLE_DEFINITIONS.md` for firmware fields and fault flags;
- `libraries/shear/docs/VARIABLE_DEFINITIONS.md` for variables used only by the shear library.

## Modeling labels

Important claims should be identified as one of:

- **EXACT / IDENTITY**
- **DEFINITION**
- **ASSUMPTION**
- **CONTROLLED APPROXIMATION**
- **EMPIRICAL INPUT**

Never silently promote an approximation into an exact statement.

## Fixed-potential rule

Active Lennard-Jones parameters must not be changed as a function of cycle count merely to imitate fatigue damage.

If one fixed potential cannot reproduce two physical quantities simultaneously, record the incompatibility as a model result. Do not hide it by switching parameter sets inside one derivation.

A new energy contribution may be introduced only when its physical origin and mathematical necessity are stated explicitly.

## Numerical-result rule

Every new numerical claim should state:

- the model and parameters;
- the numerical method;
- the relevant convergence or conservation check;
- whether the result is dimensional or nondimensional;
- whether it is a physical prediction, a null test, a calibration study, or a proof of principle.

A discrepancy much larger than the numerical convergence error must be treated as a model result, not tuned away automatically.

---

# 한국어 번역 — 저장소 작성 규칙

## GitHub Markdown 수식

저장소의 모든 Markdown 파일은 GitHub에서 정상 렌더링되는 수식 delimiter를 사용한다.

### 인라인 수식

단일 dollar 기호를 사용한다.

```markdown
$P(a,t)$
```

### 독립 수식 블록

이중 dollar 기호를 사용한다.

```markdown
$$
\partial_t P + \partial_a(Pv)=0
$$
```

Markdown에서는 `\[`와 `\]`를 사용하지 않는다.

## 모든 Markdown의 한국어 번역 의무

**research library 내부 파일까지 포함하여 저장소의 모든 `.md`에는 전체 한국어 번역이 있어야 한다.**

작성 순서는

1. 영문 기술 원문;
2. `---`;
3. `# 한국어 번역` 또는 동등한 제목;
4. 영문 기술내용 전체의 한국어 번역

으로 통일한다. 짧은 한국어 요약은 전체 번역을 대신할 수 없다.

## 활성 normal-mainline 규칙

repository root는 **수직변형 / normal-opening** 연구경로 전용이다.

새로운 활성 theory, simulation, test, result, documentation은

$$
\sigma_n(t)
\rightarrow
\{a_i(t)\}
\rightarrow
P(a,t)
\rightarrow
\text{수직 memory / secular evolution}
\rightarrow
\text{normal-opening instability}
$$

를 지원해야 한다.

normal problem 자체에서 수학적으로 필요하다는 것이 드러나지 않는 한 non-normal mechanism을 활성 mainline으로 승격하지 않는다.

## 보조 library 규칙

과거 연구나 side-path 연구는 활성 mainline과 분리하여 `libraries/` 아래 보존할 수 있다.

현재 보존된 보조 library는

- `libraries/shear/` — Rubin-chain, non-affine slip, gamma-surface, 전단 지향 simulation/test/data/report

이다.

기본 root-level normal simulation과 test는 이 library를 import하거나 의존하면 안 된다.

## 변수사전 규칙

새로운 active symbol은 도입하는 변경과 같은 변경에서 반드시 정의한다.

다음 파일을 사용한다.

- `docs/VARIABLE_DEFINITIONS_NORMAL_LJ.md` — spacing, probability, 1D normal dynamics, time scale, first-passage 변수;
- `docs/VARIABLE_DEFINITIONS_FCC_NORMAL_LJ.md` — FCC geometry, deformation gradient, lattice sum, directional elasticity, FCC calibration 변수;
- `firmware/VARIABLE_DEFINITIONS.md` — firmware field와 fault flag;
- `libraries/shear/docs/VARIABLE_DEFINITIONS.md` — shear library 전용 변수.

## 모델링 분류 라벨

중요한 주장은 다음 중 하나로 분류한다.

- **EXACT / IDENTITY**
- **DEFINITION**
- **ASSUMPTION**
- **CONTROLLED APPROXIMATION**
- **EMPIRICAL INPUT**

근사를 설명 없이 정확식으로 승격하지 않는다.

## Fixed-potential 규칙

fatigue damage를 흉내내기 위해 cycle 수에 따라 active Lennard-Jones parameter를 변경하지 않는다.

하나의 fixed potential이 두 물리량을 동시에 재현하지 못하면 그 incompatibility 자체를 model result로 기록한다. 하나의 derivation 안에서 parameter set을 바꿔 문제를 숨기지 않는다.

새 energy contribution은 물리적 기원과 수학적 필요성이 명시된 경우에만 도입한다.

## 수치결과 규칙

새로운 수치 주장은 반드시 다음을 명시한다.

- 사용한 model과 parameter;
- numerical method;
- 관련 convergence 또는 conservation check;
- dimensional result인지 nondimensional result인지;
- physical prediction, null test, calibration study, proof of principle 중 무엇인지.

numerical convergence error보다 훨씬 큰 discrepancy는 자동 tuning 대상이 아니라 model result로 취급한다.
