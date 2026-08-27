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

The default root-level normal simulations and tests must not import or depend on this library. Library material may be reused later only through an explicit comparison or coupling study.

## Variable dictionary rule

Any new active normal theory or simulation variable must be added to `docs/VARIABLE_DEFINITIONS_NORMAL_LJ.md` in the same commit.

Any new firmware field or flag must be added to `firmware/VARIABLE_DEFINITIONS.md` in the same commit.

Variables used only by the shear library belong to `libraries/shear/docs/VARIABLE_DEFINITIONS.md`.

## Modeling labels

Important claims should be identified as one of:

- **EXACT / IDENTITY**
- **DEFINITION**
- **ASSUMPTION**
- **CONTROLLED APPROXIMATION**
- **EMPIRICAL INPUT**

Never silently promote an approximation into an exact statement.

## Numerical-result rule

Every new numerical claim should state:

- the model and parameters;
- the numerical method;
- the relevant convergence or conservation check;
- whether the result is dimensional or nondimensional;
- whether it is a physical prediction, a null test, or a proof of principle.

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

기본 root-level normal simulation과 test는 이 library를 import하거나 의존하면 안 된다. 향후 명시적인 비교 또는 coupling 연구를 할 때만 필요한 부분을 다시 사용한다.

## 변수사전 규칙

새로운 활성 normal theory/simulation 변수는 같은 commit에서 `docs/VARIABLE_DEFINITIONS_NORMAL_LJ.md`에 추가한다.

새로운 firmware field 또는 flag는 같은 commit에서 `firmware/VARIABLE_DEFINITIONS.md`에 추가한다.

shear library에서만 쓰는 변수는 `libraries/shear/docs/VARIABLE_DEFINITIONS.md`에 둔다.

## 모델링 분류 라벨

중요한 주장은 다음 중 하나로 분류한다.

- **EXACT / IDENTITY**
- **DEFINITION**
- **ASSUMPTION**
- **CONTROLLED APPROXIMATION**
- **EMPIRICAL INPUT**

근사를 설명 없이 정확식으로 승격하지 않는다.

## 수치결과 규칙

새로운 수치 주장은 반드시 다음을 명시한다.

- 사용한 model과 parameter;
- numerical method;
- 관련 convergence 또는 conservation check;
- dimensional result인지 nondimensional result인지;
- physical prediction, null test, proof of principle 중 무엇인지.
