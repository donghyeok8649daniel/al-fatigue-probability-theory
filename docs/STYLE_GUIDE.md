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

**Every `.md` file in this repository, including auxiliary libraries, must contain a complete Korean translation.**

Required order:

1. English technical version;
2. horizontal rule `---`;
3. heading `# 한국어 번역` or equivalent;
4. Korean translation of the complete technical content.

A short summary is not a substitute for the full translation.

## Active-mainline rule

The repository root is reserved for the active **normal-deformation / normal-opening** mainline. Earlier shear, slip, gamma-surface, and Rubin-chain work belongs under `libraries/shear/` unless it is explicitly promoted back into the active theory.

Default root-level simulations and tests must not depend on the shear library.

## Variable dictionary rule

New active normal variables must be added to `docs/VARIABLE_DEFINITIONS_NORMAL_LJ.md` in the same commit. Shear-library variables belong to `libraries/shear/docs/VARIABLE_DEFINITIONS.md`.

## Modeling labels

Important claims should be identified as one of:

- **EXACT / IDENTITY**
- **DEFINITION**
- **ASSUMPTION**
- **CONTROLLED APPROXIMATION**
- **EMPIRICAL INPUT**

Never silently promote an approximation into an exact statement.

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

**보조 library를 포함하여 저장소의 모든 `.md` 파일에는 전체 한국어 번역이 있어야 한다.**

순서는

1. 영문 기술 원문;
2. `---`;
3. `# 한국어 번역` 또는 동등한 제목;
4. 전체 한국어 번역

으로 한다. 짧은 요약으로 전체 번역을 대신할 수 없다.

## 활성 mainline 규칙

repository root는 활성 **normal-deformation / normal-opening** mainline 전용이다. 과거 shear, slip, gamma-surface, Rubin-chain 연구는 다시 활성 이론으로 승격시키는 경우가 아니면 `libraries/shear/` 아래에 둔다.

기본 root-level simulation과 test는 shear library에 의존하면 안 된다.

## 변수사전 규칙

새로운 활성 normal 변수는 같은 commit에서 `docs/VARIABLE_DEFINITIONS_NORMAL_LJ.md`에 추가한다. shear-library 변수는 `libraries/shear/docs/VARIABLE_DEFINITIONS.md`에 둔다.

## 모델링 분류 라벨

중요한 주장은 다음 중 하나로 분류한다.

- **EXACT / IDENTITY**
- **DEFINITION**
- **ASSUMPTION**
- **CONTROLLED APPROXIMATION**
- **EMPIRICAL INPUT**

근사를 설명 없이 정확식으로 승격시키지 않는다.
