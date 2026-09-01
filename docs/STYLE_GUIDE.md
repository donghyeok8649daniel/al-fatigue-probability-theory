# Repository Style Guide

## GitHub Markdown math

Use GitHub-compatible math delimiters in every Markdown file.

Inline math:

```markdown
$P(a,t)$
```

Display math:

```markdown
$$
\partial_tP+\partial_a(Pv)=0
$$
```

Do not use `\[` and `\]` in repository Markdown.

## Mandatory bilingual Markdown rule

**Every `.md` file must contain a complete Korean translation.**

Required order:

1. complete English technical source;
2. `---`;
3. Korean translation heading;
4. complete Korean translation.

A short Korean summary is not sufficient.

## Active-mainline dimensionality rule

The active research root is **one-dimensional and normal-only**.

Active theory must be expressible in terms of 1D normal coordinates such as

$$
a_i(t),
\qquad
P(a,t),
\qquad
\mu(t),
\qquad
\mathcal E(t).
$$

Three-dimensional FCC work is preserved under `libraries/fcc_normal/`. Shear/Rubin/slip code is intentionally absent; `libraries/shear/README.md` is an inactive tombstone only. Archive code must not be imported by default active simulations or tests.

## Continuous-time rule

The fundamental state-evolution coordinate is physical time $t$.

Do not introduce fatigue cycle count as an independent state variable.

For constant frequency only, a cycle count may be reported afterward as

$$
N=ft.
$$

If the symbol $N$ is used in a finite empirical density, it must be clearly identified as a finite system/sample count rather than fatigue cycle count. Prefer another finite-size symbol such as $M$ when practical.

## Energy-feasibility rule

Any claim that stored configurational energy forces a tensile tail must state the additional compression-side constraint explicitly.

In particular, the exact safe-energy ceiling

$$
\mathcal E_{\rm safe}^{\max}(t)
$$

may be used only after declaring a physically justified lower support bound $\lambda_L(t)$ or an equivalent rigorous compression-side constraint.

Illustrative $\lambda_L$ values must never be presented as Al material constants.

## Variable dictionary rule

Any new active theory or simulation variable must be added to the appropriate variable-definition Markdown file in the same change.

Current active dictionaries include:

- `docs/VARIABLE_DEFINITIONS_NORMAL_LJ.md`;
- `docs/VARIABLE_DEFINITIONS_ENERGY_FEASIBILITY.md`;
- `docs/VARIABLE_DEFINITIONS_NORMAL_TIMESCALE.md`;
- `firmware/VARIABLE_DEFINITIONS.md`.

## Modeling labels

Important statements must be labeled as one of:

- **EXACT / IDENTITY**
- **DEFINITION**
- **ASSUMPTION**
- **CONTROLLED APPROXIMATION**
- **EMPIRICAL INPUT / PHYSICAL CONSTRAINT**

Never silently promote an assumption or controlled approximation to an exact result.

## Forbidden shortcuts

Do not:

- fit Gaussian, Weibull, or another named family to $P(a,t)$ merely for convenience;
- vary LJ parameters with time to imitate damage;
- insert an empirical damage variable and call it mechanics-derived;
- insert a retained-energy fraction without an energy balance;
- insert a fitted relaxation time or damping coefficient solely to obtain a desired fatigue life;
- use the energy ceiling without controlling reverse compression;
- convert a physical-time criterion into a cycle-life law and then treat the cycle count as fundamental.

## Numerical-result rule

Every numerical result must state:

- the exact model being evaluated;
- dimensional versus nondimensional quantities;
- numerical method;
- conservation/convergence check where applicable;
- whether any parameter is illustrative rather than a material input;
- whether the result is an exact theorem check, null test, proof of principle, or physical prediction.

## Korean code-file header rule

Every active source, simulation, test, and firmware code/header file must begin with a Korean navigation comment that states at least:

- the role of the file;
- the main functions/classes or public symbols;
- the implementation/output role when useful;
- a reminder that scientific classification belongs to the detailed docstrings and theory documents.

Use `python tools/add_korean_code_headers.py` to regenerate these headers after adding or renaming functions. Use `--check` to detect stale headers.

---

# 한국어 번역 — Repository 작성 규칙

## GitHub Markdown 수식

모든 Markdown 파일에서 GitHub-compatible math delimiter를 사용한다.

Inline math:

```markdown
$P(a,t)$
```

Display math:

```markdown
$$
\partial_tP+\partial_a(Pv)=0
$$
```

repository Markdown에서는 `\[`와 `\]`를 사용하지 않는다.

## 모든 Markdown의 한국어 번역 의무

**모든 `.md` 파일에는 완전한 한국어 번역이 있어야 한다.**

필수 순서는 다음과 같다.

1. 완전한 영문 technical source;
2. `---`;
3. 한국어 번역 heading;
4. 전체 한국어 번역.

짧은 한국어 summary는 충분하지 않다.

## Active-mainline 차원 규칙

활성 research root는 **1차원 normal-only**다.

활성 theory는

$$
a_i(t),
\qquad
P(a,t),
\qquad
\mu(t),
\qquad
\mathcal E(t)
$$

같은 1D normal coordinate로 표현할 수 있어야 한다.

3차원 FCC 연구는 `libraries/fcc_normal/`에 보존한다. shear/Rubin/slip code는 의도적으로 제외하고 `libraries/shear/README.md`만 비활성 tombstone으로 둔다. archive code를 기본 active simulation이나 test에서 import하면 안 된다.

## 연속시간 규칙

근본 state-evolution coordinate는 물리적 시간 $t$다.

fatigue cycle count를 독립적인 state variable로 도입하지 않는다.

주파수가 일정할 때만 필요하면 나중에

$$
N=ft
$$

로 cycle count를 표시할 수 있다.

finite empirical density에서 $N$을 사용할 경우 fatigue cycle count가 아니라 finite system/sample count라는 점을 명확히 해야 한다. 가능하면 finite-size symbol은 $M$처럼 다른 기호를 사용한다.

## Energy-feasibility 규칙

저장된 configurational energy가 tensile tail을 강제한다고 주장하려면 compression-side constraint를 반드시 명시해야 한다.

특히 정확한 safe-energy ceiling

$$
\mathcal E_{\rm safe}^{\max}(t)
$$

은 물리적으로 정당화된 lower support bound $\lambda_L(t)$ 또는 동등한 엄밀한 compression-side constraint를 선언한 뒤에만 사용할 수 있다.

illustrative $\lambda_L$ 값을 Al material constant로 제시하면 안 된다.

## 변수사전 규칙

새로운 활성 theory 또는 simulation variable은 같은 change에서 적절한 variable-definition Markdown에 추가한다.

현재 활성 dictionary는 다음을 포함한다.

- `docs/VARIABLE_DEFINITIONS_NORMAL_LJ.md`;
- `docs/VARIABLE_DEFINITIONS_ENERGY_FEASIBILITY.md`;
- `docs/VARIABLE_DEFINITIONS_NORMAL_TIMESCALE.md`;
- `firmware/VARIABLE_DEFINITIONS.md`.

## 모델링 분류 label

중요한 statement는 다음 중 하나로 분류한다.

- **EXACT / IDENTITY**
- **DEFINITION**
- **ASSUMPTION**
- **CONTROLLED APPROXIMATION**
- **EMPIRICAL INPUT / PHYSICAL CONSTRAINT**

assumption이나 controlled approximation을 설명 없이 exact result로 승격하지 않는다.

## 금지 shortcut

다음은 금지한다.

- 편의를 위해 $P(a,t)$에 Gaussian, Weibull 또는 다른 named family를 fitting하는 것;
- damage를 흉내내기 위해 시간에 따라 LJ parameter를 바꾸는 것;
- empirical damage variable을 넣고 mechanics-derived라고 부르는 것;
- energy balance 없이 retained-energy fraction을 넣는 것;
- 원하는 fatigue life를 만들기 위해 fitted relaxation time 또는 damping coefficient를 넣는 것;
- reverse compression을 제어하지 않고 energy ceiling을 사용하는 것;
- physical-time criterion을 cycle-life law로 바꾼 뒤 cycle count를 근본변수처럼 취급하는 것.

## 수치결과 규칙

모든 numerical result는 다음을 명시해야 한다.

- 평가한 정확한 model;
- dimensional / nondimensional quantity 구분;
- numerical method;
- 가능한 경우 conservation/convergence check;
- parameter가 material input이 아니라 illustrative 값인지 여부;
- exact theorem check, null test, proof of principle, physical prediction 중 무엇인지.

## 한국어 코드 파일 상단 안내 규칙

모든 활성 source, simulation, test, firmware code/header 파일은 맨 위에 최소한 다음 내용을 적은 한국어 탐색용 주석을 둔다.

- 파일이 하는 일;
- 주요 함수/클래스 또는 공개 symbol;
- 필요하면 구현 또는 출력 결과의 역할;
- 물리적 가정·근사·정확성 분류는 세부 docstring과 theory 문서를 따라야 한다는 주의.

함수 추가/이름변경 뒤에는 `python tools/add_korean_code_headers.py`로 헤더를 다시 생성한다. `--check` 옵션으로 오래된 헤더를 검사할 수 있다.
