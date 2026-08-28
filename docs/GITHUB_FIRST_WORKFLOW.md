# GitHub-First Research Workflow

## Status

GitHub is the **source of truth** for this project.

All project-relevant theory development, numerical work, implementation, falsification, failed approaches, variable changes, interpretation changes, and next-step decisions must be reflected in the repository as the work is performed. The chat must not become the only place where a substantive result exists.

## Required workflow

For every substantive research step:

1. **Inspect current `main` first.** Reuse the current repository state rather than relying on an older chat description.
2. **Write the derivation or decision into the repository while doing the work.** Do not postpone documentation until the end of a long reasoning chain.
3. **Update variable definitions in the same change** whenever a new symbol, state variable, diagnostic, or parameter is introduced.
4. **Implement executable theory and numerical checks** under `theory/`, `simulations/`, and `tests/` when the claim is computationally testable.
5. **Store generated evidence** under `results/data/`, `results/figures/`, and `results/reports/` when numerical work produces a reusable result.
6. **Record failed or falsified paths** rather than silently deleting them. Use `docs/FAILED_APPROACHES.md`, a milestone document, or an archived library as appropriate.
7. **Update `results/data/result_manifest.json`** whenever the active files, current result, or next target materially changes.
8. **Update the README** when the active mainline or top-level scientific interpretation changes.
9. **Run validation before claiming success.** At minimum, run the relevant unit tests and numerical convergence/conservation/null checks that apply to the change.
10. **Verify GitHub after writing.** Confirm the final `main` head and fetch the important changed file(s) before reporting that the repository has been updated.

## Commit discipline

Prefer a coherent commit for one scientific step. If connector or API limitations require several sequential commits, the final repository state must still be internally consistent and the final `main` head must be verified before reporting completion.

Do not claim that a file, result, or commit exists until it has been verified on GitHub.

## Theory-documentation mapping

Use the following mapping by default:

- exact identities and derivations: `docs/EXACT_DERIVATIONS.md` or a dedicated milestone;
- active research milestones: `docs/MILESTONE*.md`;
- assumptions and controlled approximations: `docs/ASSUMPTIONS.md` and the relevant milestone;
- failed/falsified approaches: `docs/FAILED_APPROACHES.md`;
- variables: `docs/VARIABLE_DEFINITIONS_*.md`;
- reusable theory code: `theory/`;
- experiment/simulation drivers: `simulations/`;
- regression and falsification tests: `tests/`;
- numerical tables and machine-readable summaries: `results/data/`;
- figures: `results/figures/`;
- human-readable result interpretation: `results/reports/`.

## Scientific integrity rule

Repository synchronization must not turn an assumption into a result. Every important statement must still be labeled according to the project classification system:

- **EXACT / IDENTITY**
- **DEFINITION**
- **ASSUMPTION**
- **CONTROLLED APPROXIMATION**
- **EMPIRICAL INPUT / PHYSICAL CONSTRAINT**
- **NUMERICAL RESULT / DIAGNOSTIC** when appropriate.

A failed test is a repository result and must be preserved just as carefully as a successful test.

---

# 한국어 번역 — GitHub 우선 연구 워크플로

## 상태

이 프로젝트에서 GitHub를 **source of truth**로 사용한다.

프로젝트와 관련된 모든 이론 전개, 수치계산, 구현, 반증시험, 실패한 접근, 변수 변경, 해석 변경, 다음 단계 결정은 실제 작업을 진행하면서 repository에 반영해야 한다. 중요한 결과가 chat에만 존재해서는 안 된다.

## 필수 워크플로

중요한 연구 단계를 수행할 때마다 다음을 지킨다.

1. **항상 현재 `main`을 먼저 확인한다.** 오래된 chat 설명에만 의존하지 않고 현재 repository 상태를 기준으로 작업한다.
2. **이론 유도나 결정을 진행하면서 repository에 기록한다.** 긴 reasoning이 끝난 뒤 문서화를 미루지 않는다.
3. 새로운 기호, 상태변수, diagnostic, parameter를 도입하면 **같은 change에서 변수정의 문서도 갱신한다.**
4. 계산으로 검증 가능한 주장은 가능한 경우 `theory/`, `simulations/`, `tests/`에 **실행 가능한 코드와 검증을 함께 구현한다.**
5. 재사용할 가치가 있는 수치결과가 나오면 **`results/data/`, `results/figures/`, `results/reports/`에 증거를 저장한다.**
6. 실패하거나 반증된 경로를 조용히 삭제하지 않는다. 필요에 따라 `docs/FAILED_APPROACHES.md`, milestone 문서, 또는 archive library에 보존한다.
7. active file, current result, next target이 실질적으로 바뀌면 **`results/data/result_manifest.json`을 갱신한다.**
8. active mainline 또는 최상위 과학적 해석이 바뀌면 **README를 갱신한다.**
9. 성공을 주장하기 전에 검증한다. 최소한 관련 unit test와 적용 가능한 convergence, conservation, null test를 수행한다.
10. **쓰기 작업 뒤 GitHub를 다시 확인한다.** 완료했다고 보고하기 전에 최종 `main` head와 중요한 변경파일을 다시 fetch해서 실제 반영 여부를 검증한다.

## Commit 규율

가능하면 하나의 과학적 단계는 하나의 coherent commit으로 묶는다. connector/API 제약 때문에 여러 sequential commit이 필요하더라도 최종 repository는 내부적으로 일관되어야 하며, 완료 보고 전에 최종 `main` head를 확인해야 한다.

GitHub에서 확인하지 않은 file, result, commit이 존재한다고 말하지 않는다.

## 이론-문서 매핑

기본적으로 다음 구조를 사용한다.

- exact identity와 derivation: `docs/EXACT_DERIVATIONS.md` 또는 전용 milestone;
- active research milestone: `docs/MILESTONE*.md`;
- assumption과 controlled approximation: `docs/ASSUMPTIONS.md` 및 관련 milestone;
- 실패/반증된 접근: `docs/FAILED_APPROACHES.md`;
- 변수: `docs/VARIABLE_DEFINITIONS_*.md`;
- 재사용 가능한 theory code: `theory/`;
- experiment/simulation driver: `simulations/`;
- regression 및 falsification test: `tests/`;
- numerical table과 machine-readable summary: `results/data/`;
- figure: `results/figures/`;
- 사람이 읽는 결과해석: `results/reports/`.

## 과학적 무결성 규칙

GitHub 반영을 철저히 한다고 해서 assumption을 result로 승격해서는 안 된다. 중요한 statement는 계속 다음 project classification을 따른다.

- **EXACT / IDENTITY**
- **DEFINITION**
- **ASSUMPTION**
- **CONTROLLED APPROXIMATION**
- **EMPIRICAL INPUT / PHYSICAL CONSTRAINT**
- 필요한 경우 **NUMERICAL RESULT / DIAGNOSTIC**

실패한 test도 성공한 test와 동일하게 repository에 보존해야 하는 연구결과다.
