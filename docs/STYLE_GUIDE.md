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

### Do not use

Do not use LaTeX display delimiters `\[` and `\]` in Markdown files because they may render as literal text in GitHub views.

Do not place mathematical expressions in backticks unless the intent is to show literal source code rather than rendered mathematics.

## Mandatory bilingual Markdown rule

**Every `.md` file in this repository must contain a Korean translation.**

Required order:

1. English technical version first;
2. horizontal rule `---`;
3. heading `# 한국어 번역` or an equivalent Korean title;
4. Korean translation of the complete technical content.

A short Korean summary is not a substitute for the translation when the English section contains additional technical claims, assumptions, numerical results, or limitations. Equations should be preserved unchanged between the English and Korean sections unless the mathematical content itself is being corrected.

This rule applies to README files, theory notes, assumptions, failed approaches, result records, planning documents, and future Markdown files.

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

저장소의 모든 Markdown 파일에서 수식은 GitHub에서 정상적으로 렌더링되는 delimiter를 사용해야 한다.

### 인라인 수식

단일 달러 기호를 사용한다.

```markdown
$P(a,t)$
```

### 독립된 수식 블록

이중 달러 기호를 사용한다.

```markdown
$$
\partial_t P + \partial_a(Pv)=0
$$
```

### 사용하지 말아야 할 형식

Markdown 파일에서는 LaTeX display delimiter인 `\[`와 `\]`를 사용하지 않는다. GitHub 화면에서 수식으로 렌더링되지 않고 문자 그대로 보일 수 있기 때문이다.

수식을 실제 코드 문자열로 보여주려는 경우가 아니라면 수학식을 backtick 안에 넣지 않는다.

## 모든 Markdown 파일의 한국어 번역 의무 규칙

**이 저장소의 모든 `.md` 파일에는 한국어 번역본이 반드시 포함되어야 한다.**

작성 순서는 다음과 같이 통일한다.

1. 영문 기술 원문을 먼저 작성한다.
2. 구분선 `---`를 넣는다.
3. `# 한국어 번역` 또는 이에 해당하는 한국어 제목을 넣는다.
4. 영문 기술내용 전체를 한국어로 번역한다.

영문 부분에 추가적인 기술 주장, 가정, 수치결과, 한계가 존재한다면 짧은 한국어 요약만으로 번역을 대신할 수 없다. 수학적 내용을 실제로 수정하는 경우가 아니라면 영문과 한국어 부분의 방정식은 동일하게 유지한다.

이 규칙은 README, 이론 노트, 가정 정리, 실패한 접근, 결과 기록, 연구계획 문서 및 앞으로 새로 만드는 모든 Markdown 파일에 적용한다.

## 모델링 분류 라벨

중요한 주장은 다음 중 하나로 분류해야 한다.

- **EXACT / IDENTITY**
- **DEFINITION**
- **ASSUMPTION**
- **CONTROLLED APPROXIMATION**
- **EMPIRICAL INPUT**

근사식을 아무 설명 없이 정확식으로 승격시켜서는 안 된다.
