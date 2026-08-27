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

## Bilingual research notes

For major research-planning documents, keep the English technical version first and add a Korean translation below a horizontal rule when useful. Equations should be shared exactly between the two versions unless the mathematical content itself changes.

## Modeling labels

Important claims should be identified as one of:

- **EXACT / IDENTITY**
- **DEFINITION**
- **ASSUMPTION**
- **CONTROLLED APPROXIMATION**
- **EMPIRICAL INPUT**

Never silently promote an approximation into an exact statement.
