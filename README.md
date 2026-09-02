# theory-core

Physical and mathematical source of truth for the active aluminum fatigue-probability theory.
현재 알루미늄 피로확률 이론의 물리·수학적 source of truth.

## Start here / 수식부터 보기

- [`README_EQUATION_INDEX.md`](README_EQUATION_INDEX.md) — **equation & symbol entry index / 수식·기호 진입 Index**
- [`docs/EQUATION_SUMMARY_1D_P_U_THETA.md`](docs/EQUATION_SUMMARY_1D_P_U_THETA.md) — **compact governing-equation sheet / 핵심 지배방정식 정리**
- [`docs/VARIABLE_INDEX_1D_P_U_THETA.md`](docs/VARIABLE_INDEX_1D_P_U_THETA.md) — **authoritative bilingual mathematical symbol dictionary / 기준 영·한 수학기호 사전**
- [`docs/MASTER_1D_P_U_THETA_FORMULATION.md`](docs/MASTER_1D_P_U_THETA_FORMULATION.md) — full differential derivation / 전체 미분형 유도
- [`docs/MILESTONE25_EXACT_INTEGRAL_REPRESENTATION.md`](docs/MILESTONE25_EXACT_INTEGRAL_REPRESENTATION.md) — exact push-forward, Volterra, characteristic and survival integrals / 정확한 적분해
- [`docs/CRACK_INITIATION_DEFINITION.md`](docs/CRACK_INITIATION_DEFINITION.md) — kinetic first-passage initiation / 위상공간 최초통과 균열개시

## Mandatory notation rule / 기호 정의 강제 규칙

A new mathematical symbol is not considered defined unless `docs/VARIABLE_INDEX_1D_P_U_THETA.md` is updated at the same time with:

\[
\boxed{
\text{equation definition}
+\text{English term}
+\text{Korean term}
+\text{mathematical meaning}
+\text{physical meaning}
+\text{unit/scaling}
+\text{status}
+\text{dependencies}
}
\]

새 기호는 위 항목을 Index에 동시에 추가해야 이론상 정의된 것으로 인정한다. 수식으로 정의 가능한 기호는 문장 설명만으로 끝내지 않는다.

## Active theory / 활성 이론

\[
\boxed{
\text{1D nonlinear LJ chain}
\rightarrow
\Phi^q
\rightarrow
F(\lambda,c,\tau)
\rightarrow
\{P,u,\Theta\}
\rightarrow
\{\bar a,\bar U,\text{first passage}\}
}
\]

The finite microscopic LJ system is closed. The reduced $P$–$u$–$\Theta$ equations are exact but hierarchical, and the same reduced fields also possess exact full-flow integral representations.

유한 미시 LJ 계는 닫혀 있다. 축약 $P$–$u$–$\Theta$ 식은 정확하지만 계층적이며, 동일한 축약장은 전체 미시흐름의 정확한 적분 투영으로도 표현된다.

## Branch ownership

Owned here:
- `theory/`: active mathematical mechanics and probability identities.
- `docs/`: derivations, assumptions, equation sheets, symbol definitions, and open theory problems.
- `libraries/`: theory-side reference libraries.

Not owned here: FEM/UI simulations, generated numerical results, desktop packaging, manuscript sources, or fatigue-tester hardware/firmware.

Validated theory flows downstream to `numerical-fem`, then through `integration`.
