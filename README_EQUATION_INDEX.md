# README — Equation & Symbol Index / 수식·기호 Index

This is the entry point for the active mathematical theory.  
이 파일은 현재 활성 이론의 **수식·기호·용어 색인 진입점**이다.

## 1. Authoritative documents / 기준 문서

1. [`docs/EQUATION_SUMMARY_1D_P_U_THETA.md`](docs/EQUATION_SUMMARY_1D_P_U_THETA.md) — **Equation sheet / 핵심 수식 정리본**
2. [`docs/VARIABLE_INDEX_1D_P_U_THETA.md`](docs/VARIABLE_INDEX_1D_P_U_THETA.md) — **Bilingual symbol dictionary / 영·한 기호 사전**
3. [`docs/AUXILIARY_SYMBOL_INDEX_1D.md`](docs/AUXILIARY_SYMBOL_INDEX_1D.md) — **Auxiliary symbol dictionary / 보조기호 사전**
4. [`docs/MASTER_1D_P_U_THETA_FORMULATION.md`](docs/MASTER_1D_P_U_THETA_FORMULATION.md) — **Full differential derivation / 전체 미분형 유도**
5. [`docs/MILESTONE25_EXACT_INTEGRAL_REPRESENTATION.md`](docs/MILESTONE25_EXACT_INTEGRAL_REPRESENTATION.md) — **Exact integral representations / 정확한 적분해 표현**
6. [`docs/CRACK_INITIATION_DEFINITION.md`](docs/CRACK_INITIATION_DEFINITION.md) — **Kinetic first passage / 위상공간 최초통과 균열개시**

## 2. Mandatory symbol-definition contract / 기호 정의 강제 규칙

A new mathematical symbol is not considered defined unless the symbol index is updated at the same time.  
새 수학기호를 도입할 때 기호 Index를 동시에 갱신하지 않으면 이론상 정의된 기호로 인정하지 않는다.

Every new symbol entry MUST include:

| Required field | English | 한국어 |
|---|---|---|
| Symbol | exact LaTeX symbol | 정확한 LaTeX 기호 |
| Equation definition | defining equation or operator identity | 정의식 또는 연산자 항등식 |
| English term | standard English name | 영문 명칭 |
| Korean term | Korean name | 한국어 명칭 |
| Mathematical definition | domain, conditioning, measure, operation, scalar/vector role | 정의역, 조건, 측도, 연산, 스칼라/벡터 역할 |
| Physical definition | represented physical quantity/process | 나타내는 물리량/과정 |
| Unit / scaling | SI unit or nondimensional scaling | SI 단위 또는 무차원 스케일 |
| Status | MODEL / DEFINITION / EXACT / CONDITIONAL / OPEN | 모델 / 정의 / 정확식 / 조건부 / 미완성 |
| Dependencies | symbols/equations required first | 선행 기호/식 |

A prose-only definition is insufficient when a mathematical defining relation exists.  
수식으로 정의 가능한 기호는 문장 설명만으로 정의하지 않는다.

A symbol must not silently change meaning between files. If the same glyph would represent a different mathematical object, rename or decorate the symbol explicitly.

## 3. Markdown and math-rendering rule / Markdown·수식 렌더링 규칙

- Display math uses `$$ ... $$` only.
- Inline math uses `$...$` only.
- Do not use `\[` or `\]` in Markdown math.
- Do not use `\(` or `\)` in Markdown math.
- `\operatorname{...}` and `\operatorname*{...}` are forbidden because the target renderer rejects them. Use `\mathrm{...}` instead.
- Write expectation as `\mathbb{E}`, not `\mathbb E`.
- Write indicator symbols with explicit braces, for example `\mathbf{1}`.
- Do not place display-math blocks inside Markdown tables.
- Keep a blank line before and after every display-math block.
- Use fenced code blocks only for literal code or plain-text pseudocode, not for equations that should render as math.

## 4. Active mathematical chain / 활성 수학 체계

$$
\boxed{
\text{1D nonlinear generalized-LJ chain}
\rightarrow
\Phi^q_{\tau,\tau_0}
\rightarrow
F(\lambda,c,\tau)
\rightarrow
\{P,u,\Theta\}
\rightarrow
\{\bar a,\bar U,S,F_{\rm ci}\}
}
$$

The finite LJ dynamics is closed. The one-point $P$–$u$–$\Theta$ PDEs are exact but hierarchical, and the same reduced fields also possess exact full-flow integral representations.

## 5. Governing-equation index / 지배방정식 Index

### E01 — generalized-LJ energy / generalized-LJ 에너지

$$
\boxed{
\phi(\lambda)
=
\frac{\lambda^{-m}}{m(m-n)}
-
\frac{\lambda^{-n}}{n(m-n)}
}
$$

### E02 — bulk spacing dynamics / 내부 spacing 동역학

$$
\boxed{
\ddot\lambda_i
=
\phi'(\lambda_{i+1})-2\phi'(\lambda_i)+\phi'(\lambda_{i-1})
}
$$

### E03 — empirical phase-space measure / 경험적 위상공간 측도

$$
\boxed{
F_M(\lambda,c,\tau)
=\frac1M\sum_i
\delta(\lambda-\lambda_i)\delta(c-c_i)
}
$$

### E04 — projected kinetic transport / 투영 위상공간 수송

$$
\boxed{
\partial_\tau F+\partial_\lambda(cF)+\partial_c(AF)=0
}
$$

### E05 — continuity / 연속방정식

$$
\boxed{
\partial_\tau P+\partial_\lambda(Pu)=0
}
$$

### E06 — mean-spacing-rate equation / 조건부 평균 spacing-rate 식

$$
\boxed{
D_\tau u
=\mathcal A-\frac1P\partial_\lambda(P\Theta)
}
$$

### E07 — exact $\Theta$ equation / 정확한 $\Theta$ 식

$$
\boxed{
D_\tau\Theta
+2\Theta\partial_\lambda u
+\frac1P\partial_\lambda(PC_3)
=2\Psi
}
$$

### E08 — density-shape identity / 확률밀도 형상 항등식

$$
\boxed{
\Theta\partial_\lambda\ln P
=\mathcal A-D_\tau u-\partial_\lambda\Theta
}
$$

### E09 — instantaneous integral form of $P$ / $P$의 순간 적분형

$$
\boxed{
P(\lambda,\tau)
=
\frac{\mathcal N_P(\tau)}{\Theta(\lambda,\tau)}
\exp\!\left[
\int_{\lambda_*}^{\lambda}
\frac{\mathcal A-D_\tau u}{\Theta}\,d\eta
\right]
}
$$

### E10 — exact full-flow push-forward / 전체 흐름의 정확한 push-forward

$$
\boxed{
F(\lambda,c,\tau)
=\frac1M\sum_i\int
\delta[\lambda-\Lambda_i(\tau;\Gamma_0)]
\delta[c-C_i(\tau;\Gamma_0)]
\,\mu_0(d\Gamma_0)
}
$$

### E11 — exact $\Theta$ characteristic integral / $\Theta$의 정확한 특성곡선 적분해

$$
\boxed{
\Theta(X(\tau),\tau)
=e^{-2\mathcal I_u(\tau;\alpha)}
\left[
\Theta_0(\alpha)
+\int_{\tau_0}^{\tau}
 e^{2\mathcal I_u(s;\alpha)}S_\Theta(X(s),s)\,ds
\right]
}
$$

with

$$
\boxed{
S_\Theta=2\Psi-\frac1P\partial_\lambda(PC_3)
}
$$

### E12 — G1 mean spacing / G1 평균 간격

$$
\boxed{
\bar a=a_0\int\lambda P(\lambda,\tau)\,d\lambda
}
$$

### E13 — G2 mean intrinsic configurational energy / G2 평균 고유 배치에너지

$$
\boxed{
\bar U
=U_{\rm ref}\int[\phi(\lambda)-\phi(1)]P(\lambda,\tau)\,d\lambda
}
$$

### E14 — G3 irreversible history / G3 비가역 에너지 이력

$$
\boxed{
E_{\rm hyst}(t)=\int_0^t\dot D_{\rm irr}(t')\,dt',
\qquad\dot D_{\rm irr}\ge0
}
$$

Current conservative baseline:

$$
\boxed{\dot D_{\rm irr}=0}
$$

### E15 — local stability threshold / 국소 안정성 상실점

$$
\boxed{
\phi''(\lambda_c)=0,
\qquad
\lambda_c=\left(\frac{m+1}{n+1}\right)^{1/(m-n)}
}
$$

### E16 — specimen survival integral / 시편 생존 적분식

$$
\boxed{
S_{\rm spec}(\tau)
=\int
\mathbf{1}\!\left[
\max_i\sup_{s\in[\tau_0,\tau]}\Lambda_i(s;\Gamma_0)<\lambda_c
\right]\mu_0(d\Gamma_0)
}
$$

## 6. Symbol groups / 기호 그룹

The complete entries live in `docs/VARIABLE_INDEX_1D_P_U_THETA.md` and `docs/AUXILIARY_SYMBOL_INDEX_1D.md`.

- **Scales / 스케일:** $t,\tau,t_0,m_a,a_0,E,A_0,U_{\rm ref},F_{\rm ref}$
- **Loading / 하중:** $\sigma_n,\sigma_m,\sigma_a,f,\omega^*,F_{\rm ext},q$
- **Microscopic chain / 미시 사슬:** $M,x_i,\lambda_i,a_i,c_i,m,n,\phi,V^*,T^*,E_{\rm mech}^*,\mathbf L,\mathbf G_\lambda$
- **Projected probability / 투영 확률:** $F_M,P_M,F,P,A,\mathcal A,J$
- **Moments / 모멘트:** $u,\Theta,C_3,\Psi,R_r,B_r$
- **Neighbour correlations / 이웃 상관:** $P_2^\pm,F_2^\pm,m_\pm$
- **Integral solution / 적분해:** $\Gamma,\Gamma_0,\Phi^q,\mu_0,\Lambda_i,C_i,A_i,X,\alpha,\mathcal I_u,S_\Theta,\lambda_*,\mathcal N_P$
- **Observables / 관측량:** $\bar\lambda,\bar a,\Delta\phi,\bar U,\dot D_{\rm irr},E_{\rm hyst}$
- **First passage / 최초통과:** $\lambda_c,\tau_i^c,\chi_i,F_b,P_b,S,j_{\rm esc},h,\widehat P_b,S_{\rm local},S_{\rm spec},F_{\rm ci}^{\rm local},F_{\rm ci}^{\rm spec}$

## 7. Term glossary / 핵심 용어 Index

| English term | 한국어 | Mathematical role / 수학적 역할 | Physical role / 물리적 역할 |
|---|---|---|---|
| normalized spacing | 무차원 원자층간격 | $\lambda=a/a_0$ | local normal opening/stretch |
| spacing rate | 간격 변화율 | $c=d\lambda/d\tau$ | local opening/closing rate |
| empirical measure | 경험적 측도 | finite sum of Dirac masses | distribution generated directly from represented spacings |
| phase-space density | 위상공간 밀도 | density in $(\lambda,c)$ | joint population of spacing and spacing-rate states |
| marginal density | 주변밀도 | $P=\int F\,dc$ | spacing distribution after velocity information is projected out |
| conditional mean | 조건부 평균 | $u=\mathbb{E}[c\mid\lambda]$ | mean opening/closing rate at fixed spacing |
| conditional variance | 조건부 분산 | $\Theta=\mathrm{Var}(c\mid\lambda)$ | unresolved spread of spacing rates at fixed spacing |
| conditional acceleration | 조건부 가속도 | $\mathcal A=\mathbb{E}[\ddot\lambda\mid\lambda]$ | mean microscopic acceleration at a given spacing |
| acceleration covariance | 가속도 공분산 | $\Psi=\mathrm{Cov}(c,\ddot\lambda\mid\lambda)$ | coupling between rate and acceleration fluctuations |
| moment hierarchy | 모멘트 계층 | $\partial_\tau R_r+\partial_\lambda R_{r+1}=rB_r$ | information lost under projection |
| push-forward | 푸시포워드/전방사상 측도 | image of $\mu_0$ under $\Phi^q$ | evolution of an initial microscopic population under deterministic mechanics |
| characteristic curve | 특성곡선 | $dX/ds=u(X,s)$ | trajectory in spacing space transported by mean flow |
| Volterra integral equation | 볼테라 적분방정식 | present state expressed using past-time integrals | exact time-history form of reduced balances |
| first passage | 최초통과 | first hitting of $\lambda_c$ | local mechanical initiation event |
| survival | 생존 | probability/mass not yet first-passed | intact local/specimen population |
| hazard | 위험률/개시율 | $h=-d\ln S/d\tau$ | instantaneous initiation rate conditional on survival |
| absorbing boundary | 흡수경계 | no inflow from failed side | prevents a failed trajectory from being counted as intact again |
| closure | 폐쇄/클로저 | higher statistics expressed through retained fields | extra assumption needed for an autonomous reduced solver |
| history dependence | 이력의존성 | same forcing value can have different reduced state | loading and unloading microscopic states do not retrace |
| irreversible dissipation | 비가역 소산 | $\dot D_{\rm irr}\ge0$ | energy permanently removed from recoverable mechanical storage |

## 8. Scope warning / 범위 주의

The current theory has an exact mathematical integral representation for the reduced state, but **G3 irreversible physics, the physical specimen measure $\mu_0$, laboratory fatigue time-scale bridging, and experimental validation remain open physical problems.**

현재 이론은 축약상태의 정확한 적분 표현까지 갖지만, **G3 비가역 물리, 실제 시편의 $\mu_0$, 실험실 피로 시간척도 연결, 실험검증은 아직 물리적으로 미완성**이다.
