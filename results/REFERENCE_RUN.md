# Reference Run — Rubin-chain mechanics-derived hysteresis

The following numbers were produced from the reference nondimensional case implemented in `theory/rubin_chain.py`.

## Parameters

- $M=m=1$
- $K_0=k=1$
- $F_a=0.1$
- $\omega=0.5$
- chain band edge $\omega_D=2$
- finite-chain simulation: 1200 masses, $\Delta t=0.02$, 60 periods, 5-period smooth ramp
- loop statistics: cycles 10 through 49, before a far-boundary reflection returns

## Analytic semi-infinite result

$$
Z(\omega)=0.875+0.4841229182759271i
$$

$$
Q_a=0.1
$$

$$
\phi=0.5053605102841573\ \mathrm{rad}=28.95502437185985^\circ
$$

$$
\boxed{A_H^{\mathrm{analytic}}=0.015209170034901047}
$$

## Full conservative finite-chain integration

$$
\boxed{\langle A_H^{\mathrm{numeric}}\rangle=0.015208839984912282}
$$

Cycle-to-cycle standard deviation:

$$
1.921149725428978\times10^{-7}
$$

Relative analytic-vs-numeric loop-area error:

$$
\boxed{2.1700723182610268\times10^{-5}}
$$

Final internal energy:

$$
E_{\mathrm{int}}=0.8644685639287875
$$

Integrated external work:

$$
W_{\mathrm{ext}}=0.8644577442919658
$$

Relative energy-balance error:

$$
\boxed{1.2516096816928344\times10^{-5}}
$$

## Interpretation

This run demonstrates that a nonzero loop in a **resolved coordinate** can arise from a fully conservative microscopic chain without adding a viscous damping coefficient. The loop area is energy transferred into propagating unresolved modes.

This is a Milestone-1 proof-of-principle, not yet a quantitative aluminum fatigue prediction and not yet Milestone 2 secular fatigue accumulation.
