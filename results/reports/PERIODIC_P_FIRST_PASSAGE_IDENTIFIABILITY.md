# Periodic P / first-passage identifiability audit

Status: **mathematical audit, not an aluminum calibration.**

## Question

Can a cycle-periodic spacing marginal $P(\lambda,t)$ produce cumulative fatigue first passage without a permanently drifting PDF shape?

The answer has two parts:

1. **yes**, cumulative first passage can increase while the normalized one-point shape repeats;
2. **no**, the snapshot $P$ alone does not tell us whether that accumulation occurs.

## Counterexample construction

Use $M=100$ represented labels. During each cycle, one label follows a smooth high-spacing excursion from $\lambda=1$ to $1.02$ and back, while the other 99 labels remain at $\lambda=1$. The threshold is $\lambda_c=1.01$.

At every phase of every cycle, the sorted empirical spacing set is identical regardless of which label is assigned the high excursion. Numerically,

$$
\max |P_A-P_B|_{\mathrm{sorted}}=0.
$$

Two deterministic histories are then compared.

### Same-label history

The same label takes the high excursion every cycle. It first crosses during cycle 1, and no new label ever crosses afterward.

Thus

$$
F_{\mathrm{ci}}=0.01
$$

for every cycle after the first.

### Rotating-label history

A new label takes the high excursion on each cycle. The instantaneous marginal remains exactly the same, but the cumulative first-passage fraction is

$$
F_{\mathrm{ci}}(n)=\min(n/100,1).
$$

Selected results:

| cycles | same-label $F_{\mathrm{ci}}$ | rotating-label $F_{\mathrm{ci}}$ |
|---:|---:|---:|
| 1 | 0.01 | 0.01 |
| 10 | 0.01 | 0.10 |
| 25 | 0.01 | 0.25 |
| 50 | 0.01 | 0.50 |
| 100 | 0.01 | 1.00 |

Therefore identical snapshot $P$ does not determine cumulative first passage.

## Consequence for the active formulation

The active crack-initiation document already contains the correct history-bearing object: the survivor phase-space subdensity $F_b$ and its marginal $P_b$.

The survivor mass is

$$
S(t)=\int P_b(\lambda,t)\,d\lambda,
$$

and

$$
F_{\mathrm{ci}}(t)=1-S(t).
$$

The normalized survivor shape can be periodic while its mass decays. If the integrated survivor hazard over one cycle is repeated,

$$
S_{n+1}=rS_n,
$$

so

$$
S_n=r^nS_0,
\qquad
F_{\mathrm{ci},n}=1-r^nS_0.
$$

Hence a permanently drifting normalized $P$ is not required for fatigue accumulation.

## What does not work

A stronger condition is label-wise periodicity:

$$
\lambda_i(t+T)=\lambda_i(t)
$$

for every represented label $i$.

Then exactly the same labels cross every cycle. The cumulative first-passage set is fixed after the first cycle, so repeating a quasistatic label-preserving trajectory cannot generate high-cycle accumulation.

This specifically blocks the simplest independent local-traction / reversible quasistatic picture as a fatigue mechanism.

## Extra stochastic diagnostic

For illustration only, if a surviving local domain were instead given an independent conditional crossing probability

$$
p_{\mathrm{cyc}}=10^{-5}
$$

on each cycle, then

$$
F_n=1-(1-p_{\mathrm{cyc}})^n.
$$

This gives approximately $F=0.632$ at $10^5$ cycles. This is **not adopted physics**; it only illustrates how a stationary mechanical distribution can coexist with cumulative rare-event survival when stochasticity supplies renewed trials.

## Verdict

The correct logical structure is

$$
\boxed{\text{periodic normalized shape does not rule out cumulative first passage}}
$$

but also

$$
\boxed{\text{snapshot }P\text{ alone is insufficient to compute cumulative first passage}.}
$$

The next model requirement is therefore narrower than forcing irreversible PDF drift. We need a predictive law for the **survivor escape flux** or an equivalent survivor propagator under cyclic stress. That law must explain why new surviving domains receive crossing opportunities on later cycles through physically justified mixing, thermal stochasticity, or evolving internal mechanics.
