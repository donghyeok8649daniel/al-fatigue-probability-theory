# Periodic P versus cumulative first passage

Status: **EXACT CLARIFICATION / no new physical law.**

This note separates two questions that must not be conflated:

1. whether the normalized spacing distribution changes permanently from cycle to cycle;
2. whether the cumulative crack-initiation first-passage fraction increases from cycle to cycle.

The second does not mathematically require the first. However, the snapshot marginal alone is not sufficient to determine cumulative first passage.

## 1. Existing active definitions

The active crack-initiation definition already distinguishes the nonabsorbing spacing marginal from the survivor subdensity.

Let

$$
P(\lambda,\tau)
$$

be the ordinary instantaneous spacing marginal, and let

$$
P_b(\lambda,\tau)
$$

be the survivor spacing subdensity obtained from the absorbing survivor phase-space density. Its mass is

$$
S(\tau)=\int_0^{\lambda_c}P_b(\lambda,\tau)\,d\lambda,
$$

and

$$
F_{\mathrm{ci}}^{\mathrm{local}}(\tau)=1-S(\tau).
$$

The normalized spacing density conditioned on survival is

$$
\widehat P_b(\lambda,\tau)=\frac{P_b(\lambda,\tau)}{S(\tau)}.
$$

Therefore

$$
P_b(\lambda,\tau)=S(\tau)\widehat P_b(\lambda,\tau).
$$

The escape-flux identity is

$$
\frac{dS}{d\tau}=-j_{\mathrm{esc}},
$$

and for $S>0$ the active hazard definition gives

$$
h_\tau=\frac{j_{\mathrm{esc}}}{S}
=-\frac{d}{d\tau}\ln S.
$$

No new damage variable is introduced here.

## 2. A periodic normalized shape can coexist with cumulative loss

Let the loading period be $T$. Suppose that after transients the survivor-conditioned shape is cycle-periodic,

$$
\widehat P_b(\lambda,\tau+T)
=
\widehat P_b(\lambda,\tau).
$$

This does not imply that the survivor mass is periodic.

From the hazard identity,

$$
S(\tau+T)
=
S(\tau)
\exp\left[-\int_\tau^{\tau+T}h_\tau(s)\,ds\right].
$$

Define the one-cycle survival factor

$$
r(\tau)
=
\exp\left[-\int_\tau^{\tau+T}h_\tau(s)\,ds\right].
$$

If the periodic regime also makes the integrated hazard repeat from cycle to cycle, then $r$ is constant and

$$
\boxed{
S(\tau+nT)=r^nS(\tau).
}
$$

Therefore

$$
\boxed{
F_{\mathrm{ci}}^{\mathrm{local}}(\tau+nT)
=1-r^nS(\tau).
}
$$

At the same time,

$$
\widehat P_b(\lambda,\tau+nT)
=
\widehat P_b(\lambda,\tau).
$$

Thus progressive drift of the normalized spacing-distribution shape is **not mathematically necessary** for cumulative first passage.

The history can reside in the decreasing survivor mass rather than in a drifting normalized PDF.

## 3. But the ordinary snapshot P alone does not determine first passage

This is the crucial identifiability point.

Consider $M$ represented labels. During each loading cycle, exactly one label follows a high-spacing excursion that crosses the threshold, while the other $M-1$ labels follow the same low-spacing history. At every phase of every cycle, the empirical spacing marginal is therefore identical because it contains the same one high trajectory and the same $M-1$ low trajectories.

Two deterministic histories can nevertheless have different first-passage fractions.

### History A: the same label is high every cycle

The same represented spacing crosses on every cycle. Then after the first cycle,

$$
F_{\mathrm{ci},M}^{\mathrm{local}}=\frac{1}{M}
$$

and it does not increase further.

### History B: the high label rotates from cycle to cycle

A different represented spacing crosses on each cycle. Then after $n$ cycles,

$$
F_{\mathrm{ci},M}^{\mathrm{local}}
=
\min\left(\frac{n}{M},1\right).
$$

History A and History B have the **same instantaneous marginal $P(\lambda,\tau)$ at every phase**, but different cumulative first-passage histories.

Therefore

$$
\boxed{
P(\lambda,\tau)\ \text{alone does not determine}\ F_{\mathrm{ci}}(\tau).
}
$$

The missing information is not necessarily permanent PDF drift. It is the survivor/path identity information carried by $F_b$, $P_b$, or an equivalent first-passage propagator.

## 4. Stronger result for label-wise periodic deterministic motion

Suppose every individual represented trajectory itself repeats exactly after one cycle,

$$
\lambda_i(\tau+T)=\lambda_i(\tau)
$$

for every $i$ after the periodic regime is reached.

Then the set of labels that cross $\lambda_c$ during cycle 2 is exactly the same set that crossed during cycle 1. Consequently

$$
\boxed{
F_{\mathrm{ci},M}^{\mathrm{local}}(nT)
=F_{\mathrm{ci},M}^{\mathrm{local}}(T),
\qquad n\ge1.
}
$$

Thus a perfectly reversible, label-preserving quasistatic cycle cannot produce many-cycle first-passage accumulation by repeating the same threshold test.

This is stronger than saying that the marginal $P$ is periodic. Marginal periodicity still allows hidden label rearrangement; label-wise periodicity does not.

## 5. What can produce cumulative first passage without permanent P-shape drift?

At least one of the following must exist:

- deterministic rearrangement or mixing of represented labels that is invisible to the one-point marginal $P$;
- finite-temperature stochastic transitions so that a surviving local state receives a new rare-event opportunity on later cycles;
- an evolving internal structural state that changes which labels are exposed to the threshold.

This note **does not adopt** any of these mechanisms. It only proves what information is mathematically required.

In particular, the earlier independent local-traction oscillator with a label-wise periodic quasistatic response cannot obtain high-cycle fatigue merely by reusing the same snapshot tail on every cycle.

## 6. Reduced target for fatigue initiation

If all represented domains are initially intact, the initial survivor subdensity is fixed by the initial spacing distribution,

$$
P_b(\lambda,0)
=
P_0(\lambda)I[\lambda<\lambda_c].
$$

Hence no independent empirical damage initial condition is required.

The useful reduced target can therefore be written as

$$
\boxed{
P_0+\sigma(0:t)
\longrightarrow
P_b(\lambda,t),\ S(t),\ F_{\mathrm{ci}}(t).
}
$$

The unresolved physics is the predictive survivor transport / escape flux. A closed model must determine which surviving population reaches the absorbing boundary on later cycles without reconstructing the full microscopic trajectory ordering.

## 7. Research consequence

The earlier requirement

$$
P_0+\sigma(0:t)\to P(t)
$$

remains useful for mechanical observables, but fatigue initiation should not be made to depend on permanent drift of the normalized $P$ unless the mechanics actually produces such drift.

The more precise separation is

$$
\boxed{
\text{mechanical state shape}
\quad\text{versus}\quad
\text{survivor mass / first-passage history}.
}
$$

A periodic normalized shape with nonzero survivor escape per cycle is sufficient for cumulative initiation. A label-wise periodic deterministic trajectory with no new survivor escape is not.

This result narrows the next modeling question to the origin of the survivor escape flux under laboratory cyclic loading.