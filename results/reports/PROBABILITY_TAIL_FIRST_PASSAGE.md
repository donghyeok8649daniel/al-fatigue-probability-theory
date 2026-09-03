# Probability-tail and first-passage validation

## Status

This is a deterministic numerical audit of the active conservative 1D normal generalized-LJ chain. It does **not** introduce an irreversible damage law or claim that snapshot tail mass is itself cumulative fatigue damage.

Protocol:

- 512 atoms
- dimensionless time step: 0.01
- mean stress mapping: 100 MPa / 69 GPa
- stress amplitude mapping: 100 MPa / 69 GPa
- dimensionless angular frequency: 0.02
- two-cycle smooth ramp
- six cycles total
- generalized-LJ exponents: m = 12.19, n = 6

## 1. Snapshot tail

For a fixed threshold lambda_*,

$$
Q_*(t)=\int_{\lambda_*}^{\infty}P_\lambda(\lambda,t)\,d\lambda.
$$

The exact continuity equation gives

$$
\frac{dQ_*}{d\tau}=J_\lambda(\lambda_*,\tau)
$$

when the upper-boundary flux vanishes.

For the finite empirical measure this becomes a crossing balance: the change in the number of spacings above the threshold equals upward crossings minus downward crossings.

At lambda_* = 1.004 the audit gives:

| cycle | upward crossings | downward crossings | net | cycle-end count above 1.004 |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 0 | 0 | 0 | 0 |
| 2 | 0 | 0 | 0 | 0 |
| 3 | 22 | 3 | 19 | 19 |
| 4 | 251 | 134 | 117 | 136 |
| 5 | 358 | 363 | -5 | 131 |
| 6 | 330 | 386 | -56 | 75 |

The count change equals the crossing net at every cycle.

Therefore the snapshot tail is **not monotone**. After growing, it can decrease when return crossings dominate.

## 2. Cumulative first passage

Define the first-passage fraction through a threshold by

$$
F_*(t)=\frac{1}{M}\sum_i I\left[\sup_{0\le s\le t}\lambda_i(s)\ge\lambda_*\right].
$$

Unlike the snapshot tail, this quantity cannot decrease.

For lambda_* = 1.004 the computed first-passage fractions at cycle ends are

| cycle | snapshot tail | cumulative first passage |
| ---: | ---: | ---: |
| 1 | 0 | 0 |
| 2 | 0 | 0 |
| 3 | 0.037182 | 0.043053 |
| 4 | 0.266145 | 0.491194 |
| 5 | 0.256360 | 0.704501 |
| 6 | 0.146771 | 0.743640 |

Thus the conservative driven chain supports a history-bearing cumulative threshold-crossing record even though the instantaneous tail can relax.

No spacing reaches the active normal instability threshold lambda_c = 1.1077715385524567 in this six-cycle protocol.

## 3. External work and energy

The conservative chain obeys

$$
\frac{dE_{\rm mech}^*}{d\tau}=q(\tau)\dot x_M.
$$

Across the audit the largest absolute difference between cumulative external work and mechanical-energy change is approximately

$$
9.15\times10^{-12}.
$$

Hence the tail/first-passage behavior is consistent with conservative energy accounting and must not be interpreted as a proof of positive irreversible dissipation.

## 4. Interpretation for the fatigue model

The numerical result supports the following distinction:

1. `P_lambda(lambda,t)` is a mechanically generated instantaneous state.
2. `Q_*(t)` is an instantaneous high-spacing tail and may increase or decrease.
3. `F_*(t)` is a cumulative first-passage observable and is nondecreasing.
4. Nonzero thermodynamic dissipation `D_irr` remains a separate open physical problem.
5. A registry/plasticity extension may map high spacing to a reduced registry barrier, but actual registry crossing and residual plasticity require their own dynamics and physical seed.

The dimensionless frequency omega*=0.02 is a proof-of-principle atomic-chain protocol, not a laboratory fatigue-frequency mapping.
