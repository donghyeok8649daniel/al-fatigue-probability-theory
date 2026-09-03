# Local-traction P0-to-P propagator audit

Status: **candidate reduced-model audit, not active finite-chain validation**.

## Question

Can the model satisfy

$$
P_0(\lambda)+q(0:\tau)\longrightarrow P(\lambda,\tau)
$$

without storing an ordered microscopic chain?

For the exact finite chain the answer is no in general because $P_0$ discards hidden neighbor ordering. This audit therefore tests a different, explicit reduction assumption: the local normal stress history is supplied as the generalized traction on a representative spacing coordinate.

## Candidate equations

For each initial spacing label $\lambda_0$,

$$
\frac{d\Lambda}{d\tau}=C,
$$

$$
\frac{dC}{d\tau}=q(\tau)-\phi'(\Lambda),
$$

with static preparation

$$
\Lambda(0)=\lambda_0,
$$

$$
C(0)=0.
$$

Then

$$
P(\lambda,\tau)
=\int P_0(\lambda_0)
\delta[\lambda-\Lambda(\tau;\lambda_0)]\,d\lambda_0.
$$

Thus the rate information needed at later times is carried internally by the characteristic flow even though the required user input is only the original $P_0$ plus the full loading history.

## Diagnostic protocol

- generalized LJ exponents: $m=12.19$, $n=6$
- $E=69$ GPa
- mean stress: 100 MPa
- stress amplitude: 100 MPa
- $q=\sigma/E$
- $\omega^*=0.02$
- 3 cycles
- $\Delta\tau=0.02$
- 5001 characteristic labels
- diagnostic initial support: $0.997\le\lambda_0\le1.003$
- initial rate: zero for every characteristic
- diagnostic tail threshold: $\lambda_*=1.004$

The bounded initial set is a numerical test input only. It is not a proposed aluminum $P_0$.

## Main checks

At almost the same applied $q$ value, the initial and half-cycle distributions differ strongly:

- KS distance: `0.5826834633`
- initial mean $\lambda$: `1.0`
- half-cycle mean $\lambda$: `1.0023869127`

The cycle-1 and cycle-2 distributions at the same forcing phase also differ:

- KS distance: `0.6926614677`

The maximum absolute per-characteristic mismatch between integrated external work and intrinsic mechanical-energy change over the saved states was

`9.19e-10`.

This supports internal consistency of the characteristic mechanics and shows that same-force non-retracing does not require storing atom indices in this candidate.

## Tail versus first passage

At $\lambda_*=1.004$, cycle-end values were:

| cycle | snapshot tail | cumulative first passage |
| ---: | ---: | ---: |
| 0 | 0 | 0 |
| 1 | 0 | 0.667866 |
| 2 | 0.279144 | 0.668066 |
| 3 | 0.183363 | 0.668066 |

The snapshot tail is non-monotone, whereas cumulative first passage is monotone. In cycle 1, a large population crossed the threshold and returned below it before the cycle-end snapshot, which is exactly why the two observables must remain distinct.

## What succeeded

The candidate gives a mathematically explicit map

$$
P_0+q(0:\tau)\to P(\tau)
$$

under the static-preparation and local prescribed-traction assumptions. It does so without assuming a named PDF family and without resolving ordered atomic positions.

## What is not established

This audit does not show that the candidate is an exact reduction of the finite nearest-neighbor chain. The existing counterexamples rule out such a universal exact $P_0$-only projection for arbitrary microscopic states.

The model also cannot create distributional spread from a perfectly degenerate initial condition:

$$
P_0(\lambda)=\delta(\lambda-1)
$$

with identical initial rates and identical local loading remains degenerate. A physical aluminum application therefore still needs a derived or measured non-degenerate $P_0$, or another physically justified source of local heterogeneity.

## Current decision

Keep this as a **candidate reduced constitutive path**. Do not replace the active exact finite-chain theory until the physical construction of $P_0$, the local inertia/stress coupling, and comparison against higher-fidelity or experimental results are completed.
