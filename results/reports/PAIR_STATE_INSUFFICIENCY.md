# Pair-state insufficiency audit

## Question

Does the nearest-neighbor form of the 1D generalized-LJ force imply that the reduced state can close exactly at the nearest-neighbor pair distribution?

## Result

No, not for arbitrary ordered initial chains.

For an interior spacing,

$$
\ddot\lambda_i
=\phi'(\lambda_{i+1})-2\phi'(\lambda_i)+\phi'(\lambda_{i-1}).
$$

The one-point acceleration source can be written from left/right pair marginals at one instant. However, the transport of the pair `(i,i+1)` requires both `lambda_{i-1}` and `lambda_{i+2}`. Hence pair transport requires triplet information.

The exact local hierarchy is therefore

$$
F_1\leftarrow F_2\leftarrow F_3\leftarrow\cdots.
$$

Nearest-neighbor interaction makes the hierarchy local but does not produce an exact pair closure.

## Counterexample

Use zero initial spacing rates and zero external force with

$$
A=0.99,\qquad B=1.01,\qquad C=1.02.
$$

The two initial ordered spacing states are

$$
(A,A,A,B,A,C)
$$

and

$$
(A,A,B,A,A,C).
$$

They have the same one-point spacing measure `P0`, the same one-point phase-space measure `F0`, the same first and last spacing, and exactly the same directed adjacent-pair counts:

- `(A,A)`: 2
- `(A,B)`: 1
- `(B,A)`: 1
- `(A,C)`: 1

Their triplet distributions differ.

Both states were propagated with the active generalized-LJ exponents

$$
m=12.19,\qquad n=6
$$

using the fixed-left / zero-force-right deterministic chain, velocity Verlet, and

$$
\Delta\tau=10^{-3}.
$$

At `tau = 1`, the later one-point empirical spacing distributions have

$$
D_{KS}=0.3333333333
$$

and the maximum difference between sorted spacing samples is

$$
5.4660\times10^{-3}.
$$

The spacing variances at that time are approximately

$$
3.7715\times10^{-5}
$$

and

$$
1.4192\times10^{-5}.
$$

Thus

$$
(P_0,F_0,F_{2,0})+\sigma(0:t)
\not\Rightarrow P(t)
$$

in general.

## Interpretation

This is not a failure of the mechanically generated probability idea. It identifies the exact point at which microscopic ordering enters the reduced dynamics.

The useful next target is not to assume pair independence. Instead, test a sequence of local block states `F_k` against the full deterministic chain and find the smallest `k` that predicts the observables actually needed by the fatigue theory to a predeclared tolerance:

- `P(a,t)`;
- snapshot tail mass;
- cumulative first passage;
- actual `lambda_c` first passage;
- G1 and G2 moments.

If a small `k` works over the loading regime of interest, that gives a controlled reduced model even though it is not an exact universal closure.

Reproducer: `simulations/verify_pair_state_insufficiency.py`.
