# Local correlation hierarchy and the minimal reduced state

## Purpose

The current exact deterministic map is

$$
\Gamma_0+\sigma(0:t)\longrightarrow P_a(a,t).
$$

If every microscopic coordinate in `Gamma_0` must be propagated, the method is a reduced description of lattice dynamics rather than an autonomous probability theory. The present question is therefore:

> What is the smallest reduced initial state that can determine the later one-point spacing probability under the same 1D nearest-neighbor generalized-LJ mechanics?

This note checks whether the nearest-neighbor pair distribution is sufficient.

## 1. One-point phase-space measure

Use normalized spacing and spacing rate

$$
z_i=(\lambda_i,c_i),
\qquad
c_i=\dot\lambda_i.
$$

The one-point empirical phase-space measure is

$$
F_1(z,\tau)
=\frac1M\sum_{i=1}^{M}\delta[z-z_i(\tau)].
$$

Its spacing marginal is

$$
P_\lambda(\lambda,\tau)=\int F_1(\lambda,c,\tau)\,dc.
$$

For an interior spacing,

$$
\ddot\lambda_i
=\phi'(\lambda_{i+1})
-2\phi'(\lambda_i)
+\phi'(\lambda_{i-1}).
$$

Therefore the one-point acceleration flux depends on neighboring spacings.

## 2. Directed nearest-neighbor pair measure

Define the directed adjacent-pair empirical measure

$$
F_2^+(z,z_+,\tau)
=\frac1{M-1}\sum_{i=1}^{M-1}
\delta[z-z_i(\tau)]
\delta[z_+-z_{i+1}(\tau)].
$$

A corresponding left-neighbor measure may be written as

$$
F_2^-(z,z_-,\tau)
=\frac1{M-1}\sum_{i=2}^{M}
\delta[z-z_i(\tau)]
\delta[z_--z_{i-1}(\tau)].
$$

For the one-point mean acceleration, these pair marginals are enough at one instant because the acceleration is additive in the left and right forces. In schematic form,

$$
P\mathcal A
=\int \phi'(\lambda_+)F_2^+\,dz_+
+\int \phi'(\lambda_-)F_2^-\,dz_-
-2\phi'(\lambda)P,
$$

with separate boundary terms for the first and last spacing.

Thus pair information improves the one-point description substantially.

## 3. Why the pair state is not autonomous

The exact pair transport equation has the form

$$
\partial_\tau F_2
+\partial_{\lambda_1}(c_1F_2)
+\partial_{\lambda_2}(c_2F_2)
+\partial_{c_1}\mathcal G_{2,1}
+\partial_{c_2}\mathcal G_{2,2}
=0.
$$

For the pair `(i,i+1)`, however,

$$
\ddot\lambda_i
=\phi'(\lambda_{i+1})
-2\phi'(\lambda_i)
+\phi'(\lambda_{i-1}),
$$

and

$$
\ddot\lambda_{i+1}
=\phi'(\lambda_{i+2})
-2\phi'(\lambda_{i+1})
+\phi'(\lambda_i).
$$

The pair state knows `lambda_i` and `lambda_{i+1}` but does not know the outside spacings `lambda_{i-1}` and `lambda_{i+2}`. Consequently the pair acceleration flux requires contiguous triplet statistics.

Symbolically,

$$
F_1\leftarrow F_2\leftarrow F_3.
$$

More generally, for a contiguous `k`-block empirical measure

$$
F_k(z_1,\ldots,z_k,\tau),
$$

the accelerations of the two edge members of the block depend on one state outside the block. Therefore

$$
\boxed{F_k\text{ transport generally requires }F_{k+1}.}
$$

Nearest-neighbor mechanics makes this hierarchy local, but it does not make the pair level exactly autonomous.

For an arbitrary finite chain, the exact finite hierarchy closes only when the full ordered microscopic state has been retained, or when an additional exact structural restriction is proven for the class of initial states being studied.

## 4. Explicit pair-insufficiency counterexample

Consider six normalized spacings with zero initial rates and zero external force. Let

$$
A=0.99,\qquad B=1.01,\qquad C=1.02.
$$

Construct

$$
\Lambda^{(A)}_0=(A,A,A,B,A,C),
$$

and

$$
\Lambda^{(B)}_0=(A,A,B,A,A,C).
$$

The two chains have exactly the same:

- one-point spacing multiset, hence the same `P_0`;
- one-point spacing/rate multiset, hence the same `F_0` because every rate is zero;
- directed nearest-neighbor pair multiset;
- first spacing and last spacing, so the fixed-left / force-right boundary labels also agree.

Their directed adjacent pairs are, in either case,

- `(A,A)` twice;
- `(A,B)` once;
- `(B,A)` once;
- `(A,C)` once.

Therefore their empirical directed pair distributions are identical at the initial time.

Their triplet multisets are not identical. The first chain contains `(A,A,A)` whereas the second contains different three-site contexts. Hence the pair state has discarded information that the future pair transport needs.

Using the active normalized generalized-LJ exponents

$$
m=12.19,\qquad n=6,
$$

and integrating both chains with the same deterministic fixed-left, zero-force-right dynamics gives different later one-point spacing distributions. At

$$
\tau=1,
$$

the empirical KS distance between the two spacing distributions is

$$
\boxed{D_{KS}=\frac13\approx0.333333.}
$$

The maximum difference between their sorted spacing samples is approximately

$$
\boxed{5.47\times10^{-3}.}
$$

The result is deterministic. No random force, PDF assumption, or stochastic closure is used.

## 5. Important interpretation

The counterexample establishes

$$
\boxed{(P_0,F_0,F_{2,0})+\sigma(0:t)\not\Rightarrow P(t)}
$$

for arbitrary initial ordered chains.

This does **not** mean reduced probability dynamics is impossible. It means that an exact, memoryless, finite-order closure cannot be declared merely from the nearest-neighbor form of the microscopic force.

There are now two honest routes:

1. **Exact hierarchy route:** retain enough local correlation order for the task, with the full ordered state as the finite exact endpoint.
2. **Reduced predictive route:** determine empirically and mathematically the smallest local block order `k` that predicts the observables of interest to a controlled error tolerance, without claiming exact closure.

The second route is the useful next research target because the objective is not to reproduce every microscopic coordinate. The target observables are the one-point probability, tail flux, and first-passage quantities used by G1, G2, and G4.

## 6. Next validation target

For

$$
F_k(z_i,\ldots,z_{i+k-1},t),
$$
construct a `k`-local conditional-acceleration predictor from exact chain trajectories and test its prediction error against the full chain for

$$
P(a,t),\qquad Q_*(t),\qquad F_*(t),\qquad F_{ci}(t).
$$

Sweep `k = 1,2,3,...` and report the smallest `k` that meets a predeclared error tolerance over the loading protocols of interest.

No neighbor-independence, Gaussian closure, Boltzmann law, or arbitrary Markov assumption is to be inserted merely to force a low-order closure.
