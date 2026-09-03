# Initial microscopic state to P0 and exact initial-data sufficiency

## Purpose

This note fixes the initial-state side of the active 1D probability theory. It answers two separate questions:

1. Given an initial microscopic state, how is the initial spacing probability `P0` constructed?
2. Is `P0` by itself sufficient to determine the later probability state under a prescribed stress history?

The answer is: the first map is exact and immediate; the second is not unique in general because marginal probability data discard velocities and spatial/neighbor ordering.

## 1. Active finite-chain initial state

Use the normalized spacing coordinate

$$
\lambda_i(0)=x_i(0)-x_{i-1}(0),
$$

and spacing rate

$$
c_i(0)=\dot\lambda_i(0)=\dot x_i(0)-\dot x_{i-1}(0).
$$

The left node is fixed:

$$
x_0=0,\qquad \dot x_0=0.
$$

The complete deterministic initial state may therefore be represented either by the ordered node state

$$
\Gamma_0=(x_1,\ldots,x_M,\dot x_1,\ldots,\dot x_M),
$$

or equivalently by the ordered spacing state

$$
\Lambda_0=(\lambda_1,\ldots,\lambda_M,c_1,\ldots,c_M).
$$

The equivalence follows from

$$
x_i(0)=\sum_{k=1}^{i}\lambda_k(0),
$$

$$
\dot x_i(0)=\sum_{k=1}^{i}c_k(0).
$$

Thus the ordered spacing and spacing-rate arrays are sufficient initial data for the deterministic finite chain.

## 2. Exact construction of P0 and F0

For one deterministic finite chain, the initial empirical phase-space measure is

$$
F_0(\lambda,c)
=\frac1M\sum_{i=1}^{M}
\delta[\lambda-\lambda_i(0)]
\delta[c-c_i(0)].
$$

The initial spacing marginal is

$$
P_{\lambda,0}(\lambda)
=\int F_0(\lambda,c)\,dc
=\frac1M\sum_{i=1}^{M}
\delta[\lambda-\lambda_i(0)].
$$

For the physical spacing `a=a0 lambda`,

$$
P_{a,0}(a)
=\frac1{a_0}P_{\lambda,0}\left(\frac{a}{a_0}\right)
=\frac1M\sum_{i=1}^{M}\delta[a-a_i(0)].
$$

No named probability family is introduced. `P0` is a push-forward/counting measure of the actual initial microscopic state.

## 3. Ensemble initial state

If the initial microscopic state is itself uncertain, use a full-state measure `mu_0`:

$$
\int\mu_0(d\Gamma_0)=1.
$$

Then

$$
F_0(\lambda,c)
=\frac1M\sum_i\int
\delta[\lambda-\Lambda_i(\Gamma_0)]
\delta[c-C_i(\Gamma_0)]
\mu_0(d\Gamma_0),
$$

and

$$
P_{\lambda,0}(\lambda)
=\frac1M\sum_i\int
\delta[\lambda-\Lambda_i(\Gamma_0)]
\mu_0(d\Gamma_0).
$$

The deterministic case is recovered with

$$
\mu_0=\delta_{\Gamma_0^*}.
$$

## 4. Homogeneous static initial state under a prescribed initial stress

If the chain is assumed homogeneous and at rest before cyclic loading, the stable initial stretch is obtained from the current calibration bridge

$$
q_0=\frac{\sigma_0}{E}
$$

and the static force balance

$$
\phi'(\lambda_{\rm eq})=q_0.
$$

For tensile loading below the peak cohesive force there can be more than one algebraic root. The initial equilibrium branch must therefore satisfy the local stability condition

$$
\phi''(\lambda_{\rm eq})>0.
$$

With zero initial spacing rates,

$$
P_{\lambda,0}(\lambda)=\delta(\lambda-\lambda_{\rm eq}),
$$

$$
F_0(\lambda,c)=\delta(\lambda-\lambda_{\rm eq})\delta(c).
$$

For the zero-prestress reference state this reduces to

$$
\lambda_{\rm eq}=1,
$$

$$
P_{\lambda,0}(\lambda)=\delta(\lambda-1).
$$

This delta form is not a guessed PDF. It follows from the explicit assumption of a perfectly homogeneous deterministic initial chain.

## 5. Why P0 alone is not sufficient in general

`P0` retains only the multiset of initial spacings. It discards at least:

- spacing rates;
- which spacing occupies which chain location;
- neighbor correlations;
- higher spatial correlations.

The active bulk equation is

$$
\ddot\lambda_i
=\phi'(\lambda_{i+1})
-2\phi'(\lambda_i)
+\phi'(\lambda_{i-1}).
$$

Therefore two initial chains can have the same `P0` but different accelerations and different later `P`.

Hence, without additional assumptions,

$$
P_0+\sigma(0:t)\not\Longrightarrow P(t)
$$

uniquely.

## 6. Even F0(lambda,c) is not sufficient in general

The one-point phase-space empirical measure

$$
F_0(\lambda,c)
$$

retains the multiset of spacing/rate pairs but still discards their ordered spatial arrangement.

Consider two chains formed by permuting the same list of pairs `(lambda_i,c_i)`. They have exactly the same `F0`, and consequently the same `P0`, `u0`, and `Theta0`. Nevertheless the neighbor forces generally differ because the acceleration depends on adjacent spacings.

A direct finite-chain numerical audit using the active generalized-LJ exponents `m=12.19`, `n=6`, zero external force, zero initial rates, and two permutations of the same eight initial spacings confirms this non-uniqueness. The two chains start with identical `P0` and `F0`, but their spacing distributions separate immediately; at normalized time `tau=1` the empirical two-sample KS distance is `0.375` in the audit.

This gives the strict information hierarchy

$$
\Gamma_0\;\text{or ordered }\Lambda_0
\Longrightarrow F_0
\Longrightarrow P_0,
$$

but neither reverse implication holds in general.

## 7. Analytic local proof using a test observable

For any smooth test function `g`, define

$$
M_g(t)=\frac1M\sum_i g(\lambda_i(t)).
$$

Then

$$
\frac{d^2M_g}{dt^2}
=\frac1M\sum_i
\left[
 g''(\lambda_i)c_i^2
 +g'(\lambda_i)\ddot\lambda_i
\right].
$$

If all initial rates vanish,

$$
\left.\frac{d^2M_g}{dt^2}\right|_{0}
=\frac1M\sum_i
 g'(\lambda_i(0))\ddot\lambda_i(0).
$$

The right-hand side depends on neighbor ordering through `ddot lambda_i`. Thus two initial states with the same one-point `F0` can already have different second time derivatives of the same probability moment.

In the numerical counterexample, choosing

$$
g(\lambda)=\lambda^2
$$

gives

$$
M_g''(0)\approx-3.51299\times10^{-3}
$$

for one ordering and

$$
M_g''(0)\approx-1.85653\times10^{-3}
$$

for the other, despite identical `P0` and `F0`.

## 8. What initial information is sufficient?

For the deterministic finite chain used by the current theory, an exact sufficient initial specification is

$$
\{\lambda_i(0),c_i(0)\}_{i=1}^{M}
$$

with the ordering retained, together with the prescribed loading history and boundary conditions.

Equivalently, use `Gamma_0`.

For an ensemble of possible initial chains, the exact sufficient statistical object is the full labeled initial-state measure

$$
\mu_0(d\Gamma_0).
$$

A finite collection of low-order unlabeled marginals is not guaranteed to be sufficient without an additional closure or spatial-process assumption.

## 9. Reduced-data cases where P0 can be enough

`P0` can become sufficient only after extra structure makes the missing information redundant. The most important current example is the perfectly homogeneous static state:

$$
P_0=\delta(\lambda-\lambda_{\rm eq}),
\qquad c_i(0)=0\;\forall i,
$$

plus the explicit assumption that every site has the same spacing and the ordering is homogeneous.

In that restricted case the delta marginal is shorthand for the fully specified microscopic state, not an autonomous probability closure.

Other reduced closures may be investigated later, but every such closure must state the assumption that reconstructs or removes the missing spatial correlations.

## 10. Current exact prediction map

The exact finite-chain map is therefore

$$
\Gamma_0,\;\sigma_{\rm ext}(0:t)
\longrightarrow
\Gamma(t)
\longrightarrow
P_a(a,t).
$$

For an initial ensemble,

$$
\mu_0,\;\sigma_{\rm ext}(0:t)
\longrightarrow
P_a(a,t).
$$

The research problem is no longer how to guess `P0`; given a specified microscopic initial state, `P0` is mechanically defined. The remaining reduced-model question is how much of the full labeled initial information can be discarded while preserving predictive uniqueness for the observables of interest.

## 11. Next tests

1. Repeat the permutation audit over chain length and perturbation amplitude.
2. Test whether adding nearest-neighbor pair phase-space densities materially reduces prediction ambiguity.
3. Determine whether a finite Markov-order spatial representation can close the initial-data problem for the current nearest-neighbor chain over the fatigue observables G1, G2, and G4.
4. Keep the exact full-state push-forward as the reference against which every reduced initialization is tested.
