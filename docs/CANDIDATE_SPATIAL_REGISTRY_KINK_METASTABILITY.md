# Candidate spatial registry-kink metastability

Status: **CANDIDATE / STRONGER MECHANISM RESULT — not active paper law.**

This note continues `CANDIDATE_SPATIAL_REGISTRY_KINK_FEEDBACK.md` and tests whether the spatial registry extension can actually possess a non-equivalent residual state after unloading.

## 1. Make the row repeat self-consistent with the same pair potential

Once the generalized LJ interaction is used both across rows and for same-row compatibility, the row repeat $b$ should not be chosen independently of the LJ length scale.

For an infinite one-dimensional row, the pair energy per repeat is

$$
U_{\parallel}(b)
=C_{mn}\epsilon
\left[
(\sigma_{LJ}/b)^m\zeta(m)
-(\sigma_{LJ}/b)^n\zeta(n)
\right].
$$

The stationarity condition $dU_{\parallel}/db=0$ gives

$$
\frac{\sigma_{LJ}}{b}
=
\left[
\frac{n\zeta(n)}{m\zeta(m)}
\right]^{1/(m-n)}.
$$

For $m=12.19$, $n=6$,

$$
\frac{\sigma_{LJ}}{b}=0.8942468263,
\qquad
\frac{b}{\sigma_{LJ}}=1.1182594901.
$$

This removes one arbitrary geometry choice from the spatial registry candidate.

## 2. Eliminate the fast normal coordinate quasistatically

The previous laboratory-timescale audits show that atomic normal motion is extremely fast compared with ordinary fatigue loading. Therefore, for a prescribed normal generalized force $Q_a$, define the stable relaxed registry substrate

$$
V_{\mathrm{eff}}(s;Q_a)
=
\min_{a\ \mathrm{stable}}
\left[U_0(a,s)-Q_a a\right].
$$

The spatial registry energy becomes

$$
E_{\mathrm{rk}}[\{s_j\};Q_a]
=
\sum_jV_{\mathrm{eff}}(s_j;Q_a)
+
E_{\parallel}[\{s_j\}].
$$

No damping coefficient, mobility, gradient coefficient, or empirical damage state is added.

## 3. Zero-stress metastability test

A 121-repeat row is initialized with a finite region close to the neighboring registry well. Both ends are fixed in the original well $s/b=0.5$, but **no central repeat and no slipped fraction is constrained during final minimization**.

If the spatial defect is not a true residual structural state, energy minimization should collapse the whole row back to the uniform well.

It does not.

The minimization converges to a nonuniform kink-antikink / finite slipped-patch state with

$$
\Delta E_{\mathrm{kp}}^{\mathrm{form}}
=11.38755235
$$

in the candidate dimensionless energy units.

The fixed-endpoint Hessian at that state has lowest eigenvalue approximately

$$
\lambda_{\min}\approx1.26\times10^{-4}>0,
$$

while the next eigenvalues are of order $10^1$. The tiny lowest mode is consistent with a very soft defect-position mode; no negative mode was found in this local-minimum audit.

Therefore, within the declared candidate energy,

$$
\text{uniform well}
\quad\text{and}\quad
\text{finite kink-pair state}
$$

are distinct local minima.

This is the first explicit residual post-transition state found in the project without inserting a phenomenological damage variable.

## 4. Diagnostic physical energy bridge

Matching the candidate registry-surface normal curvature to the retained normal stiffness $EA_0/a_0$ gives, for this self-consistent row-repeat choice,

$$
E_{\mathrm{unit}}\approx0.0811219\ \mathrm{eV}.
$$

Hence the metastable defect formation/local-minimum energy is approximately

$$
\Delta E_{\mathrm{kp}}^{\mathrm{form}}
\approx0.92378\ \mathrm{eV}.
$$

This is **not** the kink-pair activation barrier. The transition-state saddle can be higher and must be calculated separately before any rate/lifetime claim.

## 5. The residual core feeds directly back into normal spacing

For every optimized registry value $s_j$, recover the stable quasistatic normal spacing

$$
a_j=a_{\mathrm{eq}}(s_j;Q_a).
$$

At zero normal stress the metastable core produces

$$
\max_j\frac{a_j}{a_{\mathrm{well}}}
\approx1.143795.
$$

Thus the core contains local states whose preferred normal spacing is about **14.38% larger** than the far-field registry-well spacing on this candidate surface.

This gives the required feedback structurally rather than empirically:

$$
\text{rare registry transition}
\rightarrow
\text{metastable kink-pair state}
\rightarrow
\text{shifted }a_{\mathrm{eq}}(s_j)
\rightarrow
\text{changed local }P_a.
$$

The active normal-chain crack threshold must not yet be numerically identified with this 14.38% value because the active normal chain and this candidate registry-row energy have not been independently reconciled into one calibrated potential. The important result here is the existence and sign of the feedback.

## 6. Normal tensile stress lowers the residual-state energy

Using only the retained normal-force bridge as a diagnostic, the same local-minimum branch gives approximately:

| normal stress | $\Delta E_{\mathrm{kp}}^{\mathrm{form}}$ |
|---:|---:|
| 0 MPa | 0.92378 eV |
| 50 MPa | 0.91735 eV |
| 100 MPa | 0.91089 eV |
| 150 MPa | 0.90440 eV |
| 200 MPa | 0.89788 eV |

The direction is physically consistent with the normal-registry coupling already found: tensile opening makes the registry-core state energetically less costly.

Again, these are local-minimum energies, not saddle barriers.

## 7. Consequence for the old coherent-patch rate

The old trial

$$
\Delta G_{s,N}=N\Delta G_s
$$

is no longer the preferred barrier model. The spatial candidate now possesses an explicit residual defect state, so the correct next quantity is the transition-state barrier

$$
\Delta G_{\mathrm{kp}}(Q_a)
=
E_{\mathrm{rk}}[\text{critical kink-pair saddle};Q_a]
-
E_{\mathrm{rk}}[\text{uniform well};Q_a].
$$

Only after that saddle is known is an activated rate physically interpretable.

## 8. What this changes in the research status

Before this audit, the project had a rare-event clock candidate but no non-equivalent final state.

After this audit, the candidate chain is

$$
P_0
\rightarrow
P_a(t)
\rightarrow
\Delta G_s(a)
\rightarrow
\text{spatial registry transition}
\rightarrow
\text{metastable kink-pair state}
\rightarrow
P_a^{\mathrm{residual}}.
$$

So the missing **post-transition structural feedback** is no longer merely hypothetical within the candidate model.

The remaining hard step is the transition pathway/rate, not the existence of a residual state.

## 9. Remaining falsification tests

This candidate must still be rejected or revised if any of the following fail:

1. a minimum-energy-path calculation finds no physically accessible kink-pair saddle under the intended tensile stress range;
2. full coupled $a_j,s_j$ dynamics destroys the metastable state that appears under quasistatic normal elimination;
3. the required activation barrier remains so large that experimental fatigue-cycle probabilities are negligible even near the intended stress range;
4. reconciliation with the active normal interaction removes the large normal-opening feedback;
5. experimental single-crystal Al behavior contradicts the predicted stress/temperature trends.

Characteristic area/volume and specimen multiplicity remain intentionally deferred to later calibration.

## 10. Current verdict

The spatial registry route survives this test.

Within a self-consistent reduced pair-potential embedding, a finite slipped patch can remain as a metastable kink-antikink state after the external normal stress is removed, and that residual core changes the locally preferred normal spacing strongly.

Therefore this route currently provides all three pieces that the single-$s$ model lacked:

- a non-equivalent residual structural state;
- a mechanical feedback into normal spacing;
- compatibility with a rare activated event rather than an artificially slow natural frequency.

It is still a candidate until $\Delta G_{\mathrm{kp}}$ and the cycle-resolved probability propagator are computed.
