# Candidate spatial registry-kink feedback

Status: **CANDIDATE / PROMISING EXTENSION — not active governing law.**

This note keeps the current research constraints:

- no FCC geometry;
- no arbitrary damage variable;
- no fitted damping added only to create a fatigue clock;
- no characteristic area/volume calibration at this stage;
- the active normal-only paper remains unchanged.

## 1. Why the single collective registry coordinate was insufficient

The retained registry surface is periodic,

$$
U_0(a,s+b)=U_0(a,s).
$$

Therefore a uniform completed shift from $s_0$ to $s_0+b$ lands in an energetically equivalent state. A single collective coordinate $s$ can therefore provide a barrier and a first-passage event, but it cannot by itself create a non-equivalent residual structure that changes the later normal-spacing distribution.

The missing information is spatial incompatibility: a local part of a row can slip while its neighbors have not yet slipped.

## 2. Promote registry from one number to a discrete local field

Let $j$ index repeated units along the registry direction and define the moving-row coordinate

$$
x_j=jb+s_j.
$$

The cross-row energy of repeat $j$ is still the already derived periodic surface

$$
U_0(a_j,s_j).
$$

The same moving repeats must also remain compatible with one another along the row. Without inserting a phenomenological gradient coefficient, retain their direct same-row pair interaction:

$$
E_{\parallel}
=
\sum_{j<k}
\left[
v_{m,n}
\left(
\left|(k-j)b+s_k-s_j\right|
\right)
-v_{m,n}((k-j)b)
\right].
$$

Hence the candidate intrinsic energy is

$$
E_{\mathrm{rk}}
=
\sum_jU_0(a_j,s_j)+E_{\parallel}.
$$

This is a reduced-row extension of the existing interaction mechanics. It is not a 3D crystal model.

## 3. Uniform slip remains exactly equivalent

For a uniform lattice-period shift,

$$
s_j\rightarrow s_j+b\qquad\text{for every }j,
$$

all differences $s_k-s_j$ are unchanged, while $U_0$ is periodic. Therefore

$$
E_{\mathrm{rk}}[\{s_j+b\}]
=
E_{\mathrm{rk}}[\{s_j\}].
$$

This is required: a global relabeling by one perfect registry period must not be falsely counted as damage.

## 4. A local slip is not equivalent

Write

$$
s_j=s_0+z_jb+\tilde s_j,
\qquad z_j\in\mathbb Z.
$$

If $z_j$ changes only in a finite region, then neighboring repeats at the region boundary have different registry histories. Their same-row separations contain $s_{j+1}-s_j$, so $E_{\parallel}$ changes.

Define the discrete topological/kink indicator

$$
q_j^{\mathrm{k}}=z_{j+1}-z_j.
$$

A finite shifted patch inside an unshifted row contains two boundaries: a kink and an antikink. An abrupt one-period step is generally not a low-energy state because it strongly distorts same-row separations. The mechanically allowed transition therefore spreads over a finite core in which some $s_j$ lie between neighboring wells.

This is the first non-arbitrary post-transition structural state found in the current reduced mechanics.

## 5. Why this can feed back into the normal spacing distribution

The local intrinsic normal generalized force is

$$
Q_{a,j}^{\mathrm{int}}
=-\partial_aU_0(a_j,s_j).
$$

At the exact well centers,

$$
\partial_aU_0(a,s_0+b)
=
\partial_aU_0(a,s_0),
$$

so far from a kink the adjacent wells remain normal-force equivalent.

Inside a kink core, however, $s_j$ is between the well centers. There is no symmetry requiring

$$
\partial_aU_0(a,s_j)
=
\partial_aU_0(a,s_0).
$$

Therefore a spatial registry core can change the locally preferred normal spacing even though the two far-field registry wells are equivalent.

The candidate feedback chain is thus

$$
P_a
\rightarrow
\Delta G_s(a)
\rightarrow
\text{kink-pair nucleation}
\rightarrow
\{s_j\}_{\mathrm{core}}
\rightarrow
\partial_aU_0(a,s_j)
\rightarrow
P_a^{\mathrm{next}}.
$$

No scalar empirical damage variable is introduced in this chain.

## 6. Direct-sum diagnostic of normal-registry coupling

Using the current registry diagnostic surface

$$
m=12.19,\qquad n=6,\qquad b=\sigma_{LJ}=1,
$$

with the stable registry well at $s_0/b=0.5$, the converged direct sum gives the zero-normal-traction registry-well equilibrium

$$
a_{0,r}/b\approx0.9910707.
$$

For fixed registry position $s$, define the stable conditional normal equilibrium by

$$
\partial_aU_0(a_{\mathrm{eq}},s)=0,
\qquad
\partial_a^2U_0(a_{\mathrm{eq}},s)>0.
$$

The diagnostic values are approximately:

| $s/b$ | $[U_0(a_{0,r},s)-U_0(a_{0,r},s_0)]$ | $a_{\mathrm{eq}}/b$ | opening relative to $a_{0,r}$ |
|---:|---:|---:|---:|
| 0.50 | 0 | 0.991071 | 0 |
| 0.40 | 0.128386 | 1.007995 | 1.71% |
| 0.30 | 0.480899 | 1.040134 | 4.95% |
| 0.25 | 0.710900 | 1.055051 | 6.46% |
| 0.20 | 0.950409 | 1.067692 | 7.73% |
| 0.10 | 1.358597 | 1.084831 | 9.46% |
| 0.00 | 1.521402 | 1.090567 | 10.04% |

Thus the registry saddle region strongly prefers a larger normal separation on this reduced surface. The saddle value corresponds to roughly a 10% normal opening relative to the registry-well equilibrium.

This is a dimensionless mechanism diagnostic. It must not yet be identified quantitatively with an Al dislocation core or with the active normal-chain crack threshold because the normal-only chain and the registry surface are not yet one independently calibrated atomistic potential.

Nevertheless the sign and strength of the coupling are important: the missing post-transition state can, in principle, perturb $P_a$ through the already existing $U_0(a,s)$ rather than through an invented feedback law.

## 7. Correction to the earlier coherent-patch Arrhenius diagnostic

The previous conditional rate used

$$
\Delta G_{s,N}=N\Delta G_s
$$

for an entire patch translating coherently. That remains a useful sensitivity toy, but it is no longer the preferred physical transition path once $s_j$ is spatially resolved.

A local extended system can instead cross by nucleating a kink-antikink pair and then moving those cores. The physically relevant activation barrier is therefore

$$
\Delta G_{\mathrm{kp}}
=
E_{\mathrm{rk}}[\text{critical kink-pair saddle}]
-
E_{\mathrm{rk}}[\text{initial well state}],
$$

not an assumed $N\Delta G_s$.

The next rate calculation must use $\Delta G_{\mathrm{kp}}$ obtained from the spatial energy landscape. No value is adopted here.

## 8. Relation to known mechanics

This structure is closely related to the physical content of Frenkel-Kontorova and Peierls-Nabarro descriptions: a periodic misfit energy competes with elastic compatibility, producing finite-width disregistry cores and kink-type defects. This literature supports the mechanism class, but does not validate the numerical parameters of the present reduced Al model.

Useful references:

- G. Schoeck, *The Peierls model: Progress and limitations*, Materials Science and Engineering A 400-401 (2005) 7-17, DOI: 10.1016/j.msea.2005.03.050.
- A. P. Sutton, *Multiscale models of dislocations*, in Physics of Elasticity and Crystal Defects, Oxford University Press (2024), DOI: 10.1093/oso/9780198908081.003.0007.
- S. P. Fitzgerald, *Kink pair production and dislocation motion*, Scientific Reports 6, 39708 (2016), DOI: 10.1038/srep39708.

## 9. What is still not solved

This candidate does **not** yet prove irreversible fatigue evolution. The remaining tasks are:

1. compute the minimum-energy kink and critical kink-pair saddle from $E_{\mathrm{rk}}$;
2. verify that the core is metastable or long-lived after unloading when a finite-temperature relaxation model is declared;
3. couple $a_j$ and $s_j$ consistently and test whether a nucleated core produces a persistent measurable change in $P_a$;
4. only after that derive an activated rate and cycle-to-cycle $P_a$ propagator;
5. leave characteristic area/volume calibration for later experimental scaling.

The candidate should be rejected if the spatially resolved minimum-energy transition relaxes immediately back with no residual structural state and no persistent change in the normal-spacing marginal.

## 10. Current verdict

The single collective-registry model could provide a barrier but not a non-equivalent final state. Spatially resolving registry produces a mechanically justified candidate for that missing state: a kink/antikink core created by local incompatibility.

The key new point is

$$
\text{uniform registry shift is equivalent, but a local registry shift is not.}
$$

Because the kink core occupies intermediate registry values where $\partial_aU_0$ differs strongly from the well value, this route can potentially close the missing feedback from a rare plastic event back into $P_a(a,t)$.

It remains a candidate until the spatial saddle, residual state, and $P_a$ feedback are solved numerically.
