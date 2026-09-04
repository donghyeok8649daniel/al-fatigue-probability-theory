# Candidate spatial registry-kink feedback

Status: **CANDIDATE / PARTIAL RESULT — not active governing law.**

This note supersedes the earlier optimistic metastability wording. The spatial registry extension remains useful, but the latest unbounded relaxation and NEB audit does **not** establish a strongly trapped residual plastic state under ideal pure-normal loading.

Research constraints retained here:

- no FCC geometry;
- no arbitrary scalar damage variable;
- no fitted damping introduced only to create a fatigue clock;
- no characteristic area/volume calibration at this stage;
- the active normal-only paper remains unchanged.

## 1. Spatial registry field

The single collective registry coordinate is periodic,

$$
U_0(a,s+b)=U_0(a,s),
$$

so a uniform completed shift by one period is energetically equivalent. To represent local incompatibility, let $j$ index repeated units and define

$$
x_j=jb+s_j.
$$

Keep the cross-row energy $U_0(a_j,s_j)$ and the direct same-row compatibility energy

$$
E_{\parallel}
=
\sum_{j<k}
\left[
v_{m,n}\left(\left|(k-j)b+s_k-s_j\right|\right)
-v_{m,n}((k-j)b)
\right].
$$

The candidate spatial energy is

$$
E_{\mathrm{rk}}
=
\sum_jU_0(a_j,s_j)+E_{\parallel}.
$$

No phenomenological gradient coefficient is inserted.

## 2. Uniform shift versus local shift

For

$$
s_j\rightarrow s_j+b\qquad\text{for every }j,
$$

all $s_k-s_j$ are unchanged and $U_0$ is periodic, so the total intrinsic energy is unchanged.

For a local shift, write

$$
s_j=s_0+z_jb+\tilde s_j,
\qquad z_j\in\mathbb Z,
$$

and define

$$
q_j^{\mathrm{k}}=z_{j+1}-z_j.
$$

A finite shifted patch contains a kink and an antikink. Their cores contain intermediate registry values and therefore are not equivalent to either far-field well.

## 3. Mechanical coupling back to normal spacing

The local intrinsic normal generalized force is

$$
Q_{a,j}^{\mathrm{int}}
=-\partial_aU_0(a_j,s_j).
$$

At two equivalent well centers,

$$
\partial_aU_0(a,s_0+b)=\partial_aU_0(a,s_0),
$$

but inside a kink core there is no such equality. Thus a kink core can temporarily change the locally preferred normal spacing without inventing a feedback law.

The mechanically allowed coupling chain remains

$$
P_a
\rightarrow
\text{registry rare event}
\rightarrow
\{s_j\}_{\mathrm{core}}
\rightarrow
\partial_aU_0(a,s_j)
\rightarrow
\text{local normal-spacing perturbation}.
$$

What is no longer justified is automatically appending a permanent $P_a^{\mathrm{next}}$ to that chain.

## 4. Self-consistent row-spacing correction

The newer audit no longer sets $b=\sigma_{LJ}$ by convenience. For an infinite same-row chain at its equilibrium spacing,

$$
\frac{\sigma_{LJ}}{b}
=
\left[
\frac{n\zeta(n)}{m\zeta(m)}
\right]^{1/(m-n)}.
$$

For

$$
m=12.19,\qquad n=6,
$$

the current value is

$$
\frac{\sigma_{LJ}}{b}\approx0.8942468263.
$$

On this corrected surface the zero-normal-traction registry-well normal equilibrium is

$$
a_{\mathrm{well}}/b\approx0.8582179181.
$$

This replaces the older $b=\sigma_{LJ}$ numerical diagnostic in this candidate route.

## 5. Correction to the first metastability claim

The earlier spatial audit used L-BFGS-B with registry bounds. That calculation retained a narrow nonuniform profile and classified it as metastable.

Repeating the calculation without those artificial registry bounds changes the result:

- initial kink-pair patches of width 16, 18 and 20 repeats relax back to the intact state;
- sufficiently wider patches, around 21 repeats and above in the present 121-site system, can remain in lattice-pinned local minima;
- the smallest Hessian eigenvalues of those retained states are only of order $10^{-4}$ in the present dimensionless units.

Therefore the correct interpretation is

$$
\boxed{
\text{wide kink pairs can be lattice-pinned, but the trapping mode is extremely shallow.}
}
$$

The older statement that the first relaxed narrow pair itself demonstrated robust metastability is withdrawn.

## 6. MEP / NEB audit

The first robust wide-pair state used in the current zero-stress MEP audit has a formation energy relative to the relaxed intact state of approximately

$$
0.923785\ \mathrm{eV}
$$

under the retained diagnostic physical energy bridge.

A 27-image FIRE/CI-NEB band does not resolve a clean interior climbing-image saddle above this very shallow product state. The band approaches the product almost monotonically at the available resolution.

A dense scan of the final NEB segment finds only a very small additional corrugation of order

$$
10^{-5}\ \mathrm{eV}.
$$

That micro-eV-scale number is path- and resolution-sensitive, so it is **not** promoted as a converged value of $\Delta G_{\mathrm{kp}}$.

The robust energetic statement is only that the forward rare-event scale is dominated by roughly

$$
\boxed{
\Delta G_{\mathrm{kp}}^{\mathrm{forward}}\sim0.924\ \mathrm{eV}
}
$$

for the present candidate zero-stress surface, while the post-formation separation/migration trapping is extremely weak.

## 7. Normal-tension sensitivity

The dominant pair-formation cost decreases under normal tensile loading in the current candidate calculation:

| normal stress | pair-formation cost |
|---:|---:|
| 0 MPa | 0.923785 eV |
| 50 MPa | 0.917359 eV |
| 100 MPa | 0.910900 eV |
| 150 MPa | 0.904408 eV |
| 200 MPa | 0.897881 eV |

Thus the previously established sign survives:

$$
a\uparrow
\quad\Longrightarrow\quad
\text{registry transition becomes energetically easier}.
$$

These are candidate formation energies, not independently calibrated Al dislocation barriers.

## 8. Conditional activated-rate interpretation

If the finite-temperature activated-rate assumptions are provisionally retained, an energy scale around $0.9$ eV can indeed separate a THz attempt scale from a many-cycle local-event probability.

For the diagnostic case 300 K, 20 Hz, and $100\pm100$ MPa normal loading, using the stress-dependent formation-energy curve as if it were the rate barrier gives a one-cycle hazard of approximately

$$
3.74\times10^{-5}.
$$

This is only a mechanism-timescale diagnostic. It is not a fatigue-life prediction and does not include characteristic area/volume scaling.

## 9. Current physical verdict

The spatial registry route now has two different conclusions.

First, it gives a mechanically derived rare transient state:

$$
\boxed{
P_a\rightarrow\text{rare kink-core event}\rightarrow\text{temporary local normal opening}
}
$$

which may still be relevant to a first-passage crack-initiation calculation.

Second, it does **not** yet give the desired persistent feedback

$$
P_a(a,0)+\sigma(0:t)\rightarrow P_a(a,t)
$$

through stored kink damage, because the post-formation separation/migration mode is too weakly trapped in the ideal homogeneous pure-normal model.

Therefore this candidate is not promoted as the missing progressive $P_a$ evolution law.

A durable post-transition state would require additional independently justified mechanics, such as a real pinning defect, boundary incompatibility, or another non-equivalent structural constraint. None is inserted here.

## 10. Relation to the active theory

The active normal-only theory remains unchanged. The safe result retained from this extension is:

1. normal opening reduces a registry rare-event energy scale;
2. a spatial registry core strongly couples back to local normal spacing while it exists;
3. rare-event first passage can occur on many-cycle timescales without slow atomic inertia;
4. permanent plastic memory has **not** been established in the present ideal pure-normal model.

The next active-theory decision should therefore be made from the original requirement

$$
P_0(a)+\sigma(0:t)\longrightarrow P(a,t),
$$

not by assuming that kink storage has already solved it.
