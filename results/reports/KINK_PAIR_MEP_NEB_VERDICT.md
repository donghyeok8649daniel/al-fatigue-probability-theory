# Kink-pair MEP / NEB verdict

Status: **candidate mechanism audit only**.

## Question

Can the spatial registry field produce a physically non-equivalent post-transition state with a finite activation barrier, without inserting an empirical damage variable or an artificial slow inertia?

## Important correction

The previous spatial-registry note treated a narrow relaxed kink-pair profile as metastable. Repeating the minimization **without the artificial registry bounds used by the earlier L-BFGS-B diagnostic** changes that conclusion.

- initial patch widths 16, 18 and 20 repeats return to the intact basin;
- wider patches, beginning around 21 repeats in the present 121-site audit, can remain in extremely shallow lattice-pinned kink-pair minima;
- the lowest Hessian curvature of those retained states is of order $10^{-4}$ in the present dimensionless units, i.e. much softer than the ordinary registry/core modes.

Therefore the existence of a nonuniform residual configuration is real only in a much weaker sense than previously stated: it is a **very shallow Peierls-pinned state**, not a strongly trapped plastic product.

## Zero-stress energy scale

For the current self-consistent row spacing choice

$$
\frac{\sigma_{LJ}}{b}=0.8942468263,
$$

and the retained normal-stiffness energy bridge, the first robust wide-pair product used in the MEP audit has a formation energy

$$
E_{\rm pair}-E_{\rm intact}\approx 0.923785\ {\rm eV}.
$$

A 27-image FIRE/CI-NEB relaxation did **not** resolve a stable interior climbing image above that product state. The band approaches the shallow product nearly monotonically. A dense scan of the last band segment finds only a very small additional corrugation, about

$$
\sim1.8\times10^{-5}\ {\rm eV},
$$

but this micro-eV-scale value is sensitive to the path and interpolation resolution and is **not promoted as a converged kink-pair saddle barrier**.

The robust statement is therefore

$$
\boxed{
\Delta G_{kp}^{\rm forward}\ \text{is dominated by an energy scale of about }0.924\ {\rm eV}
}
$$

within this candidate surface, while the final Peierls trapping corrugation is extremely small.

## Normal-tension sensitivity

Using the same row model and relaxing the normal coordinate quasistatically, the retained pair-formation cost decreases as normal tensile stress increases:

| normal stress | formation cost |
|---:|---:|
| 0 MPa | 0.923785 eV |
| 50 MPa | 0.917359 eV |
| 100 MPa | 0.910900 eV |
| 150 MPa | 0.904408 eV |
| 200 MPa | 0.897881 eV |

The sign is consistent with the earlier observation that normal opening lowers registry corrugation.

## Conditional cycle-rate diagnostic

If, and only if, the finite-temperature activated-rate assumptions from the previous candidate note are retained, and if the dominant formation energy is provisionally used as the rate barrier, then for 300 K, 20 Hz and $100\pm100$ MPa normal loading the current diagnostic gives

$$
H_{\rm cycle}\approx3.74\times10^{-5},
$$

so the local transition probability per cycle is also approximately

$$
3.74\times10^{-5}.
$$

The corresponding conditional median is about

$$
1.85\times10^4\ \text{cycles}.
$$

This is **not a fatigue-life prediction**. No characteristic area/volume, independent-domain count, experimental calibration or validated transition-state prefactor has been used.

## Main physical verdict

The calculation separates two issues:

1. **Rare formation is plausible.** A $\sim0.9$ eV collective registry event can naturally separate a THz attempt scale from a many-cycle local-event probability.
2. **Persistent plastic memory is not established.** The retained pair has an extremely soft separation/migration mode and the NEB audit does not reveal a substantial reverse trapping barrier under ideal pure-normal loading.

Thus the spatial registry candidate currently supports

$$
P_a\to\text{rare kink-core first passage}\to\text{temporary local normal opening}
$$

better than it supports

$$
P_a\to\text{permanent kink storage}\to P_a^{\rm next}.
$$

The latter should **not** be inserted into the active theory without additional physical structure such as a real pinning defect, boundary incompatibility, or another independently justified non-equivalent state.

## Consequence for the main research direction

The active requirement remains

$$
P_0(a)+\sigma(0:t)\longrightarrow P(a,t)
$$

without resolving every atom.

This kink audit does not yet close that progressive evolution law. It does, however, provide a potentially useful mechanically derived **rare transient state** that can couple to the normal spacing and therefore may enter a first-passage crack-initiation calculation.

The active normal-only paper should remain unchanged until this distinction is resolved.
