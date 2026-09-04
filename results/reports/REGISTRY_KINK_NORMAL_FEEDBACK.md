# Registry-kink normal-feedback audit

**Classification:** dimensionless direct-sum mechanism diagnostic; candidate only.

The calculation uses the retained multilayer registry surface with

$$
m=12.19,\qquad n=6,\qquad b=\sigma_{LJ}=1,
$$

and the stable registry well at $s/b=0.5$.

The stable zero-normal-traction well root is

$$
a_{0,r}/b=0.9910707232.
$$

For each fixed registry position $s$, the stable conditional normal equilibrium is obtained from

$$
\partial_aU_0(a_{\mathrm{eq}},s)=0,
\qquad
\partial_a^2U_0(a_{\mathrm{eq}},s)>0.
$$

| $s/b$ | energy excess at $a_{0,r}$ | $\partial_aU_0(a_{0,r},s)$ | stable $a_{\mathrm{eq}}/b$ | relative opening |
|---:|---:|---:|---:|---:|
| 0.50 | 0.000000 | ~0 | 0.991071 | 0.000% |
| 0.40 | 0.128386 | -2.091306 | 1.007995 | 1.708% |
| 0.30 | 0.480899 | -7.893606 | 1.040134 | 4.951% |
| 0.25 | 0.710900 | -11.723411 | 1.055051 | 6.456% |
| 0.20 | 0.950409 | -15.745595 | 1.067692 | 7.731% |
| 0.10 | 1.358597 | -22.674485 | 1.084831 | 9.460% |
| 0.00 | 1.521402 | -25.462384 | 1.090567 | 10.039% |

## Result

The registry-saddle side of the core strongly changes the preferred normal separation. At $s/b=0$, the stable normal equilibrium on this reduced surface is about **10.04% larger** than at the registry well.

This establishes the missing *mechanical feedback possibility*:

$$
\text{spatial registry core}
\rightarrow
\partial_aU_0(a,s)\text{ changes}
\rightarrow
\text{local normal equilibrium changes}
\rightarrow
P_a\text{ can change}.
$$

The far-field wells remain exactly equivalent because $U_0(a,s+b)=U_0(a,s)$. Therefore a uniform lattice-period shift is not counted as damage; only a spatially nonuniform core can create a non-equivalent local structural state.

## What this does not prove

This audit does not yet compute a kink profile, a critical kink-pair saddle, a kink-pair activation barrier, an irreversible residual state, a characteristic area/volume, or an Al fatigue life.

The earlier coherent-patch barrier $N\Delta G_s$ is therefore not promoted. The next calculation is the minimum-energy spatial kink and the critical kink-pair saddle from the same pair-interaction energy.
