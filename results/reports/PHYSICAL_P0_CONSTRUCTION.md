# Physical P0 construction audit

Status: **candidate initialization study; not a calibrated aluminum P0 and not a promotion of the local-traction model.**

## Result

The strict P0-only candidate is internally consistent only if the initial probability is interpreted as a **slow structural/coarse-grained spacing distribution** rather than as the full instantaneous thermal spacing distribution.

Recommended construction:

$$
P_0^{\mathrm{str}}(a)
=\int_\Omega w(x)\,\delta[a-a_0^{\mathrm{str}}(x)]\,dx.
$$

Equivalently, if a residual structural microstrain map is available,

$$
\lambda_0(x)=1+\epsilon_0^{\mathrm{str}}(x),
$$

$$
P_{0,\lambda}^{\mathrm{str}}(\lambda)
=\int_\Omega w(x)\,\delta[\lambda-\lambda_0(x)]\,dx.
$$

This directly converts measured/computed specimen heterogeneity into P0 without fitting a named PDF family.

## Why the thermal instantaneous P0 is not the same input

Near lambda=1 the retained normal calibration gives

$$
K_a=\frac{EA_0}{a_{\mathrm{ref}}}=14.5431151764\ \mathrm{N/m}.
$$

If one explicitly assumes a classical canonical single harmonic spacing coordinate,

$$
\mathrm{Var}(a-a_{\mathrm{ref}})=\frac{k_BT}{K_a}.
$$

The resulting diagnostic positional standard deviations are approximately:

| T | sigma_a | sigma_lambda | naive upper tail above current lambda_c |
|---:|---:|---:|---:|
| 80 K | 8.715 pm | 0.03044 | 0.000200 |
| 293 K | 16.678 pm | 0.05826 | 0.03217 |
| 300 K | 16.876 pm | 0.05895 | 0.03376 |

The current normal-curvature threshold is

$$
\lambda_c=1.1077715386,
$$

which is only about 30.852 pm above the retained reference spacing. Therefore a naive instantaneous classical harmonic positional marginal would already place a non-negligible fraction above the current threshold at room temperature.

**This is not a crack-probability prediction.** It is a falsification warning against combining the current instantaneous thermal coordinate and current G4 threshold without accounting for thermal phase-space statistics, relative-displacement correlations, phonons/quantum effects, coarse-graining and the physical meaning of the threshold.

Moreover, a genuine thermal preparation has nonzero initial rate statistics. Hence generally

$$
F_0^{\mathrm{th}}(a,c)\ne P_0^{\mathrm{th}}(a)\delta(c),
$$

so an instantaneous thermal P0 does not satisfy the strict P0-only initialization used by the candidate local-traction propagator.

## Experimental route

The preferred experimental route is spatial residual-strain / d-spacing mapping, followed by an empirical push-forward to P0. This avoids interpreting raw diffraction peak width as a probability distribution.

Residual-strain mapping with synchrotron X-ray diffraction has demonstrated strain accuracy better than 1e-4 in materials including Al. Al single-crystal synchrotron studies also show measurable subgrain structure and residual strain. The 1e-4 scale used in the JSON summary is included only as a measurement/resolution demonstration and is **not** adopted as a universal P0 width.

If line-profile inversion is attempted instead, the profile must first be corrected/separated for instrumental broadening, coherent-domain size, overlap, mosaic/orientation spread, anisotropy and other non-strain contributions. Only under an explicit strain-only incoherent-mixture interpretation can the Bragg-law Jacobian be used to push an angular profile to a spacing distribution.

## Literature checks

- M. Kresch et al., "Phonons in aluminum at high temperatures studied by inelastic neutron scattering," Physical Review B 77, 024301 (2008), DOI 10.1103/PhysRevB.77.024301. Al phonon DOS measured from 10 to 775 K and temperature-dependent spectral changes reported.
- T. A. Herring, Microscopy 62 Suppl. 1, S99-S106 (2013). Thermal diffuse scattering analysis reported an Al atomic displacement scale of order 12 pm perpendicular to Bragg planes.
- A. M. Korsunsky et al., Journal of Synchrotron Radiation 9, 77-81 (2002), DOI 10.1107/S0909049502001905. Residual-strain mapping in materials including Al with reported strain accuracy better than 1e-4.
- T. Okada et al., "Tensile deformation and recrystallization of aluminum single crystals with sub-grained structures studied by synchrotron X-ray radiation," Mechanical Engineering Journal (2020), article 19-00634. As-grown Al single crystals were reported to contain subgrains; a residual-strain level of about 1e-4 was reported for the studied deformed specimen.
- P. Scardi et al., "Size-strain separation in diffraction line profile analysis," Journal of Applied Crystallography 51 (2018). Instrumental-profile and overlap treatment are required before specimen size/strain broadening is interpreted.

## Decision

For now:

1. keep the local-traction propagator **candidate-only**;
2. define its compatible P0 as structural/coarse-grained P0 built by spatial push-forward;
3. do not insert a Boltzmann/Gaussian thermal P0 merely to obtain width;
4. do not interpret the naive 300 K tail as fatigue initiation;
5. next test whether the candidate has any physically defensible laboratory-cycle evolution after the atomic-to-laboratory time-scale issue is imposed.

If that time-scale test fails, return to the exact finite-chain/correlation-hierarchy checkpoint as previously agreed.
