# Local-traction P0 propagator — laboratory-timescale verdict

## Verdict

The candidate does satisfy the formal reduced input/output requirement

$$
P_0+\sigma(0:t)\rightarrow P(t),
$$

but in its current **conservative atomic-inertia form it does not provide a physically defensible mechanism for progressive subcritical fatigue accumulation at laboratory frequencies**.

Accordingly it remains candidate-only and the main theory path returns to the exact finite-chain / correlation-hierarchy checkpoint.

## Why the first candidate audit looked strongly history dependent

The first diagnostic used a nonzero initial mean traction corresponding to 100 MPa while the initial spacing support was centered around lambda=1 and all initial rates were zero. Since phi'(1)=0, the central characteristic had nonzero initial acceleration

$$
\ddot\lambda(0)\approx1.4493\times10^{-3}.
$$

Across the initial support [0.997, 1.003], the initial acceleration ranged from approximately -1.4570e-3 to 4.5464e-3. Thus the initial state was not mechanically relaxed and launched atomic-scale transients.

This does not invalidate the characteristic mathematics, but it invalidates using that strong non-retracing signal as evidence of laboratory fatigue accumulation.

## Frequency mapping

With the retained normal time scale

$$
t_0=5.55046\times10^{-14}\ \mathrm{s},
$$

physical frequency maps as

$$
\omega^*=2\pi f t_0.
$$

| f | omega* | linear inertial correction for kappa=1 |
|---:|---:|---:|
| 1 Hz | 3.4875e-13 | 1.2162e-25 |
| 20 Hz | 6.9749e-12 | 4.8649e-23 |
| 100 Hz | 3.4875e-11 | 1.2162e-21 |
| 1000 Hz | 3.4875e-10 | 1.2162e-19 |

The previous diagnostic value omega*=0.02 maps to approximately 5.735e10 Hz = 57.35 GHz.

## Analytic slow-loading check

For a mechanically equilibrated characteristic, linearization gives

$$
\ddot y+\kappa y=q_a\sin(\omega^*\tau).
$$

Away from resonance,

$$
y_p=\frac{q_a}{\kappa-(\omega^*)^2}\sin(\omega^*\tau).
$$

The response has no dissipative phase lag, and

$$
\oint q\,dy_p=0.
$$

The relative departure from the quasistatic amplitude is

$$
\frac{(\omega^*)^2}{\kappa-(\omega^*)^2}.
$$

At 20 Hz with kappa near unity this is about 4.865e-23. The local atomic coordinate is therefore effectively quasistatic under ordinary fatigue frequencies.

Nonlinear dynamics can differ near resonance, a separatrix, or direct instability, but that does not supply a generic progressive low-frequency subcritical-fatigue mechanism.

## Consequence

The physical P0 work remains useful:

$$
\text{prepared residual spacing/microstrain field}\rightarrow P_0.
$$

But the propagation mechanism must still explain why the structural probability state changes cycle after cycle at laboratory frequencies. The present local conservative oscillator does not.

No arbitrary viscosity, diffusion coefficient, white-noise bath, or fitted damage variable should be inserted to force that behavior.

The candidate can be reopened only if a slow internal mechanism is derived physically or if a multiscale reduction shows how collective dynamics generate the slow probability evolution.
