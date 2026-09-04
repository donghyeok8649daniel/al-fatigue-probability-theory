# Local-traction P0 propagator: laboratory-time-scale verdict

Status: **CANDIDATE VERDICT — the present conservative local-traction model is not promoted as a laboratory fatigue-accumulation mechanism.**

This note follows the agreed decision rule: try the P0-only local-traction reduction, and if its physical cycle evolution fails, return to the exact finite-chain / correlation-hierarchy checkpoint rather than adding arbitrary damping, diffusion, or a fitted lifetime law.

## 1. A preparation inconsistency in the first candidate audit

The first local-traction audit used

$$
q(0)=\frac{100\ \mathrm{MPa}}{69\ \mathrm{GPa}}\approx1.4493\times10^{-3}
$$

while the diagnostic initial spacing support was centered at lambda=1 and the initial rates were set to zero.

Because

$$
\phi'(1)=0,
$$

the central characteristic had the initial acceleration

$$
\ddot\lambda(0)=q(0)-\phi'(1)\approx1.4493\times10^{-3},
$$

not zero.

Across the diagnostic interval 0.997 <= lambda_0 <= 1.003, the initial acceleration ranged approximately from

$$
-1.4570\times10^{-3}
$$

to

$$
4.5464\times10^{-3}.
$$

Therefore the strong history dependence in that audit cannot be interpreted as established fatigue accumulation. It contains a deliberately unrelaxed atomic-scale startup transient.

## 2. Structural P0 requires mechanical support of the initial heterogeneity

If P0 is interpreted as the slow structural distribution defined in `CANDIDATE_PHYSICAL_P0_CONSTRUCTION.md`, every initial spacing must be mechanically compatible with the prepared state.

One explicit candidate realization is to associate each initial spacing label lambda_0 with a residual normalized traction

$$
q_0^{\mathrm{res}}(\lambda_0)
=\phi'(\lambda_0)-q_{\mathrm{ext}}(0).
$$

Then the characteristic equation becomes

$$
\frac{d\lambda}{d\tau}=c,
$$

$$
\frac{dc}{d\tau}
=q_{\mathrm{ext}}(\tau)+q_0^{\mathrm{res}}(\lambda_0)-\phi'(\lambda),
$$

with

$$
\lambda(0)=\lambda_0,
\qquad
c(0)=0.
$$

This gives

$$
\left.\frac{dc}{d\tau}\right|_{\tau=0}=0.
$$

No extra probability distribution is needed because the residual traction is a deterministic function of the initial spacing label under this particular construction. Other structural-support mechanisms would need their own derivation.

## 3. Linearized slow-loading limit

Linearize one mechanically prepared characteristic about its stable equilibrium. Let y be the small spacing perturbation and let

$$
\kappa=\phi''(\lambda_0)>0.
$$

For a small sinusoidal traction perturbation,

$$
\Delta q(\tau)=q_a\sin(\omega^*\tau),
$$

the conservative linearized equation is

$$
\ddot y+\kappa y=q_a\sin(\omega^*\tau).
$$

Away from resonance, the periodic particular solution is

$$
y_p(\tau)
=\frac{q_a}{\kappa-(\omega^*)^2}
\sin(\omega^*\tau).
$$

There is no dissipative phase lag. The cycle work of the periodic particular solution is exactly

$$
\oint \Delta q\,dy_p=0.
$$

The quasistatic response is

$$
y_{\mathrm{qs}}=\frac{q_a}{\kappa}\sin(\omega^*\tau),
$$

so the relative inertial correction is

$$
\frac{y_p-y_{\mathrm{qs}}}{y_{\mathrm{qs}}}
=\frac{(\omega^*)^2}{\kappa-(\omega^*)^2}.
$$

Near lambda=1, kappa=1. Therefore at very small omega* the conservative local coordinate becomes quasistatic and reversible to order (omega*)^2.

## 4. Laboratory frequency mapped to the retained atomic time scale

Using the retained normal-chain calibration

$$
t_0\approx5.55046\times10^{-14}\ \mathrm{s},
$$

a physical loading frequency f maps to

$$
\omega^*=2\pi f t_0.
$$

Representative values are

- 1 Hz: omega* ~= 3.4875e-13
- 20 Hz: omega* ~= 6.9749e-12
- 100 Hz: omega* ~= 3.4875e-11
- 1000 Hz: omega* ~= 3.4875e-10

At 20 Hz and kappa=1, the relative inertial correction is only

$$
\frac{(\omega^*)^2}{1-(\omega^*)^2}
\approx4.865\times10^{-23}.
$$

By contrast, the first candidate audit used

$$
\omega^*=0.02,
$$

which corresponds, under the same retained time scale, to approximately

$$
f\approx5.735\times10^{10}\ \mathrm{Hz}
$$

or 57.35 GHz.

Therefore that audit is a high-frequency mechanism test, not a laboratory-fatigue-frequency validation.

## 5. Nonlinear implication

The exact nonlinear conservative coordinate can depart from the linear result near resonance, a separatrix, or a direct instability. However, for laboratory loading that is extremely slow compared with the retained atomic time scale and remains on a stable branch, the local coordinate approaches adiabatic/quasistatic following.

If it begins from a mechanically prepared equilibrium and the applied traction returns to the same value without crossing a direct instability, the present model contains no identified mechanism that produces a systematic irreversible cycle-to-cycle drift of the structural P0 population.

A direct threshold crossing under the peak load could still occur, but that is a load-threshold event rather than a progressive subcritical fatigue-accumulation mechanism.

## 6. Verdict

The local-traction P0 propagator succeeds at the mathematical input-output requirement

$$
P_0+\sigma(0:t)\longrightarrow P(t)
$$

under its stated preparation assumptions.

But **in its present conservative atomic-inertia form it fails the stronger physical requirement of explaining progressive laboratory-frequency fatigue accumulation below direct instability.**

The earlier strong non-retracing audit is not sufficient evidence because its initial state was not mechanically equilibrated and its dimensionless loading frequency corresponds to tens of GHz under the retained time calibration.

Therefore:

1. do not promote the candidate local-traction model to the active paper theory;
2. retain the physical P0 push-forward construction as useful initialization work;
3. return the main theory path to the exact finite-chain / correlation-hierarchy checkpoint;
4. treat the atomic-to-laboratory cycle-evolution mechanism as an unresolved central physics problem;
5. do not repair the failure by inserting arbitrary viscosity, diffusion, white noise, or empirical damage.

## 7. What would reopen this candidate

The candidate may be reconsidered only if a physically derived slow internal mechanism is added independently of the desired fatigue fit, or if a rigorous multiscale reduction shows how collective microscopic modes generate a slow probability evolution from the supplied stress field.
