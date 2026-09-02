# Registry inertia and fatigue-frequency scale audit

## Scope

This numerical audit uses no FCC geometry.  It combines the active reduced
row/layer energy with the existing normal calibration to determine whether the
registry inertia can legitimately be tuned to produce a fatigue-frequency
parametric mode.

The current row kernel is counted per upper atom/repeat, so the finite kinetic
coordinate is one reference repeat (or a finite coherent patch) translating by
the physical registry displacement $s$ relative to its background.

## 1. Normalized curvature ratio

At the normalized equilibrium used in the prior spatial-chain calculations,

$$
a_0^*=0.9919601754,\qquad s_0/b=0.5,
$$

with $m=12$, $n=6$, $b=\sigma_{LJ}=\epsilon_{LJ}=1$, the direct double-sum
curvatures converge to

$$
U_{aa}=106.7616293,
$$

$$
U_{ss}=25.71792262,
$$

and therefore

$$
\boxed{r_K=U_{ss}/U_{aa}=0.2408910654.}
$$

The convergence from $(k_{max},p_{max})=(20,50)$ through $(200,500)$ is saved
in `curvature_convergence.csv`.

## 2. Kinetic interpretation of $s$

For reference-repeat mass $m_r$ and participating background mass $M_b$, the
relative coordinate $s=y_r-y_b$ has

$$
\boxed{\mu_s=\frac{m_rM_b}{m_r+M_b}.}
$$

Hence the simplest current-coordinate embedding gives

$$
0<\rho_\mu\equiv\frac{\mu_s}{m_r}\le1,
$$

with $\rho_\mu=1$ for a fixed/heavy background and $1/2$ for equal moving
repeat masses.

A coherent $N$-repeat patch does not lower the natural frequency because both
mass and energy curvature scale by $N$.

## 3. Physical frequency from the existing normal calibration

Use the retained normal calibration

$$
a_0=2.8627442948\times10^{-10}\,\mathrm{m},
$$

$$
E=69\,\mathrm{GPa},
$$

$$
A_0=6.0338\times10^{-20}\,\mathrm{m^2},
$$

$$
t_0=5.55046\times10^{-14}\,\mathrm{s}.
$$

The implied repeat mass is

$$
m_r=t_0^2EA_0/a_0
\approx4.48039\times10^{-26}\,\mathrm{kg}.
$$

Matching the active normal curvature to the normal stiffness $EA_0/a_0$ gives

$$
K_a\approx14.5431\,\mathrm{N/m},
$$

$$
K_s=K_a r_K\approx3.50331\,\mathrm{N/m}.
$$

Thus

$$
\boxed{
f_s=\frac{1}{2\pi t_0}\sqrt{\frac{r_K}{\rho_\mu}}.}
$$

For the largest natural inertia, $\rho_\mu=1$,

$$
\boxed{f_s\approx1.40735\times10^{12}\,\mathrm{Hz}.}
$$

For equal moving masses, $\rho_\mu=1/2$,

$$
f_s\approx1.99029\times10^{12}\,\mathrm{Hz}.
$$

The current registry coordinate is therefore an atomic/THz mode, not a direct
fatigue-frequency mode.

## 4. Principal parametric-resonance requirement

The small-modulation principal condition is approximately

$$
\omega_{load}=2\omega_s.
$$

Solving for the inertia ratio,

$$
\boxed{
\rho_{\mu,req}=\frac{r_K}{(\pi f_{load}t_0)^2}.
}
$$

At 20 Hz this gives

$$
\boxed{\rho_{\mu,req}\approx1.98063\times10^{22}.}
$$

Keeping the one-repeat stiffness, this corresponds to

$$
\boxed{\mu_{s,req}\approx8.87398\times10^{-4}\,\mathrm{kg}.}
$$

This is incompatible with the per-repeat microscopic coordinate.  Enlarging a
coherent patch cannot supply this slow mode consistently because its registry
energy grows with its mass.

## 5. Correction to the previous free-$\mu_s$ scan

The previous normalized variational scan found strong amplification near
$\mu_s^*\sim8\times10^2$ for a deliberately fast numerical drive
$\Omega^*=0.35$.  That result remains mathematically correct as a transfer
matrix diagnostic, but the present kinetic derivation shows that such a large
inertia is outside the natural range of the current physical registry
coordinate.

Therefore the claim must be narrowed to:

> time-dependent normal opening modulates registry stiffness, but the direct
> principal parametric route is not available at laboratory fatigue frequencies
> for the present per-repeat registry coordinate.

## 6. Consequence

For ideal exact symmetry with no physical seed,

$$
P(a,s,t)=P(a,t)\delta(s-s_0)
$$

remains the safe pure-normal baseline.

The next theory problem is not to tune $\mu_s$.  It is to identify a physically
justified symmetry-breaking source or a genuinely slow internal structural
state.  Random Gaussian noise or an arbitrary registry drive is not required
and should not be inserted merely to activate $s$.
