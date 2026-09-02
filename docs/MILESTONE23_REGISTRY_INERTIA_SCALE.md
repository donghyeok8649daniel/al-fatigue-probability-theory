# Milestone 23 — Registry inertia from the current reduced coordinate

## 0. Scope

This milestone keeps the current reduced row/layer model only.  No FCC
geometry, EAM/DFT surface, stochastic forcing, damping, Boltzmann law, or
Fokker--Planck closure is introduced.

The active intrinsic state energy is

$$
U_0(a,s)=\sum_{k\ge1}\sum_{p\in\mathbb Z}
v_{m,n}\!\left(\sqrt{k^2a^2+(pb+s)^2}\right).
$$

The unresolved quantity left by Milestone 22 was the registry inertia
$\mu_s$ in

$$
\mu_s\ddot s+\partial_sU_0(a,s)=0.
$$

The question is whether $\mu_s$ is actually a freely adjustable reduced
parameter once the meaning and counting of $U_0$ are respected.

---

## 1. Counting convention fixes the kinetic object

The audited two-row kernel

$$
W(d,s)=\sum_{p\in\mathbb Z}
v_{m,n}\!\left(\sqrt{d^2+(pb+s)^2}\right)
$$

is the cross-row interaction energy **per upper atom/repeat**.  The active
multilayer sum preserves that same reference-repeat counting: one reference
repeat interacts with the rows at $a,2a,3a,\ldots$.

This gives a finite kinetic interpretation that is consistent with the energy
counting:

> $s$ is the physical tangential displacement of one reference repeat (or a
> finite coherent reference patch) relative to its surrounding background.

The same $s$ appears in the interaction with every background row because the
same reference repeat has been translated relative to all of them.  It is not
necessary to interpret infinitely many background layers as all moving by the
same $s$; that interpretation was correctly rejected in Milestone 19 because
its kinetic mass diverges.

---

## 2. Exact reduced mass of the relative registry coordinate

Let $y_r$ be the tangential coordinate of the reference repeat of mass $m_r$,
and $y_b$ the coordinate of the participating background with effective mass
$M_b$. Define

$$
\boxed{s=y_r-y_b.}
$$

With center-of-mass coordinate removed, the exact two-body kinetic energy is

$$
T_s=\frac12\mu_s\dot s^2,
$$

where

$$
\boxed{
\mu_s=\frac{m_rM_b}{m_r+M_b}.
}
$$

Therefore

$$
\boxed{0<\frac{\mu_s}{m_r}<1}
$$

for finite positive $M_b$, and the fixed/heavy-background limit is

$$
\boxed{\mu_s\to m_r.}
$$

For two equal moving repeats,

$$
\boxed{\mu_s=\frac12m_r.}
$$

Hence, under the present physical-length definition of $s$, an inertia such as
$\mu_s\gg m_r$ is not produced by this local relative-coordinate embedding.

Classification: exact kinematics after the reference-repeat/background
embedding is declared.

---

## 3. Finite coherent patch does not generate a slow registry mode

Suppose a coherent patch contains $N$ identical repeat units and moves with one
common registry coordinate $s$.

The mass is extensive,

$$
M_N=N m_r,
$$

but the registry energy is also extensive,

$$
U_N(a,s)=N U_0(a,s).
$$

Thus the small-oscillation stiffness is

$$
K_{s,N}=N U_{ss},
$$

and

$$
\boxed{
\omega_{s,N}^2
=\frac{N U_{ss}}{N m_r}
=\frac{U_{ss}}{m_r}.
}
$$

Therefore simply increasing the coherent patch size cannot create a slow
fatigue-frequency registry mode.  Multiplying the inertia by $N$ while leaving
the one-repeat $U_0$ unchanged would mix incompatible extensive quantities.

---

## 4. Calibration-compatible frequency formula

Let the normalized active $U_0$ have equilibrium curvature ratio

$$
\boxed{
r_K=\frac{U_{ss}(a_0,s_0)}{U_{aa}(a_0,s_0)}.
}
$$

Use the existing normal calibration

$$
K_a=\frac{EA_0}{a_0},
$$

and

$$
\boxed{
t_0=\sqrt{\frac{m_r a_0}{EA_0}}.}
$$

If

$$
\rho_\mu=\frac{\mu_s}{m_r},
$$

then the physical registry curvature is

$$
K_s=K_a r_K,
$$

so

$$
\omega_s^2=\frac{K_s}{\mu_s}
=\frac{1}{t_0^2}\frac{r_K}{\rho_\mu}.
$$

Therefore

$$
\boxed{
\omega_s t_0=\sqrt{\frac{r_K}{\rho_\mu}},
\qquad
f_s=\frac{1}{2\pi t_0}\sqrt{\frac{r_K}{\rho_\mu}}.
}
$$

No numerical value of the row repeat $b$ is needed in this curvature-ratio
mapping.

The step $K_a=EA_0/a_0$ is a calibration bridge between the current reduced
multilayer energy and the existing normal stiffness scale.  It is not a new
crystal-geometry assumption.

---

## 5. Current normalized curvature ratio

For the same normalized diagnostic used in Milestones 21--22,

$$
m=12,\quad n=6,\quad b=\sigma_{LJ}=\epsilon_{LJ}=1,
$$

at

$$
a_0=0.9919601754,\qquad s_0/b=0.5,
$$

a converged direct double-sum calculation gives approximately

$$
\boxed{U_{aa}(a_0,s_0)=106.7616293}
$$

and

$$
\boxed{U_{ss}(a_0,s_0)=25.7179226.}
$$

Thus

$$
\boxed{r_K=0.2408910654.}
$$

The mixed curvature $U_{as}$ is zero at the symmetric registry point to
numerical precision, as required by the symmetry analysis.

---

## 6. Physical registry frequency under the existing Al normal calibration

The already retained normal calibration is

$$
a_0^{\rm phys}=2.8627442948\times10^{-10}\ {\rm m},
$$

$$
E=69\ {\rm GPa},
$$

$$
A_0=6.0338\times10^{-20}\ {\rm m^2},
$$

with

$$
\boxed{t_0=5.55046\times10^{-14}\ {\rm s}.}
$$

The implied reference-repeat mass is the Al atomic mass scale,

$$
m_r\approx4.48039\times10^{-26}\ {\rm kg}.
$$

For the fixed/heavy-background limit $\rho_\mu=1$,

$$
\boxed{
f_s\approx1.40735\times10^{12}\ {\rm Hz}.}
$$

For two equal participating repeat masses, $\rho_\mu=1/2$, so the frequency is
even larger by $\sqrt2$.

Hence the natural registry coordinate remains an atomic/THz-scale mode under
the present reduced-energy and normal-stiffness calibration.

---

## 7. Re-evaluation of the Milestone-22 parametric-resonance scan

Milestone 22 scanned a formally free normalized $\mu_s$ and found strong
five-cycle variational amplification around

$$
\mu_s^*\sim8\times10^2
$$

for the deliberately fast numerical loading $\Omega^*=0.35$.

The present kinetic derivation changes the interpretation of that result.
Under the natural per-repeat relative coordinate,

$$
\boxed{0<\rho_\mu\le1,}
$$

not $\rho_\mu\sim800$.

Therefore the previously observed resonance band lies outside the natural
inertial range of the current coordinate definition.  It remains a mathematical
mechanism diagnostic, but it must not be presented as a physically available
registry resonance without redefining the kinetic object and its energy
consistently.

---

## 8. Laboratory fatigue frequency makes the scale conflict much stronger

For principal small-modulation parametric resonance,

$$
\omega_{\rm load}\approx2\omega_s.
$$

Solving this condition for the required inertia ratio gives

$$
\boxed{
\rho_{\mu,{\rm req}}
=\frac{r_K}{(\pi f_{\rm load}t_0)^2}.
}
$$

At $f_{\rm load}=20\ {\rm Hz}$,

$$
\boxed{
\rho_{\mu,{\rm req}}\approx1.98\times10^{22}.
}
$$

Equivalently, while retaining a one-repeat registry stiffness, this would
correspond to

$$
\mu_{s,{\rm req}}\approx8.87\times10^{-4}\ {\rm kg},
$$

which is incompatible with the microscopic per-repeat coordinate.

A coherent patch cannot repair this mismatch because its energy and mass both
scale with the number of repeats and the frequency does not decrease.

Therefore the direct principal parametric-resonance explanation for laboratory
fatigue-frequency registry activation is rejected under the current reduced
coordinate.

This does not mathematically exclude all high-order, nonlinear, defect-mediated,
or thermally seeded mechanisms.  It specifically rejects the attempt to obtain
a 20-Hz-scale registry mode by assigning a huge free inertia to the existing
one-repeat $s$ coordinate.

---

## 9. Consequence for the active mainline

For the ideal symmetric, zero-temperature, defect-free reduced model under pure
normal loading,

$$
\boxed{P(a,s,t)=P(a,t)\,\delta(s-s_0)}
$$

remains the physically defensible baseline.

The normal density $P(a,t)$ is mechanically generated by spatial wave diversity.
The registry coordinate becomes active only when a physically justified
symmetry-breaking seed or additional slow structural degree of freedom is
identified.

Therefore the next target is **not** to tune $\mu_s$.  It is to identify a
physical source of registry symmetry breaking or a different slow internal
state whose time scale follows from the mechanics.
