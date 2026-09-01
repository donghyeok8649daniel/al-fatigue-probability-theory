# Milestone 19 — Assumption-explicit coupled mechanics for $(a,s)$

## 0. Purpose

The active probability state is

$$
\mathbf q=(a,s),
$$

with intrinsic local lattice energy

$$
U_0(a,s)
=\sum_{k\ge1}\sum_{p\in\mathbb Z}
v_{m,n}\!\left(\sqrt{k^2a^2+(pb+s)^2}\right).
$$

Milestone 18 established the exact smooth-moment identity

$$
\boldsymbol\Theta\nabla\ln P
=\boldsymbol{\mathcal A}-D_t\mathbf u-\nabla\cdot\boldsymbol\Theta.
$$

The remaining question is how mechanics supplies the conditional acceleration
$\boldsymbol{\mathcal A}$.  This milestone derives the most that can be obtained
without silently assuming a mass, damping law, Fokker--Planck equation, thermal
bath, or arbitrary spatial spring coupling.

---

## 1. First assumption ledger

The following statements have different logical status and must not be mixed.

### Already fixed by the active reduced theory

1. The state coordinates are collective normal opening/spacing $a$ and unwrapped collective registry $s$.
2. The intrinsic state energy is the multiplicity-free $U_0(a,s)$ above.
3. External work is not inserted into $U_0$.
4. $P(a,s,t)$ is a state density, not a prescribed Gaussian, Weibull, or Boltzmann family.

### Not fixed by $U_0(a,s)$

1. The inertial mass associated with changing $a$ or $s$.
2. Whether the reduced coordinates are overdamped, underdamped, or quasistatic.
3. Spatial coupling between distinct representative patches.
4. A dissipative/irreversible mechanism.
5. The finite representative volume/area whose mass and energy are paired with $U_0$.

Therefore none of these may be inserted merely to obtain a broad PDF or a
hysteresis loop.

---

## 2. Finite-coordinate embedding

Choose a **finite** representative atomic/layer set with microscopic positions

$$
\mathbf R_A=\mathbf R_A(a,s),
$$

where $A$ indexes atoms or rigid represented objects and $m_A$ are their masses.
This is a coordinate embedding, not a probability assumption.

The velocity of object $A$ is

$$
\dot{\mathbf R}_A
=\frac{\partial\mathbf R_A}{\partial q_i}\dot q_i.
$$

Hence its exact kinetic energy on this two-coordinate manifold is

$$
T=\frac12\sum_A m_A|\dot{\mathbf R}_A|^2
=\frac12G_{ij}(\mathbf q)\dot q_i\dot q_j,
$$

with the generalized mass metric

$$
\boxed{
G_{ij}(\mathbf q)
=\sum_A m_A
\frac{\partial\mathbf R_A}{\partial q_i}\cdot
\frac{\partial\mathbf R_A}{\partial q_j}.
}
$$

A diagonal $G$ is therefore **not** assumed.  If the chosen geometry produces
$G_{as}\neq0$, it must be retained.

Classification: exact kinematics **within the chosen finite coordinate
embedding**.

---

## 3. External generalized force from virtual work

For microscopic external forces $\mathbf F_A^{\rm ext}$,

$$
\delta W_{\rm ext}
=\sum_A\mathbf F_A^{\rm ext}\cdot\delta\mathbf R_A
=Q_i\,\delta q_i,
$$

so

$$
\boxed{
Q_i
=\sum_A\mathbf F_A^{\rm ext}\cdot
\frac{\partial\mathbf R_A}{\partial q_i}.
}
$$

Thus $Q_a$ and $Q_s$ must ultimately come from the declared loading geometry.
Expressions such as $Q_a=A_0\sigma$ or a Schmid-projected $Q_s$ are geometry
mappings that require the representative area/orientation definition; they are
not identities supplied by $U_0$ itself.

---

## 4. Euler--Lagrange mechanics

Use

$$
L=T-U_0(a,s).
$$

With generalized external force $Q_i$, Euler--Lagrange gives

$$
\frac{d}{dt}\frac{\partial L}{\partial\dot q_i}
-\frac{\partial L}{\partial q_i}=Q_i.
$$

Define the first-kind metric connection

$$
\Gamma_{i,jk}
=\frac12\left(
\partial_jG_{ik}+\partial_kG_{ij}-\partial_iG_{jk}
\right).
$$

Then

$$
\boxed{
G_{ij}\ddot q_j
+\Gamma_{i,jk}\dot q_j\dot q_k
+\partial_iU_0
=Q_i.
}
$$

This is the minimal conservative coupled $(a,s)$ mechanics on the chosen
coordinate manifold.

If the embedding gives a constant metric,

$$
\boxed{
\ddot{\mathbf q}
=\mathbf G^{-1}\left(\mathbf Q-\nabla U_0\right).
}
$$

No damping term has been added.

---

## 5. Conditional acceleration for the probability-shape equation

Define

$$
\mathbf u(\mathbf q,t)
=\mathbb E[\dot{\mathbf q}\mid\mathbf q],
$$

and

$$
\boldsymbol\Theta(\mathbf q,t)
=\operatorname{Cov}(\dot{\mathbf q}\mid\mathbf q).
$$

Then

$$
\mathbb E[\dot q_j\dot q_k\mid\mathbf q]
=u_ju_k+\Theta_{jk}.
$$

If, at fixed $\mathbf q,t$, the generalized force is deterministic and no
unresolved force has been omitted, conditional averaging of the reduced
Euler--Lagrange equation gives

$$
\boxed{
\mathcal A_i
=(G^{-1})_{i\ell}
\left[
Q_\ell-\partial_\ell U_0
-\Gamma_{\ell,jk}
\left(u_ju_k+\Theta_{jk}\right)
\right].
}
$$

For a constant metric,

$$
\boxed{
\boldsymbol{\mathcal A}
=\mathbf G^{-1}(\mathbf Q-\nabla U_0).
}
$$

This is the direct mechanics-to-probability bridge that was missing in
Milestone 18.

Substitution into the exact shape identity gives, for a constant metric,

$$
\boxed{
\boldsymbol\Theta\nabla\ln P
=
\mathbf G^{-1}(\mathbf Q-\nabla U_0)
-D_t\mathbf u
-\nabla\cdot\boldsymbol\Theta.
}
$$

This still does **not** constitute a closed probability evolution law:
$\mathbf u$ and $\boldsymbol\Theta$ still require dynamics/ensemble information.

---

## 6. Energy consistency

The reduced conservative mechanical energy is

$$
E_{\rm mech}
=\frac12\dot{\mathbf q}^{T}\mathbf G\dot{\mathbf q}+U_0(a,s).
$$

For a time-independent coordinate embedding and intrinsic potential,
Euler--Lagrange gives

$$
\boxed{
\frac{dE_{\rm mech}}{dt}
=\mathbf Q\cdot\dot{\mathbf q}.
}
$$

At zero external power, the conservative reduced subsystem conserves
$E_{\rm mech}$.

This must be distinguished from governing relation G2,

$$
\boxed{
\bar U(t)=\iint\Delta U_0(a,s)P(a,s,t)\,da\,ds,
}
$$

which is the ensemble mean **intrinsic lattice energy**, not total kinetic plus
potential energy.

Therefore a numerical $(a,s)$ mechanics implementation must check both:

1. trajectory-level mechanical work/energy balance;
2. density-level G2 mean intrinsic energy.

---

## 7. Important obstruction: $U_0$ alone does not define an inertial mass

The current $U_0(a,s)$ is a multiplicity-free local state energy obtained by
summing the interaction of one reference layer with layers at $a,2a,3a,\ldots$
and using the same collective registry parameter $s$ in every term.

It is tempting to interpret all of those infinitely many layers as literal
moving masses with coordinates proportional to $ka$ and the same $s$.  That is
not legitimate.  Such a naive infinite moving-layer embedding produces mass
sums of the form

$$
\sum_{k\ge1}m_k k^2
$$

for the $a$ coordinate and

$$
\sum_{k\ge1}m_k
$$

for the $s$ coordinate, which do not define a finite local inertial metric.

Therefore:

$$
\boxed{
U_0(a,s)\ \text{defines the intrinsic energy landscape, but not a unique
finite kinetic metric.}
}
$$

A finite representative moving set must be specified before inertial dynamics
is claimed quantitatively.

This is not a failure of the energy model; it is a separation between local
state energy and the kinetic embedding of that local state.

---

## 8. Consequence for the PDF

A single deterministic two-coordinate trajectory with identical initial
conditions for every represented state gives

$$
P(a,s,t)=\delta[a-a(t)]\delta[s-s(t)],
$$

and therefore

$$
\boldsymbol\Theta=0.
$$

The divided $\boldsymbol\Theta^{-1}$ shape formula is then inapplicable.
Consequently a broad physical $P$ requires a real source of state diversity,
such as spatially distinct coupled cells, physically specified microscopic
initial variation, or other unresolved degrees of freedom.  That source must
be derived or measured; it must not be inserted only to create a desired PDF.

---

## 9. Plasticity and irreversibility

Keep $s$ unwrapped,

$$
s=s_0+zb+\tilde s.
$$

Crossing from one registry well to another is a slip event and changes the
well-index population.  However, in the **purely conservative** mechanics above,
a barrier crossing alone is not sufficient to prove permanent plastic strain
or hysteretic dissipation.  The system can retain kinetic energy or recross the
barrier.

Therefore the physically strong criterion remains residual well-index change
after unloading **and an explicitly justified relaxation/irreversible
mechanism**.

G3,

$$
E_{\rm hyst}=\int\dot D_{\rm irr}\,dt,
$$

is intentionally not derived from the conservative equations in this
milestone.

---

## 10. What is now closed and what remains open

### Derived now

- generalized mass metric from a finite coordinate embedding;
- exact reduced Euler--Lagrange coupling of $a$ and $s$;
- direct expression for conditional acceleration $\boldsymbol{\mathcal A}$;
- mechanics-to-$\Theta$-density bridge;
- conservative work/energy identity;
- proof that $U_0$ alone does not specify a unique finite inertial metric.

### Still physically unresolved

1. The finite representative embedding consistent with the energy calibration area $A_0$.
2. Spatial coupling between multiple represented $(a,s)$ states.
3. The physical origin/evolution of $\boldsymbol\Theta$ rather than merely measuring it.
4. The irreversible mechanism required for G3 and permanent plasticity.
5. Quantitative Al orientation/area/mass mapping.

These are the next derivation targets.  No Smoluchowski, Langevin noise, damping,
or ad hoc nearest-patch spring should be promoted to the mainline until its
assumptions are checked against the intended physical regime.
