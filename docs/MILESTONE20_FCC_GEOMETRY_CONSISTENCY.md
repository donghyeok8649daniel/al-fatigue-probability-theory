# Milestone 20 — FCC geometry consistency audit for $(a,s)$

## 0. Why this audit is necessary

Milestone 19 showed that a finite coordinate embedding is required before an
inertial metric can be assigned to the reduced state

$$
\mathbf q=(a,s).
$$

The next natural idea is to identify $a$ with normal opening of an FCC slip
plane and $s$ with translation along a perfect FCC slip direction.  Doing that
forces a precise virtual-work calculation.  The result exposes an important
geometric inconsistency in the older reduced driving pair

$$
Q_a=A_0\sigma,\qquad Q_s=A_0M\sigma
$$

when $M\neq0$ and both coordinates are interpreted as components of one
literal plane-relative displacement.

This milestone records the correction before that reduced pair is used in the
new mechanics-generated probability model.

---

## 1. Declared FCC slip-plane basis

For a perfect FCC slip system,

$$
\mathbf n\in\{111\},\qquad
\mathbf d\in\langle110\rangle,\qquad
\mathbf n\cdot\mathbf d=0,
$$

where $\mathbf n$ is the unit slip-plane normal and $\mathbf d$ is the unit
in-plane slip direction.

For a cubic lattice parameter $a_{\rm lat}$,

$$
\boxed{d_{111}=\frac{a_{\rm lat}}{\sqrt3}}
$$

and the perfect $a_{\rm lat}/2\langle110\rangle$ translation magnitude is

$$
\boxed{b_{\rm perf}=\frac{a_{\rm lat}}{\sqrt2}}.
$$

The triangular $(111)$ surface primitive-cell area is

$$
\boxed{A_{111}=\frac{\sqrt3}{4}a_{\rm lat}^2}.
$$

These are crystallographic geometry identities.  They do **not** by themselves
identify the unresolved fatigue energy area $A_0$ with $A_{111}$.

---

## 2. Literal plane-relative interpretation of $(a,s)$

Suppose one finite plane patch has relative displacement

$$
\mathbf r(a,s)=a\mathbf n+s\mathbf d.
$$

Then

$$
\delta\mathbf r=\mathbf n\,\delta a+\mathbf d\,\delta s.
$$

This gives a clean physical meaning to the coordinates:

- $a$: opening normal to the declared slip plane;
- $s$: translation along the declared slip direction.

It is also exactly the type of orthogonal normal/registry geometry suggested
by distances of the form $\sqrt{d^2+(pb+s)^2}$ in the current row model.

---

## 3. Generalized forces from virtual work

Let the uniaxial tensile axis be the unit vector $\boldsymbol\ell$ and the
Cauchy stress tensor be

$$
\boldsymbol\sigma=\sigma\,\boldsymbol\ell\otimes\boldsymbol\ell.
$$

Traction on the declared plane is

$$
\mathbf t=\boldsymbol\sigma\mathbf n
=\sigma(\boldsymbol\ell\cdot\mathbf n)\boldsymbol\ell.
$$

For plane-patch area $A_0$,

$$
\delta W=A_0\mathbf t\cdot\delta\mathbf r
=Q_a\,\delta a+Q_s\,\delta s.
$$

Therefore

$$
\boxed{
Q_a=A_0\sigma(\boldsymbol\ell\cdot\mathbf n)^2
}
$$

and

$$
\boxed{
Q_s=A_0\sigma
(\boldsymbol\ell\cdot\mathbf n)
(\boldsymbol\ell\cdot\mathbf d).
}
$$

Define the signed Schmid factor

$$
\boxed{
M=(\boldsymbol\ell\cdot\mathbf n)
(\boldsymbol\ell\cdot\mathbf d).
}
$$

Then the second equation is exactly

$$
Q_s=A_0M\sigma.
$$

Thus the old $Q_s$ expression is geometrically consistent with slip-plane
virtual work, but the old $Q_a=A_0\sigma$ expression is **not** the corresponding
normal generalized force unless the plane normal is parallel to the loading
axis.

---

## 4. The incompatibility is exact

If

$$
Q_a=A_0\sigma
$$

is to be recovered from the literal plane-opening formula, then

$$
|\boldsymbol\ell\cdot\mathbf n|=1,
$$

so the slip-plane normal is parallel or antiparallel to the tensile axis.
Because the slip direction lies in that plane,

$$
\mathbf n\cdot\mathbf d=0,
$$

which immediately gives

$$
\boldsymbol\ell\cdot\mathbf d=0
$$

and hence

$$
\boxed{M=0.}
$$

Therefore

$$
\boxed{
Q_a=A_0\sigma\ \text{and}\ Q_s=A_0M\sigma\ (M\neq0)
}
$$

cannot both be obtained from one literal FCC plane-relative coordinate pair
$(a,s)$.

This is a geometry result, not a modeling preference.

---

## 5. Example: [001] tension on one FCC slip system

Take

$$
\boldsymbol\ell=[001],\qquad
(111)[10\bar1].
$$

Then

$$
\boldsymbol\ell\cdot\mathbf n=\frac1{\sqrt3},
\qquad
\boldsymbol\ell\cdot\mathbf d=-\frac1{\sqrt2}.
$$

Thus

$$
\boxed{
\frac{Q_a}{A_0\sigma}=\frac13
}
$$

and

$$
\boxed{
M=-\frac1{\sqrt6}\approx-0.408248.
}
$$

The same applied uniaxial stress therefore produces both a normal traction and
a resolved shear traction on the slip plane, but the normal coordinate is not
driven by the full axial force $A_0\sigma$.

---

## 6. Finite two-rigid-plane kinetic metric

A finite kinematic embedding can be made without an arbitrary effective mass.
Let the upper and lower rigid plane patches have masses $M_+$ and $M_-$ and fix
their common center of mass.  Their relative coordinate is again

$$
\mathbf r=a\mathbf n+s\mathbf d.
$$

The reduced mass is

$$
\mu=\frac{M_+M_-}{M_++M_-}.
$$

Therefore

$$
T=\frac12\mu|\dot{\mathbf r}|^2
$$

and

$$
\boxed{
\mathbf G
=\mu
\begin{bmatrix}
1 & \mathbf n\cdot\mathbf d\\
\mathbf n\cdot\mathbf d & 1
\end{bmatrix}.
}
$$

For a valid slip system, $\mathbf n\cdot\mathbf d=0$, so

$$
\boxed{
\mathbf G=\mu\mathbf I.
}
$$

For equal patch masses $M_+=M_-=M_{\rm patch}$,

$$
\mu=\frac{M_{\rm patch}}2.
$$

This is a finite, exactly derived inertial metric for the declared two-rigid-
plane embedding.  It is **not yet** the kinetic metric of the current $U_0$
because the energy-to-patch-area mapping remains unresolved.

---

## 7. Critical energy-geometry mismatch that remains

The current intrinsic energy is

$$
U_0(a,s)
=\sum_{k\ge1}\sum_{p\in\mathbb Z}
v_{m,n}\!\left(\sqrt{k^2a^2+(pb+s)^2}\right).
$$

This is an ideal row/layer geometry with one normal coordinate and one row
translation coordinate.  It is not yet a demonstrated full three-dimensional
FCC $(111)$ stacking energy per surface primitive cell.

In particular, real FCC close-packed planes have three-dimensional triangular
in-plane geometry and ABC stacking offsets.  Those features are not established
by the present $(k,p)$ row sum.

Therefore the following identification is **not allowed yet**:

$$
A_0\stackrel{?}{=}A_{111}
$$

followed by treating current $U_0$ as the exact energy of that atomistic FCC
plane patch.

The geometry module intentionally computes $A_{111}$ but does not assign it to
$A_0$.

---

## 8. Two physically distinct routes from here

### Route A — slip-plane coupled state

Interpret

$$
a=\text{slip-plane normal opening},\qquad
s=\text{slip translation}.
$$

Then use the exact generalized-force pair

$$
Q_a=A_0\sigma(\ell\cdot n)^2,
\qquad
Q_s=A_0M\sigma.
$$

This route is geometrically coherent for the two-coordinate plane-relative
mechanics, but requires upgrading/validating $U_0(a,s)$ against actual FCC
stacking.

### Route B — loading-axis spacing plus internal slip

Keep

$$
a=\text{spacing/stretch along the tensile axis}
$$

and retain full axial drive $Q_a=A_0\sigma$.

Then $s$ is no longer simply the orthogonal tangential coordinate of the same
plane pair.  A physically literal energy must be re-derived from FCC lattice
vectors under axial deformation plus slip; the simple
$\sqrt{k^2a^2+(pb+s)^2}$ geometry cannot be assumed to remain exact.

Because the project has repeatedly defined the macroscopic problem as 1D
uniaxial tension with an internal slip coordinate, **Route B is conceptually
closest to the original research goal**, but it needs the harder FCC
energy-geometry derivation.

No route is silently selected in this milestone.

---

## 9. Consequence for the probability model

The exact moment identity

$$
\boldsymbol\Theta\nabla\ln P
=\boldsymbol{\mathcal A}-D_t\mathbf u-\nabla\cdot\boldsymbol\Theta
$$

is unaffected by this audit.  What changes is the mechanics used to calculate
$\boldsymbol{\mathcal A}$.

The probability framework therefore survives; the current task is to make the
coordinate-energy-load map physically self-consistent before using it for a
quantitative Al prediction.

---

## 10. Status

Implemented in `theory/fcc_slip_kinematics.py`:

- FCC $d_{111}$;
- perfect $\langle110\rangle$ translation magnitude;
- $(111)$ primitive area;
- validation of one `{111}<110>` slip system;
- signed Schmid factor;
- exact uniaxial plane-patch generalized forces;
- finite two-rigid-plane reduced mass metric.

The next derivation should decide between Route A and Route B **before** a
quantitative FCC/Al $(a,s)$ mechanics simulation is claimed.
