# Canonical Poisson/Bessel infinite-lattice potential

## Geometry

The lower row is infinite,

\[
(x_n^-,y_n^-)=(nb,0),\qquad n\in\mathbb Z,
\]

and one upper cell is at

\[
(x^+,y^+)=(b/2+s,a).
\]

Therefore the pair distance is

\[
r_n(a,s)=\sqrt{a^2+[(n+1/2)b-s]^2}.
\]

For the Lennard--Jones pair potential

\[
\phi_{\rm LJ}(r)=4\epsilon\left[(\sigma/r)^{12}-(\sigma/r)^6\right],
\]

the exact infinite lower-row interaction is

\[
W(a,s)=4\epsilon\left[\sigma^{12}S_6(a,s)-\sigma^6S_3(a,s)\right],
\]

where

\[
S_p(a,s)=\sum_{n=-\infty}^{\infty}
\left[a^2+((n+1/2)b-s)^2\right]^{-p}.
\]

## Poisson-summed form

With the Fourier-transform convention used by Poisson summation,

\[
S_p(a,s)=A_p(a)+\sum_{m=1}^{\infty}B_{p,m}(a)
\cos\left[2\pi m\left(\frac12-\frac{s}{b}\right)\right],
\]

with

\[
A_p(a)=\frac{\sqrt\pi\,\Gamma(p-1/2)}{b\,\Gamma(p)}a^{1-2p},
\]

and

\[
B_{p,m}(a)=
\frac{4\sqrt\pi}{b\Gamma(p)}
\left(\frac{\pi m}{ab}\right)^{p-1/2}
K_{p-1/2}\left(\frac{2\pi ma}{b}\right).
\]

The half-cell phase is not optional. It represents the staggered geometry in which `s=0` places the upper row halfway between lower-row lattice sites. Equivalently,

\[
\cos\left[2\pi m\left(\frac12-\frac{s}{b}\right)\right]
=(-1)^m\cos\left(\frac{2\pi ms}{b}\right).
\]

For the LJ powers,

\[
A_6(a)=\frac{63\pi}{256b}a^{-11},\qquad
A_3(a)=\frac{3\pi}{8b}a^{-5}.
\]

Thus the registry-averaged normal cohesive part is

\[
W_0(a)=4\epsilon\left[
\frac{63\pi\sigma^{12}}{256b\,a^{11}}
-\frac{3\pi\sigma^6}{8b\,a^5}
\right].
\]

The nonzero reciprocal-lattice modes describe registry corrugation. Since

\[
K_\nu(z)\sim\sqrt{\frac{\pi}{2z}}e^{-z}\quad(z\to\infty),
\]

the corrugation decays approximately as

\[
B_{p,m}(a)\propto e^{-2\pi ma/b}.
\]

This gives the direct physical coupling between normal opening and configurational/slip barriers: increasing `a` exponentially weakens the registry corrugation.

## Exact first derivatives

Define

\[
\theta_m=2\pi m\left(\frac12-\frac{s}{b}\right),
\qquad
z_m=\frac{2\pi ma}{b}.
\]

Then

\[
\frac{\partial S_p}{\partial s}
=
\sum_{m=1}^{\infty}
\frac{2\pi m}{b}B_{p,m}(a)\sin\theta_m.
\]

For the normal derivative,

\[
\frac{\partial S_p}{\partial a}
=A_p'(a)+\sum_{m=1}^{\infty}B_{p,m}'(a)\cos\theta_m,
\]

where

\[
A_p'(a)=\frac{1-2p}{a}A_p(a),
\]

and `B'` is evaluated analytically using

\[
\frac{dK_\nu(z)}{dz}
=-\frac12\left[K_{\nu-1}(z)+K_{\nu+1}(z)\right].
\]

Therefore

\[
W_a=4\epsilon\left[\sigma^{12}S_{6,a}-\sigma^6S_{3,a}\right],
\]

\[
W_s=4\epsilon\left[\sigma^{12}S_{6,s}-\sigma^6S_{3,s}\right].
\]

These are the conservative microscopic forces used by the probability dynamics:

\[
F_a=-W_a,\qquad F_s=-W_s.
\]

After external work is included in `G_N`, the same gradients enter the deterministic probability PDE,

\[
\partial_tP_N=\nabla_{\mathbf q}\cdot
\left[\mathbf M\left(P_N\nabla_{\mathbf q}G_N+k_BT\nabla_{\mathbf q}P_N\right)\right].
\]

Thus the Bessel representation is not a separate dynamical law. It is the analytic infinite-lattice energy kernel whose gradients generate the drift field of the Smoluchowski equation.

## N-cell many-body energy

For an N-cell representative correlated region,

\[
\mathcal U_N(\mathbf a,\mathbf s)
=
\sum_{i=1}^{N}W(a_i,s_i)
+
\sum_{i<j}\phi_{\rm LJ}(R_{ij}),
\]

with

\[
R_{ij}=\sqrt{[(j-i)b+s_j-s_i]^2+(a_j-a_i)^2}.
\]

The first term contains each upper cell's complete interaction with the infinite lower row. The second term contains interactions among upper cells. This is not a double count of the old normal-chain and two-row potentials; the analytic normal-chain formula is retained only as a restricted benchmark/projection.

## Numerical implementation

`lattice_bessel.py` evaluates the reciprocal-lattice series until a declared numerical tolerance is reached. A finite `max_modes` exists only as a safety ceiling. `two_row_lj_direct_reference()` retains a large finite real-space sum exclusively for verification tests.

The production model in `model.py` uses the Bessel kernel for the lower-row interaction; the historical `lower_images` parameter remains only for API compatibility and no longer controls production lower-row physics.
