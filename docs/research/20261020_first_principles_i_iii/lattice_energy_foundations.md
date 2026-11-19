# 01 — Lattice-energy foundations

Date: 20 October 2026. Scope: Part I of the critical-research directive.
The primary manuscript is the actual Papeeria **Revision 11**, not the older
paper on the solver branch. The source and implementation domains are recorded
in `source_provenance.json`. Everything numerical here is dimensionless and
illustrative; no aluminium parameter or fatigue life is inferred.

## 1. The microscopic object and its forces

For a finite, declared atomistic system, let

$$
H(R,P;t)=\tfrac12 P^T\mathsf m^{-1}P+V_{\rm pair}(R)+V_{\rm env}(R)
                                      +V_{\rm load}(R,t),
\qquad
V_{\rm pair}=\sum_{i<j}\phi(r_{ij}),
$$
$$
\phi(r)=4\epsilon_{\rm LJ}[(\sigma_{\rm LJ}/r)^{12}
                          -(\sigma_{\rm LJ}/r)^6].
$$

Here $\phi$ is a **pair** energy, not a configuration energy. Boundary
constraints, periodic-image conventions, pinned remote atoms and any
thermostat have to be specified separately. A prescribed displacement and a
dead applied traction are different experiments. No universal external-work
term follows just from choosing LJ interactions.

Writing $d_{ij}=r_i-r_j$, the pair contribution is

$$
\dot r_i=p_i/m_i,\qquad
\dot p_i=24\epsilon_{\rm LJ}\sum_{j\ne i}
 [2\sigma_{\rm LJ}^{12}r_{ij}^{-14}-\sigma_{\rm LJ}^{6}r_{ij}^{-8}]d_{ij}
 -\nabla_i(V_{\rm env}+V_{\rm load}).
$$

For derivative audits, one bond has Cartesian stiffness

$$
D_{ij}=\frac{\phi'(r)}r I+
  [\phi''(r)-\phi'(r)/r]\hat d\hat d^T,
\quad
\phi''(r)=24\epsilon_{\rm LJ}[26\sigma_{\rm LJ}^{12}r^{-14}
                                      -7\sigma_{\rm LJ}^{6}r^{-8}].
$$

Its contributions to the configuration Hessian are $+D_{ij}$ in the two
diagonal blocks and $-D_{ij}$ in the off-diagonal blocks. Energy gradients
are resisting derivatives; the internal force is their negative.

## 2. One atom and an infinite row

Put a selected upper atom at $(s,a)$ and lower sites at
$((j+1/2)b,0)$. Then, for $a,b>0$,

$$
S_p(a,s)=\sum_{j\in\mathbb Z}[a^2+((j+1/2)b-s)^2]^{-p},
\qquad W_{\rm row}=4\epsilon_{\rm LJ}(\sigma_{\rm LJ}^{12}S_6
                                             -\sigma_{\rm LJ}^{6}S_3).
$$

This counts each interaction with that selected atom once. It is neither the
whole infinite-row energy nor a crystal's cohesive energy. For two identical
periodic rows of length $N_b b$, their cross interaction is
$N_b W_{\rm row}$, with the appropriate relative registry and with edges
or periodic images treated consistently. The total for two infinite rows is
infinite; the finite object is the energy per repeat or per length.

### Independent Poisson derivation

Use the Fourier convention $\widehat f(k)=\int_{\mathbb R}f(x)e^{-ikx}dx$.
For $f_p(x)=(a^2+x^2)^{-p}$, insert its Gamma integral and perform the
Gaussian transform:

$$
\widehat f_p(k)=\frac{\sqrt\pi}{\Gamma(p)}
 \int_0^\infty t^{p-3/2}e^{-a^2t-k^2/(4t)}dt.
$$

The integral representation of $K_\nu$ gives

$$
\widehat f_p(k)=\frac{2\sqrt\pi}{\Gamma(p)}
  \left(\frac{|k|}{2a}\right)^{p-1/2}K_{p-1/2}(a|k|),\quad k\ne0,
$$
$$
\widehat f_p(0)=\frac{\sqrt\pi\Gamma(p-1/2)}{\Gamma(p)}a^{1-2p}.
$$

The integral identity used here is [NIST DLMF 10.32.10](https://dlmf.nist.gov/10.32.E10).
Applying Poisson summation at the **shifted** row gives

$$
S_p=\frac{\sqrt\pi\Gamma(p-1/2)}{b\Gamma(p)a^{2p-1}}
 +\frac{4\sqrt\pi}{b\Gamma(p)}\sum_{n=1}^\infty
 \left(\frac{\pi n}{ab}\right)^{p-1/2}
 K_{p-1/2}(2\pi na/b)
 \cos[2\pi n(1/2-s/b)].                                      \tag{L1}
$$

Thus LJ uses $K_{5/2}$ and $K_{11/2}$, including every nonzero mode.
The half-cell phase is essential. The zero reciprocal mode is the lateral
average of this row kernel; it is **not** a separately identified physical
bulk energy. Nonzero modes are registry corrugation, not automatically a
surface excess. Bulk/surface distinctions require the real-space cut count.

### Domains, convergence and derivatives

The direct series is absolutely convergent for $p>1/2$. On each compact
set $a\in[a_-,a_+]$, $a_->0$, and one registry period, it and its first
two differentiated series converge uniformly: the differentiated terms have
summable inverse-power bounds independent of the point in that compact set.
The reciprocal terms and each fixed derivative have a polynomial-in-$n$
factor times $e^{-2\pi n a_-/b}$, so their series converge uniformly as
well. Termwise force/Hessian differentiation is consequently valid there.
This is not a uniform statement as the noncollision lower bound disappears.
At $a=0$ and an aligned row site the interaction is singular.

With $x_j=(j+1/2)b-s$, direct differentiation gives

$$
S_{p,a}=-2paS_{p+1},\quad
S_{p,s}=2p\sum_jx_j(a^2+x_j^2)^{-p-1},
$$
$$
S_{p,aa}=-2pS_{p+1}+4p(p+1)a^2S_{p+2},
$$
$$
S_{p,as}=-4p(p+1)a\sum_jx_j(a^2+x_j^2)^{-p-2},
\quad
S_{p,ss}=-2pS_{p+1}+4p(p+1)\sum_jx_j^2(a^2+x_j^2)^{-p-2}.       \tag{L2}
$$

Multiply the respective $p=6,3$ channels by the LJ coefficients, and negate
the gradient for the force. Differentiating (L1) produces the same expressions.
For $|s|\le b/2$, a deliberately conservative direct-value tail bound is
$2b^{-2p}(N-1)^{1-2p}/(2p-1)$ after indices $-N,\ldots,N$.
It bounds a value channel, not every differentiated or LJ-subtracted channel
without the corresponding coefficients and derivative bounds.

## 3. An explicitly declared cleavage geometry

The following is a reconstructible **2D zigzag reference**, not a claim that
Rev.11 has selected a real aluminium cleavage plane. Define the uncut sites

$$
r_{j\ell}=(jb+\ell b/2,\ell h),\quad h>0.
$$

Keep $\ell\le0$ fixed and translate $\ell\ge1$ by $(-s,\delta a)$.
For a lower/upper pair of layers with $k=\ell_+-\ell_-$, the separation is
$kh+\delta a$ and lateral phase $kb/2-s$. The cross kernel per repeat is

$$
w_k(\delta a,s)=W_{\rm row}(kh+\delta a,s+(1-k)b/2).
$$

Odd and even layers therefore have different reference registries. Reusing
the same half-staggered row at every layer would silently change the crystal.
For each $k$, the pairs crossing the cut are
$(\ell_-,\ell_+)=(0,k),(-1,k-1),\ldots,(1-k,1)$, exactly $k$ pairs.
Hence the rigid excess per lateral repeat is

$$
W_\Gamma(\delta a,s)=\sum_{k=1}^\infty k
 \{W_{\rm row}(kh+\delta a,s+(1-k)b/2)
                  -W_{\rm row}(kh,(1-k)b/2)\}.                \tag{L3}
$$

Cross pairs are counted once: no additional $1/2$ or factor two applies.
For finite half-stacks with $L,U$ rows the count is instead

$$
c_k=\max\{0,\min(0,U-k)-\max(1-L,1-k)+1\}.
$$

Its sum is $LU$, and $c_k=k$ only before the finite edges are felt.

### Cleavage is not homogeneous dilation

In this cut path $\partial_{\delta a}(kh+\delta a)=1$. Therefore
$W_{\Gamma,\delta a}=\sum_k k\,w_{k,d}$ and
$W_{\Gamma,\delta a\delta a}=\sum_k k\,w_{k,dd}$.
The extra chain-rule factor $k$ of a different path $d_k=ka$ is absent.
The fixed reference values have zero derivatives; their gradients are not
subtracted. Slip and mixed derivatives follow the same declared phase map.

The row power kernel decays as $d^{1-2p}$. Its absolute weighted unseparated
sum is summable for $p>3/2$, which covers both LJ channels and suffices for
the separation endpoint. On compact $\delta a>-h$ sets, derivative sums
converge uniformly. A finite-opening difference has an improved algebraic
tail, but this improvement must not be confused with an endpoint bound.
The zero-mode layer sum is retained explicitly, for $r=2p-1$,

$$
\sum_{k\ge1}k[(kh+\delta a)^{-r}-(kh)^{-r}]
=h^{-r}\{\zeta(r-1,1+x)-x\zeta(r,1+x)-\zeta(r-1)\},\quad x=\delta a/h.
$$

This is an accelerator for the same infinite energy, not a zero-mode-only
model; the nonzero registry modes remain in (L3).

Dominated convergence, using the absolute weighted baseline bound for
$\delta a\ge0$, gives

$$
W_{\rm sep}=-\sum_{k\ge1}k\,W_{\rm row}(kh,(1-k)b/2).
$$

For fixed $b$, the bulk per-atom energy is the half-counted in-row sum plus
$\sum_{k\ge1}w(kh,\Delta_k^0)$. Its normal derivative is the cut traction
at zero opening. Thus a matching zero-normal-stress reference can have
$W_{\Gamma,\delta a}(0,0)=0$; it is not the old single-row minimum by fiat.
The resisting traction is $\partial W_\Gamma/\partial\delta a$ divided by
the **declared** cut-cell area, if an area exists. A 2D row example alone
does not define a 3D area or a surface energy in J/m².

## 4. A three-dimensional plane needs a two-dimensional transform

For a planar Bravais lattice $\Lambda$ with primitive area $A_\parallel$,
define $S_p^{(2)}(d,\Delta)=\sum_{R\in\Lambda}(d^2+|R+\Delta|^2)^{-p}$.
The same Gamma/Gaussian calculation in two integration dimensions gives

$$
\widehat f_p^{(2)}(G)=\frac{2\pi}{\Gamma(p)}
 \left(\frac{|G|}{2d}\right)^{p-1}K_{p-1}(d|G|),\quad G\ne0,
\qquad \widehat f_p^{(2)}(0)=\frac{\pi}{p-1}d^{2-2p},
$$
$$
S_p^{(2)}=A_\parallel^{-1}\sum_{G\in\Lambda^*}
                      \widehat f_p^{(2)}(G)e^{iG\cdot\Delta}.   \tag{L4}
$$

Here $p>1$, $d>0$; LJ uses $K_2,K_5$, not $K_{5/2},K_{11/2}$.
The reciprocal area prefactor and the entire zero mode are required.
The planar absolute cut baseline $\sum_k k\,S_p^{(2)}(kh,\Delta_k)$
converges for $p>2$; LJ satisfies that condition.

For a declared FCC(111) reference one may take triangular vectors
$a_1=b(1,0),a_2=b(1/2,\sqrt3/2)$, $A_\parallel=\sqrt3b^2/2$,
$b=a_{\rm lat}/\sqrt2$, $h=a_{\rm lat}/\sqrt3$, and ABC shift
$\tau=(a_1+a_2)/3$. The actual cut then uses the planar kernel at
$(kh+\delta a,k\tau+u)$ with count $k$. A scalar slip path and its full
translation period must also be declared. This is exactly where (L4), rather
than the 1D row transform, is needed.

## 5. Rigid, relaxed and nonlinear embedding energies

If separation creates equivalent surfaces in the **same** geometry, boundary
conditions and relaxation convention, then $W_{\rm sep}/A_\Gamma=2\gamma$.
Rigid cleavage and relaxed surface formation are distinct calculations.
At each fixed constraint the relaxed potential can be no greater than the
admissible rigid configuration's potential; this does not order two excess
curves after different relaxed reference energies have been subtracted.

For EAM, first form each site's complete density, then apply $F$:
$\Delta V_{\rm emb}=\sum_i[F(\rho_i(q))-F(\rho_i^{\rm ref})]$.
The site densities include all adopted neighbors. One cannot replace this by
pair-crossing multiplicities times $F$ of an individual layer density.
That extension requires its own counting, gauge, derivatives and convergence.

## 6. Executed checks and their limits

The isolated code in `research/first_principles_20261020/` imports no production
force routine for its direct sums. The saved `results/run_001/verification.json`
contains values, gradients, Hessians and truncations. In the row comparison
the maximum scaled channel error is $9.84\times10^{-15}$. A separate single
production-row comparison differs by $7.11\times10^{-15}$ in energy/gradient.
The planar comparison records its slower real-space tail and the
energy-channel analytic bound; it does not pretend to have the row's error.

For the new finite $4\times5$ half-stack count, 20 row pairs give
$0.4273971979607903$ by direct summation and
$0.42739719796078934$ by counted Bessel summation. The cut derivative is
$3.941084252594483$; the erroneous homogeneous chain factor gives
$5.116347430416690$. These are not alternative physical calibrations.
For the illustrative stress-free stack $h=0.761436142462363$, the cut
separation limit is $2.563620429191983$; the finite-opening approach and
derivative differences are recorded separately. A zero difference between
two double-precision truncations is not a proof of a zero omitted tail.

**Conclusion:** the Poisson--Bessel energy is retained and independently
reconstructed. What remains unestablished is its identification with the
finite-temperature effective free energy and kinetics of the actual material.
