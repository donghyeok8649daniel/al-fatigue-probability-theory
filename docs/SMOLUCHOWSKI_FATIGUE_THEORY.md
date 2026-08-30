# 1D normal-tensile Smoluchowski fatigue theory

## Exact energies and interpretation

For the generalized pair potential

$$v(r)=\varepsilon_{\rm LJ}[(\sigma_{\rm LJ}/r)^m-(\sigma_{\rm LJ}/r)^n],
\qquad m>n>1,$$

the number of pairs separated by $ka$ in an $N$-atom uniform chain is $N-k$:

$$U_N=\sum_{k=1}^{N-1}(N-k)v(ka).$$

An interior atom sees both directions, but every unordered pair is shared by
two atoms. Hence

$$U_\infty=\lim_{N\to\infty}U_N/N
=\tfrac12\sum_{k\ne0}v(|k|a)=\sum_{k\ge1}v(ka)$$

and therefore, exactly and without a physical cutoff,

$$U_\infty(a)=\varepsilon_{\rm LJ}\{\zeta(m)(\sigma_{\rm LJ}/a)^m
-\zeta(n)(\sigma_{\rm LJ}/a)^n\}.$$

This is J per atom (or representative lattice cell), not the total energy of
an $N$-atom chain. The prefactor is the generalized-LJ convention used here.
The common 12--6 form $4\epsilon[(\sigma/r)^{12}-(\sigma/r)^6]$ uses a
different energy parameter: for that special shape only,
$\varepsilon_{\rm LJ}=4\epsilon$.

This $U_\infty$ is a homogeneous equation of state: every spacing changes to
$a$. For one gap $g$ between otherwise undeformed half-chains of spacing
$a_0$, collect cross-gap pairs by $\ell=i+j$; there are $\ell+1$ of them.
With $q=g/a_0$,

$$\sum_{\ell\ge0}(\ell+1)(q+\ell)^{-p}
=\sum_{\ell\ge0}[(q+\ell)+(1-q)](q+\ell)^{-p}
=\zeta(p-1,q)+(1-q)\zeta(p,q)=S_p(q),$$

so

$$U_{\rm gap}=\varepsilon_{\rm LJ}\left[(\sigma_{\rm LJ}/a_0)^mS_m(q)
-(\sigma_{\rm LJ}/a_0)^nS_n(q)\right].$$

Thus $U_\infty$ is exact for homogeneous spacing; $U_{\rm gap}$ is the direct
candidate for a local opening. Using $U_\infty(a)$ for a local random spacing
is a **local-homogeneous/mean-field assumption**. The active kinetic
demonstration retains that assumption for continuity, and never mixes the two
energies in one potential.

## Calibration, normalization, and instability

$U_\infty'(a_0)=0$ gives

$$a_0=\sigma_{\rm LJ}\left[{m\zeta(m)\over n\zeta(n)}\right]^{1/(m-n)}.$$

Writing $\lambda=a/a_0$ and using the equilibrium identity gives exactly

$$U_\infty=E_0\phi(\lambda),\qquad
\phi={\lambda^{-m}\over m(m-n)}-{\lambda^{-n}\over n(m-n)},$$

where (an arbitrary common additive constant $C$ may be restored)

$$E_0=\varepsilon_{\rm LJ}m(m-n)\zeta(m)(\sigma_{\rm LJ}/a_0)^m
=\varepsilon_{\rm LJ}n(m-n)\zeta(n)(\sigma_{\rm LJ}/a_0)^n.$$

Consequently $\phi'(1)=0$, $\phi''(1)=1$, and
$\varepsilon_{\rm LJ}=E_0/[m(m-n)\zeta(m)(\sigma_{\rm LJ}/a_0)^m]$.
The cell calibration is $E_0=EA_0a_0$ when
$E=(a_0/A_0)U''(a_0)$. Neither this relation nor the present theory fixes
$A_0$ without an orientation/coarse-graining definition.

Tangent stability is lost at

$$a_c=\sigma_{\rm LJ}\left[{m(m+1)\zeta(m)\over
n(n+1)\zeta(n)}\right]^{1/(m-n)},\quad
\lambda_c=\left({m+1\over n+1}\right)^{1/(m-n)}=1.1077715386\ldots.$$

## Controlled ensemble and four working observables

With $F=\sigma A_0$ (Pa = J m$^{-3}$, m$^2$, hence N),
$\Phi=U-F(a-a_0)$. For tensile $F$, $\Phi\to-\infty$ as $a\to\infty$;
even at $F=0$, $U\to0$. A global one-spacing Gibbs density on $(0,\infty)$
is therefore invalid. The present solver uses a declared finite intact domain
and its normalized metastable Gibbs density only as an initial/conditional
state. A fixed-total-length chain is an alternative canonical ensemble but
couples all spacings and is not the one-coordinate kinetic model.

The required observables are

$$\int Pda=1,\quad \bar a=\int aPda,\quad
\bar U=\int[U(a)-U(a_0)]Pda,$$
$$w_h(t)=\int_0^t\sigma\dot{\bar\lambda}\,ds,\qquad
H_k=\oint_k\sigma\,d\bar\lambda.$$

$U$ is a state function; $H_k$ is path-dependent work density and is never
added to $U$.

## Langevin reduction, current, and hysteresis

Assume only an isothermal eliminated bath, separation between fast velocity
relaxation and slow spacing evolution, constant friction $\gamma$, and white
noise on the resolved timescale. The underdamped bath equation loses its
inertial term after velocity relaxation. Fluctuation--dissipation then fixes
the noise amplitude (it is not fitted):

$$\gamma\dot a=-\partial_a\Phi+\sqrt{2\gamma k_BT}\,\xi,
\quad\langle\xi(t)\xi(t')\rangle=\delta(t-t').$$

The Itô forward equation (constant mobility makes Itô and Stratonovich
identical) is

$$\partial_tP=-\partial_aJ,\qquad
J=-M_a(\Phi_aP+k_BT P_a),\quad M_a=1/\gamma,
\quad D_a=k_BT/\gamma.$$

| quantity | meaning | SI unit |
|---|---|---|
| $a$ | layer spacing | m |
| $P$ | spacing density | m$^{-1}$ |
| $\Phi$ | cell effective energy | J |
| $\Phi_a$ | generalized force | N |
| $\gamma$ | spacing friction | N s m$^{-1}$ |
| $M_a$ | mobility | m N$^{-1}$ s$^{-1}$ |
| $D_a$ | spacing diffusivity | m$^2$ s$^{-1}$ |
| $J$ | probability current | s$^{-1}$ |

For $\sigma=\sigma_m+\sigma_a\sin(2\pi ft)$, the operator has relaxation
time $\tau_r$. At finite $f\tau_r$, $J$ cannot transport $P$ to the
instantaneous conditional equilibrium before the load changes. Loading and
unloading therefore have different distributions at the same stress. This
phase lag—not the diffusion term alone—is the reversible Markov hysteresis
mechanism. Very fast forcing moves almost no probability, intermediate rates
can maximize the loop, and as $f\tau_r\to0$ the conditional state is recovered
and $H_k\to0$.

For $\mathcal F[P]=\int P[\Phi+k_BT\ln P]da$, integration by parts yields

$$\dot{\mathcal F}=\int P\partial_t\Phi\,da
-\int{J^2\over M_aP}\,da-[\mu J]_{a_L}^{a_R},
\quad\mu=\Phi+k_BT(\ln P+1).$$

The nonnegative second integral is rate-dependent entropy production times
temperature. Boundary work/free-energy transport must be retained for escape.

## Escape, survival, hazard, and precursor

A reflecting Markov model may reach a periodic state and need not accumulate
damage. For $0<F<F_c$, $U'(a)=F$ has a stable root below $a_c$ and an unstable
barrier $a_b(F)>a_c$; they merge at $F_c=U'(a_c)$. The minimal irreversible
extension evolves unnormalised intact density $\rho$ with an absorbing (or
physically calibrated radiation) boundary:

$$S=\int_\Omega\rho da,\quad P_{\rm init}=1-S,\quad
h=J_{\rm out}/S.$$

The conditional precursor

$$p_{\rm tail}=S^{-1}\int_{a_c}^{a_b(t)}\rho da$$

is distinct from reflecting tail mass, outgoing flux, and cumulative
initiation. The code uses a declared fixed absorbing coordinate so a moving
barrier is not silently approximated; $a_b(t)$ is evaluated as a diagnostic.

## Exact moment balances

For $Q_q=\int_L^R q(a)\rho da$,

$$\dot Q_q=-[qJ]_L^R+\int_L^Rq'Jda.$$

With $D=M_ak_BT$ this gives

$$\dot S=J_L-J_R,$$
$$\dot Q_a=-[aJ]_L^R-M_a\int\Phi_a\rho da-D[\rho]_L^R,$$
$$\dot Q_{a^2}=-[a^2J]_L^R-2M_a\int a\Phi_a\rho da
-2D[a\rho]_L^R+2DS,$$
$$\dot Q_U=-[UJ]_L^R-M_a\int U_a\Phi_a\rho da
-D[U_a\rho]_L^R+D\int U_{aa}\rho da.$$

For any conditional moment $\langle q\rangle_c=Q_q/S$,
$\dot{\langle q\rangle}_c=(\dot Q_q-\langle q\rangle_c\dot S)/S$.
Thus
$\dot{\mathrm{Var}}_c=\dot{\langle a^2\rangle}_c
-2\langle a\rangle_c\dot{\langle a\rangle}_c$.
Reflecting boundaries set both currents to zero; absorbing boundaries do not.
Because $\Phi_a$ contains nonlinear inverse powers, these equations generate
unclosed higher/negative moments. No Gaussian or Taylor closure is active.

## Parameters and limitations

Physical inputs are $m,n,\varepsilon_{\rm LJ},\sigma_{\rm LJ},a_0,E,A_0,T$,
$M_a$ (or $\gamma$), $\sigma_m,\sigma_a,f$, and an escape-boundary/radiation
parameter. Mass is needed only for an inertial comparison. Correlation scales
$l_c,A_c,V_c$ remain symbolic and are not $A_0$. Numerical inputs—domain,
cells, timestep, tolerance, output interval and any visualization mesh—are
kept separate. FEM element count is not a count of independent samples.

References: NIST DLMF sections 25.2 and 25.11 (Riemann/Hurwitz zeta);
H. Risken, *The Fokker--Planck Equation*, 2nd ed.; C. W. Gardiner,
*Handbook of Stochastic Methods*; J. S. Chang and G. Cooper,
J. Comput. Phys. **6** (1970) 1--16; D. L. Scharfetter and H. K. Gummel,
IEEE Trans. Electron Devices **16** (1969) 64--77.
