# Plasticity and strain audit

## 1. Scope and active case

This audit changes no mobility, temperature, barrier, stress map, or energy
parameter. The representative load is the existing nonlinear mechanism probe
with $E=69$ GPa, mean stress 900 MPa, amplitude 2000 MPa, model frequency 25,
$M_a=1$, $M_s=0.05$, $k_BT=0.02$, and $\chi=0.2$. All dynamics below are in
model time, not seconds or hertz.

The desktop path before this audit constructed
`solver_v1.model.TwoRowLJ` directly. It did not use either analytic LJ--EAM
hybrid. The default is retained for backward compatibility, but the UI now
labels it explicitly and offers three distinguishable energy surfaces.

## 2. Canonical state and well indexing

The state is $q=(a,s)$ and

$$
s=bn+\xi,\qquad n\in\mathbb Z,\qquad -\frac b2\leq\xi<\frac b2.
$$

The half-open well is

$$
\Omega_n=\left[(n-\tfrac12)b,(n+\tfrac12)b\right).
$$

The implementation uses

$$
n=\left\lfloor\frac{s+b/2}{b}\right\rfloor,qquad \xi=s-bn,
$$

including for negative $s$ and exact half-well boundaries. `TwoRowLJ` now
delegates to the shared `registry_decomposition` implementation, eliminating a
duplicate convention.

## 3. Surviving-ensemble strain decomposition

Let

$$
S(t)=\int_{\Omega_{\mathrm{intact}}}P(a,s,t)\,da\,ds
$$

and

$$
\langle X\rangle_S=\frac{1}{S(t)}
\int_{\Omega_{\mathrm{intact}}}X(a,s)P(a,s,t)\,da\,ds.
$$

Production observables use this conditional survivor density after opening
absorption. They are

$$
\epsilon_a=\left\langle\frac{a-a_0}{a_0}\right\rangle_S,
$$

$$
\epsilon_\xi=\frac{\chi}{a_0}\langle\xi\rangle_S,qquad
\epsilon_p=\frac{\chi b}{a_0}\langle n\rangle_S,
$$

$$
\epsilon_{\mathrm{total}}=\epsilon_a+\epsilon_\xi+\epsilon_p.
$$

Scaling the intact density by an arbitrary survivor fraction leaves all four
conditional strains unchanged. Opened mass is not included or renormalized
back. Across the refinement suite the largest decomposition residual was
$6.16\times10^{-18}$.

For the finest aligned cyclic run, the ranges were

| Quantity | Minimum | Maximum |
|---|---:|---:|
| $\epsilon_{\mathrm{total}}$ | 0.00963985 | 0.0389568 |
| $\epsilon_a$ | 0.00763415 | 0.0369122 |
| $\epsilon_\xi$ | 0.00200131 | 0.00205030 |
| $\epsilon_p$ | 0 | $2.61235\times10^{-12}$ |

Thus the plastic curve is not rescaled or clipped; it is about ten orders of
magnitude below the total strain in this run.

## 4. Well occupation

Absolute intact well populations are

$$
P_n(t)=\int_{\Omega_n}\int P(a,s,t)\,da\,ds.
$$

Every represented $s$ cell belongs to exactly one well, so

$$
\sum_nP_n=M_{\mathrm{raw}}
$$

up to FV quadrature and conservative-solve error. For the finest 3-well run at
the end of three cycles,

$$
P_{-1}\simeq1.35\times10^{-16},\quad
P_{+1}\simeq1.00752\times10^{-11},
$$

and the remaining survivor mass is in $P_0$. The outside-central-well mass was
$1.00754\times10^{-11}$.

## 5. Exact interwell flux and well balance

Positive $\mathcal J_{n+1/2}$ points toward increasing $s$. With opening loss
$\dot A_n$,

$$
\dot P_n=\mathcal J_{n-1/2}-\mathcal J_{n+1/2}-\dot A_n.
$$

Scharfetter--Gummel one-way activities are retained separately:

$$
j_{\mathrm{net}}=j_+-j_-,qquad
j_{\mathrm{gross}}=j_++j_-.
$$

They are never replaced by $|j_{\mathrm{net}}|$. In the finest cyclic run,

| Integrated quantity | Value |
|---|---:|
| forward $\sum\int j_+dt$ | $1.02930\times10^{-11}$ |
| backward $\sum\int j_-dt$ | $2.17627\times10^{-13}$ |
| signed net | $1.00754\times10^{-11}$ |
| gross | $1.05107\times10^{-11}$ |

The maximum well-population balance residual was $3.33\times10^{-15}$ and the
registry-moment balance residual was $5.08\times10^{-26}$.

## 6. Plastic flow and selective opening

For the unnormalized registry moment $M_n=\sum_n nP_n$,

$$
\dot M_n=\sum_n\mathcal J_{n+1/2}-\sum_n n\dot A_n.
$$

The conditional plastic-strain rate is

$$
\dot\epsilon_p=\frac{\chi b}{a_0S}\sum_n\mathcal J_{n+1/2}
+\frac{\chi b}{a_0S}
\left[-\sum_n n\dot A_n+\langle n\rangle_S\sum_n\dot A_n\right].
$$

The first term is model plastic flow. The second is selective opening and is
not plastic flow. The finest run gave maximum absolute rates

$$
|\dot\epsilon_{p,\mathrm{flow}}|=3.60\times10^{-11},qquad
|\dot\epsilon_{p,\mathrm{open}}|=1.24\times10^{-14}.
$$

The absorbed registry moment was $3.22\times10^{-16}$. Selective opening does
not explain the observed signed registry transfer.

## 7. Spatial/time refinement

The original 31-cell $s$ grid misses each half-well interface by 0.03226. The
42- and 60-cell grids align them to roundoff. Selected results are:

| Grid | $dt$ | Integrator | max $|\epsilon_p|$ | net transfer | alignment error |
|---|---:|---|---:|---:|---:|
| 21x31 | $10^{-3}$ | explicit | $1.4685\times10^{-12}$ | $5.6642\times10^{-12}$ | 0.03226 |
| 21x31 | $10^{-3}$ | implicit | $1.4679\times10^{-12}$ | $5.6622\times10^{-12}$ | 0.03226 |
| 21x42 | $10^{-3}$ | implicit | $3.0589\times10^{-12}$ | $1.1799\times10^{-11}$ | $1.1\times10^{-16}$ |
| 31x42 | $5\times10^{-4}$ | implicit | $2.5872\times10^{-12}$ | $9.9787\times10^{-12}$ | $1.1\times10^{-16}$ |
| 31x60 | $5\times10^{-4}$ | implicit | $2.6450\times10^{-12}$ | $1.0202\times10^{-11}$ | $1.1\times10^{-16}$ |
| 41x60 | $2.5\times10^{-4}$ | implicit | $2.6123\times10^{-12}$ | $1.0075\times10^{-11}$ | $1.1\times10^{-16}$ |

Across all aligned B/E/F/G runs, the maximum-to-reference difference is about
17%; the two finest 60-cell registry grids F/G differ by about 1.2%. The
conservative observed $\epsilon_p$ resolution floor from the full aligned set
is $4.47\times10^{-13}$, giving a signal/floor ratio about 5.8. The correct
classification for this particular cyclic case is **numerically resolved but
extremely small model interwell flow**. The unaligned UI-preview value is not
independently reliable.

Opening absorbed mass ranged from $6.38\times10^{-12}$ to
$4.07\times10^{-14}$ across aligned grids and is not convergence-resolved.
Crack-probability and plasticity resolution therefore remain separate.

## 8. Unload and zero-stress hold

After one cycle, a ramp to zero, and a 2.0-model-time hold:

| Grid | $dt$ | $\epsilon_p$ at unload | $\epsilon_p$ after hold |
|---|---:|---:|---:|
| 21x42 | 0.002 | $1.037\times10^{-12}$ | $1.994\times10^{-11}$ |
| 31x60 | 0.001 | $9.114\times10^{-13}$ | $1.584\times10^{-11}$ |

The registry shift persists and continues evolving, but the two results differ
by about 20%, while opening loss and remaining intrawell relaxation are also
grid-sensitive. Consequently this is a **persistent residual candidate, not a
converged residual plastic deformation**. The minimal periodic model also has
no hardening law; long-time zero-load spreading between equivalent wells must
not be equated with experimental Al residual plasticity.

## 9. Barrier and event-ordering audit

At the 2900 MPa tensile peak, the unchanged TwoRowLJ mapping gives
$f^*=3.62681$. The bound-basin diagnostics are

$$
\Delta G_s^\dagger=0.0229361,
$$

$$
\Delta G_{\mathrm{open}}^\dagger=0.220088
\quad\text{at the configurational minimum},
$$

and $0.0997005$ at the configurational saddle. The configurational barrier is
lower, but the period 0.04 is much shorter than the pristine
$\tau_s\simeq0.397$; the existing diagnostic is
$\omega\tau_s\simeq62.4$. Only about $10^{-11}$ probability crosses registry
interfaces in three cycles.

At this load, tiny aligned-grid interwell flow is resolved before any
convergence-resolved opening loss. This is a model-time mechanism result for
the dimensionless reference LJ surface, not calibrated Al event ordering.

## 10. First-cycle transient

For the finest aligned run:

| Cycle | mean $\epsilon$ | amplitude | mean $\epsilon_a$ | mean $\epsilon_\xi$ | $\Delta\epsilon_p$ | gross activity | absorbed mass |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.025495 | 0.014292 | 0.023470 | 0.002025 | $9.00\times10^{-13}$ | $3.53\times10^{-12}$ | $8.79\times10^{-15}$ |
| 2 | 0.022028 | 0.013031 | 0.020003 | 0.002026 | $8.70\times10^{-13}$ | $3.50\times10^{-12}$ | $1.71\times10^{-14}$ |
| 3 | 0.021813 | 0.012954 | 0.019787 | 0.002026 | $8.42\times10^{-13}$ | $3.47\times10^{-12}$ | $1.48\times10^{-14}$ |

The mean and amplitude approach a periodic response while intrawell strain is
nearly periodic. This is consistent primarily with relaxation of the initial
Gibbs survivor distribution. The opening loss is tiny, does not show a dominant
first-cycle decrement, and is grid-unresolved; vulnerable-tail depletion is
not the main demonstrated cause in this run.

## 11. Active UI energy models

The selector is explicit:

1. `two_row_lj_reference`: `TwoRowLJ`, dimensionless reference, default;
2. `analytic_lj_eam_hypothetical`: `AnalyticLJEAM`, uncalibrated sensitivity set;
3. `al_target_best_feasible_hybrid`: `AnalyticLJEAM`, stored best-feasible and
   practically non-identifiable Al-target fit.

Every result records the ID, Python class, parameter source, calibration
status, energy/geometry scale, $a_0$, $b$, $\chi$, $k_BT$, and $\kappa$. The
best-feasible parameters were not refitted. At $\chi=0.2$, its existing stored
surface gives $a_0=0.8164957222$ and $\kappa=13.85505145$.

## 12. Physical-time status

The dimensionalization remains

$$
M_i^*=\frac{t_0E_0}{L_0^2}M_{i,\mathrm{phys}},qquad
t_{\mathrm{phys}}=t_0t^*,qquad
f_{\mathrm{Hz}}=f_{\mathrm{model}}/t_0.
$$

$M_{a,\mathrm{phys}}$, $M_{s,\mathrm{phys}}$, and $t_0$ remain unavailable.
The UI physical-time option stays disabled and all uncalibrated result metadata
uses `model time` and `cycles / model time`. No export-ready metadata reports
seconds or hertz in this state.

## 13. Final classification and limitations

The representative current-model run is:

- elastic-dominated;
- elastic plus intrawell anelastic/configurational response;
- numerically resolved, extremely small model plastic flow during cycling;
- unresolved with respect to residual plasticity after unload/hold;
- opening first passage below the demonstrated spatial-convergence level.

The equations and plotting scale are correct. The previous UI ambiguity was
the undisclosed use of the original TwoRowLJ energy model. None of these
model-time results is a quantitative pure-Al or physical-Hz prediction.
