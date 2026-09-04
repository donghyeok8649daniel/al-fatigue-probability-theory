# Activated Registry Rare-Event Diagnostic

## Classification

**CONDITIONAL CANDIDATE — PARTIAL SUCCESS, NOT ACTIVE THEORY**

This audit tests whether the already-defined registry barrier can supply a slow many-cycle event clock even though all conservative normal/registry modes remain atomic or elastic-fast relative to laboratory fatigue loading.

The diagnostic adds a finite-temperature activated-rate hypothesis. It does **not** add a fitted damping coefficient or mobility, but it does require a thermal bath, fast intrawell equilibration, rare escape, and a coherent patch of $N$ repeats.

## Current direct-sum calibration

For the current $m=12.19$, $n=6$ multilayer registry surface with $s_0/b=1/2$:

$$
a_{0,r}/b=0.9910707144,
$$

$$
U_{aa}=108.3424,\qquad U_{ss}=26.1868,\qquad U_{ss}/U_{aa}=0.2417043.
$$

Matching the normal curvature to the retained Al normal calibration gives a diagnostic energy scale

$$
E_r\approx0.0699044\ {m eV}
$$

per dimensionless energy unit and a reference registry barrier

$$
\Delta G_s(a_{0,r})\approx0.10635\ {m eV/repeat}.
$$

The harmonic registry attempt frequency is

$$
\nu_s\approx1.41\times10^{12}\ {m Hz}.
$$

## Conditional coherent-patch rate

For a coherent patch of $N$ repeats,

$$
\Delta G_{s,N}=N\Delta G_s,
$$

and the conditional diagnostic rate is

$$
k_N(a,T)=\nu_s\exp[-N\Delta G_s(a)/(k_BT)].
$$

$N$ is **not calibrated** here. It is only scanned to determine whether an exponential barrier factor can in principle bridge the THz attempt scale and fatigue-cycle scale.

## 300 K, 20 Hz, 100 +/- 100 MPa diagnostic

The normal coordinate is treated quasistatically through the stable root $\phi'(\lambda)=\sigma/E$. Over the 0--200 MPa cycle the mapped registry barrier changes only from about 0.10635 eV/repeat to 0.10120 eV/repeat.

| $N$ | hazard/cycle | survival/cycle | median cycles from this conditional rate |
|---:|---:|---:|---:|
| 7 | $4.956\times10^{-2}$ | 0.951646 | $1.40\times10^1$ |
| 8 | $9.265\times10^{-4}$ | 0.999074 | $7.48\times10^2$ |
| 9 | $1.739\times10^{-5}$ | 0.9999826 | $3.99\times10^4$ |
| 10 | $3.275\times10^{-7}$ | 0.99999967 | $2.12\times10^6$ |
| 11 | $6.192\times10^{-9}$ | 0.999999994 | $1.12\times10^8$ |

This is the first tested mechanism in the current project that can create a many-cycle rare-event clock **without forcing a conservative natural frequency down to laboratory Hz**. The local survival can remain extremely close to one per cycle, which is compatible with later spatial multiplicity / characteristic-domain scaling.

The numerical values are not a calibration of $N$ or fatigue life.

## Why this is only a partial success

The ideal registry potential is periodic:

$$
U_0(a,s+b)=U_0(a,s).
$$

A transition from well $z$ to $z+1$ therefore changes the unwrapped registry label but lands in an energetically equivalent ideal well. In the current ideal model, this does not automatically alter the normal interaction branch.

Therefore the mechanism can supply

- a slow registry-transition / first-passage clock;
- a high local survival probability per cycle;
- a natural exponential sensitivity to a coherent-event size $N$;

but it does **not** yet supply progressive evolution of the full normal spacing marginal $P_a(a,t)$.

Also, the same activated formula predicts some zero-load hopping. Equivalent-well hopping must not be called fatigue damage unless a completed transition is tied to a physical irreversible post-transition defect state.

## Verdict

Retain as a candidate rare-event/plastic first-passage mechanism. Do not promote it to the active paper or G3 yet.

The next question is whether the current mechanics can generate a **non-equivalent post-transition state** that feeds back on the normal energy and therefore on $P_a(a,t)$. If that cannot be derived, the active theory should keep irreversible plastic feedback open rather than inventing one.
