# Final reduced closure: general structural P0 verification

Status: **verification of the active reduced hypothesis; not an Al fatigue-life calibration.**

## Purpose

This audit tests the final reduced mapping

$$
P_0(\lambda)+\sigma(0:t)
\to
P_b(\lambda,t),\ S(t),\ F_{\rm ci}(t)
$$

for a nontrivial structural initial distribution without assigning a named PDF family. The synthetic $P_0$ is a seven-point empirical distribution used only to exercise the equations.

The diagnostic inputs are $m=12.19$, $n=6$, $E=69$ GPa, $T=300$ K, mean/amplitude stress $100/100$ MPa, $f=20$ Hz, and $A_c/A_0=50$. The area ratio is a sensitivity value, not a calibration.

## Mechanical closed-cycle check

For each initial structural spacing $\lambda_0$, the quasistatic map is

$$
\phi'[\Lambda(\lambda_0,t)]
=
\phi'(\lambda_0)+q(t)-q_{\rm ref}.
$$

Because the stress returns to the reference phase after one cycle, every represented spacing must return to its own $\lambda_0$ if the stable branch is retained. Numerically, the maximum closed-cycle error is

$$
2.22\times10^{-16}.
$$

Thus the mechanical part is reversible to numerical precision. No cycle-by-cycle permanent spacing drift was inserted.

## One-cycle survivor loss

The seven synthetic $P_0$ labels have integrated one-cycle hazards ranging from approximately

$$
3.87\times10^{-7}
\quad\text{to}\quad
1.28\times10^{-5}.
$$

The weighted one-cycle survival is

$$
S_1=0.99999591754,
$$

so the local first-passage fraction in one cycle is

$$
1-S_1=4.08246\times10^{-6}.
$$

This is the intended high-local-survival regime: the normalized structural mechanics can close reversibly while a very small survivor mass is removed by thermal first passage.

## Long-cycle survival and selection

Using the exact characteristic survival factor for each $\lambda_0$ gives

| cycles | survival $S_N$ | first-passage fraction $1-S_N$ | conditional mean initial spacing |
|---:|---:|---:|---:|
| 1 | 0.9999959 | $4.08\times10^{-6}$ | 1.0016320 |
| $10^4$ | 0.960556 | 0.039444 | 1.0016136 |
| $10^5$ | 0.699360 | 0.300640 | 1.0014719 |
| $3\times10^5$ | 0.415339 | 0.584661 | 1.0012677 |
| $10^6$ | 0.134255 | 0.865745 | 1.0009337 |

The conditional survivor population shifts toward the initially smaller/lower-hazard spacings. This is **selection of survivors**, not an imposed permanent deformation of every surviving spacing.

## Frequency signature of the strict quasistatic limit

At the reference structural spacing with the same stress waveform and $A_c/A_0=50$,

| frequency | integrated hazard/cycle | median cycles | median time |
|---:|---:|---:|---:|
| 1 Hz | $4.596\times10^{-5}$ | $1.508\times10^4$ | 15082 s |
| 10 Hz | $4.596\times10^{-6}$ | $1.508\times10^5$ | 15082 s |
| 20 Hz | $2.298\times10^{-6}$ | $3.016\times10^5$ | 15082 s |
| 100 Hz | $4.596\times10^{-7}$ | $1.508\times10^6$ | 15082 s |

Thus

$$
f\mathcal H_c=\text{constant}
$$

for a purely phase-controlled quasistatic trajectory, and the median failure time in seconds is frequency-independent while cycles to initiation scale linearly with $f$.

This is not presented as an established property of real Al fatigue. It is a strong falsifiable signature of the current normal-instability/fast-equilibration hypothesis. If dedicated experiments instead show a genuinely cycle-controlled life over the regime where the mechanical response remains quasistatic, an additional slow structural state is required.

## Temperature signature

The same diagnostic gives strong activation sensitivity:

| temperature | integrated hazard/cycle | median cycles |
|---:|---:|---:|
| 250 K | $1.16\times10^{-9}$ | $5.97\times10^8$ |
| 275 K | $7.27\times10^{-8}$ | $9.53\times10^6$ |
| 300 K | $2.30\times10^{-6}$ | $3.02\times10^5$ |
| 325 K | $4.29\times10^{-5}$ | $1.62\times10^4$ |
| 350 K | $5.28\times10^{-4}$ | $1.31\times10^3$ |
| 400 K | $3.16\times10^{-2}$ | 21.95 |

Again, these numbers depend exponentially on the uncalibrated characteristic area and are not life predictions. Their value is that the model now exposes a sharp experimental falsification test rather than hiding frequency or temperature dependence in a fitted fatigue law.

## Final interpretation

This audit supports the internal mathematical consistency of the final reduced closure:

$$
\boxed{
\text{reversible structural transport}
+
\text{rare absorbing first passage}
\Rightarrow
\text{cumulative local survival loss}
}
$$

without requiring a permanently drifting normalized $P$ or an empirical scalar damage variable.

What remains unresolved is physical identification of $A_c$, measurement/construction of the actual structural $P_0$, and experimental validation of the operational $\lambda_c$ boundary and the predicted temperature/frequency signatures.
