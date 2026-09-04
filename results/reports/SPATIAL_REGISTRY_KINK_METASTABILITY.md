# Spatial registry-kink metastability audit

**Classification:** candidate-only, self-consistent dimensionless mechanism diagnostic.

## 1. Row repeat is not fitted

Once the same generalized LJ pair potential is used for compatibility along the moving row, the repeat $b$ should not be left inconsistent with that potential. For an infinite 1D row,

$$
U_{\parallel}(b)
=C_{mn}\epsilon
\left[
(\sigma/b)^m\zeta(m)
-(\sigma/b)^n\zeta(n)
\right].
$$

Stationarity gives

$$
\frac{\sigma}{b}
=
\left[
\frac{n\zeta(n)}{m\zeta(m)}
\right]^{1/(m-n)}.
$$

For $m=12.19$, $n=6$,

$$
\sigma/b=0.8942468263,
\qquad
b/\sigma=1.1182594901.
$$

Thus the row-compatibility scale used in this audit is derived from the same pair potential rather than fitted.

## 2. Fast normal coordinate is eliminated quasistatically

At laboratory fatigue frequency, normal atomic motion is adiabatic. Therefore, for each registry position $s$ and normal generalized force $Q_a$, the local substrate used in the spatial registry problem is

$$
V_{\mathrm{eff}}(s;Q_a)
=
\min_{a\ \mathrm{stable}}
\left[U_0(a,s)-Q_a a\right].
$$

This is not an added phenomenological potential; it is the stable normal minimum of the existing $U_0(a,s)$ surface.

## 3. Spatial energy

For moving-row positions $x_j=jb+s_j$, the candidate energy is

$$
E_{\mathrm{rk}}
=
\sum_jV_{\mathrm{eff}}(s_j;Q_a)
+
\sum_{j<k}
\left[
v_{m,n}(|(k-j)b+s_k-s_j|)
-v_{m,n}((k-j)b)
\right].
$$

A 121-repeat row was initialized with a finite shifted patch and both endpoints fixed in the original registry well $s/b=0.5$. No central site or slipped fraction was constrained during the final minimization.

## 4. Metastable residual state

The final minimization did **not** collapse to the uniform well. It retained a finite kink-antikink / slipped-patch state.

At zero applied normal stress:

- dimensionless formation/local-minimum energy: $11.38755235$;
- diagnostic physical energy bridge: $0.92377947$ eV;
- lowest fixed-endpoint Hessian eigenvalue: $1.26\times10^{-4}>0$;
- next Hessian eigenvalue is of order $10^1$, so the near-zero mode is consistent with a soft defect-translation mode rather than an obvious negative instability.

This is evidence of a **local minimum**, not a transition-state calculation. The activation barrier remains unknown.

## 5. Residual normal-spacing feedback

The same relaxed core gives a local stable normal-spacing field $a_j=a_{\mathrm{eq}}(s_j)$.

Relative to the registry-well equilibrium, the maximum local preferred opening is

$$
\max_j\frac{a_j}{a_{\mathrm{well}}}
=1.143795.
$$

Thus the residual core produces a local normal-equilibrium opening of about **14.38%** on this candidate surface.

This is the first calculation in the project that simultaneously gives:

1. equivalent far-field registry wells;
2. a non-equivalent residual local structural state after a finite slip patch forms;
3. a direct mechanical change of the normal spacing field without inserting an empirical damage variable.

## 6. Tensile-stress sensitivity of the residual-state energy

Using the retained $EA_0/a_0$ stiffness bridge only as a diagnostic conversion, the local-minimum energy changes as follows:

| normal stress | residual-state energy |
|---:|---:|
| 0 MPa | 0.92378 eV |
| 50 MPa | 0.91735 eV |
| 100 MPa | 0.91089 eV |
| 150 MPa | 0.90440 eV |
| 200 MPa | 0.89788 eV |

Tension lowers the residual-state energy in the expected direction because registry-core states prefer larger normal separation.

These values are **formation/local-minimum energies, not activation barriers**. They must not yet be inserted into an Arrhenius fatigue-life formula.

## 7. Consequence for the model architecture

The candidate chain is now

$$
P_a(a,0)
\xrightarrow{\sigma(0:t)}
P_a(a,t)
\rightarrow
\text{registry-barrier accessibility}
\rightarrow
\text{kink-pair nucleation}
\rightarrow
\text{metastable }\{s_j\}
\rightarrow
\text{shifted local }a_{\mathrm{eq}}
\rightarrow
P_a^{\mathrm{next}}.
$$

The missing step is the actual transition probability. That requires the **critical kink-pair saddle / minimum-energy path**, not the old coherent $N\Delta G_s$ approximation.

## 8. What is still open

- critical kink-pair activation barrier $\Delta G_{\mathrm{kp}}(\sigma,T)$;
- cycle-resolved transition rate derived from that barrier;
- whether the residual core survives when the full coupled $a_j,s_j$ dynamics, rather than quasistatic elimination of $a_j$, is used;
- mapping to the active normal-chain crack threshold without mixing incompatible calibrations;
- characteristic area/volume and specimen scale-up, intentionally deferred.

The paper mainline is not changed by this diagnostic yet.
