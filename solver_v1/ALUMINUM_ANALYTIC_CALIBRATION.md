# Pure-Al target calibration and identifiability audit

## 1. Scope and status

This audit asks whether the analytic LJ--EAM hybrid can reproduce a small set
of pure-Al static atomistic targets after an explicit FCC-to-two-row mapping.
It does not fit fatigue life, mobility, model time, yield stress, or specimen
statistics. The LJ pair functional form remains fixed by design. The outcome
is a **best-feasible static fit**, not a unique calibrated Al potential,
because its local inverse problem is severely ill-conditioned.

## 2. Traceable reference data

The selected fit set uses the 0 K Mishin-1999 Al EAM geometry, cohesion, and
elastic constants, the Lu-2000 LDA-DFT direct $a/4\langle101\rangle$ GSF
point, and the Wei-2016 DFT clean Al(111) surface energy. This mixture is
deliberate and visible; it is not described as one internally consistent DFT
dataset. Temperatures, paths, sources, and held-out flags are recorded in
`data/aluminum_reference_targets.csv` and `ALUMINUM_CALIBRATION_SOURCES.md`.

The five reduced fit targets are

$$
\left(rac{a_0}{b},E_{\mathrm{coh}},W_{,aa},
\Delta W_s^\dagger,W_{\mathrm{sep}}\right)
=
(0.816496581,3.36,18.9690584,0.110825653,0.93980153),
$$

where energies are eV per atom or per reduced cell as labeled in the result
tables. The normalization scales are respectively 0.5%, 2%, 5%, 10%, and 10%
of the target. They are transparent model-discrepancy scales, not claimed
experimental standard deviations.

## 3. Reduced crystallographic mapping

For an ideal FCC(111) primitive surface cell,

$$
b=\frac{a_{\mathrm{lat}}}{\sqrt{2}},\qquad
h_{111}=\frac{a_{\mathrm{lat}}}{\sqrt{3}},\qquad
A_{\mathrm{atomic\,cell}}=\frac{\sqrt{3}}{4}a_{\mathrm{lat}}^2.
$$

At $a_{\mathrm{lat}}=4.05$ angstrom this gives $b=2.86378246$
angstrom, $h_{111}=2.33826859$ angstrom, and
$A_{\mathrm{atomic\,cell}}=7.10249084$ angstrom squared. Consequently,

$$
\frac{a_0}{b}=\sqrt{\frac{2}{3}}.
$$

The microscopic area converts an interfacial energy only:

$$
E_{\mathrm{cell}}[\mathrm{eV}]
=\frac{\gamma[\mathrm{J/m^2}],A_{\mathrm{atomic\,cell}}[\mathrm{m^2}]}
{1.602176634\times10^{-19}}.
$$

$A_{\mathrm{atomic\,cell}}$ is not the statistical correlation area $A_c$.
$A_c$ never enters the energy, stress mapping, GSF normalization, or this
calibration.

The symmetric model saddle $s=b/2$ is the direct $a/4\langle101\rangle$
displacement when the period is the full $a/2\langle101\rangle$ Burgers
translation. It therefore maps directly to 0.250 J/m squared from Lu et al.
The lower-energy $a/6\langle112\rangle$ path is held out because it is a
different path. The intrinsic fault and vacancy energy are not currently
mappable: the reduced path has no independent intrinsic-fault minimum and no
vacancy degree of freedom.

The fixed-lateral longitudinal curvature uses

$$
C_{111}=\frac{C_{11}+2C_{12}+4C_{44}}{3}=122\ \mathrm{GPa},
$$

$$
W_{,aa}^{\mathrm{target}}
=C_{111}\frac{A_{\mathrm{atomic\,cell}}b^2}{h_{111}}.
$$

Bulk cohesion is only an effective/coarse-grained target: its 3D coordination
is not identical to a two-row cell. The opening target assumes two equivalent
clean surfaces, $W_{\mathrm{sep}}=2\gamma_{111}A_{\mathrm{atomic\,cell}}$.

## 4. Parameterization and exact density gauge

The raw square-root family begins with

$$
\theta=(\epsilon_{\mathrm{LJ}},\sigma_{\mathrm{LJ}},
\beta_\varrho,r_e,\varrho_0,A,\varrho_{\mathrm{ref}}),
$$

$$
f_\varrho(r)=\varrho_0
\exp[-\beta_\varrho(r/r_e-1)],\qquad
F_0(\varrho)=-A\sqrt{\varrho/\varrho_{\mathrm{ref}}}.
$$

Writing $b=1$, the energy depends on the raw density scales only through

$$
d=\frac{\beta_\varrho b}{r_e},\qquad
\widehat A=A\sqrt{\frac{\varrho_0e^{\beta_\varrho}}
{\varrho_{\mathrm{ref}}}}.
$$

Thus five raw density/embedding quantities do not define five independent
directions. The exact gauge is fixed with $b=r_e=\varrho_0=1$ and
$\varrho_{\mathrm{ref}}=\varrho(a=\sqrt{2/3},s=0)$. The square-root fit then
has four effective parameters

$$
(\epsilon_{\mathrm{LJ}},\sigma_{\mathrm{LJ}}/b,d,\widehat A).
$$

## 5. Nondimensional residual and staged fit

Every fit minimizes

$$
L(\theta)=\sum_k r_k^2,\qquad
r_k=\frac{y_k(\theta)-y_k^{\mathrm{target}}}{s_k},
$$

in logarithmic positive parameters with deterministic bounded trust-region
least squares. Infeasible proposals lacking an equilibrium, a positive
definite pristine Hessian, or positive registry/opening barriers are rejected;
the energy or Hessian is never clipped.

The row order in the sensitivity analysis implements the stages: geometry and
cohesion; add curvature; add registry barrier; add separation. The stage
singular values are stored with the joint matrix in the machine-readable
results.

## 6. Square-root embedding result

The best reproducible square-root result is

$$
(\epsilon_{\mathrm{LJ}},\sigma_{\mathrm{LJ}}/b,d,\widehat A)
=(3.30092\times10^{-4},1.35,1.32265,3.15754).
$$

Its normalized residuals are

$$
(0.507,-4.667,-1.011,2.072,10.190),
$$

and $L=131.19$. The LJ energy nearly vanishes and $\sigma/b$ reaches its upper
bound; the work of separation is 1.89748 eV/cell instead of 0.93980. This is a
quantified structural failure, not a usable Al fit.

## 7. Minimum analytic extension

Only after that failure, one analytic term is added:

$$
F(x)=-A\sqrt{x}+B(x-1),\qquad x=\varrho/\varrho_{\mathrm{ref}},
$$

$$
F'(\varrho)=
-\frac{A}{2\varrho_{\mathrm{ref}}\sqrt{x}}
+\frac{B}{\varrho_{\mathrm{ref}}},qquad
F''(\varrho)=
\frac{A}{4\varrho_{\mathrm{ref}}^2x^{3/2}}.
$$

In unrestricted EAM, a density-linear term can be moved into the pair term by
an EAM gauge transformation. Here it is separately named only under the
declared convention that the pair functional form is LJ. A tested positive
quadratic term did not improve the baseline optimum and was not retained.

## 8. Best-feasible extended fit

The deterministic solution is

$$
(\epsilon_{\mathrm{LJ}},\sigma_{\mathrm{LJ}}/b,d,A,B)
=(0.155132673,0.852468869,1.695130172,6.15126590,3.10771759)\ \mathrm{eV},
$$

where the final two entries have energy units. It predicts

$$
\frac{a_0}{b}=0.816495722,quad
E_{\mathrm{coh}}=3.35985487\ \mathrm{eV/atom},quad
W_{,aa}=18.9862261\ \mathrm{eV},
$$

$$
\Delta W_s^\dagger=0.110372478\ \mathrm{eV/cell},\quad
W_{\mathrm{sep}}=0.941236070\ \mathrm{eV/cell}.
$$

The pristine minimum is stable: $\lambda_{\min}(H_0)=6.38824$ eV. With
$\chi=0$, its reduced relaxed axial scale is
$\kappa_{\mathrm{axial}}=15.50217$ eV. This is recomputed for this surface;
the historical LJ-only value remains unchanged.

## 9. Identifiability result

The normalized five-by-five Jacobian with respect to log parameters has
singular values

$$
(267.822, 156.773, 26.1984, 4.43142, 1.02810\times10^{-6}),
$$

so its numerical rank is five under machine-epsilon rank counting but its
condition number is $2.61\times10^8$. The near-null right-singular direction
is approximately

$$
(-0.2051, 0.000316, 0.97867, 0.00510, -0.01074).
$$

It is dominated by correlated changes in LJ energy and density decay. The
very good five-target residual does not identify a unique potential. More
independent gamma-surface and traction-separation points are needed before
confidence intervals have physical meaning. The complete matrix and singular
values are in `results/aluminum_calibration`.

Moving $\pm0.25$ in log-parameter space along this near-null direction gives
fits with $L=0.101$ and 0.173, still small compared with the square-root
failure. Yet at $f^*=0.3$ their opening barriers at the configurational saddle
span 0.183--0.358 eV while the configurational barriers span
0.0187--0.0219 eV. These representative sets are saved rather than hiding the
large extrapolative uncertainty.

## 10. Held-out validation

The direct-path prediction 0.24898 J/m squared is 11.2% above the held-out
0.224 J/m squared lowest $\langle112\rangle$ path; this is also a path
difference, not only error. The predicted surface energy 1.06162 J/m squared
is 22.0% above the held-out Mishin-EAM 0.870 value. The predicted equivalent
$C_{111}=122.11$ GPa is 2.30% above the held-out 100 K Mishin-MD value 119.37
GPa. Intrinsic fault energy and vacancy formation remain predictions that the
present state space cannot make, rather than being silently inserted into the
loss.

## 11. Static barrier hierarchy

No orientation is invented as part of fitting. For a limited mechanism probe,
the explicitly hypothetical maximum-Schmid projection $\chi=0.5$ is used. At
$k_BT=0.02$ eV, the finite-temperature bound-basin configurational barrier
and the normal opening barrier at that saddle are:

| $f^*$ | $\Delta\mathcal F_s^\dagger$ (eV) | $\Delta G_{\mathrm{open}}^\dagger$ (eV) | ratio |
|---:|---:|---:|---:|
| 0.05 | 0.09006 | 0.67443 | 0.1335 |
| 0.10 | 0.07112 | 0.56454 | 0.1260 |
| 0.20 | 0.04606 | 0.39449 | 0.1168 |
| 0.30 | 0.01873 | 0.27566 | 0.0679 |
| 0.40 | 0.00258 | 0.18854 | 0.0137 |

A distinct configurational saddle is lost between $f^*=0.4$ and 0.5. The
minimum registry-dependent opening threshold is 0.82088 and the principal
registry value is 0.93226. Thus a static slip-favored interval exists for this
projection. Barrier ordering alone does not prove event ordering.

For loading normal to the same (111) plane,
$(\mathbf e\cdot\mathbf m)(\mathbf e\cdot\mathbf n)=0$ and $\chi=0$; such a
load supplies no direct resolved shear bias in this single-slip reduction.

## 12. Limited deterministic PDE check

With $\chi=0.5$, $M_a=1$, $M_s=0.05$, $k_BT=0.02$ eV, a coarse but resolved
21-by-42 implicit calculation held at $f^*=0.45$ for 0.5 model time moves about
$3.4\times10^{-4}$ of intact mass outside the central well and gives
$\epsilon_p\approx8.7\times10^{-5}$. Refining 15-by-30, 21-by-42, and
27-by-54 changes $\epsilon_p$ from $8.82$, $8.74$, to $8.59$ times
$10^{-5}$, while outer-well populations converge more slowly.

After return to zero force, a short 0.5-model-time hold did not erase the
registry moment; it changed from $8.71\times10^{-5}$ to
$8.94\times10^{-5}$. That is not a residual-plasticity validation: the hold
is short relative to uncalibrated kinetics and the result shows sensitivity
to normal-grid/opening treatment (one run accumulated $1.27\times10^{-8}$
opening mass). The honest dynamic classification is **coupled/unresolved**,
not calibrated slip-first Al plasticity.

## 13. Three-model comparison

`model_comparison.csv` separates the original dimensionless `TwoRowLJ`, the
previous uncalibrated hybrid sensitivity set, the failed square-root fit, and
the best-feasible extension. The original LJ reference retains
$a_0=0.7713438268704838$. Its $\chi=0$ frozen-normal scale is
$94.7109672665$, while the separately verified default $\chi=0.2$ relaxed
scale remains $86.29296488740997$ exactly; calibration does not modify either
value. The comparison CSV uses $\chi=0$ throughout so the listed scales are
comparable.

## 14. Static energy versus probability dynamics

This task changes no probability equation. Substitution of the candidate
surface into

$$
\partial_tP=\nabla\cdot[\mathbf M(P\nabla G+k_BT\nabla P)]
$$

retains deterministic probability evolution, configurational well flux, and
opening absorbed-mass bookkeeping. $M_s$ and model-time frequency remain
uncalibrated. Static energy calibration cannot turn model time into seconds or
predict fatigue cycles.

## 15. Analytic and numerical status

The LJ and exponential-density infinite lattice sums and their derivatives are
analytic reciprocal series. The embedding functions and chain-rule Hessians
are analytic. Tolerance-truncated Bessel series, equilibrium/saddle roots,
bounded least squares, sensitivity differences, and PDE evolution are
controlled numerical operations. None is renamed a closed-form Al potential.

## 16. Required future data and limitations

The two-row model omits full FCC coordination, multiple slip systems, the
intrinsic-fault branch, lateral relaxation, defects, and a complete
traction-separation curve. A defensible next inverse problem needs a mutually
consistent 0 K atomistic dataset containing several GSF-path points, several
opening-curve points, and independent environment-sensitive observables. That
may show that even the extended analytic family is structurally insufficient.

Still uncalibrated are coordinate mobilities, physical seconds, experimental
fatigue life, the statistical correlation area $A_c$, specimen aggregation,
and experimental validation. No one of those can be inferred from this static
five-target fit.

## 17. Reproducibility

The mapping, deterministic calibration, sensitivity matrix, and tests are in
`aluminum_calibration.py` and `test_aluminum_calibration.py`. Machine-readable
targets, mappings, parameter sets, residuals, singular values, barrier table,
and model comparison are under `results/aluminum_calibration/`.
