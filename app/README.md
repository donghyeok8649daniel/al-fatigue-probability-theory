# Probability-PDE desktop application

Launch from the repository root:

```bash
py -3 -m pip install -r solver_v1/requirements.txt
py -3 -m app.desktop_ui
```

For an import and canonical-calibration smoke test that does not open a window:

```bash
py -3 -m app.desktop_ui --smoke
```

## Canonical data path

`app.solver_adapter` is the only solver-facing path used by the desktop UI:

```text
physical axial stress
  -> reduced stress sigma/E
  -> f* = kappa_axial sigma/E
  -> N=1 direct probability PDE
  -> strain, probability, and numerical-diagnostic fields
```

The relaxed calibration is

$$
\mathbf c = [1,\chi]^T,
\qquad
\kappa_{\mathrm{axial}}
= \frac{a_0}{\mathbf c^T H_0^{-1}\mathbf c}.
$$

The frozen-normal value $a_0W_{aa}$ is diagnostic only. The UI does not use raw
$\sigma/E$ or the frozen-normal scale as its generalized force.

The initial condition is the conditional principal-well Gibbs density at
$\sigma(t=0)$. For the current sinusoidal input this is the entered mean stress,
so a nonzero mean load does not create an artificial zero-preload jump.

## Result fields

The strain view displays four distinct PDE results:

- normal opening;
- intrawell registry;
- well-index plastic;
- total axial strain, equal to the sum of the other three.

Plastic strain is not forced to be nonzero. It appears only when the solved PDE
carries probability into a nonzero configurational well index.

A nonzero floating-point value is not, by itself, evidence of plasticity. The
solver also reports $P_n(t)$, signed and gross Scharfetter--Gummel fluxes across
$s=(n+1/2)b$, and a discrete well-population balance residual. A residual
plastic-strain claim additionally requires persistence after unloading and a
zero-stress hold. The convergence audit is documented in
[`solver_v1/PROBABILITY_AND_PLASTICITY_RESOLUTION.md`](../solver_v1/PROBABILITY_AND_PLASTICITY_RESOLUTION.md).

Probability views use the PDE outputs `survival`, `initiation_probability`, and
`first_passage_flux`. At the UI boundary the physical fields are defined by

```text
initiation_probability = cumulative_absorbed_mass
survival_probability   = 1 - cumulative_absorbed_mass
```

The raw intact mass, `1 - raw_intact_mass`, mass-balance residual, and
negative-mass correction remain explicitly labeled numerical diagnostics. A
roundoff-scale raw mass discrepancy is never presented as physical crack
initiation.

## Local and specimen probability

The local PDE probability and specimen aggregation are separate layers. The
optional UI inputs are the characteristic correlation area $A_c$ and the
effective stressed area $A_{\mathrm{stressed}}$, both in mm$^2$. They define

$$
N_{\mathrm{eff}}=\frac{A_{\mathrm{stressed}}}{A_c}.
$$

Under the explicitly declared independent-equivalent-region approximation,

$$
S_{\mathrm{spec}}=S_{\mathrm{local}}^{N_{\mathrm{eff}}},
\qquad
P_{\mathrm{init,spec}}=1-S_{\mathrm{spec}}.
$$

The implementation uses `log1p`/`expm1` arithmetic. $A_c$ is an external
statistical calibration input; no physical default is supplied. It never
enters `f* = kappa_axial sigma/E` and never scales strain or well activity.
Changing either area only reruns this cheap post-processing, not the PDE.

Specimen probability is withheld unless the local absorbed-mass signal has a
grid/time/integrator convergence certificate. The mathematical extrapolation
is retained as a separately named diagnostic, but an unresolved local signal
is never area-amplified and displayed as resolved physical probability.

## Language

The desktop defaults to Korean and can switch between **한국어** and **English**
from the header without restarting. `app.i18n` owns stable translation keys and
both locale dictionaries. Switching language redraws widget and plot text only:
entered values, solver settings, selected result, completed arrays, and saved
zoom/pan limits are retained, and the PDE is not rerun.

The calculation runs on a background worker thread. Completed PDE records are
sent to Tk through a queue so the GUI event loop remains responsive. Plot views
support cursor-centred wheel zoom, Shift+wheel X-only zoom, Ctrl+wheel Y-only
zoom, left-drag pan, double-click/Home reset, and per-result-field view memory.

The default frequency control is **cycles per model-time unit**, not hertz. The result
metadata and Solve summary expose local-pristine $\tau_{\mathrm{fast}}$,
$\tau_{\mathrm{slow}}$, $\omega\tau_{\mathrm{fast}}$,
$\omega\tau_{\mathrm{slow}}$, and small-signal transfer diagnostics. These are
diagnostics of the current finite-mobility model, not exact nonlinear response
values. See
[`solver_v1/DYNAMICS_FAST_A_DIAGNOSTICS.md`](../solver_v1/DYNAMICS_FAST_A_DIAGNOSTICS.md).

The default 21 by 31 explicit grid is labeled **Preview** because it does not
resolve the small Case-A normal displacement. **Resolved** selects an 81 by 91
grid and the same SG operator with backward-Euler time integration. The latter
is slower, remains a numerical-resolution mode rather than a different physical
model, and still requires ordinary grid/time convergence checks for publication
results.

The **Time basis** selector retains this model-time mode. Physical seconds and
hertz appear only when a validated kinetic calibration supplies both
reduced-coordinate mobilities and the common time scale. The committed Al
kinetic file is deliberately uncalibrated, so it cannot enable fake seconds.

The specimen panel separately displays local PDE probability, an explicitly
uncertified mathematical independent-region extrapolation, and certified
physical specimen probability. The last requires the background **Run
convergence check** workflow; Preview never certifies it. Completed solves also
include a per-cycle absorbed-mass/barrier/flux/survival table.

Editable load presets provide small-signal, moderate mechanism-probe,
all-compressive control, and extreme mechanism-stress-test inputs. They do not
clamp or reinterpret entered stresses. The UI displays $\sigma_{\min}$,
$\sigma_{\max}$, and their signed $\sigma/E$ values and marks multi-GPa inputs
as mechanism stress tests rather than calibrated pure-Al fatigue loads.

## Scientific scope

This is a dimensionless mechanism solver, not a quantitatively calibrated
pure-aluminum fatigue-life predictor. Solver time is not laboratory seconds.
The N=3 tensor-train code remains an initial Gibbs compression prototype and
does not perform time-dependent 6D Smoluchowski evolution.
