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

Probability views use the PDE outputs `survival`, `initiation_probability`, and
`first_passage_flux` directly. The mass-diagnostics view also shows intact mass,
cumulative absorbed mass, mass-balance residual, and negative-mass correction.

The calculation runs on a background worker thread. Completed PDE records are
sent to Tk through a queue so the GUI event loop remains responsive. Plot views
support cursor-centred wheel zoom, Shift+wheel X-only zoom, Ctrl+wheel Y-only
zoom, left-drag pan, double-click/Home reset, and per-result-field view memory.

The frequency control is **cycles per model-time unit**, not hertz. The result
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

## Scientific scope

This is a dimensionless mechanism solver, not a quantitatively calibrated
pure-aluminum fatigue-life predictor. Solver time is not laboratory seconds.
The N=3 tensor-train code remains an initial Gibbs compression prototype and
does not perform time-dependent 6D Smoluchowski evolution.
