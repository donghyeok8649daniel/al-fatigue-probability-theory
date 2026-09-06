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

## Scientific scope

This is a dimensionless mechanism solver, not a quantitatively calibrated
pure-aluminum fatigue-life predictor. Solver time is not laboratory seconds.
The N=3 tensor-train code remains an initial Gibbs compression prototype and
does not perform time-dependent 6D Smoluchowski evolution.
