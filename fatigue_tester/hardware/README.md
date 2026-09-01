# Fatigue Tester Hardware

The physical tester is intentionally separated from the probability solver.

## Hardware responsibility

MCU/tester side:

- cyclic axial loading;
- load-cell force feedback;
- displacement/strain sensing;
- DCPD and temperature acquisition;
- E-stop, travel limits and safe shutdown;
- deterministic timestamps and telemetry.

PC side:

- `P(a,s,t)` probability evolution;
- the four governing-equation outputs;
- plastic-slip and hysteresis analysis;
- fatigue crack-initiation inference.

## Sizing

For specimen area `A`, mean stress `sigma_m` and amplitude `sigma_a`, use

\[
F_{req}=A\max(|\sigma_m+\sigma_a|,|\sigma_m-\sigma_a|).
\]

Actuator, load cell, grips and frame must all be sized above the resulting force with a separately justified engineering margin. For a sinusoidal displacement amplitude `x_a`,

\[
v_{pk}\approx 2\pi f x_a.
\]

Static force rating alone is therefore not enough for actuator selection.

## BOM

- `bom/fatigue_tester_bom.csv` — Git-editable BOM and current low-cost price/link basis.
- The generated `.xlsx` version contains automatic purchase-quantity and budget formulas.

Prices are estimates for planning. Re-check product option, stock, shipping and compatibility before purchase. Mains-voltage servo wiring and the hardwired E-stop/power-interruption path require appropriate electrical safety practice.
