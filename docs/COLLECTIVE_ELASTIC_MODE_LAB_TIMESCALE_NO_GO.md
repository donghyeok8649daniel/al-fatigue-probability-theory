# Collective normal elastic modes do not provide the laboratory fatigue clock

Status: **ACTIVE CONSTRAINT / NO-GO RESULT**

Purpose: test the proposed intermediate route between the exact finite chain and the rejected independent local oscillator while preserving the hard requirement

$$
P_0(a)+\sigma(0:t)\longrightarrow P(a,t).
$$

## 1. Linearized collective normal mechanics

For the active normalized generalized-LJ chain, linearization about the stable homogeneous reference state gives the usual nearest-neighbor acoustic equation. In the long-wavelength limit,

$$
\eta_{\tau\tau}\simeq\eta_{\xi\xi}.
$$

The retained physical calibration is

$$
a_0=2.8627442948\times10^{-10}\ {\rm m},
$$

$$
t_0=5.55046\times10^{-14}\ {\rm s}.
$$

Therefore the corresponding acoustic speed is

$$
c_a=\frac{a_0}{t_0}\simeq5.158\times10^3\ {\rm m/s}.
$$

For a fixed-left / traction-right segment of physical length $L$, the lowest longitudinal elastic mode is, to leading quarter-wave order,

$$
f_1\simeq\frac{c_a}{4L}.
$$

This estimate concerns the collective **elastic** normal mode only. It does not include plasticity, defect motion, thermal activation, damping, or registry transition.

## 2. Laboratory-frequency scaling

At 20 Hz the segment length required for the lowest elastic mode itself to occur at the loading frequency is

$$
L_{20}\simeq\frac{c_a}{4(20\ {\rm Hz})}\simeq64.47\ {\rm m}.
$$

Equivalently, the represented number of active spacings would be approximately

$$
M_{20}\simeq2.25\times10^{11}.
$$

Representative values are:

- $L=1\ \mu$m: $f_1\simeq1.289\times10^9$ Hz;
- $L=0.1$ mm: $f_1\simeq1.289\times10^7$ Hz;
- $L=1$ mm: $f_1\simeq1.289\times10^6$ Hz;
- $L=10$ mm: $f_1\simeq1.289\times10^5$ Hz;
- $L=100$ mm: $f_1\simeq1.289\times10^4$ Hz.

Hence any small local or characteristic region lies even deeper in the adiabatic regime than a full laboratory specimen.

This result is independent of the later choice of characteristic area or volume. It does not calibrate those quantities.

## 3. Connection to the previous quasistatic result

The existing quasistatic numerical protocol already showed that, for the conservative homogeneous normal chain on the stable branch, decreasing $\omega M$ collapses the residual spacing variance toward the homogeneous quasistatic state.

The present physical scaling shows that ordinary fatigue loading frequencies place realistic small domains at extremely small loading-frequency / elastic-mode-frequency ratios. Therefore collective conservative normal elasticity does not rescue the missing slow fatigue timescale.

## 4. Consequence for the P0-to-P requirement

The failure of the collective elastic-mode route does **not** change the required model interface:

$$
P_0(a)+\sigma(0:t)\longrightarrow P(a,t).
$$

It constrains the generator of that map.

The generator cannot be only reversible conservative normal elasticity if the model is intended to produce slow cycle-by-cycle fatigue evolution under ordinary laboratory frequencies.

A physically distinct slow or rare internal mechanism is required. Candidate mechanisms may later include a justified plastic/registry transition, defect evolution, or a thermally activated escape process, but none is promoted here.

## 5. Return point

Because the collective elastic-mode route fails its laboratory-timescale test, return to the exact finite-chain / correlation-hierarchy checkpoint for the normal mechanics and keep the following as separate established pieces:

1. $P_0$ can be constructed from a physically measured or computed structural spacing/microstrain field without assuming a named PDF family.
2. $P_0$ alone is not an exact sufficient state for the full deterministic finite chain.
3. The independent local-traction oscillator gives a formal $P_0\to P(t)$ propagator but fails the laboratory-timescale test.
4. Collective conservative normal elastic modes also fail the laboratory-timescale test for realistic small domains.
5. Characteristic length/area/volume and specimen-scale survival multiplication remain later calibration tasks and are not used to repair the present local-evolution problem.
