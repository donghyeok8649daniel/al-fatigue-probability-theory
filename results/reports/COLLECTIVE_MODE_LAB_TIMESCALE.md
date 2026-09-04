# Collective elastic-mode laboratory-timescale audit

## Question

Can the active 1D normal chain retain collective mechanics and, by its long-wavelength elastic modes alone, produce the slow laboratory-time evolution required for

$$
P_0(a)+\sigma(0:t)\longrightarrow P(a,t)
$$

without an additional slow internal mechanism?

## Linearized acoustic scaling

Near the stable reference state $\lambda=1$, the normalized nearest-neighbor chain has the long-wavelength acoustic limit

$$
\eta_{\tau\tau}\simeq\eta_{\xi\xi}.
$$

Using the retained calibration

$$
a_0=2.8627442948\times10^{-10}\ {\rm m},
$$

$$
t_0=5.55046\times10^{-14}\ {\rm s},
$$

the corresponding physical acoustic speed is

$$
c_a=\frac{a_0}{t_0}\simeq5.158\times10^3\ {\rm m/s}.
$$

For a fixed-left / traction-right 1D segment of length $L$, the lowest longitudinal elastic mode is approximately the quarter-wave mode

$$
f_1\simeq\frac{c_a}{4L}.
$$

This estimate is sufficient for the present order-of-magnitude falsification test.

## Required size for laboratory frequencies

To move the lowest elastic mode down to 20 Hz requires

$$
L\simeq\frac{c_a}{4(20\ {\rm Hz})}\simeq64.47\ {\rm m},
$$

or roughly

$$
M\simeq2.25\times10^{11}
$$

represented spacings.

By comparison:

| Segment length | Lowest fixed-free longitudinal mode |
|---:|---:|
| $1\ \mu$m | $1.289\times10^9$ Hz |
| $0.1$ mm | $1.289\times10^7$ Hz |
| $1$ mm | $1.289\times10^6$ Hz |
| $10$ mm | $1.289\times10^5$ Hz |
| $100$ mm | $1.289\times10^4$ Hz |

A characteristic region smaller than the specimen is therefore even farther from ordinary 1--100 Hz fatigue frequencies.

For a 10 mm segment at 20 Hz, the squared frequency ratio is only

$$
\left(\frac{20}{f_1}\right)^2\simeq2.41\times10^{-8}.
$$

For a 1 mm region it is approximately $2.41\times10^{-10}$, and for a $1\ \mu$m region approximately $2.41\times10^{-16}$.

## Relation to the existing quasistatic diagnosis

The existing quasistatic protocol already established that, on the stable homogeneous branch, reducing the protocol parameter $\omega M$ drives the conservative normal chain toward the homogeneous quasistatic state and collapses the residual spacing variance. The present physical scaling shows that realistic laboratory fatigue tests lie extremely deep in that adiabatic regime for any small local region.

Therefore moving from an independent atomic oscillator to a collective **elastic** normal mode does not solve the laboratory-time-scale problem.

## Verdict

**NO-GO for collective normal elasticity as the missing fatigue clock.**

Collective elastic modes remain important for transmitting stress and for the spatial stress field, but they do not by themselves produce a slow cycle-by-cycle probability evolution at laboratory fatigue frequencies. A separate slow or rare internal mechanism is still required.

This result does **not** determine a characteristic length, area, or volume. Those remain later calibration/scale-up quantities.

The modeling requirement remains

$$
P_0(a)+\sigma(0:t)\longrightarrow P(a,t),
$$

but the operator that produces irreversible or cumulative evolution cannot be supplied by conservative normal elastic dynamics alone.
