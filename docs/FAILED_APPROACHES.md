# Failed / Rejected Approaches

Keep failed approaches here instead of silently deleting them. A failed model can still be useful as a null test.

## Reversible one-coordinate LJ model

Model:

\[
\sigma(t) \leftrightarrow a(t),\qquad U(a)\text{ single-valued and conservative.}
\]

Result:

\[
\oint \sigma\,d\epsilon = 0.
\]

Reason for rejection as a fatigue model: no internal irreversible state evolution. It remains useful as a reversible baseline/unit test.

## Prescribed Weibull-like spacing density

Earlier work prescribed a Weibull-type density, later multiplied by an oscillatory factor. This is not accepted as a foundational evolution law because the distribution was imposed rather than derived from mechanics.

## Arbitrary stochastic kernel / Kramers rates

Transition kernels, barrier-crossing rates, or damping constants must not be introduced solely to create hysteresis or fatigue accumulation. Such models may be used only after a derivation or controlled coarse-graining argument.

## Instantaneous tail probability as crack probability

\[
Q_c(t)=\int_{a_c}^{\infty}P(a,t)\,da
\]

is an instantaneous unstable fraction, not automatically a cumulative crack-initiation probability. A first-passage or absorbing-boundary formulation is required for initiation probability.
