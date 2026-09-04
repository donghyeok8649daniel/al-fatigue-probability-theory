# Peak-hazard asymptotic verification

Status: **verification of a derived approximation to the active reduced closure; not an Al calibration.**

## Derived formula

For a sinusoidal reduced traction

$$
q(t)=q_m+q_a\sin(2\pi f t),
$$

and a rare event dominated by the tensile maximum, the one-cycle hazard has the leading approximation

$$
\mathcal H_c
\sim
\frac{\nu_p}{f\sqrt{2\pi Bq_a(\lambda_c-\lambda_p)}}
\exp(-Bg_p),
$$

where

$$
B=\frac{EA_ca_0}{k_BT},
$$

$\lambda_p$ is the stable spacing at peak traction, $g_p=\Delta\psi_c(\lambda_p)$, and $\nu_p=\sqrt{\phi''(\lambda_p)}/(2\pi t_0)$.

The derivation uses the exact identity

$$
\frac{d\Delta\psi_c}{dq}=-(\lambda_c-\lambda_s)
$$

and a local phase expansion around the tensile maximum. No global Taylor approximation of the retained potential is used.

## Direct quadrature comparison

At mean stress 100 MPa, 20 Hz and 300 K:

| $A_c/A_0$ | amplitude | direct $\mathcal H_c$ | peak approximation | relative error |
|---:|---:|---:|---:|---:|
| 40 | 50 MPa | $2.925\times10^{-3}$ | $2.486\times10^{-3}$ | -15.0% |
| 40 | 100 MPa | $4.675\times10^{-3}$ | $4.207\times10^{-3}$ | -10.0% |
| 40 | 150 MPa | $8.721\times10^{-3}$ | $8.167\times10^{-3}$ | -6.36% |
| 50 | 50 MPa | $1.169\times10^{-6}$ | $1.005\times10^{-6}$ | -14.0% |
| 50 | 100 MPa | $2.298\times10^{-6}$ | $2.118\times10^{-6}$ | -7.83% |
| 50 | 150 MPa | $5.369\times10^{-6}$ | $5.111\times10^{-6}$ | -4.80% |
| 60 | 50 MPa | $4.746\times10^{-10}$ | $4.148\times10^{-10}$ | -12.6% |
| 60 | 100 MPa | $1.161\times10^{-9}$ | $1.088\times10^{-9}$ | -6.26% |
| 60 | 150 MPa | $3.395\times10^{-9}$ | $3.265\times10^{-9}$ | -3.83% |

The leading peak formula is therefore already useful for order-of-magnitude and parameter-sensitivity analysis in the audited rare-event regime. Direct quadrature remains the authoritative numerical evaluation.

## Temperature-slope inversion

The peak formula implies

$$
\ln\left(\frac{f\mathcal H_c}{\sqrt T}\right)
\simeq \text{constant}
-\frac{EA_ca_0g_p}{k_B}\frac1T.
$$

A controlled numerical test generated direct one-cycle hazards from 260 to 340 K using an injected sensitivity value

$$
A_c/A_0=50.
$$

A linear fit of the transformed quantity against $1/T$ gave slope

$$
m_T\approx-11253.05\ \mathrm{K},
$$

while the peak formula predicts

$$
-11223.66\ \mathrm{K}.
$$

Using

$$
A_c
\simeq
-\frac{k_Bm_T}{Ea_0g_p}
$$

recovers

$$
\boxed{A_c/A_0\approx50.13}
$$

from the synthetic data, an error of about 0.26% relative to the injected value.

## Meaning

This result does **not** justify choosing $A_c/A_0=50$. It shows that, if the normal-instability hypothesis is correct and a local first-passage hazard can be measured, a temperature sweep supplies an identifiable route to the characteristic cohesive area instead of treating it as an invisible S-N fitting knob.

Specimen-level S-N data cannot be inserted directly into this inversion until the specimen correlation area/volume and local-to-specimen survival mapping are separately established.
