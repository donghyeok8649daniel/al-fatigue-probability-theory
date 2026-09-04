# Peak-hazard asymptotics and characteristic-area identifiability

Status: **DERIVED CONSEQUENCE OF THE ACTIVE REDUCED CLOSURE.**  
This note does not introduce a new fatigue law. It derives a high-barrier approximation to the one-cycle hazard already defined by `FINAL_REDUCED_P0_THERMAL_FIRST_PASSAGE_CLOSURE.md`.

## 1. Starting point

For one structural label, the active reduced rate is

$$
k_c(\lambda,T;A_c)
=\nu_s(\lambda)
\exp[-B\,g(\lambda)],
$$

where, locally in this derivation,

$$
B=\frac{EA_ca_0}{k_BT},
$$

$$
g(\lambda)=\Delta\psi_c(\lambda),
$$

and

$$
\nu_s(\lambda)=\frac{\sqrt{\phi''(\lambda)}}{2\pi t_0}.
$$

The one-cycle integrated hazard is

$$
\mathcal H_c
=\int_0^{T_f}k_c[\Lambda(t),T;A_c]dt.
$$

The following derivation assumes the same conditions as the active rate law: quasistatic stable-branch mechanics, fast intrawell thermal re-equilibration, and a rare-event regime. In addition, the cycle integral is assumed to be dominated by one smooth tensile maximum.

## 2. Exact stress derivative of the operational energy climb

Let

$$
q=\phi'(\lambda_s)
$$

be the reduced traction corresponding to the stable spacing $\lambda_s$. The effective-potential climb to the operational boundary is

$$
g(q)
=\left[\phi(\lambda_c)-q\lambda_c\right]
-\left[\phi(\lambda_s)-q\lambda_s\right].
$$

Because $\lambda_s$ is a stationary point of the effective potential,

$$
\phi'(\lambda_s)-q=0.
$$

Differentiating $g$ with respect to $q$ therefore gives exactly

$$
\boxed{
\frac{dg}{dq}
=-(\lambda_c-\lambda_s).
}
$$

Thus tensile traction always lowers the operational energy climb while the stable branch exists.

A second derivative gives

$$
\boxed{
\frac{d^2g}{dq^2}
=\frac{d\lambda_s}{dq}
=\frac{1}{\phi''(\lambda_s)}>0.
}
$$

These identities are consequences of the retained potential and do not require a Taylor approximation of the full potential. The local expansion below is only an asymptotic expansion of the already defined cycle integral near its dominant phase.

## 3. Sinusoidal tensile peak

Consider

$$
q(t)=q_m+q_a\sin(2\pi f t),
\qquad q_a>0.
$$

Let the unique tensile maximum occur at $t=t_p$, and define

$$
q_p=q_m+q_a,
$$

$$
\lambda_p=\lambda_s(q_p),
$$

$$
g_p=g(\lambda_p),
$$

$$
\nu_p=\nu_s(\lambda_p).
$$

Near the maximum,

$$
q(t)
=q_p-\frac12q_a(2\pi f)^2(t-t_p)^2+O[(t-t_p)^4].
$$

Using the exact derivative $dg/dq=-(\lambda_c-\lambda_s)$ at the peak,

$$
g[q(t)]
=g_p
+\frac12(\lambda_c-\lambda_p)q_a(2\pi f)^2(t-t_p)^2
+O[(t-t_p)^4].
$$

In the high-barrier limit, the phase interval contributing to the integral is narrow, so $\nu_s$ may be replaced by $\nu_p$ to leading order. Gaussian integration then yields

$$
\boxed{
\mathcal H_c
\sim
\frac{\nu_p}{f\sqrt{2\pi Bq_a(\lambda_c-\lambda_p)}}
\exp(-Bg_p).
}
$$

This is the peak-hazard asymptotic formula.

## 4. Derived cycle-life form

For a narrow structural initial state in the periodic regime,

$$
S_N\simeq\exp(-N\mathcal H_c).
$$

The cycle count corresponding to cumulative probability $p$ is therefore

$$
N_p
=\frac{-\ln(1-p)}{\mathcal H_c}.
$$

Using the peak approximation,

$$
\boxed{
N_p
\sim
\frac{-\ln(1-p)\,f}{\nu_p}
\sqrt{2\pi Bq_a(\lambda_c-\lambda_p)}
\exp(Bg_p).
}
$$

For the median,

$$
N_{50}
\sim
\frac{(\ln2)f}{\nu_p}
\sqrt{2\pi Bq_a(\lambda_c-\lambda_p)}
\exp(Bg_p).
$$

This has an S-N-like exponential stress sensitivity, but no Basquin exponent or empirical life distribution has been inserted. Stress enters through the mechanically derived $\lambda_p$ and $g_p$.

## 5. Frequency signature

The leading formula contains the explicit factor

$$
\mathcal H_c\propto\frac1f.
$$

Hence

$$
N_p\propto f
$$

and the physical time to a fixed probability,

$$
t_p=\frac{N_p}{f},
$$

is frequency-independent to leading order inside the strict quasistatic/fast-equilibration regime.

This is a strong falsification signature. It is not assumed to be a universal property of real Al fatigue.

## 6. Temperature transform and Ac inversion

Since

$$
B=\frac{EA_ca_0}{k_BT}=\frac{C}{T},
\qquad
C=\frac{EA_ca_0}{k_B},
$$

the peak formula can be written

$$
\mathcal H_c
\sim
\frac{\nu_p}{f}
\sqrt{\frac{T}{2\pi Cq_a(\lambda_c-\lambda_p)}}
\exp\left(-\frac{Cg_p}{T}\right).
$$

Therefore

$$
\boxed{
\ln\left(\frac{f\mathcal H_c}{\sqrt{T}}\right)
\sim
\ln\left[
\frac{\nu_p}{\sqrt{2\pi Cq_a(\lambda_c-\lambda_p)}}
\right]
-\frac{Cg_p}{T}.
}
$$

At fixed stress waveform, a plot of

$$
\ln\left(\frac{f\mathcal H_c}{\sqrt{T}}\right)
$$

against $1/T$ should therefore be approximately linear in the peak-dominated rare-event regime.

If a **local** hazard is independently measurable and the retained normal mechanism is valid, the slope $m_T$ gives

$$
m_T
\simeq
-\frac{EA_ca_0g_p}{k_B},
$$

so that

$$
\boxed{
A_c
\simeq
-\frac{k_Bm_T}{Ea_0g_p}.
}
$$

This is an identifiability relation, not a present calibration. Specimen-scale S-N data cannot be inserted here without first resolving the local-to-specimen correlation scaling.

## 7. Numerical check

For the current diagnostic parameters

- $m=12.19$, $n=6$;
- $E=69$ GPa;
- $100\pm100$ MPa sinusoidal loading;
- $f=20$ Hz;
- $T=300$ K;
- sensitivity value $A_c/A_0=50$;

the direct one-cycle quadrature gives approximately

$$
\mathcal H_c^{\mathrm{direct}}
=2.2979\times10^{-6},
$$

while the peak asymptotic formula gives

$$
\mathcal H_c^{\mathrm{peak}}
=2.1180\times10^{-6}.
$$

The leading-order error is about $-7.8\%$ at this point. Across the diagnostic sweep $A_c/A_0\in\{40,50,60\}$ and stress amplitudes 50--150 MPa, the leading peak approximation remains within roughly 4--15% of direct quadrature.

A second audit generates exact one-cycle hazards from 260--340 K for a synthetic $A_c/A_0=50$ case and applies the transformed-temperature slope above. The inferred area ratio is approximately

$$
A_c/A_0\approx50.13,
$$

showing that the inversion recovers the injected sensitivity value to about 0.3% in that controlled numerical test.

## 8. Interpretation

The final reduced model now has three useful levels:

$$
\boxed{
\text{full cycle quadrature}
\quad\to\quad
\text{peak asymptotic life formula}
\quad\to\quad
\text{temperature-slope identifiability of }A_c.
}
$$

The asymptotic formula is not used when the event is not rare, when several comparable peaks contribute, when the stable branch approaches deterministic loss, or when the fast-equilibration assumption fails.

The main experimental value of this result is that the reduced hypothesis predicts a specific joint dependence on stress, temperature, and frequency that can be falsified before specimen-scale probability products are introduced.
