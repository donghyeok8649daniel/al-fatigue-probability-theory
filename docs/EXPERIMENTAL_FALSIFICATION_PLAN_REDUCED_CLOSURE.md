# Experimental falsification plan for the reduced P0-to-survival closure

Status: **VALIDATION PLAN — not calibration results.**

The final reduced equation is now mathematically closed under its declared assumptions. The next stage is not to add another state variable by default. It is to test the signatures that would distinguish the present normal-instability/thermal-first-passage hypothesis from a missing slow plastic or defect state.

## 1. Quantities that must be kept separate

The validation program must not mix the following scales:

1. **structural initial state** $P_0(\lambda)$;
2. **local first-passage survival** $S(t)$ and $F_{\mathrm{ci}}(t)$;
3. **characteristic cohesive event area** $A_c$ appearing in the activation energy;
4. **specimen-scale correlation area/volume**, which is required only when local survival is converted to specimen survival.

The present plan addresses items 1--3 first. Specimen-scale multiplicity is deliberately postponed.

## 2. Structural P0 before fatigue

The reduced law treats $P_0$ as structural/prestress spacing information, not as instantaneous thermal displacement. Candidate measurements include high-resolution diffraction or another independently justified residual microstrain mapping technique.

For each specimen/preparation state, record a reference loading phase and construct

$$
P_0(\lambda)
$$

from measured structural spacing or strain data without fitting a named PDF unless the measurement itself justifies that reduction.

The first falsification question is whether specimens with measurably different $P_0$ but otherwise matched loading conditions show the ordering of local initiation hazard predicted by the reduced model.

## 3. Frequency test: the strongest low-cost discriminator

In the strict quasistatic and fast-intrawell-equilibration regime, the active closure predicts

$$
\mathcal H_c\propto\frac1f,
$$

so a fixed local initiation probability should occur at approximately constant **physical time**, while the corresponding number of cycles scales approximately as

$$
N_p\propto f.
$$

A direct test should hold temperature, mean stress, stress amplitude, waveform, surface preparation, and reference $P_0$ as fixed as possible while changing frequency over a range in which macroscopic heating and rate-dependent plasticity remain negligible.

### Pass pattern

- $N_p/f$ approximately constant over the tested quasistatic range;
- no additional cycle-count dependence beyond the predicted phase-time exposure.

### Fail pattern

If $N_p$ is approximately frequency-independent while physical time changes strongly with $f$, the current thermal-only normal first-passage clock is missing a cycle-evolving structural state. In that case the next candidate should be a physically derived slip/dislocation/defect state, not an arbitrary modification of $P$.

## 4. Temperature test and Ac identifiability

At fixed stress waveform, the peak-dominated rare-event asymptotic predicts

$$
\ln\left(\frac{f\mathcal H_c}{\sqrt T}\right)
\approx
\mathrm{const}
-\frac{EA_ca_0g_p}{k_B}\frac1T.
$$

Therefore the transformed local hazard should be approximately linear in $1/T$ over a regime where the same mechanism remains active.

If the measured local hazard supports that linearity, the slope $m_T$ gives the characteristic cohesive area through

$$
A_c
\approx
-\frac{k_Bm_T}{Ea_0g_p}.
$$

This is the preferred future route to $A_c$ because it uses a predicted temperature slope rather than silently choosing $A_c$ to match one S-N point.

### Consistency test

Repeat the temperature-slope inversion at more than one mean/amplitude stress pair. A physically meaningful approximately constant $A_c$ should not change arbitrarily with the selected test condition inside the same mechanism regime.

If the inferred $A_c$ varies strongly and systematically, at least one of the following is wrong or incomplete:

- the operational boundary $\lambda_c$;
- the fast-equilibration/TST approximation;
- the assumption of one coherent area scale;
- the normal-instability mechanism itself.

## 5. Mean-stress and amplitude matrix

After an independently constrained $A_c$ and measured $P_0$ are available, run a small matrix of mean stress and stress amplitude.

The model predicts stress sensitivity through the mechanically derived stable spacing and energy climb,

$$
\frac{d\Delta\psi_c}{dq}=-(\lambda_c-\lambda_s)<0.
$$

Thus increased tensile traction lowers the local energy climb without introducing a fitted mean-stress correction law.

The model should predict the ranking and approximate magnitude of local hazard changes across this matrix using the same $A_c$.

## 6. Direct check of the operational boundary

The present initiation event is defined by first passage to

$$
\phi''(\lambda_c)=0.
$$

This is a model-based dividing surface, not yet an experimentally proven crack-initiation coordinate. A higher-fidelity calculation or a sufficiently resolved experiment should test whether reaching this local normal condition corresponds to irreversible initiation rather than frequent harmless recrossing.

If the operational boundary fails, it should be revised from physics or higher-fidelity evidence. The failure should not be hidden by changing $A_c$ alone.

## 7. Recommended order of experiments

1. Establish repeatable structural $P_0$ measurement for the chosen specimen preparation.
2. Run a frequency sweep at one moderate stress state and fixed temperature.
3. If the frequency signature is compatible, run a temperature sweep and infer a provisional $A_c$ from the transformed slope.
4. Repeat the temperature slope at a second stress state to test whether the inferred $A_c$ is stable.
5. With the same $A_c$, predict a mean-stress/amplitude matrix before testing it.
6. Only after local behavior is credible, determine specimen-scale correlation area/volume and local-to-specimen scaling.

This order minimizes the risk of fitting multiple unknown scales simultaneously.

## 8. Decision logic

### Outcome A: frequency and temperature signatures both pass

Continue with independent $P_0$ characterization, $A_c$ identification, and higher-fidelity validation of $\lambda_c$. The normal-instability submodel remains viable.

### Outcome B: temperature activation passes but frequency test is cycle-controlled

The rare-event barrier may still be relevant, but the assumption of immediate intrawell renewal is incomplete. Introduce a slow state only if its mechanics can be derived; registry/dislocation evolution becomes the leading research target.

### Outcome C: frequency signature passes but transformed temperature law fails

Revisit the TST/rare-event approximation, coherent-area model, or operational boundary before altering the structural transport law.

### Outcome D: both signatures fail

Reject the active reduced laboratory closure as the dominant fatigue-initiation mechanism for the tested regime. Retain the exact finite-chain probability projection as a reference framework, but move the physical mainline to a justified plastic/dislocation mechanism.

## 9. What is not allowed in the validation stage

- choosing $A_c$ separately for every stress level;
- replacing measured $P_0$ with a convenient named PDF without evidence;
- using FEM element area as $A_c$ or as a count of independent failure opportunities;
- multiplying local probabilities into specimen probabilities before a correlation scale is established;
- adding an empirical damage variable solely because one of the above tests fails.

The model should be allowed to fail cleanly. That is the purpose of the present reduced closure: its unknowns and falsification signatures are explicit.
