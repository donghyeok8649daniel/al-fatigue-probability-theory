# Secular Slip Reference Run

This file records the current nondimensional proof-of-principle reference values for the conservative nonlinear slip-bath model.

## Parameters

$$
M=m=k=k_c=b=1,
$$

$$
\Delta_\gamma=0.1,\qquad \omega=0.2,
$$

with a smooth two-cycle force ramp and no viscous damping.

The reference direct integration used 800 bath sites and velocity-Verlet integration with

$$
\Delta t=0.01.
$$

## Cycle-end states

### $F_a=0.34$

Last six cycle-end values:

$$
-0.0240,\,-0.0239,\,-0.0240,\,-0.0239,\,-0.0240,\,-0.0240.
$$

Interpretation: bounded intrawell response; no secular cycle drift.

### $F_a=0.40$

Last six cycle-end values:

$$
-1.9650,\,-1.9648,\,-1.9650,\,-1.9649,\,-1.9650,\,-1.9649.
$$

Interpretation: finite transient interwell relocation followed by a periodic state.

### $F_a=0.50$

Last six cycle-end values:

$$
-5.8529,\,-6.8542,\,-7.8523,\,-8.8538,\,-9.8519,\,-10.8534.
$$

The late-cycle increment is approximately

$$
\boxed{\Delta s_{\rm cycle}\approx-1.00.}
$$

Representative late-cycle work is approximately

$$
\boxed{\oint F\,ds\approx2.994}
$$

in the nondimensional units of this model.

## Energy balance

For the $F_a=0.50$ reference run, the relative error in

$$
E_{\rm int}(t)-E_{\rm int}(0)=\int_0^tF\dot s\,dt
$$

was approximately

$$
\boxed{1.8\times10^{-7}}.
$$

## Spacing-like distribution diagnostic

For the finite bath, define local relative-displacement samples

$$
q_0=s-u_1,
\qquad
q_j=u_{j+1}-u_j.
$$

In the $F_a=0.50$ run, the variance over these samples increased from approximately

$$
1.88\times10^{-3}
$$

at cycle 2 to

$$
3.53\times10^{-2}
$$

at cycle 11.

This demonstrates redistribution of deformation into unresolved lattice modes. **It must not yet be interpreted as a thermodynamic-limit fatigue prediction for $P(a,t)$**, because the statistic depends on the finite observation domain and contains propagating phonon strain as well as structural change.

## Interpretation boundary

The strongest valid claim from this reference run is:

$$
\boxed{
\text{conservative microscopic dynamics}
\;\Rightarrow\;
\text{hysteresis + inter-basin cycle-state change is possible}
}
$$

for a nonlinear periodic non-affine coordinate coupled to a lattice bath.

The run does **not** establish an Al S–N curve, a crack-initiation life, or a low-stress fatigue threshold.
