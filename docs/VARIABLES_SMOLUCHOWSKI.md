# Variables and units

| Symbol | Definition | SI unit | Kind/status |
|---|---|---:|---|
| $m,n$ | generalized-LJ exponents ($12.19,6$ current calibration) | 1 | physical |
| $\varepsilon_{LJ}$ | generalized-potential prefactor (not standard LJ well depth) | J | physical, unresolved coarse graining |
| $\sigma_{LJ}$ | pair-potential length | m | physical |
| $a_0$ | equilibrium homogeneous spacing | m | derived from LJ parameters |
| $E_0$ | representative-cell energy scale | J | $EA_0a_0$ |
| $E$ | loading-axis Young modulus | Pa | physical |
| $A_0$ | mechanical representative layer area | m$^2$ | physical, not $A_c$ |
| $T$ | bath temperature | K | physical |
| $\gamma$ | spacing friction | N s m$^{-1}$ | physical, unresolved |
| $M_a$ | spacing mobility $1/\gamma$ | m N$^{-1}$ s$^{-1}$ | physical |
| $D_a$ | spacing diffusivity $k_BT/\gamma$ | m$^2$ s$^{-1}$ | derived, not fitted independently |
| $t_r=a_0^2/(M_aE_0)=\gamma a_0^2/E_0$ | spacing relaxation time used for reduced time | s | derived from unresolved mobility and calibrated energy/length |
| $\beta=E_0/(k_BT)$ | reduced inverse temperature | 1 | derived from physical inputs |
| $f_n=Fa_0/E_0=\sigma/E_{[hkl]}$ | reduced normal force when $E_0=E_{[hkl]}A_0a_0$ | 1 | derived; $A_0$ cancels |
| $T^*=T_c/t_r=1/(f_{\rm physical}t_r)$ | reduced loading period | 1 | derived; requires physical mobility |
| $\sigma_m,\sigma_a,f$ | mean stress, amplitude, frequency | Pa, Pa, Hz | loading |
| $F=\sigma A_0$ | representative normal force | N | derived |
| $P$ | normalized reflecting/conditional density | m$^{-1}$ | state |
| $\rho$ | unnormalized intact density | m$^{-1}$ | state |
| $J$ | spacing-space probability current | s$^{-1}$ | state |
| $S$ | intact survival probability | 1 | state |
| $h$ | initiation hazard | s$^{-1}$ | state |
| $T_c$ | physical loading period | s | loading |
| $\mathcal U(t_2,t_1)$ | linear intact-density evolution operator | 1 | derived operator |
| $\mathcal K=\mathcal U(T_c,0)$ | one-cycle absorbing evolution operator | 1 | derived operator |
| $q_0$ | unit-mass principal conditional intact density | m$^{-1}$ | derived state |
| $w_0$ | left Perron survival-propensity weight, $\int w_0q_0da=1$ | 1 | derived adjoint diagnostic |
| $r$ | principal multiplier; asymptotic survival ratio per cycle | 1 | derived, not fitted |
| $r_2$ | second cycle-operator eigenvalue; $|r_2|/r$ controls leading transient | 1 | derived numerical spectrum |
| $1-r$ | principal-mode escape probability per cycle | 1 | derived, not a damage parameter |
| $\mathcal H_c=-\log r$ | integrated hazard over one cycle | 1 | derived |
| $-\log(r)/T_c$ | mean physical hazard rate over a cycle | s$^{-1}$ | derived |
| $N_p=\log(1-p)/\log r$ | principal-mode cycle count to cumulative probability $p$ | cycle | derived; predictive only after calibration |
| $K_{\rm init}$ | first-passage cycle index; geometric with parameter $1-r$ in the principal mode | cycle | derived cycle-binned random variable |
| $l_c,A_c,V_c$ | correlation scales | m, m$^2$, m$^3$ | unresolved; distinct from $A_0$ |
| domain/grid/$\Delta t$ | finite-volume controls | mixed | numerical, not fitted physics |
| FEM mesh density | visualization/discretization control | 1 or m$^{-d}$ | numerical; not sample count |
