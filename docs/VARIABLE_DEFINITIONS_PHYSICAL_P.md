# Variable Definitions — Physical Statistical-Mechanical $P$

## Classification labels

- **EXACT / IDENTITY** — exact under the stated reduced Hamiltonian/calibration or stated equilibrium ensemble.
- **DEFINITION** — mathematical or thermodynamic definition.
- **PHYSICAL INPUT** — must be fixed independently from physical coarse graining or experiment, not fitted to obtain a desired fatigue life.
- **CONTROLLED APPROXIMATION** — valid only under the explicitly stated reservoir, quasistatic, or metastable time-scale assumptions.
- **DIMENSIONLESS MODEL DIAGNOSTIC** — result of the calibrated reduced model, not automatically a material prediction.

## Variables

| Symbol | Definition | Meaning | Unit | Classification |
|---|---|---|---|---|
| $a$ | physical layer spacing | normal separation of represented layers | m | DEFINITION |
| $a_0$ | equilibrium layer spacing | spacing used to normalize $a$ | m | PHYSICAL INPUT |
| $\lambda$ | $a/a_0$ | normalized layer spacing | dimensionless | DEFINITION |
| $A_0$ | representative layer-patch area | area associated with one coarse-grained layer interaction | m$^2$ | PHYSICAL INPUT |
| $E$ | Young's modulus used in the 1D normal calibration | normal elastic modulus | Pa | PHYSICAL INPUT |
| $E_0$ | $EA_0a_0$ | physical energy scale multiplying $\phi$ | J | EXACT / IDENTITY |
| $T$ | thermodynamic temperature | heat-bath temperature when canonical equilibrium is invoked | K | PHYSICAL INPUT |
| $k_B$ | Boltzmann constant | thermodynamic conversion constant | J/K | PHYSICAL CONSTANT |
| $\chi$ | $E_0/(k_BT)$ | inverse reduced temperature | dimensionless | EXACT / IDENTITY |
| $\sigma$ | applied normal stress | physical normal stress on represented patch | Pa | PHYSICAL INPUT |
| $F$ | $A_0\sigma$ | physical normal force on represented patch | N | DEFINITION |
| $f$ | $Fa_0/E_0=\sigma/E$ | dimensionless tensile force | dimensionless | EXACT / IDENTITY |
| $\phi(\lambda)$ | calibrated generalized-LJ reduced energy | full nonlinear layer potential | dimensionless | DEFINITION |
| $w_f(\lambda)$ | $\phi(\lambda)-f\lambda$ | force-biased reduced potential | dimensionless | DEFINITION |
| $f_c$ | $\max_{\lambda>0}\phi'(\lambda)$ | idealized critical dimensionless tensile force | dimensionless | EXACT / IDENTITY |
| $\lambda_c$ | $\phi''(\lambda_c)=0$ | tangent-instability spacing | dimensionless | EXACT / IDENTITY |
| $\lambda_s(f)$ | stable root of $\phi'(\lambda)=f$ | intact force-biased local minimum | dimensionless | EXACT / IDENTITY |
| $\lambda_b(f)$ | unstable root of $\phi'(\lambda)=f$ | opening barrier location for $0<f<f_c$ | dimensionless | EXACT / IDENTITY |
| $\Delta w(f)$ | $w_f(\lambda_b)-w_f(\lambda_s)$ | dimensionless metastable barrier | dimensionless | EXACT / IDENTITY |
| $E_{\rm tot}$ | total Hamiltonian energy | fixed energy in microcanonical ensemble | J | PHYSICAL INPUT |
| $L$ | $\sum_i\lambda_i$ | total normalized chain length | dimensionless | DEFINITION / CONSTRAINT |
| $d$ | number of independent quadratic momentum DOFs | exponent source after momentum integration | count | DEFINITION |
| $Z_M(L,\chi)$ | fixed-length canonical configurational partition function | equilibrium normalization at fixed $L,T$ | integration measure unit | DEFINITION |
| $f_{\rm th}$ | $-(1/\chi)\partial_L\ln Z_M$ | thermodynamic force conjugate to total normalized length | dimensionless | DEFINITION |
| $P_M^{\rm mc}$ | microcanonical one-spacing marginal | isolated equilibrium spacing distribution | inverse stretch | EXACT ENSEMBLE FORM |
| $P_M(\lambda\mid L,\chi)$ | $e^{-\chi\phi(\lambda)}Z_{M-1}(L-\lambda)/Z_M(L)$ | canonical fixed-length one-spacing marginal | inverse stretch | EXACT ENSEMBLE FORM |
| $P_{\rm ms}$ | Gibbs density conditioned on $0<\lambda<\lambda_b$ | intact-basin metastable local-equilibrium density | inverse stretch | CONTROLLED APPROXIMATION |
| $Q_c^{\rm ms}$ | $\int_{\lambda_c}^{\lambda_b}P_{\rm ms}d\lambda$ | instantaneous population beyond tangent-instability spacing but inside basin cutoff | dimensionless | CONTROLLED APPROXIMATION |

## Central physical forms

The calibration fixes

$$
\boxed{E_0=EA_0a_0},
$$

and therefore

$$
\boxed{\chi=\frac{EA_0a_0}{k_BT}}.
$$

For zero-temperature homogeneous quasistatic force control,

$$
\boxed{
P_{T=0,\mathrm{qs}}(\lambda\mid f)
=\delta[\lambda-\lambda_s(f)].
}
$$

For canonical fixed total normalized length,

$$
\boxed{
P_M(\lambda\mid L,\chi)
=
\frac{e^{-\chi\phi(\lambda)}Z_{M-1}(L-\lambda,\chi)}{Z_M(L,\chi)}.
}
$$

For a subcritical tensile metastable basin under the local-equilibrium assumption,

$$
\boxed{
P_{\rm ms}(\lambda\mid f,\chi)
\propto
\exp\{-\chi[\phi(\lambda)-f\lambda]\}
\mathbf 1_{0<\lambda<\lambda_b(f)}.
}
$$

A full-domain tensile Gibbs density is not normalizable for $f>0$ because $\phi(\lambda)-f\lambda\to-\infty$ as $\lambda\to\infty$.

---

# 한국어 번역 — 물리적 통계역학 $P$ 변수정의

## 분류 라벨

- **EXACT / IDENTITY** — 명시한 reduced Hamiltonian/calibration 또는 명시한 equilibrium ensemble 아래 정확.
- **DEFINITION** — 수학적 또는 thermodynamic 정의.
- **PHYSICAL INPUT** — 물리적 coarse graining 또는 실험으로 독립 결정해야 하며 원하는 fatigue life에 맞춰 fitting하면 안 됨.
- **CONTROLLED APPROXIMATION** — 명시한 reservoir, quasistatic, metastable time-scale assumption 아래에서만 유효.
- **DIMENSIONLESS MODEL DIAGNOSTIC** — calibration된 reduced model 결과이며 자동으로 material prediction은 아님.

## 변수

| 기호 | 정의 | 의미 | 단위 | 분류 |
|---|---|---|---|---|
| $a$ | physical layer spacing | represented layer의 normal separation | m | DEFINITION |
| $a_0$ | equilibrium layer spacing | $a$를 normalize하는 spacing | m | PHYSICAL INPUT |
| $\lambda$ | $a/a_0$ | normalized layer spacing | 무차원 | DEFINITION |
| $A_0$ | representative layer-patch area | coarse-grained layer interaction 하나에 대응하는 면적 | m$^2$ | PHYSICAL INPUT |
| $E$ | 1D normal calibration에 사용하는 Young's modulus | normal elastic modulus | Pa | PHYSICAL INPUT |
| $E_0$ | $EA_0a_0$ | $\phi$에 곱해지는 physical energy scale | J | EXACT / IDENTITY |
| $T$ | thermodynamic temperature | canonical equilibrium을 쓸 때의 heat-bath temperature | K | PHYSICAL INPUT |
| $k_B$ | Boltzmann constant | thermodynamic conversion constant | J/K | PHYSICAL CONSTANT |
| $\chi$ | $E_0/(k_BT)$ | inverse reduced temperature | 무차원 | EXACT / IDENTITY |
| $\sigma$ | applied normal stress | represented patch에 작용하는 physical normal stress | Pa | PHYSICAL INPUT |
| $F$ | $A_0\sigma$ | represented patch의 physical normal force | N | DEFINITION |
| $f$ | $Fa_0/E_0=\sigma/E$ | dimensionless tensile force | 무차원 | EXACT / IDENTITY |
| $\phi(\lambda)$ | calibrated generalized-LJ reduced energy | full nonlinear layer potential | 무차원 | DEFINITION |
| $w_f(\lambda)$ | $\phi(\lambda)-f\lambda$ | force-biased reduced potential | 무차원 | DEFINITION |
| $f_c$ | $\max_{\lambda>0}\phi'(\lambda)$ | idealized critical dimensionless tensile force | 무차원 | EXACT / IDENTITY |
| $\lambda_c$ | $\phi''(\lambda_c)=0$ | tangent-instability spacing | 무차원 | EXACT / IDENTITY |
| $\lambda_s(f)$ | $\phi'(\lambda)=f$의 stable root | intact force-biased local minimum | 무차원 | EXACT / IDENTITY |
| $\lambda_b(f)$ | $\phi'(\lambda)=f$의 unstable root | $0<f<f_c$에서 opening barrier 위치 | 무차원 | EXACT / IDENTITY |
| $\Delta w(f)$ | $w_f(\lambda_b)-w_f(\lambda_s)$ | dimensionless metastable barrier | 무차원 | EXACT / IDENTITY |
| $E_{\rm tot}$ | total Hamiltonian energy | microcanonical ensemble의 fixed energy | J | PHYSICAL INPUT |
| $L$ | $\sum_i\lambda_i$ | total normalized chain length | 무차원 | DEFINITION / CONSTRAINT |
| $d$ | independent quadratic momentum DOF 수 | momentum integration 뒤 exponent를 결정 | count | DEFINITION |
| $Z_M(L,\chi)$ | fixed-length canonical configurational partition function | fixed $L,T$ equilibrium normalization | integration measure 단위 | DEFINITION |
| $f_{\rm th}$ | $-(1/\chi)\partial_L\ln Z_M$ | total normalized length에 conjugate한 thermodynamic force | 무차원 | DEFINITION |
| $P_M^{\rm mc}$ | microcanonical one-spacing marginal | isolated equilibrium spacing distribution | inverse stretch | EXACT ENSEMBLE FORM |
| $P_M(\lambda\mid L,\chi)$ | $e^{-\chi\phi(\lambda)}Z_{M-1}(L-\lambda)/Z_M(L)$ | canonical fixed-length one-spacing marginal | inverse stretch | EXACT ENSEMBLE FORM |
| $P_{\rm ms}$ | $0<\lambda<\lambda_b$에 조건부인 Gibbs density | intact-basin metastable local-equilibrium density | inverse stretch | CONTROLLED APPROXIMATION |
| $Q_c^{\rm ms}$ | $\int_{\lambda_c}^{\lambda_b}P_{\rm ms}d\lambda$ | basin cutoff 안에서 tangent-instability spacing을 넘은 instantaneous population | 무차원 | CONTROLLED APPROXIMATION |

## 핵심 물리 함수형

calibration은

$$
\boxed{E_0=EA_0a_0}
$$

를 강제하고 따라서

$$
\boxed{\chi=\frac{EA_0a_0}{k_BT}}
$$

이다.

zero-temperature homogeneous quasistatic force control에서는

$$
\boxed{
P_{T=0,\mathrm{qs}}(\lambda\mid f)
=\delta[\lambda-\lambda_s(f)]
}
$$

이다.

canonical fixed total normalized length에서는

$$
\boxed{
P_M(\lambda\mid L,\chi)
=
\frac{e^{-\chi\phi(\lambda)}Z_{M-1}(L-\lambda,\chi)}{Z_M(L,\chi)}
}
$$

이다.

subcritical tensile metastable basin에서 local-equilibrium assumption을 쓰면

$$
\boxed{
P_{\rm ms}(\lambda\mid f,\chi)
\propto
\exp\{-\chi[\phi(\lambda)-f\lambda]\}
\mathbf 1_{0<\lambda<\lambda_b(f)}
}
$$

이다.

$f>0$에서는 $\phi(\lambda)-f\lambda\to-\infty$이므로 full-domain tensile Gibbs density는 normalizable하지 않다.
