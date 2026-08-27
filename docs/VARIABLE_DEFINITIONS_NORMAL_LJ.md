# Variable Definitions — Active Normal-LJ Mainline

This document defines the variables used by the active normal-deformation theory, the normal generalized-Lennard-Jones chain simulation, the probability formulation, and the normal-opening crack-initiation formulation.

## Classification labels

- **EXACT / IDENTITY** — exact under the explicitly stated microscopic model.
- **DEFINITION** — chosen mathematical definition.
- **ASSUMPTION** — modeling assumption.
- **CONTROLLED APPROXIMATION** — simplification that must be checked.
- **EMPIRICAL INPUT** — measured or externally supplied material quantity.

## 1. Microscopic normal geometry

| Symbol | Definition | Physical meaning | Unit | Classification |
|---|---|---|---|---|
| $i,j,k$ | integer indices | atom, bond, or neighbor index | dimensionless | DEFINITION |
| $N$ | finite number of sampled spacings/atoms when used in $P_N$ | finite-system size | dimensionless | DEFINITION |
| $x_i(t)$ | atomic coordinate along the loading axis | position of atom $i$ in the reduced 1D normal chain | m physically; dimensionless in normalized simulation | DEFINITION |
| $a_i(t)$ | $x_{i+1}(t)-x_i(t)$ | local normal interatomic spacing | m | DEFINITION |
| $a$ | continuous spacing-space coordinate | argument of $P(a,t)$ | m | DEFINITION |
| $a_0$ | reference/equilibrium spacing | normal spacing of the reference lattice | m | EMPIRICAL INPUT or calibration output |
| $\lambda_i$ | $a_i/a_0$ | local normalized normal stretch | dimensionless | DEFINITION |
| $\lambda$ | $a/a_0$ | continuous normalized spacing coordinate | dimensionless | DEFINITION |
| $\lambda_c$ | solution of $\phi''(\lambda_c)=0$ | idealized local normal stability-loss stretch | dimensionless | EXACT under the stated normalized LJ model |
| $a_c$ | $a_0\lambda_c$ | critical physical normal spacing corresponding to $\lambda_c$ | m | DEFINITION |
| $R_i^{(k)}$ | $a_i+\cdots+a_{i+k-1}$ | distance from atom $i$ to its $k$-th neighbor in the reduced geometry | m | EXACT / IDENTITY under the chain geometry |
| $r$ | continuous pair-distance coordinate | argument of $v(r)$ or $P_k(r,t)$ | m | DEFINITION |

## 2. Normal-spacing probability variables

For finite $N$,

$$
P_N(a,t)=\frac1N\sum_{i=1}^{N}\delta\!\left(a-a_i(t)\right).
$$

The central thermodynamic-limit state is

$$
\boxed{P(a,t)=\lim_{N\to\infty}P_N(a,t)}
$$

when the limit exists.

| Symbol | Definition | Physical meaning | Unit | Classification |
|---|---|---|---|---|
| $P_N(a,t)$ | finite empirical spacing density | finite-system distribution of local normal spacings | m$^{-1}$ | DEFINITION |
| $P(a,t)$ | thermodynamic-limit spacing density | main distribution-valued structural state | m$^{-1}$ | DEFINITION |
| $P_k(r,t)$ | density of $R_i^{(k)}$ | exact $k$-th-neighbor distance distribution | m$^{-1}$ | DEFINITION |
| $P^{*k}$ | $k$-fold convolution of $P$ | approximation to $P_k$ only under adjacent-spacing independence | m$^{-1}$ | CONTROLLED APPROXIMATION when used as $P_k\approx P^{*k}$ |
| $\bar a$ | $\int aP(a,t)\,da$ | mean normal spacing | m | DEFINITION |
| $\operatorname{Var}(a)$ | $\int(a-\bar a)^2P(a,t)\,da$ | variance of normal spacing | m$^2$ | DEFINITION |
| $Q_c(t)$ | $\int_{a_c}^{\infty}P(a,t)\,da$ | instantaneous probability mass above the candidate normal instability spacing | dimensionless | DEFINITION |
| $\delta(\cdot)$ | Dirac delta | distribution used to define empirical state densities | inverse unit of its argument | DEFINITION |

Normalization is

$$
\boxed{\int P(a,t)\,da=1.}
$$

## 3. Exact spacing-space kinematics

| Symbol | Definition | Physical meaning | Unit | Classification |
|---|---|---|---|---|
| $v(a,t)$ | $\langle\dot a_i\mid a_i=a\rangle$ | conditional transport velocity in spacing space | m/s | DEFINITION |
| $c$ | $\dot a$ | local spacing velocity in phase space | m/s | DEFINITION |
| $F(a,c,t)$ | joint density of spacing and spacing velocity | phase-space lift of $P$ | s/m$^2$ | DEFINITION |
| $A(a,c,t)$ | $\langle\ddot a_i\mid a_i=a,\dot a_i=c\rangle$ | conditional normal-spacing acceleration | m/s$^2$ | DEFINITION |

Exact kinematic continuity equation:

$$
\boxed{\partial_tP+\partial_a(Pv)=0.}
$$

Exact phase-space transport form under the stated projected variables:

$$
\boxed{\partial_tF+\partial_a(cF)+\partial_c(AF)=0.}
$$

Moment identities include

$$
\frac{d}{dt}\langle a^n\rangle=n\langle a^{n-1}v\rangle,
$$

and

$$
\boxed{\frac{d}{dt}\operatorname{Var}(a)=2\operatorname{Cov}(a,v).}
$$

## 4. Normal loading variables

| Symbol | Definition | Physical meaning | Unit | Classification |
|---|---|---|---|---|
| $t$ | time | physical or nondimensional simulation time | s physically | DEFINITION |
| $\sigma_n(t)$ | applied normal stress | externally imposed axial/normal cyclic stress | Pa | DEFINITION / experimental input history |
| $\sigma_m$ | mean normal stress | cycle mean | Pa | DEFINITION |
| $\sigma_a$ | normal stress amplitude | cyclic stress amplitude | Pa | DEFINITION |
| $\epsilon_n(t)$ | normal strain | strain work-conjugate to $\sigma_n$ under the chosen convention | dimensionless | DEFINITION |
| $f$ | cyclic frequency | cycles per second | Hz | DEFINITION |
| $\omega$ | $2\pi f$ | angular frequency | rad/s | DEFINITION |
| $T$ | $2\pi/\omega$ | loading period | s | DEFINITION |
| $N_{\rm cyc}$ | completed loading-cycle count | fatigue-cycle index | dimensionless | DEFINITION |
| $F_{\rm ext}(t)$ | externally applied normal force in a finite chain/tester | force conjugate to end displacement | N physically | DEFINITION |
| $F_m$ | mean normal force | mean of $F_{\rm ext}$ | N | DEFINITION |
| $F_a$ | normal-force amplitude | amplitude of $F_{\rm ext}$ | N physically | DEFINITION |

Sinusoidal normal stress:

$$
\boxed{\sigma_n(t)=\sigma_m+\sigma_a\sin(\omega t).}
$$

## 5. Generalized Lennard-Jones variables

The active physical pair-potential baseline is

$$
v(r)=\varepsilon_{\rm LJ}
\left[
\left(\frac{\sigma_{\rm LJ}}{r}\right)^m
-
\left(\frac{\sigma_{\rm LJ}}{r}\right)^n
\right].
$$

| Symbol | Definition | Physical meaning | Unit | Classification |
|---|---|---|---|---|
| $v(r)$ | generalized LJ pair energy | microscopic pair interaction energy | J | ASSUMPTION when adopted as material potential |
| $\varepsilon_{\rm LJ}$ | LJ energy scale | pair-energy parameter | J or eV | calibration input/output |
| $\sigma_{\rm LJ}$ | LJ length scale | characteristic distance parameter | m | calibration input/output |
| $m$ | repulsive exponent | short-range repulsive shape parameter | dimensionless | calibration output/currently 12.19 |
| $n$ | attractive exponent | attractive shape parameter | dimensionless | fixed/currently 6 in the active baseline |
| $\phi(\lambda)$ | normalized generalized-LJ energy used by the current 1D code | dimensionless pair energy | dimensionless | DEFINITION |
| $\phi'(\lambda)$ | $d\phi/d\lambda$ | normalized tensile force coordinate | dimensionless | EXACT derivative |
| $\phi''(\lambda)$ | $d^2\phi/d\lambda^2$ | normalized tangent stiffness | dimensionless | EXACT derivative |

The active normalized form is

$$
\boxed{
\phi(\lambda)=
\frac{\lambda^{-m}}{m(m-n)}
-
\frac{\lambda^{-n}}{n(m-n)}
}
$$

with

$$
\phi'(1)=0,
\qquad
\phi''(1)=1.
$$

Current exponents:

$$
m=12.19,
\qquad
n=6.
$$

The idealized normal stability condition is

$$
\boxed{\phi''(\lambda_c)=0}
$$

and therefore

$$
\boxed{
\lambda_c=
\left(\frac{m+1}{n+1}\right)^{1/(m-n)}.
}
$$

For the current exponents,

$$
\lambda_c\approx1.1077715386.
$$

## 6. Chain-model variables used in `theory/normal_lj_chain.py`

| Code field / symbol | Meaning | Unit in current code | Classification |
|---|---|---|---|
| `repulsive_exponent` / $m$ | generalized-LJ repulsive exponent | dimensionless | model parameter |
| `attractive_exponent` / $n$ | generalized-LJ attractive exponent | dimensionless | model parameter |
| `mean_force` | mean dimensionless end force | dimensionless | input |
| `force_amplitude` | amplitude of dimensionless end force | dimensionless | input |
| `omega` / $\omega^*$ in the normalized simulation | dimensionless angular loading frequency | dimensionless | input |
| `ramp_cycles` | number of cycles used for smooth loading ramp | cycles | numerical protocol parameter |
| `atoms` | number of atoms in finite chain | dimensionless count | CONTROLLED APPROXIMATION parameter |
| `dt` | velocity-Verlet time step | dimensionless time | numerical parameter |
| `cycles` | maximum number of cycles integrated | cycles | numerical parameter |
| `record_stride` | integration steps between stored samples | steps | data-storage parameter |
| `runaway_spacing` | spacing at which integration is stopped after large separation | dimensionless stretch | numerical stop only, **not** crack criterion |
| `period` | $2\pi/\omega$ in normalized simulation | dimensionless time | DEFINITION |
| `time` | stored simulation times | dimensionless time | output |
| `force` | stored external dimensionless force | dimensionless | output |
| `max_spacing` | maximum instantaneous $\lambda_i$ | dimensionless | diagnostic output |
| `internal_energy` | kinetic plus normalized LJ configurational energy | dimensionless | output |
| `external_work` | integrated external power | dimensionless | output |
| `cycle_mean_spacing` | mean $\lambda_i$ at cycle endpoints | dimensionless | output |
| `cycle_variance_spacing` | variance of $\lambda_i$ at cycle endpoints | dimensionless$^2$ | output |
| `cycle_max_spacing` | maximum $\lambda_i$ at cycle endpoints | dimensionless | output |
| `cycle_min_spacing` | minimum $\lambda_i$ at cycle endpoints | dimensionless | output |
| `cycle_snapshots` | arrays of all $\lambda_i$ at selected cycle endpoints | dimensionless | output |
| `first_instability` | first event for which $\max_i\lambda_i\ge\lambda_c$ | event record | diagnostic based on stated LJ criterion |
| `energy_balance_relative_error` | $|\Delta E_{\rm int}-W_{\rm ext}|/|W_{\rm ext}|$ with numerical protection near zero | dimensionless | numerical verification metric |

The current finite-chain potential is

$$
\boxed{V=\sum_i\phi(\lambda_i).}
$$

The leftmost atom is fixed and the prescribed normal force is applied to the rightmost atom.

## 7. Stress/force mapping variables

The active normalized-chain mapping uses

$$
f^*=\frac{\sigma_n}{E}
$$

under the current normalization.

| Symbol | Definition | Meaning | Unit | Classification |
|---|---|---|---|---|
| $f^*$ | $\sigma_n/E$ | dimensionless normal force/stress coordinate used by the reduced chain | dimensionless | DEFINITION under current mapping |
| $E$ | reference Young's modulus | stress scale | Pa | EMPIRICAL INPUT; current reference 69 GPa |
| $A_0$ | effective reference area for 1D-to-force mapping | cross-sectional conversion scale | m$^2$ | model calibration quantity |

The current 100 MPa reference amplitude is

$$
f_a^*=\frac{100\,\mathrm{MPa}}{69\,\mathrm{GPa}}
\approx1.44927536\times10^{-3}.
$$

## 8. Energy hierarchy variables

| Symbol | Definition | Meaning | Unit | Classification |
|---|---|---|---|---|
| $V$ | sum of pair energies in finite chain | finite-chain configurational energy | J physically or dimensionless after normalization | DEFINITION |
| $\mathcal U(t)$ | $\sum_k\int v(r)P_k(r,t)\,dr$ | distribution-level pair-potential energy per atom under the chosen counting convention | J | EXACT under the stated pair-potential hierarchy |
| $E_{\rm int}$ | kinetic plus configurational energy excluding external loading potential | internal energy | J or normalized | DEFINITION |
| $W_{\rm ext}$ | time integral of external power | external mechanical work | J or normalized | DEFINITION |
| $A_H$ | $\oint\sigma_n\,d\epsilon_n$ | normal stress-strain hysteresis work per volume per cycle | J/m$^3$ | DEFINITION |

For the exact pair-distance hierarchy,

$$
\boxed{
\mathcal U(t)=\sum_{k=1}^{\infty}\int_0^\infty v(r)P_k(r,t)\,dr.
}
$$

## 9. Time-scale variables

The current 1D mechanical scale is

$$
\boxed{
t_0=\sqrt{\frac{m_{\rm Al}a_0}{EA_0}}.
}
$$

| Symbol | Definition | Meaning | Unit | Classification |
|---|---|---|---|---|
| $t_0$ | equation above | atomic mechanical time scale in the current normalization | s | DEFINITION under the scaling |
| $m_{\rm Al}$ | atomic mass of Al | inertia scale | kg | EMPIRICAL INPUT |
| $\omega^*$ | $\omega t_0=2\pi f t_0$ | dimensionless angular frequency | dimensionless | DEFINITION |
| $f_{\rm phys}$ | $\omega^*/(2\pi t_0)$ | physical frequency corresponding to normalized $\omega^*$ | Hz | DEFINITION |

A laboratory-frequency statement must explicitly map $f$ and $\omega^*$ through $t_0$.

## 10. Normal-opening initiation variables

| Symbol | Definition | Meaning | Unit | Classification |
|---|---|---|---|---|
| $\tau_c$ | first time a stated mechanical normal-instability event occurs | first-passage initiation time | s | DEFINITION |
| $S(t)$ | $\Pr(\tau_c>t)$ | survival probability | dimensionless | DEFINITION |
| $F_{\rm ci}(t)$ | $\Pr(\tau_c\le t)=1-S(t)$ | cumulative crack-initiation probability | dimensionless | DEFINITION |
| $h(t)$ | $-\dot S/S$ when defined | initiation hazard | s$^{-1}$ | DEFINITION |

The active theory explicitly distinguishes the instantaneous tail quantity $Q_c(t)$ from the cumulative first-passage quantity $F_{\rm ci}(t)$.

---

# 한국어 번역 — 활성 Normal-LJ 변수정의

이 문서는 활성 수직변형 이론, normal generalized-Lennard-Jones chain simulation, 확률분포 정식화, normal-opening crack-initiation 정식화에 사용되는 변수를 정의한다.

## 분류 라벨

- **EXACT / IDENTITY** — 명시된 미시모델 아래 정확히 성립.
- **DEFINITION** — 선택한 수학적 정의.
- **ASSUMPTION** — 모델링 가정.
- **CONTROLLED APPROXIMATION** — 검증이 필요한 축약.
- **EMPIRICAL INPUT** — 측정되거나 외부에서 공급되는 재료량.

## 1. 미시 수직기하

| 기호 | 정의 | 물리적 의미 | 단위 | 분류 |
|---|---|---|---|---|
| $i,j,k$ | 정수 index | atom, bond 또는 neighbor index | 무차원 | DEFINITION |
| $N$ | $P_N$에서 사용하는 유한 spacing/atom 수 | finite-system size | 무차원 | DEFINITION |
| $x_i(t)$ | loading axis 방향 atomic coordinate | 축약 1D normal chain의 atom $i$ 위치 | 물리적으로 m; normalized simulation에서는 무차원 | DEFINITION |
| $a_i(t)$ | $x_{i+1}-x_i$ | 국부 수직 원자간격 | m | DEFINITION |
| $a$ | 연속 spacing-space coordinate | $P(a,t)$의 argument | m | DEFINITION |
| $a_0$ | reference/equilibrium spacing | 기준 lattice의 수직 원자간격 | m | EMPIRICAL INPUT 또는 calibration output |
| $\lambda_i$ | $a_i/a_0$ | 국부 normalized normal stretch | 무차원 | DEFINITION |
| $\lambda$ | $a/a_0$ | 연속 normalized spacing coordinate | 무차원 | DEFINITION |
| $\lambda_c$ | $\phi''(\lambda_c)=0$의 해 | 이상화된 국부 normal stability-loss stretch | 무차원 | stated LJ model 아래 EXACT |
| $a_c$ | $a_0\lambda_c$ | $\lambda_c$에 대응하는 physical critical normal spacing | m | DEFINITION |
| $R_i^{(k)}$ | $a_i+\cdots+a_{i+k-1}$ | reduced geometry에서 $i$ atom부터 $k$번째 neighbor까지 거리 | m | chain geometry 아래 EXACT / IDENTITY |
| $r$ | 연속 pair-distance coordinate | $v(r)$ 또는 $P_k(r,t)$의 argument | m | DEFINITION |

## 2. 수직 원자간격 확률변수

유한 $N$에서

$$
P_N(a,t)=\frac1N\sum_{i=1}^{N}\delta\!\left(a-a_i(t)\right)
$$

를 정의한다.

중심 thermodynamic-limit state는

$$
\boxed{P(a,t)=\lim_{N\to\infty}P_N(a,t)}
$$

이다.

| 기호 | 정의 | 물리적 의미 | 단위 | 분류 |
|---|---|---|---|---|
| $P_N(a,t)$ | finite empirical spacing density | 유한계의 국부 normal-spacing distribution | m$^{-1}$ | DEFINITION |
| $P(a,t)$ | thermodynamic-limit spacing density | 중심 distribution-valued structural state | m$^{-1}$ | DEFINITION |
| $P_k(r,t)$ | $R_i^{(k)}$의 density | 정확한 $k$-th-neighbor distance distribution | m$^{-1}$ | DEFINITION |
| $P^{*k}$ | $P$의 $k$-fold convolution | adjacent-spacing independence일 때만 $P_k$에 대한 approximation | m$^{-1}$ | 사용 시 CONTROLLED APPROXIMATION |
| $\bar a$ | $\int aP(a,t)\,da$ | 평균 수직 원자간격 | m | DEFINITION |
| $\operatorname{Var}(a)$ | $\int(a-\bar a)^2P(a,t)\,da$ | normal-spacing variance | m$^2$ | DEFINITION |
| $Q_c(t)$ | $\int_{a_c}^{\infty}P(a,t)\,da$ | candidate normal-instability spacing 위의 순간 probability mass | 무차원 | DEFINITION |
| $\delta(\cdot)$ | Dirac delta | empirical state density 정의에 쓰는 distribution | argument 역단위 | DEFINITION |

정규화는

$$
\boxed{\int P(a,t)\,da=1}
$$

이다.

## 3. 정확한 spacing-space 운동학

| 기호 | 정의 | 물리적 의미 | 단위 | 분류 |
|---|---|---|---|---|
| $v(a,t)$ | $\langle\dot a_i\mid a_i=a\rangle$ | spacing space의 conditional transport velocity | m/s | DEFINITION |
| $c$ | $\dot a$ | phase space의 local spacing velocity | m/s | DEFINITION |
| $F(a,c,t)$ | spacing과 spacing velocity의 joint density | $P$의 phase-space lift | s/m$^2$ | DEFINITION |
| $A(a,c,t)$ | $\langle\ddot a_i\mid a_i=a,\dot a_i=c\rangle$ | conditional normal-spacing acceleration | m/s$^2$ | DEFINITION |

정확한 kinematic continuity equation은

$$
\boxed{\partial_tP+\partial_a(Pv)=0}
$$

이다.

phase-space transport form은

$$
\boxed{\partial_tF+\partial_a(cF)+\partial_c(AF)=0}
$$

이다.

moment identity는

$$
\frac{d}{dt}\langle a^n\rangle=n\langle a^{n-1}v\rangle
$$

이고,

$$
\boxed{\frac{d}{dt}\operatorname{Var}(a)=2\operatorname{Cov}(a,v)}
$$

이다.

## 4. 수직하중 변수

| 기호 | 정의 | 물리적 의미 | 단위 | 분류 |
|---|---|---|---|---|
| $t$ | 시간 | physical 또는 nondimensional simulation time | 물리적으로 s | DEFINITION |
| $\sigma_n(t)$ | applied normal stress | loading axis 방향 외부 반복응력 | Pa | DEFINITION / experimental input history |
| $\sigma_m$ | mean normal stress | cycle mean | Pa | DEFINITION |
| $\sigma_a$ | normal stress amplitude | 반복응력 amplitude | Pa | DEFINITION |
| $\epsilon_n(t)$ | normal strain | chosen convention에서 $\sigma_n$과 work-conjugate strain | 무차원 | DEFINITION |
| $f$ | cyclic frequency | 초당 cycle 수 | Hz | DEFINITION |
| $\omega$ | $2\pi f$ | angular frequency | rad/s | DEFINITION |
| $T$ | $2\pi/\omega$ | loading period | s | DEFINITION |
| $N_{\rm cyc}$ | 완료한 loading cycle 수 | fatigue-cycle index | 무차원 | DEFINITION |
| $F_{\rm ext}(t)$ | finite chain/tester의 external normal force | end displacement와 conjugate인 force | 물리적으로 N | DEFINITION |
| $F_m$ | mean normal force | $F_{\rm ext}$ 평균 | N | DEFINITION |
| $F_a$ | normal-force amplitude | $F_{\rm ext}$ amplitude | N | DEFINITION |

사인 수직응력은

$$
\boxed{\sigma_n(t)=\sigma_m+\sigma_a\sin(\omega t)}
$$

이다.

## 5. Generalized Lennard-Jones 변수

활성 physical pair-potential baseline은

$$
v(r)=\varepsilon_{\rm LJ}
\left[
\left(\frac{\sigma_{\rm LJ}}{r}\right)^m
-
\left(\frac{\sigma_{\rm LJ}}{r}\right)^n
\right]
$$

이다.

| 기호 | 정의 | 물리적 의미 | 단위 | 분류 |
|---|---|---|---|---|
| $v(r)$ | generalized LJ pair energy | microscopic pair interaction energy | J | material potential로 채택할 때 ASSUMPTION |
| $\varepsilon_{\rm LJ}$ | LJ energy scale | pair-energy parameter | J 또는 eV | calibration input/output |
| $\sigma_{\rm LJ}$ | LJ length scale | characteristic distance parameter | m | calibration input/output |
| $m$ | repulsive exponent | short-range repulsive shape parameter | 무차원 | calibration output; 현재 12.19 |
| $n$ | attractive exponent | attractive shape parameter | 무차원 | 현재 baseline에서 6 |
| $\phi(\lambda)$ | 현재 1D code에서 쓰는 normalized generalized-LJ energy | dimensionless pair energy | 무차원 | DEFINITION |
| $\phi'(\lambda)$ | $d\phi/d\lambda$ | normalized tensile force coordinate | 무차원 | EXACT derivative |
| $\phi''(\lambda)$ | $d^2\phi/d\lambda^2$ | normalized tangent stiffness | 무차원 | EXACT derivative |

활성 normalized form은

$$
\boxed{
\phi(\lambda)=
\frac{\lambda^{-m}}{m(m-n)}
-
\frac{\lambda^{-n}}{n(m-n)}
}
$$

이고,

$$
\phi'(1)=0,
\qquad
\phi''(1)=1
$$

을 만족한다.

현재 exponent는

$$
m=12.19,
\qquad
n=6
$$

이다.

이상화된 normal stability condition은

$$
\boxed{\phi''(\lambda_c)=0}
$$

이고,

$$
\boxed{
\lambda_c=\left(\frac{m+1}{n+1}\right)^{1/(m-n)}
}
$$

이다.

현재 exponent에서

$$
\lambda_c\approx1.1077715386
$$

이다.

## 6. `theory/normal_lj_chain.py`의 chain-model 변수

| code field / 기호 | 의미 | 현재 code 단위 | 분류 |
|---|---|---|---|
| `repulsive_exponent` / $m$ | generalized-LJ repulsive exponent | 무차원 | model parameter |
| `attractive_exponent` / $n$ | generalized-LJ attractive exponent | 무차원 | model parameter |
| `mean_force` | mean dimensionless end force | 무차원 | input |
| `force_amplitude` | dimensionless end-force amplitude | 무차원 | input |
| `omega` / normalized $\omega^*$ | dimensionless angular loading frequency | 무차원 | input |
| `ramp_cycles` | smooth loading ramp에 쓰는 cycle 수 | cycle | numerical protocol parameter |
| `atoms` | finite chain atom 수 | count | CONTROLLED APPROXIMATION parameter |
| `dt` | velocity-Verlet time step | dimensionless time | numerical parameter |
| `cycles` | 최대 integration cycle 수 | cycle | numerical parameter |
| `record_stride` | 저장 sample 사이 integration step 수 | step | data-storage parameter |
| `runaway_spacing` | 큰 separation 이후 계산을 중지하는 spacing | dimensionless stretch | numerical stop이며 **crack criterion이 아님** |
| `period` | normalized simulation에서 $2\pi/\omega$ | dimensionless time | DEFINITION |
| `time` | 저장된 simulation time | dimensionless time | output |
| `force` | 저장된 external dimensionless force | 무차원 | output |
| `max_spacing` | instantaneous $\max_i\lambda_i$ | 무차원 | diagnostic output |
| `internal_energy` | kinetic + normalized LJ configurational energy | 무차원 | output |
| `external_work` | integrated external power | 무차원 | output |
| `cycle_mean_spacing` | cycle endpoint 평균 $\lambda_i$ | 무차원 | output |
| `cycle_variance_spacing` | cycle endpoint $\lambda_i$ variance | 무차원$^2$ | output |
| `cycle_max_spacing` | cycle endpoint 최대 $\lambda_i$ | 무차원 | output |
| `cycle_min_spacing` | cycle endpoint 최소 $\lambda_i$ | 무차원 | output |
| `cycle_snapshots` | 선택 cycle endpoint의 전체 $\lambda_i$ array | 무차원 | output |
| `first_instability` | 처음 $\max_i\lambda_i\ge\lambda_c$가 된 event | event record | stated LJ criterion 기반 diagnostic |
| `energy_balance_relative_error` | $|\Delta E_{\rm int}-W_{\rm ext}|/|W_{\rm ext}|$에 near-zero numerical protection 적용 | 무차원 | numerical verification metric |

현재 finite-chain potential은

$$
\boxed{V=\sum_i\phi(\lambda_i)}
$$

이다.

왼쪽 atom은 고정하고 오른쪽 끝 atom에 prescribed normal force를 가한다.

## 7. Stress/force mapping 변수

활성 normalized-chain mapping은

$$
f^*=\frac{\sigma_n}{E}
$$

이다.

| 기호 | 정의 | 의미 | 단위 | 분류 |
|---|---|---|---|---|
| $f^*$ | $\sigma_n/E$ | reduced chain에서 쓰는 dimensionless normal force/stress coordinate | 무차원 | current mapping 아래 DEFINITION |
| $E$ | reference Young's modulus | stress scale | Pa | EMPIRICAL INPUT; 현재 reference 69 GPa |
| $A_0$ | effective reference area | 1D force/stress conversion scale | m$^2$ | model calibration quantity |

현재 100 MPa reference amplitude는

$$
f_a^*=\frac{100\,\mathrm{MPa}}{69\,\mathrm{GPa}}
\approx1.44927536\times10^{-3}
$$

이다.

## 8. 에너지 hierarchy 변수

| 기호 | 정의 | 의미 | 단위 | 분류 |
|---|---|---|---|---|
| $V$ | finite chain의 pair-energy sum | finite-chain configurational energy | 물리적으로 J 또는 normalized | DEFINITION |
| $\mathcal U(t)$ | $\sum_k\int v(r)P_k(r,t)\,dr$ | chosen counting convention에서 distribution-level pair-potential energy per atom | J | stated pair-potential hierarchy 아래 EXACT |
| $E_{\rm int}$ | external loading potential을 제외한 kinetic + configurational energy | internal energy | J 또는 normalized | DEFINITION |
| $W_{\rm ext}$ | external power 시간적분 | external mechanical work | J 또는 normalized | DEFINITION |
| $A_H$ | $\oint\sigma_n\,d\epsilon_n$ | cycle당 normal stress-strain hysteresis work per volume | J/m$^3$ | DEFINITION |

정확한 pair-distance hierarchy는

$$
\boxed{\mathcal U(t)=\sum_{k=1}^{\infty}\int_0^\infty v(r)P_k(r,t)\,dr}
$$

이다.

## 9. 시간척도 변수

현재 1D mechanical scale은

$$
\boxed{t_0=\sqrt{\frac{m_{\rm Al}a_0}{EA_0}}}
$$

이다.

| 기호 | 정의 | 의미 | 단위 | 분류 |
|---|---|---|---|---|
| $t_0$ | 위 식 | current normalization의 atomic mechanical time scale | s | scaling 아래 DEFINITION |
| $m_{\rm Al}$ | aluminum atomic mass | inertia scale | kg | EMPIRICAL INPUT |
| $\omega^*$ | $\omega t_0=2\pi f t_0$ | dimensionless angular frequency | 무차원 | DEFINITION |
| $f_{\rm phys}$ | $\omega^*/(2\pi t_0)$ | normalized $\omega^*$에 대응하는 physical frequency | Hz | DEFINITION |

실험 frequency를 주장할 때는 반드시 $t_0$를 통해 $f$와 $\omega^*$를 명시적으로 연결한다.

## 10. Normal-opening initiation 변수

| 기호 | 정의 | 의미 | 단위 | 분류 |
|---|---|---|---|---|
| $\tau_c$ | stated mechanical normal-instability event가 처음 발생한 시간 | first-passage initiation time | s | DEFINITION |
| $S(t)$ | $\Pr(\tau_c>t)$ | survival probability | 무차원 | DEFINITION |
| $F_{\rm ci}(t)$ | $\Pr(\tau_c\le t)=1-S(t)$ | cumulative crack-initiation probability | 무차원 | DEFINITION |
| $h(t)$ | $-\dot S/S$ when defined | initiation hazard | s$^{-1}$ | DEFINITION |

활성 이론은 순간 tail quantity $Q_c(t)$와 cumulative first-passage quantity $F_{\rm ci}(t)$를 명확히 구분한다.
