# Milestone 3 — FCC normal-deformation generalized-LJ lattice sum

## Status

**Mechanics validation / calibration study. Not a fatigue-life prediction.**

The active mainline is normal deformation only. This milestone removes the largest geometric simplification in the current 1D nearest-neighbor chain by evaluating the generalized Lennard-Jones pair potential on a three-dimensional FCC Bravais lattice.

The purpose is not to add a new phenomenological mechanism. It is to determine, with as few hidden assumptions as possible, which normal mechanical properties can be produced by the same fixed pair potential and which cannot.

## 1. FCC lattice geometry

Let $a_{\rm lat}$ be the conventional FCC lattice constant. Nonzero lattice vectors can be represented as

$$
\mathbf R=\frac{a_{\rm lat}}{2}(i,j,k),
$$

where $i,j,k\in\mathbb Z$ and $i+j+k$ is even.

The reference atomic volume is

$$
\boxed{\Omega_0=\frac{a_{\rm lat}^3}{4}.}
$$

For a homogeneous deformation gradient $\mathbf F$, the pair-potential energy per atom is

$$
\boxed{
U(\mathbf F)
=\frac12\sum_{\mathbf R\neq0}
v\!\left(|\mathbf F\mathbf R|\right).
}
$$

The factor $1/2$ prevents double counting of pairs.

### Classification

- **EXACT / IDENTITY under the stated central-pair model:** the infinite FCC lattice sum.
- **CONTROLLED APPROXIMATION:** truncating the infinite lattice sum at a finite spherical cutoff.
- **ASSUMPTION:** representing Al by a fixed central generalized-LJ pair potential.

## 2. Generalized Lennard-Jones interaction

The active pair law is

$$
\boxed{
v(r)=\varepsilon_{\rm LJ}
\left[
\left(\frac{\sigma_{\rm LJ}}{r}\right)^m
-
\left(\frac{\sigma_{\rm LJ}}{r}\right)^n
\right],
}
$$

with the inherited normal-mainline exponents

$$
m=12.19,
\qquad
n=6.
$$

The parameters remain fixed during cycling. No damage-dependent potential is permitted.

For a chosen $(m,n)$ and a chosen reference lattice constant, the zero-pressure condition determines the dimensionless ratio

$$
q=\frac{\sigma_{\rm LJ}}{a_{\rm lat}}.
$$

With the converged FCC lattice sum used here,

$$
\boxed{\sigma_{\rm LJ}\approx2.62721\ \text{Å}.}
$$

## 3. [001] normal loading with free transverse relaxation

For a high-symmetry [001] normal test, use

$$
\boxed{
\mathbf F
=\operatorname{diag}(\lambda_t,\lambda_t,\lambda_n).
}
$$

Here $\lambda_n$ is the imposed normal stretch and $\lambda_t$ is the common transverse stretch.

A traction-free homogeneous transverse boundary condition is represented by minimizing the energy with respect to $\lambda_t$ at fixed $\lambda_n$:

$$
\boxed{
\frac{\partial U}{\partial\lambda_t}=0.
}
$$

The axial engineering / nominal stress used in the numerical study is

$$
\boxed{
P_{33}
=\frac1{\Omega_0}
\frac{\partial U}{\partial\lambda_n}
}
$$

with $\lambda_t$ relaxed at the same $\lambda_n$.

This calculation has no damping, no slip variable, no plasticity law, and no fatigue-damage variable.

## 4. Directional experimental reference

The external Al benchmark used for this validation is the experimental-property table reproduced in the 2013 Journal of Applied Crystallography paper *The interpretation of polycrystalline coherent inelastic neutron scattering from aluminium*, DOI `10.1107/S0021889813023728`.

Values used here are

$$
a_{\rm lat}=4.0495\ \text{Å},
$$

$$
E_{\rm coh}=3.43\ \text{eV/atom},
$$

and

$$
C_{11}=107\ \text{GPa},\qquad
C_{12}=61\ \text{GPa},\qquad
C_{44}=29\ \text{GPa}.
$$

For a cubic crystal loaded along [001] with free transverse contraction,

$$
\boxed{
E_{[001]}
=\frac{(C_{11}-C_{12})(C_{11}+2C_{12})}
{C_{11}+C_{12}}
}
$$

and

$$
\boxed{
\nu_{[001]}=\frac{C_{12}}{C_{11}+C_{12}}.
}
$$

Therefore the directional reference values are

$$
\boxed{E_{[001]}\approx62.7024\ \text{GPa}},
$$

$$
\boxed{\nu_{[001]}\approx0.363095.}
$$

This corrects an earlier temporary comparison that used the tabulated isotropic/aggregate Young modulus directly as if it were the [001] single-crystal modulus.

## 5. Result A — fitting the cohesive energy makes the normal lattice far too stiff

If $a_{\rm lat}$ and the experimental cohesive energy $3.43$ eV/atom are imposed, the required LJ energy scale is

$$
\varepsilon_{\rm LJ}\approx1.56683\ \text{eV}.
$$

The resulting relaxed [001] Young modulus is

$$
\boxed{E_{[001]}^{\rm LJ}\approx220.466\ \text{GPa}.}
$$

This is about $3.52$ times the directional reference modulus.

The calculated Poisson ratio is

$$
\nu_{[001]}^{\rm LJ}\approx0.363410.
$$

Thus the geometric normal contraction ratio is already close to the experimental [001] value, but the energy scale required by cohesion makes the lattice much too stiff.

## 6. Result B — fitting the normal elastic modulus gives excellent normal elastic constants

Instead calibrate only the LJ energy scale to the experimental [001] modulus while keeping

- the FCC lattice geometry,
- $m=12.19$,
- $n=6$,
- the zero-pressure lattice spacing,
- and the same fixed pair law.

The required energy scale is

$$
\boxed{\varepsilon_{\rm LJ}\approx0.445621\ \text{eV}.}
$$

The predicted directional modulus and Poisson ratio are

$$
E_{[001]}^{\rm LJ}\approx62.7024\ \text{GPa},
$$

$$
\nu_{[001]}^{\rm LJ}\approx0.363410.
$$

The corresponding cubic elastic constants are

$$
\boxed{C_{11}^{\rm LJ}\approx107.169\ \text{GPa}},
$$

$$
\boxed{C_{12}^{\rm LJ}\approx61.180\ \text{GPa}}.
$$

These are strikingly close to the external normal-elastic reference values $107$ GPa and $61$ GPa.

This is an important positive result for the normal-only research direction: the inherited $(m,n)=(12.19,6)$ generalized-LJ shape reproduces the small-strain [001] normal elastic sector very well after one energy-scale calibration.

## 7. Exact central-pair limitation — Cauchy relation

For a cubic crystal at zero pressure with a central pair potential, the Cauchy relation requires

$$
\boxed{C_{12}=C_{44}.}
$$

The numerical FCC calculation recovers this to the finite-difference/lattice-cutoff error:

$$
C_{12}^{\rm LJ}\approx61.17959\ \text{GPa},
$$

$$
C_{44}^{\rm LJ}\approx61.17953\ \text{GPa}.
$$

The difference is only about

$$
0.060\ \text{MPa}.
$$

Real Al instead has approximately

$$
C_{12}\approx61\ \text{GPa},
\qquad
C_{44}\approx29\ \text{GPa}.
$$

Therefore a fixed central pair potential cannot be a complete three-dimensional quantitative Al potential.

This is not a numerical failure. It is a mathematical structural limitation of the model class. The same limitation is explicitly discussed in DOI `10.1107/S0021889813023728`.

For the present project, shear physics remains archived in `libraries/shear/`; the active question is whether the pair model remains sufficiently accurate for the **normal sector** before a mathematically required many-body correction is introduced.

## 8. Result C — normal ideal-strength prediction without fitting the peak

With the energy scale fitted only to $E_{[001]}$, the relaxed FCC stress-strain curve has a maximum engineering stress of approximately

$$
\boxed{\sigma_{\rm ideal}^{\rm LJ}\approx9.045\ \text{GPa}}
$$

at an engineering strain of about

$$
\boxed{\epsilon_n\approx0.25.}
$$

The ideal tensile peak was **not** used in the calibration.

For comparison, a first-principles [001] calculation for pure Al reported an ideal tensile strength of approximately $10.63$ GPa in Metals 12, 2143 (2022), DOI `10.3390/met12122143`.

The agreement in strength scale is useful but must not be overstated: the methods, temperature assumptions, relaxation protocols, and peak strains differ.

## 9. Result D — stiffness and cohesive energy cannot both be correct with this fixed pair law

The price of the excellent normal elastic calibration is a predicted cohesive energy of only

$$
\boxed{E_{\rm coh}^{\rm LJ}\approx0.9755\ \text{eV/atom}},
$$

compared with the experimental scale

$$
3.43\ \text{eV/atom}.
$$

Conversely, fitting cohesion produces $E_{[001]}\approx220.5$ GPa.

Therefore the current pair law has two very different possible interpretations:

1. **effective normal-mechanics potential** — calibrate the tangent normal mechanics and obtain realistic $C_{11}$, $C_{12}$ and an order-correct ideal normal strength;
2. **thermodynamic cohesive potential** — calibrate the separation energy, but then the normal elastic stiffness becomes badly wrong.

The same fixed pair potential cannot currently serve both roles quantitatively.

## 10. Consequence for the probability fatigue theory

This matters directly for the next step.

The current normal probability theory wants to use

$$
P(a,t)
$$

and ultimately first-passage or thermal-ensemble reasoning for normal opening.

If a thermal probability contains a Boltzmann factor such as

$$
\exp\!\left(-\frac{\Delta U}{k_B T}\right),
$$

then the absolute energy barrier matters exponentially. A potential whose cohesive energy is wrong by a factor of roughly $3.5$ cannot be used for a quantitative thermal escape rate merely because its small-strain normal stiffness is accurate.

Therefore **thermal activation will not be added yet**.

The next theoretical problem is now sharply defined:

$$
\boxed{
\text{retain the successful LJ normal mechanics}
\quad+\quad
\text{derive the minimum physically justified energy correction needed for cohesion}
}
$$

without introducing a fitted fatigue law or cycle-dependent LJ parameters.

A many-body electronic/cohesive contribution is a candidate only because the central-pair model has now been mathematically falsified as a complete 3D energy model; it is not introduced merely to improve a curve fit.

## 11. Numerical convergence

The lattice sum was repeated with spherical cutoffs of 5, 6, 8, 10, 12, and 15 conventional lattice constants. For the cohesive-energy calibration, the predicted $E_{[001]}$ changes from about $220.617$ GPa at cutoff 5 to $220.466$ GPa at cutoff 15.

The normal-fit ideal-strength result changes by only a few MPa over the same range.

Thus the stiffness/cohesion conflict is much larger than the finite-cutoff numerical uncertainty.

## 12. Files

- `theory/fcc_normal_lj.py`
- `simulations/run_fcc_normal_lj.py`
- `tests/test_fcc_normal_lj.py`
- `results/data/fcc_normal_lj_summary.json`
- `results/data/fcc_normal_lj_stress_strain.csv`
- `results/figures/fcc_lj_001_stress_strain.svg`
- `results/figures/fcc_lj_elastic_constants.svg`
- `results/figures/fcc_lj_cohesion_conflict.svg`
- `results/figures/fcc_lj_cutoff_convergence.svg`

---

# 한국어 번역 — FCC 수직변형 generalized-LJ lattice sum

## 상태

**역학 검증 / calibration 연구다. 피로수명 예측은 아직 아니다.**

활성 mainline은 수직변형만 다룬다. 이번 단계에서는 기존 1D 최근접이웃 chain의 가장 큰 기하학적 근사를 줄이기 위해 generalized Lennard-Jones pair potential을 3차원 FCC Bravais lattice 전체에 직접 합산한다.

목적은 새로운 경험적 메커니즘을 추가하는 것이 아니다. 같은 고정 pair potential 하나가 어떤 수직 물성을 동시에 만들 수 있고 어떤 물성에서는 수학적으로 실패하는지를 확인하는 것이다.

## 1. FCC lattice geometry

conventional FCC lattice constant를 $a_{\rm lat}$라 한다. 0이 아닌 lattice vector는

$$
\mathbf R=\frac{a_{\rm lat}}{2}(i,j,k)
$$

로 표현할 수 있고 $i,j,k\in\mathbb Z$, $i+j+k$는 짝수다.

원자 하나당 기준부피는

$$
\boxed{\Omega_0=\frac{a_{\rm lat}^3}{4}}
$$

이다.

균질 deformation gradient $\mathbf F$에 대해 원자 하나당 pair-potential energy는

$$
\boxed{
U(\mathbf F)
=\frac12\sum_{\mathbf R\neq0}
v\!\left(|\mathbf F\mathbf R|\right)
}
$$

이다.

$1/2$는 pair double counting을 제거한다.

### 분류

- **EXACT / IDENTITY under the stated central-pair model:** 무한 FCC lattice sum.
- **CONTROLLED APPROXIMATION:** 무한 합을 유한 spherical cutoff로 자르는 것.
- **ASSUMPTION:** Al을 고정 central generalized-LJ pair potential로 나타내는 것.

## 2. Generalized Lennard-Jones interaction

활성 pair law는

$$
\boxed{
v(r)=\varepsilon_{\rm LJ}
\left[
\left(\frac{\sigma_{\rm LJ}}{r}\right)^m
-
\left(\frac{\sigma_{\rm LJ}}{r}\right)^n
\right]
}
$$

이고 기존 normal-mainline에서 사용한 exponent

$$
m=12.19,\qquad n=6
$$

을 그대로 사용한다.

parameter는 cycle에 따라 변하지 않는다. damage-dependent potential은 허용하지 않는다.

주어진 $(m,n)$과 기준 lattice constant에 대해 zero-pressure condition이

$$
q=\frac{\sigma_{\rm LJ}}{a_{\rm lat}}
$$

를 결정한다.

이번 converged FCC lattice sum에서는

$$
\boxed{\sigma_{\rm LJ}\approx2.62721\ \text{Å}}
$$

가 나온다.

## 3. [001] 수직인장과 자유 횡수축

고대칭 [001] normal test에 대해

$$
\boxed{
\mathbf F=\operatorname{diag}(\lambda_t,\lambda_t,\lambda_n)
}
$$

를 사용한다.

$\lambda_n$은 imposed normal stretch이고 $\lambda_t$는 두 횡방향의 공통 stretch다.

균질한 횡방향 traction-free condition은 고정 $\lambda_n$에서 $\lambda_t$에 대해 에너지를 최소화하여

$$
\boxed{\frac{\partial U}{\partial\lambda_t}=0}
$$

으로 둔다.

수치계산에 사용한 axial engineering / nominal stress는

$$
\boxed{
P_{33}
=\frac1{\Omega_0}
\frac{\partial U}{\partial\lambda_n}
}
$$

이고 같은 $\lambda_n$에서 $\lambda_t$를 relaxation한다.

여기에는 damping, slip variable, plasticity law, fatigue-damage variable이 없다.

## 4. 방향성 실험 reference

외부 Al benchmark는 2013 Journal of Applied Crystallography 논문 *The interpretation of polycrystalline coherent inelastic neutron scattering from aluminium*, DOI `10.1107/S0021889813023728`에 정리된 experimental-property table을 사용했다.

사용값은

$$
a_{\rm lat}=4.0495\ \text{Å},
$$

$$
E_{\rm coh}=3.43\ \text{eV/atom},
$$

$$
C_{11}=107\ \text{GPa},\qquad
C_{12}=61\ \text{GPa},\qquad
C_{44}=29\ \text{GPa}
$$

이다.

cubic crystal을 [001] 방향으로 인장하고 횡수축이 자유라면

$$
\boxed{
E_{[001]}
=\frac{(C_{11}-C_{12})(C_{11}+2C_{12})}
{C_{11}+C_{12}}
}
$$

이고

$$
\boxed{\nu_{[001]}=\frac{C_{12}}{C_{11}+C_{12}}}
$$

이다.

따라서 방향성 reference는

$$
\boxed{E_{[001]}\approx62.7024\ \text{GPa}},
$$

$$
\boxed{\nu_{[001]}\approx0.363095}
$$

가 된다.

이 값은 이전 중간계산에서 표의 aggregate Young modulus를 [001] single-crystal modulus처럼 직접 사용했던 부분을 수정한 것이다.

## 5. 결과 A — cohesive energy를 맞추면 normal lattice가 지나치게 단단해진다

$a_{\rm lat}$와 experimental cohesive energy $3.43$ eV/atom을 동시에 강제하면 필요한 LJ energy scale은

$$
\varepsilon_{\rm LJ}\approx1.56683\ \text{eV}
$$

이다.

이때 relaxed [001] Young modulus는

$$
\boxed{E_{[001]}^{\rm LJ}\approx220.466\ \text{GPa}}
$$

가 된다.

방향성 reference modulus의 약 $3.52$배다.

계산된 Poisson ratio는

$$
\nu_{[001]}^{\rm LJ}\approx0.363410
$$

이다.

즉 수직변형의 기하학적 횡수축 비율은 이미 실험 [001] 값과 가깝지만 cohesion을 맞추는 energy scale이 lattice stiffness를 너무 크게 만든다.

## 6. 결과 B — normal elastic modulus를 맞추면 수직 elastic constant는 매우 잘 맞는다

이번에는 FCC geometry, $m=12.19$, $n=6$, zero-pressure lattice spacing, 고정 pair law를 그대로 둔 상태에서 energy scale 하나만 experimental $E_{[001]}$에 맞춘다.

필요한 energy scale은

$$
\boxed{\varepsilon_{\rm LJ}\approx0.445621\ \text{eV}}
$$

이다.

예측된 방향성 modulus와 Poisson ratio는

$$
E_{[001]}^{\rm LJ}\approx62.7024\ \text{GPa},
$$

$$
\nu_{[001]}^{\rm LJ}\approx0.363410
$$

이다.

이에 해당하는 cubic elastic constant는

$$
\boxed{C_{11}^{\rm LJ}\approx107.169\ \text{GPa}},
$$

$$
\boxed{C_{12}^{\rm LJ}\approx61.180\ \text{GPa}}
$$

이다.

외부 normal-elastic reference $107$ GPa, $61$ GPa와 매우 가깝다.

이것은 normal-only 연구방향에 상당히 중요한 긍정적 결과다. 기존 $(m,n)=(12.19,6)$ generalized-LJ shape는 energy scale 하나를 calibration하면 small-strain [001] normal elastic sector를 매우 잘 재현한다.

## 7. Central-pair의 정확한 한계 — Cauchy relation

zero pressure의 cubic crystal에서 central pair potential은

$$
\boxed{C_{12}=C_{44}}
$$

라는 Cauchy relation을 만족해야 한다.

수치 FCC 계산도 finite-difference/lattice-cutoff 오차 범위에서 이를 그대로 복원한다.

$$
C_{12}^{\rm LJ}\approx61.17959\ \text{GPa},
$$

$$
C_{44}^{\rm LJ}\approx61.17953\ \text{GPa}.
$$

차이는 약

$$
0.060\ \text{MPa}
$$

뿐이다.

반면 실제 Al은 대략

$$
C_{12}\approx61\ \text{GPa},\qquad
C_{44}\approx29\ \text{GPa}
$$

다.

따라서 고정 central pair potential은 3차원 Al 전체를 정량적으로 완전하게 표현할 수 없다.

이것은 numerical failure가 아니라 model class 자체의 수학적 구조 제한이다. DOI `10.1107/S0021889813023728`에서도 같은 Cauchy violation 문제가 명시적으로 논의된다.

현재 프로젝트에서는 shear physics는 `libraries/shear/`에 보존하고, active mainline에서는 **normal sector에서 pair model이 어디까지 정확한가**를 우선 사용한다. many-body correction은 수학적으로 필요하다고 확인되는 지점에서만 도입한다.

## 8. 결과 C — peak를 fitting하지 않은 normal ideal-strength 예측

$E_{[001]}$만 맞춘 energy scale을 사용하면 relaxed FCC stress-strain curve의 최대 engineering stress는 약

$$
\boxed{\sigma_{\rm ideal}^{\rm LJ}\approx9.045\ \text{GPa}}
$$

이고 engineering strain은 약

$$
\boxed{\epsilon_n\approx0.25}
$$

이다.

ideal tensile peak는 calibration에 사용하지 않았다.

비교를 위해 Metals 12, 2143 (2022), DOI `10.3390/met12122143`의 first-principles [001] 계산은 pure Al의 ideal tensile strength를 약 $10.63$ GPa로 보고했다.

strength scale이 비슷한 것은 유용하지만 과대해석하면 안 된다. 계산방법, 온도조건, relaxation protocol, peak strain이 서로 다르다.

## 9. 결과 D — stiffness와 cohesive energy는 이 고정 pair law 하나로 동시에 맞지 않는다

normal elasticity를 잘 맞추는 대신 예측 cohesive energy는

$$
\boxed{E_{\rm coh}^{\rm LJ}\approx0.9755\ \text{eV/atom}}
$$

밖에 되지 않는다.

experimental scale

$$
3.43\ \text{eV/atom}
$$

보다 크게 작다.

반대로 cohesion을 맞추면 $E_{[001]}\approx220.5$ GPa가 된다.

따라서 현재 pair law는 두 가지 역할을 동시에 정량적으로 수행할 수 없다.

1. **effective normal-mechanics potential** — tangent normal mechanics를 맞추면 현실적인 $C_{11}$, $C_{12}$와 order-correct ideal normal strength를 얻는다.
2. **thermodynamic cohesive potential** — separation energy를 맞추면 normal elastic stiffness가 크게 틀어진다.

## 10. 확률 피로이론에 주는 의미

이 결과는 다음 단계와 직접 연결된다.

현재 normal probability theory는

$$
P(a,t)
$$

에서 시작해 결국 normal-opening first passage나 thermal ensemble을 다루려 한다.

thermal probability에

$$
\exp\!\left(-\frac{\Delta U}{k_BT}\right)
$$

같은 항이 들어가면 absolute energy barrier가 지수적으로 중요하다. cohesive energy가 약 $3.5$배 틀린 potential을 small-strain stiffness가 잘 맞는다는 이유만으로 quantitative thermal escape rate에 그대로 쓰면 안 된다.

따라서 **thermal activation은 아직 넣지 않는다.**

다음 이론문제는 이제 아주 명확하다.

$$
\boxed{
\text{성공적인 LJ normal mechanics는 유지}
\quad+\quad
\text{cohesion을 복원하는 최소한의 물리적 energy correction을 유도}
}
$$

이다.

many-body electronic/cohesive contribution은 pair model이 complete 3D energy model로 수학적으로 부족하다는 것이 확인됐기 때문에 후보가 되는 것이지, curve fitting을 잘하기 위해 임의로 추가하는 것이 아니다.

## 11. 수치 convergence

lattice sum을 conventional lattice constant 기준 cutoff 5, 6, 8, 10, 12, 15에서 반복했다. cohesive-energy calibration에서 $E_{[001]}$은 cutoff 5의 약 $220.617$ GPa에서 cutoff 15의 약 $220.466$ GPa로 수렴한다.

normal-fit ideal-strength도 같은 범위에서 수 MPa 정도만 변한다.

따라서 stiffness/cohesion conflict는 finite-cutoff numerical uncertainty보다 압도적으로 크다.

## 12. 파일

- `theory/fcc_normal_lj.py`
- `simulations/run_fcc_normal_lj.py`
- `tests/test_fcc_normal_lj.py`
- `results/data/fcc_normal_lj_summary.json`
- `results/data/fcc_normal_lj_stress_strain.csv`
- `results/figures/fcc_lj_001_stress_strain.svg`
- `results/figures/fcc_lj_elastic_constants.svg`
- `results/figures/fcc_lj_cohesion_conflict.svg`
- `results/figures/fcc_lj_cutoff_convergence.svg`
