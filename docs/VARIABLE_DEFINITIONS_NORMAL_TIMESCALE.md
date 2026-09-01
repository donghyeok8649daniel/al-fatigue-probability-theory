# Variable Definitions — Normal-LJ time-scale analysis

This file defines variables introduced in the active normal-only time-scale analysis.

| Symbol | Definition | Physical meaning | Unit | Classification |
|---|---|---|---|---|
| $t_0$ | $\sqrt{M_{\rm Al}a_0/(EA_0)}$ | Atomic time scale of the normalized normal chain | s | DEFINITION under stated normalization |
| $c_0$ | $a_0/t_0$ | Lattice-wave speed scale implied by the normalized chain | m/s | DERIVED QUANTITY |
| $L$ | Number of moving atoms in the fixed-free linearized chain | Chain size used in modal analysis | dimensionless | DEFINITION |
| $q_j$ | $(2j-1)\pi/(2L+1)$ | Allowed wave number of mode $j$ | rad | EXACT for linearized fixed-free chain |
| $\omega_j^*$ | $2\sqrt{k^*}\sin(q_j/2)$ | Dimensionless angular frequency of normal mode $j$ | dimensionless | EXACT for linearized chain |
| $\omega_{\min}^*$ | Lowest value of $\omega_j^*$ | Slowest conservative normal mode | dimensionless | DERIVED QUANTITY |
| $f_{\rm loc}$ | $\sqrt{\phi''(\lambda)}/(2\pi t_0)$ | Local small-oscillation frequency at stretch $\lambda$ | Hz | CONTROLLED APPROXIMATION because of local linearization |
| $\lambda_{100}$ | Stable root of $\phi'(\lambda)=100\,\mathrm{MPa}/E$ | Homogeneous normal stretch at 100 MPa in the 1D mapping | dimensionless | DERIVED QUANTITY |
| $\phi'''(\lambda_c)$ | Third derivative of normalized LJ energy at instability | Slope of tangent stiffness near the LJ inflection | dimensionless | EXACT under normalized LJ model |
| $\Delta\lambda_{20}$ | Linearized estimate of $\lambda_c-\lambda$ required for a 20 Hz local mode | Critical-softening proximity diagnostic | dimensionless | CONTROLLED APPROXIMATION |

The 20 Hz target frequency is an experimental loading scale, not a microscopic parameter inserted into the potential.

---

# 한국어 번역 — Normal-LJ 시간척도 변수정의

이 문서는 활성 normal-only 시간척도 분석에서 새로 도입된 변수를 정의한다.

| 기호 | 정의 | 물리적 의미 | 단위 | 분류 |
|---|---|---|---|---|
| $t_0$ | $\sqrt{M_{\rm Al}a_0/(EA_0)}$ | 정규화된 normal chain의 원자 시간척도 | s | 명시한 normalization 아래 DEFINITION |
| $c_0$ | $a_0/t_0$ | 정규화된 chain이 암시하는 lattice-wave speed scale | m/s | DERIVED QUANTITY |
| $L$ | fixed-free 선형화 chain의 moving atom 수 | modal analysis에 사용하는 chain 크기 | 무차원 | DEFINITION |
| $q_j$ | $(2j-1)\pi/(2L+1)$ | mode $j$의 허용 wave number | rad | 선형화 fixed-free chain에서 EXACT |
| $\omega_j^*$ | $2\sqrt{k^*}\sin(q_j/2)$ | normal mode $j$의 무차원 각주파수 | 무차원 | 선형화 chain에서 EXACT |
| $\omega_{\min}^*$ | $\omega_j^*$의 최솟값 | 가장 느린 conservative normal mode | 무차원 | DERIVED QUANTITY |
| $f_{\rm loc}$ | $\sqrt{\phi''(\lambda)}/(2\pi t_0)$ | stretch $\lambda$에서 local small-oscillation frequency | Hz | local linearization 때문에 CONTROLLED APPROXIMATION |
| $\lambda_{100}$ | $\phi'(\lambda)=100\,\mathrm{MPa}/E$의 stable root | 1D mapping에서 100 MPa homogeneous normal stretch | 무차원 | DERIVED QUANTITY |
| $\phi'''(\lambda_c)$ | instability에서 normalized LJ energy의 3차 미분 | LJ inflection 부근 tangent stiffness 변화율 | 무차원 | normalized LJ model에서 EXACT |
| $\Delta\lambda_{20}$ | 20 Hz local mode에 필요한 $\lambda_c-\lambda$의 선형화 추정값 | critical-softening proximity 진단값 | 무차원 | CONTROLLED APPROXIMATION |

20 Hz는 실험적 loading scale이며 potential에 삽입한 microscopic parameter가 아니다.
